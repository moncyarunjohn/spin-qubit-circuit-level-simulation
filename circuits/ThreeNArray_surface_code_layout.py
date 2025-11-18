from typing import Set, List, Dict
import math
from collections import defaultdict

def generate_rotated_surface_code_circuit_layout(distance: int = None,
                                                 x_distance: int = None,
                                                 z_distance: int = None):
    """
    Generate the layout of a rotated surface code with the given distance.

    Args:
        distance (int): Distance of the surface code.

    Returns:
        Tuple: A tuple containing the following elements (x_observable, z_observable, data_coords, x_measure_coords, 
        z_measure_coords, q2p, p2q, data_qubits, x_measurement_qubits, measurement_qubits, cnot_targets, measure_coord_to_order,
        data_coord_to_order, z_order):
            - x_observable (List[complex]): List of complex numbers representing the x-observable qubits.
            - z_observable (List[complex]): List of complex numbers representing the z-observable qubits.
            - data_coords (Set[complex]): Set of complex numbers representing the data qubits.
            - x_measure_coords (Set[complex]): Set of complex numbers representing the x-measurement qubits.
            - z_measure_coords (Set[complex]): Set of complex numbers representing the z-measurement qubits.
            - data_snake (List[complex]): List of complex numbers representing the data qubits in the snake order.
            - check_snake (List[complex]): List of complex numbers representing the check qubits in the snake order.
            - x_measure_index (List[int]): List of indices for the x-measurement qubits.
            - z_measure_index (List[int]): List of indices for the z-measurement qubits.
            - X_Map (Dict[list[int]]): Dictionary with key:the difference in position of the entangling qubits (data and X check) and value is the control and target qubits arranged in alternating order.
            - Z_Map (Dict[list[int]]): Dictionary with key:the difference in position of the entangling qubits (data and Z check) and value is the control and target qubits arranged in alternating order.
    """
    if distance is not None and x_distance is None and z_distance is None:
        # check if distance is at least 3
        if distance < 3:
            raise ValueError("Distance must be at least 3.")
        # x and z distance treated similarly for now.
        x_distance = distance
        z_distance = distance
    elif x_distance is not None and z_distance is not None and distance is None:
        if x_distance < 3 or z_distance < 3:
            raise ValueError("x_distance and z_distance must be at least 3.")
    else:
        raise ValueError("Exactly one of distance or (x_distance and z_distance) must be provided.")
        

    # Place data qubits
    data_coords: Set[complex] = set()
    x_observable: List[complex] = []
    z_observable: List[complex] = []
    for x in [i + 0.5 for i in range(z_distance)]:
        for y in [i + 0.5 for i in range(x_distance)]: # two different usages of for loops seen here
            q = x * 2 + y * 2 * 1j
            data_coords.add(q)
            if y == 0.5:
                x_observable.append(q) # z-observable qubits are along the y=0.5 line? 
                ### Have switched x and z
            if x == 0.5:
                z_observable.append(q) # x-observable qubits are along the x=0.5 line?
                ### Have switched x and z

    # Place measurement qubits.
    x_measure_coords: Set[complex] = set()
    z_measure_coords: Set[complex] = set()
    for x in range(z_distance + 1):
        for y in range(x_distance + 1):
            q = x * 2 + y * 2j
            on_boundary_1 = x == 0 or x == z_distance
            on_boundary_2 = y == 0 or y == x_distance
            parity = (int(round(x)) % 2) != (int(round(y)) % 2)
            if on_boundary_2 and parity:
                ### have switched boundary_1 and boundary_2
                continue
            if on_boundary_1 and not parity:
                ### have switched boundary_1 and boundary_2
                continue
            if parity:
                x_measure_coords.add(q)
                ### Think this should remain
            else:
                z_measure_coords.add(q)
                ### Think this should remain as well


    data_snake=sorted(sorted(data_coords, key=lambda c: (-1)**(((c.imag - c.real) % 4)/2) * c.real,reverse=True),key=lambda c: c.real-c.imag)
    # i=0
    # for c in data_snake:
    #     i+=1
    #     data_snake_map.append((c,i))

    X_observable_index: List[int] = []
    Z_observable_index: List[int] = []

    for c in x_observable:
        X_observable_index.append(data_snake.index(c))

    for c in z_observable:
        Z_observable_index.append(data_snake.index(c))

    ## sorted preserves the order of equal elements. So we can sort by one criterion and then by another.
    ## Here, we are sorting first by the north east line (which incidently alternates to opposite between every adjacent northwest line.
    ## hence the (-1)**((c.imag- c.real)%2) factor.)
    ## Reverse because the biggest northeast line (which is also a diagonal) has to be in opposite/descending direction/order.
    ## Then we sort the northwest lines, which is c.real-c.imag. So, these lines (northeast directional) are going from the northwest corner to the southeast corner.
    ##
    check_snake=sorted(sorted(z_measure_coords | x_measure_coords, key=lambda c: (-1)**(((c.imag - c.real) % 4)/2) * c.real),key = lambda c: c.real-c.imag)
    ## Same algorithm for the measurement qubits. Except that the directions of the northeast lines are opposite to those of the data qubits (while their order is he same: from northwest to southeast).

    # i=0
    # for c in check_snake:
    #     i+=1
    #     check_snake_map.append((c,i))

    measurement_qubits:List[int] = list(range(x_distance*z_distance - 1))
    data_qubits:List[int] = list(range(x_distance*z_distance - 1,2*x_distance*z_distance - 1))

    x_measure_index: List[int] = []
    z_measure_index: List[int] = []

    for c in x_measure_coords:
        x_measure_index.append(check_snake.index(c))

    for c in z_measure_coords:
        z_measure_index.append(check_snake.index(c))

    # X_Map: List[tuple[int,int,int]] = []
    X_Map_2N = defaultdict(list)
    X_Map = defaultdict(list)
    ## This map is from X check qubit index to data qubit index.

    for c in x_measure_coords:
        on_boundary_11 = c.real == 0
        on_boundary_12 = c.real == 2*z_distance
        if on_boundary_11:
            X_Map_2N[check_snake.index(c)-data_snake.index(c+1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1+1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1+1j))
            X_Map_2N[check_snake.index(c)-data_snake.index(c+1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1-1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1-1j))
        elif on_boundary_12:
            X_Map_2N[check_snake.index(c)-data_snake.index(c-1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1+1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1+1j))
            X_Map_2N[check_snake.index(c)-data_snake.index(c-1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1-1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1-1j))
        else:
            X_Map_2N[check_snake.index(c)-data_snake.index(c+1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1+1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1+1j))
            X_Map_2N[check_snake.index(c)-data_snake.index(c-1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1+1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1+1j))
            X_Map_2N[check_snake.index(c)-data_snake.index(c+1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1-1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1-1j))
            X_Map_2N[check_snake.index(c)-data_snake.index(c-1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1-1j)))
            X_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1-1j))


    Z_Map_2N = defaultdict(list)
    ## This map is from Z check qubit index to data qubit index.
    Z_Map = defaultdict(list)

    for c in z_measure_coords:
        on_boundary_21 = c.imag == 0
        on_boundary_22 = c.imag == 2*x_distance
        if on_boundary_21:
            Z_Map_2N[check_snake.index(c)-data_snake.index(c+1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1+1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1+1j))
            Z_Map_2N[check_snake.index(c)-data_snake.index(c-1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1+1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1+1j))
        elif on_boundary_22:
            Z_Map_2N[check_snake.index(c)-data_snake.index(c+1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1-1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1-1j))
            Z_Map_2N[check_snake.index(c)-data_snake.index(c-1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1-1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1-1j))
        else:
            Z_Map_2N[check_snake.index(c)-data_snake.index(c+1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1+1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1+1j))
            Z_Map_2N[check_snake.index(c)-data_snake.index(c-1+1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1+1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1+1j))
            Z_Map_2N[check_snake.index(c)-data_snake.index(c+1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c+1-1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c+1-1j))
            Z_Map_2N[check_snake.index(c)-data_snake.index(c-1-1j)].append((check_snake.index(c),z_distance*x_distance-1+data_snake.index(c-1-1j)))
            Z_Map[check_snake.index(c)].append(z_distance*x_distance-1+data_snake.index(c-1-1j))



    Map = defaultdict(list)

    for k in X_Map.keys() | Z_Map.keys():   # union of keys
        Map[k] = X_Map[k] + Z_Map[k]         # union of sets

    # Map = sorted(Map.items())
    
    # Map_2N = defaultdict(list)

    # for k in X_Map_2N.keys() | Z_Map_2N.keys():   # union of keys
    #     Map_2N[k] = X_Map_2N[k] + Z_Map_2N[k]         # union of sets
    Map_2N:List[int] = []
    Map_2N = list(range(min(min(X_Map_2N),min(Z_Map_2N)),1+max(max(X_Map_2N),max(Z_Map_2N))))
    # Map_2N = sorted(Map_2N.items())
    # Map = X_Map | Z_Map ## Apparently, this doesn't work, since they overwrite when there are common keys
    # Map_2N = X_Map_2N | Z_Map_2N ## Apparently, this doesn't work, since they overwrite when there are common keys

    return (x_distance,
            z_distance,
            X_observable_index,
            Z_observable_index,
            # data_coords,
            # x_measure_coords,
            # z_measure_coords,
            data_snake,
            check_snake,
            measurement_qubits,
            data_qubits,
            x_measure_index,
            z_measure_index,
            # X_Map,
            Map,
            X_Map_2N,
            # Z_Map,
            Z_Map_2N,
            Map_2N)