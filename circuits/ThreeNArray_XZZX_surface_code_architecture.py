# from collections import defaultdict
import stim
from typing import Set, List, Dict, ClassVar
from dataclasses import dataclass
import math
from collections import defaultdict
from itertools import islice
import sys
sys.path.append('../spin_qubit_architecture_circuits')
from circuits.ThreeNArray_surface_code_layout import generate_rotated_surface_code_circuit_layout # type: ignore

"""
Code for simulating a rotated CSS surface code over the spin qubit architecture using Arrays

"""


def append_anti_basis_error(circuit: stim.Circuit, targets: List[int], p: float, basis: str) -> None: #Adds the errors of the basis (X error for Z-basis) to target qubits appended to the circuit
    if p > 0:
        if basis == "X":
            circuit.append("Z_ERROR", targets, p)
        else:
            circuit.append("X_ERROR", targets, p) # What does this do for basis ZX?


@dataclass
class CircuitGenParametersCSS():
    """
    Parameters for generating a rotated CSS surface code circuit over the architecture.
    
    Idling errors (biased):
    The probability of the idling error is scaled according the subsequent operation being done.
        parameter : before_round_data_bias_probability tuple (rate,bias)

        p = p_x + p_y + p_z and η = p_z / (p_x + p_y)
        
        therefore:
        p_x = p_y = p / [2*(1+η)]
        p_z = p*η / [1+η]
    
    Hadamard gates:
        parameter : after_clifford1_depolarization

        Single qubit depolarizing errors
        

    Two Qubit Gate errors (CNOT,CZ):
        parameter : after_clifford2_depolarization

        Two qubit depolarizing errors

    State preparation error:
        parameter : after_reset_flip_probability

        flip state to orthogonal state

    Measurement error:
        parameter : before_measure_flip_probability

        Flip measurement outcome.
    
    Swap error:   Swaps were used in this architecture? instead of shuttle? - Because of MEC. For shuttling, different architechture is needed.
        parameter : pswap_depolarization
        parameter : nswaps

        nswap depolarizing errors in qubits that are swapped to get closer. 
    """
    
    rounds: int
    distance: int = None
    x_distance: int = None
    z_distance: int = None
    after_clifford1_depolarization: float = 0 # this will be for the single qubit gates
    before_round_data_bias_probability: tuple = (0, 0) # (p, eta) this is the idling (we relate it with the subsequent operations)
    before_measure_flip_probability: float = 0 # # this is for  the measurement errors
    after_reset_flip_probability: float = 0 # # this is for the reset errors
    exclude_other_basis_detectors: bool = False
    after_clifford2_depolarization: float = 0 # this is for the two qubit gates
    pswap_depolarization: float = 0 ## This has not been altered anywhere else ## error for shuttling the qubits ### pswap_depolarization: float = 0 # error for swapping the qubits
    nswaps: tuple = (0,0) # (Ny,Nx) number of swaps in each direction to get closer ### nswaps: tuple = (0,0) # (Ny,Nx) number of swaps in each direction to get closer







def create_rotated_surface_code_CSS_architecture(params: CircuitGenParametersCSS,
                                is_memory_H: bool = False,
                                *, 
                                exclude_other_basis_detectors: bool = False,
                                ) -> stim.Circuit:
    
    pshuttle_depolarization = params.pswap_depolarization

    if params.rounds < 2:
        raise ValueError("Need rounds >= 2")
    if params.distance is not None and params.distance < 2:
        raise ValueError("Need a distance >= 2")    
    (x_distance,
    z_distance,
    X_observable_index,
    Z_observable_index,
    data_snake, ## Is not used at all?
    check_snake, ## Probably only used for visualization
    measurement_qubits,
    data_qubits,
    x_measure_index,
    z_measure_index,
    # X_Map, ## Probably not needed. Can be combined with Z_Map and sent as Map.
    Map, ## Used only in the last round. Really needed?
    X_Map_2N,
    # Z_Map,## Probably not needed. Can be combined with Z_Map and sent as Map.
    Z_Map_2N,
    Positions_2N ## This can be substituted by the 'range (-d,d)'. DONE
    ) = generate_rotated_surface_code_circuit_layout(params.distance, params.x_distance, params.z_distance)

    chosen_basis_observable = Z_observable_index if is_memory_H else X_observable_index ### Done. But caution.This should be edited. But the x and z needed to be swapped according the other xzzx file?
    chosen_basis_measure_index = z_measure_index if is_memory_H else x_measure_index ### This too should be edited. Has been.

    
    Map_2N = defaultdict(list)

    for k, v in X_Map_2N.items():
        Map_2N[k] += v

    for k, v in Z_Map_2N.items():
        Map_2N[k] += v
    
    
    if is_memory_H: 
        data_qubits_x = data_qubits[::2]
        data_qubits_z = data_qubits[1::2]
    else:
        data_qubits_x = data_qubits[1::2]
        data_qubits_z = data_qubits[::2]


    # backandforth = True




    #                                      000000   0       0   000000  0        000000000         11
    #                                     0          0     0   0        0        0                  1
    #                                     0           0   0    0        0        0                  1
    #                                     0            0 0     0        0        0                  1
    #                                     0             0      0        0        0000000            1
    #                                     0             0      0        0        0                  1
    #                                     0             0      0        0        0                  1
    #                                     0             0      0        0        0                  1
    #                                      000000       0       000000  0000000  000000000         111

    #####--CYCLE--###################################################
    # Build the repeated actions that make up the surface code cycle
    ## We add the cycle actions to the head, body and tail. This is the part which is repeated.
    cycle_actions_1 = stim.Circuit()

    # backandforth = not backandforth

    # Hadamard gates to check qubits
    cycle_actions_1.append("TICK", []) 
    cycle_actions_1.append("H", measurement_qubits)
    if params.after_clifford1_depolarization > 0:
        cycle_actions_1.append("DEPOLARIZE1", measurement_qubits, params.after_clifford1_depolarization)
    # Biased channel on the data qubits as per idling in this time step
    ## Same as reset above, but now the time step is the Hadamard gate
    if params.before_round_data_bias_probability[0] > 0:
        # We consider that for noisier operations related with longer times the idling should be higher
        idling_fact = params.after_clifford1_depolarization / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        cycle_actions_1.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])



    # Entangling gates
    ## This is probably where we need to alter.
    ## We need to add the shuttling errors before and after the CNOTs
    ## We need to add the idling errors to the qubits that are not involved in the CNOTs
    ### Plan 1: Do the (x,y)-coordinates of all the qubits in the RSC. This already exists in q2p.
    ### Plan10: Find the figure/graph that Ruben wanted. Adjacency graph. Algebraic, not hardcoded (except it being for RSC). For any d.
    ### Plan11: Using (x,y), plot the snake pattern and affix indices.
    ### Plan12: With the (x,y), create map which shows data to check indices.
    ### PLan13: Find starting and ending positions of the check qubits = (+d,1-d) [Assuming check and data align exactly, exept for the last data qubit]
    ### Plan14: Use maps of the data, check indices to find the CNOT pairs. (Although there should be a formula for the figure. Not worth it though)

    for position in sorted(Positions_2N):#, reverse=backandforth): ## Careful with this. backandforth do not run inside a loop. We are only defining a type of circuit.
        cycle_actions_1.append("TICK", [])
        
        for k in measurement_qubits:
            cycle_actions_1.append("QUBIT_COORDS", [k], [position+k, 0])
        
        if position == sorted(Positions_2N)[0]: # We do not need to do shuttling before the first CNOTs ### DOES THIS WORK IN THE FORTH (DIRECTION IN EVERY SECOND ROUND)?
            pass
        else:
            if pshuttle_depolarization>0:
                cycle_actions_1.append("DEPOLARIZE1", measurement_qubits,pshuttle_depolarization)
                
            if params.before_round_data_bias_probability[0] > 0:
                # We consider that for noisier operations related with longer times the idling should be higher
                idling_fact = pshuttle_depolarization / params.after_clifford2_depolarization
                p = idling_fact * params.before_round_data_bias_probability[0]
                eta = params.before_round_data_bias_probability[1]
                p_x = p/(2*(1+eta))
                p_y = p/(2*(1+eta))
                p_z = p*eta / (1+eta)
                cycle_actions_1.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])
        
        # if position in X_Map_2N:
        #     for pair in X_Map_2N[position]:
        #         cycle_actions_1.append("CNOT", pair)
        #         if params.after_clifford2_depolarization > 0:
        #             cycle_actions_1.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)
        # if position in Z_Map_2N:
        #     for pair in Z_Map_2N[position]:
        #         cycle_actions_1.append("CZ", pair)
        #         if params.after_clifford2_depolarization > 0:
        #             cycle_actions_1.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)
        # # cycle_actions.append("DEPOLARIZE1", sorted_ReMap(position),pshuttle_depolarization)

        for pair in Map_2N[position]:
            if check_snake[pair[0]]-data_snake[pair[1]-x_distance*z_distance+1] in ((1+1j), (-1-1j)):
                cycle_actions_1.append("CNOT", pair)
                if params.after_clifford2_depolarization > 0:
                    cycle_actions_1.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)
            else:
                cycle_actions_1.append("CZ", pair)
                if params.after_clifford2_depolarization > 0:
                    cycle_actions_1.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)

        if params.before_round_data_bias_probability[0] > 0:
            p = params.before_round_data_bias_probability[0]
            eta = params.before_round_data_bias_probability[1]
            p_x = p/(2*(1+eta))
            p_y = p/(2*(1+eta))
            p_z = p*eta / (1+eta)
            # We select the data_qubits and measurement_qubits that do not belong to the targets
            cycle_actions_1.append("PAULI_CHANNEL_1", list((set(data_qubits) | set(measurement_qubits))-(set().union(*X_Map_2N[position]) | set().union(*Z_Map_2N[position]))), [p_x, p_y, p_z])


    # Hadamard gates to check qubits
    cycle_actions_1.append("TICK", [])    
    cycle_actions_1.append("H", measurement_qubits)
    if params.after_clifford1_depolarization > 0:
        cycle_actions_1.append("DEPOLARIZE1", measurement_qubits, params.after_clifford1_depolarization)
    # Biased channel on the data qubits as per idling in this time step

    if params.before_round_data_bias_probability[0] > 0:
        # rescale idling as 1 qubit gate
        idling_fact = params.after_clifford1_depolarization / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0] 
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        cycle_actions_1.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])
    
    # Measure the check qubits
    cycle_actions_1.append("TICK", [])    
    if params.before_measure_flip_probability > 0:
        append_anti_basis_error(cycle_actions_1, measurement_qubits, params.before_measure_flip_probability, basis="Z")
    cycle_actions_1.append("M" + "Z", measurement_qubits)





    #                                      000000   0       0   000000  0        000000000          2222222
    #                                     0          0     0   0        0        0                 2       2
    #                                     0           0   0    0        0        0                         22 
    #                                     0            0 0     0        0        0                        22
    #                                     0             0      0        0        0000000                 22
    #                                     0             0      0        0        0                      2
    #                                     0             0      0        0        0                    22
    #                                     0             0      0        0        0                  22
    #                                      000000       0       000000  0000000  000000000         2222222222

    #####--CYCLE--###################################################
    # Build the repeated actions that make up the surface code cycle
    ## We add the cycle actions to the head, body and tail. This is the part which is repeated.
    cycle_actions_2 = stim.Circuit()

    # backandforth = not backandforth

    # Hadamard gates to check qubits
    cycle_actions_2.append("TICK", []) 
    cycle_actions_2.append("H", measurement_qubits)
    if params.after_clifford1_depolarization > 0:
        cycle_actions_2.append("DEPOLARIZE1", measurement_qubits, params.after_clifford1_depolarization)
    # Biased channel on the data qubits as per idling in this time step
    ## Same as reset above, but now the time step is the Hadamard gate
    if params.before_round_data_bias_probability[0] > 0:
        # We consider that for noisier operations related with longer times the idling should be higher
        idling_fact = params.after_clifford1_depolarization / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        cycle_actions_2.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])



    # Entangling gates
    ## This is probably where we need to alter.
    ## We need to add the shuttling errors before and after the CNOTs
    ## We need to add the idling errors to the qubits that are not involved in the CNOTs
    ### Plan 1: Do the (x,y)-coordinates of all the qubits in the RSC. This already exists in q2p.
    ### Plan10: Find the figure/graph that Ruben wanted. Adjacency graph. Algebraic, not hardcoded (except it being for RSC). For any d.
    ### Plan11: Using (x,y), plot the snake pattern and affix indices.
    ### Plan12: With the (x,y), create map which shows data to check indices.
    ### PLan13: Find starting and ending positions of the check qubits = (+d,1-d) [Assuming check and data align exactly, exept for the last data qubit]
    ### Plan14: Use maps of the data, check indices to find the CNOT pairs. (Although there should be a formula for the figure. Not worth it though)

    for position in sorted(Positions_2N,reverse=True):#, reverse=backandforth): ## Careful with this. backandforth do not run inside a loop. We are only defining a type of circuit.
        cycle_actions_2.append("TICK", [])

        for k in measurement_qubits:
            cycle_actions_2.append("QUBIT_COORDS", [k], [position+k, 0])

        if position == sorted(Positions_2N,reverse=True)[0]: # We do not need to do shuttling before the first CNOTs ### SHOULD THIS BE MAX INSTEAD? BCAUSE IT'S FORTH AND NOT BACK?
            pass
        else:
            if pshuttle_depolarization>0:
                cycle_actions_2.append("DEPOLARIZE1", measurement_qubits,pshuttle_depolarization)
                
            if params.before_round_data_bias_probability[0] > 0:
                # We consider that for noisier operations related with longer times the idling should be higher
                idling_fact = pshuttle_depolarization / params.after_clifford2_depolarization
                p = idling_fact * params.before_round_data_bias_probability[0]
                eta = params.before_round_data_bias_probability[1]
                p_x = p/(2*(1+eta))
                p_y = p/(2*(1+eta))
                p_z = p*eta / (1+eta)
                cycle_actions_2.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])
        
        # if position in X_Map_2N:
        #     for pair in X_Map_2N[position]:
        #         cycle_actions_2.append("CNOT", pair)
        #         if params.after_clifford2_depolarization > 0:
        #             cycle_actions_2.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)
        # if position in Z_Map_2N:
        #     for pair in Z_Map_2N[position]:
        #         cycle_actions_2.append("CZ", pair)
        #         if params.after_clifford2_depolarization > 0:
        #             cycle_actions_2.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)
        # # cycle_actions.append("DEPOLARIZE1", sorted_ReMap(position),pshuttle_depolarization)
        for pair in Map_2N[position]:
            if check_snake[pair[0]]-data_snake[pair[1]-x_distance*z_distance+1] in ((1+1j), (-1-1j)):
                cycle_actions_2.append("CNOT", pair)
                if params.after_clifford2_depolarization > 0:
                    cycle_actions_2.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)
            else:
                cycle_actions_2.append("CZ", pair)
                if params.after_clifford2_depolarization > 0:
                    cycle_actions_2.append("DEPOLARIZE2", pair, params.after_clifford2_depolarization)
        if params.before_round_data_bias_probability[0] > 0:
            p = params.before_round_data_bias_probability[0]
            eta = params.before_round_data_bias_probability[1]
            p_x = p/(2*(1+eta))
            p_y = p/(2*(1+eta))
            p_z = p*eta / (1+eta)
            # We select the data_qubits and measurement_qubits that do not belong to the targets
            cycle_actions_2.append("PAULI_CHANNEL_1", list((set(data_qubits) | set(measurement_qubits))-(set().union(*X_Map_2N[position]) | set().union(*Z_Map_2N[position]))), [p_x, p_y, p_z])


    # Hadamard gates to check qubits
    cycle_actions_2.append("TICK", [])    
    cycle_actions_2.append("H", measurement_qubits)
    if params.after_clifford1_depolarization > 0:
        cycle_actions_2.append("DEPOLARIZE1", measurement_qubits, params.after_clifford1_depolarization)
    # Biased channel on the data qubits as per idling in this time step

    if params.before_round_data_bias_probability[0] > 0:
        # rescale idling as 1 qubit gate
        idling_fact = params.after_clifford1_depolarization / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0] 
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        cycle_actions_2.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])
    
    # Measure the check qubits
    cycle_actions_2.append("TICK", [])    
    if params.before_measure_flip_probability > 0:
        append_anti_basis_error(cycle_actions_2, measurement_qubits, params.before_measure_flip_probability, basis="Z")
    cycle_actions_2.append("M" + "Z", measurement_qubits)










    #                                     0     0  0000000  00000   000000 
    #                                     0     0  0        0   0   0     0
    #                                     0     0  0        0   0   0     0
    #                                     0000000  00000    0   0   0     0
    #                                     0     0  0        00000   0     0
    #                                     0     0  0        0   0   0     0
    #                                     0     0  0000000  0   0   000000

    ####--HEAD--####################################################
    # Build the start of the circuit, getting a state that's ready to cycle
    # In particular, the first cycle has different detectoors and so has to be handled special.
    ### Why though? Are we not doing the same as in the body?

    ## This is for visualization. We'll fix this later
    head = stim.Circuit()
    for k in data_qubits:
        head.append("QUBIT_COORDS", [k], [k-x_distance*z_distance + 1, 1])
    

    # Reset the data qubits
    head.append("TICK", []) 
    head.append("RX", data_qubits_x)
    head.append("R", data_qubits_z)
    # head.append("R" + "ZX"[is_memory_x], data_qubits) # What exactly is a "ZX". There is no definition of it in the function (append_anti_basis_error) defined in the beginning of this script.
    ### Why reset in ZX basis? Because we want to prepare the data qubits in the logical 0 or + state? 
    ### How do you reset in logical basis for a random code? This was discussed with Reza in the beginning. No clear answer yet.
    ### In the other file for XZZX code, it seems that the alternate data qubits are reset in states 0 and states + respectively.

    if params.after_reset_flip_probability > 0:
        append_anti_basis_error(head, data_qubits_x, params.after_reset_flip_probability, "X")
        append_anti_basis_error(head, data_qubits_z, params.after_reset_flip_probability, "Z")

    # if params.after_reset_flip_probability > 0:
    #     append_anti_basis_error(head, data_qubits, params.after_reset_flip_probability, "ZX"[is_memory_x]) ## the function doesn't have a ZX specific action, so why?
    #### Pay attention here. What to do?


    # We have compiled the rotated planar code with CZ for the Z checks and every Hadamard gate in state prep
    # has been made to form a |+> state for the check, here explicitly with the Hadamard gate
    ## The check qubits were not present until now? OR did the reset makes it trivial?
    head.append("R" + "Z", measurement_qubits)
    if params.after_reset_flip_probability > 0:
        append_anti_basis_error(head, measurement_qubits, params.after_reset_flip_probability, "Z")

    head += cycle_actions_1

    # Biased channel on the data qubits as per idling in this time step
    if params.before_round_data_bias_probability[0] > 0:
        # rescale idling as measurement
        idling_fact = params.before_measure_flip_probability / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        head.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])
    ## how is this done? How is head and cycle_actions combined?
    ## Just the detectoor is added here. Maybe it's different from cycle_actions. Cycle_actions do not have detectoors.
    ## From my memory, the data qubit initializaion and the detectoors are differnt from cycle_actions.
    ## Why wasn't the rest of the parts not added to head. Like head = dataqubit_initialization + cycle_actions + head_detectoors

    # for measure in sorted(chosen_basis_measure_coords, key=lambda c: (c.real, c.imag)):
    #     head.append(
    #         "DETECTOR",
    #         [stim.target_rec(-len(measurement_qubits) + measure_coord_to_order[measure])],
    #         [measure.real, measure.imag, 0.0]
    #     )

    for measure in chosen_basis_measure_index:
        head.append(
            "DETECTOR",
            [stim.target_rec(-len(measurement_qubits) + measure)],
            [check_snake[measure].real, check_snake[measure].imag, 0.0] # What is this visualization?
            # If we are visualizing 2N array, then the detectoors might need to keep moving. Currently the detectoors are in RSC
        )








    #                                     00000    0000    00000   0     0     11
    #                                     0    0  0    0   0    0   0   0       1
    #                                     0    0  0    0   0    0    0 0        1
    #                                     00000   0    0   0    0     0         1
    #                                     0    0  0    0   0    0     0         1
    #                                     0    0  0    0   0    0     0         1
    #                                     00000    0000    00000      0        111


    ####--BODY--####################################################
    # Build the repeated body of the circuit, including the detectoors comparing to previous cycles.
    body_1 = stim.Circuit()

    body_1.append("TICK", [])
    body_1.append("R" + "Z", measurement_qubits) # Reset all measurement qubits to Z measurement basis (0)
    ## R probably means reset
    ## Note that only measurement qubits are reset. Data qubits have data (but do they get error corrected every round? No, right?)
    if params.after_reset_flip_probability > 0:
        append_anti_basis_error(body_1, measurement_qubits, params.after_reset_flip_probability, "Z")

    # Biased channel on the data qubits as per idling in this time step
    ## What time step? The one in which the checks are being reset?
    if params.before_round_data_bias_probability[0] > 0:
        # We consider that for noisier operations related with longer times the idling should be higher
        idling_fact = params.after_reset_flip_probability / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        body_1.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])

    body_1 += cycle_actions_2

    # Biased channel on the data qubits as per idling in this time step
    if params.before_round_data_bias_probability[0] > 0:
        # rescale idling as measurement
        idling_fact = params.before_measure_flip_probability / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        body_1.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])

    m = len(measurement_qubits)
    body_1.append("SHIFT_COORDS", [], [0.0, 0.0, 1.0])
    for m_index in measurement_qubits:
        # m_coord = check_snake[m_index]
        k = len(measurement_qubits) - m_index
        if not exclude_other_basis_detectors or m_index in chosen_basis_measure_index:
            body_1.append(
                "DETECTOR",
                [stim.target_rec(-k), stim.target_rec(-k - m)],
                [check_snake[m_index].real, check_snake[m_index].imag, 0.0] # Again, this visualizaion is in RSC, not 2N array
            )
    ## This looks pretty short. Just the indices are updated. Then the corresponding detectoors are defined.








    #                                     00000    0000    00000   0     0       222222
    #                                     0    0  0    0   0    0   0   0       2      2
    #                                     0    0  0    0   0    0    0 0              22
    #                                     00000   0    0   0    0     0              2
    #                                     0    0  0    0   0    0     0            22
    #                                     0    0  0    0   0    0     0          22
    #                                     00000    0000    00000      0         22222222


    ####--BODY--####################################################
    # Build the repeated body of the circuit, including the detectoors comparing to previous cycles.
    body_2 = stim.Circuit()

    body_2.append("TICK", [])
    body_2.append("R" + "Z", measurement_qubits) # Reset all measurement qubits to Z measurement basis (0)
    ## R probably means reset
    ## Note that only measurement qubits are reset. Data qubits have data (but do they get error corrected every round? No, right?)
    if params.after_reset_flip_probability > 0:
        append_anti_basis_error(body_2, measurement_qubits, params.after_reset_flip_probability, "Z")

    # Biased channel on the data qubits as per idling in this time step
    ## What time step? The one in which the checks are being reset?
    if params.before_round_data_bias_probability[0] > 0:
        # We consider that for noisier operations related with longer times the idling should be higher
        idling_fact = params.after_reset_flip_probability / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        body_2.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])

    body_2 += cycle_actions_1

    # Biased channel on the data qubits as per idling in this time step
    if params.before_round_data_bias_probability[0] > 0:
        # rescale idling as measurement
        idling_fact = params.before_measure_flip_probability / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        body_2.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])

    m = len(measurement_qubits)
    body_2.append("SHIFT_COORDS", [], [0.0, 0.0, 1.0])
    for m_index in measurement_qubits:
        # m_coord = check_snake[m_index]
        k = len(measurement_qubits) - m_index
        if not exclude_other_basis_detectors or m_index in chosen_basis_measure_index:
            body_2.append(
                "DETECTOR",
                [stim.target_rec(-k), stim.target_rec(-k - m)],
                [check_snake[m_index].real, check_snake[m_index].imag, 0.0] # Again, this visualizaion is in RSC, not 2N array
            )
    ## This looks pretty short. Just the indices are updated. Then the corresponding detectoors are defined.






    #                                     TTTTTTT  0000000  IIIII  L        
    #                                        T     0     0    I    L        
    #                                        T     0     0    I    L        
    #                                        T     0000000    I    L        
    #                                        T     0     0    I    L        
    #                                        T     0     0    I    L        
    #                                        T     0     0  IIIII  LLLLL    


    ####--TAIL--####################################################
    # Build the end of the circuit, getting out of the cycle state and terminating.
    # In particular, the data measurements create detectoors that have to be handled special.
    # Also, the tail is responsible for identifying the logical observable.
    tail = stim.Circuit()
     # Reset the checks
    tail.append("TICK", []) 
    tail.append("R" + "Z", measurement_qubits)
    if params.after_reset_flip_probability > 0:
        append_anti_basis_error(tail, measurement_qubits, params.after_reset_flip_probability, "Z")
    # Biased channel on the data qubits as per idling in this time step
    
    if params.before_round_data_bias_probability[0] > 0:
        idling_fact = params.after_reset_flip_probability / params.after_clifford2_depolarization
        p = idling_fact * params.before_round_data_bias_probability[0]
        eta = params.before_round_data_bias_probability[1]
        p_x = p/(2*(1+eta))
        p_y = p/(2*(1+eta))
        p_z = p*eta / (1+eta)
        tail.append("PAULI_CHANNEL_1", data_qubits, [p_x, p_y, p_z])
    
    if (params.rounds)%2 == 0:
        tail+= cycle_actions_2
    else:
        tail+= cycle_actions_1

    
    # Detectoors of checks last round
    m = len(measurement_qubits)
    tail.append("SHIFT_COORDS", [], [0.0, 0.0, 1.0])
    for m_index in measurement_qubits:
        # m_coord = q2p[m_index]
        k = len(measurement_qubits) - m_index
        if not exclude_other_basis_detectors or m_index in chosen_basis_measure_index:
            tail.append(
                "DETECTOR",
                [stim.target_rec(-k), stim.target_rec(-k - m)],
                [check_snake[m_index].real, check_snake[m_index].imag, 0.0] # Visualization in RSC
            )

 #   # Measure the data qubits
 #   if params.before_measure_flip_probability > 0:
 #       append_anti_basis_error(tail, data_qubits, params.before_measure_flip_probability, "ZX"[is_memory_x])
 #   tail.append("M" + "ZX"[is_memory_x], data_qubits)
 #   ### Pay attention here. What to do?

 #    Measure the data qubits. Similar to the reset, every other qubit is measured in the X basis and the rest in the Z basis.
    if params.before_measure_flip_probability > 0:
        append_anti_basis_error(tail, data_qubits_x, params.before_measure_flip_probability, "X")
        append_anti_basis_error(tail, data_qubits_z, params.before_measure_flip_probability, "Z")  
 #   for q in data_qubits_z:
 #       # this keeps the order of measuring the data qubits consistent with the order of the data qubits
 #       tail.append("M" + "XZ"[q in data_qubits_z], [q])
 #   for q in data_qubits_x:
 #       # this keeps the order of measuring the data qubits consistent with the order of the data qubits
 #       tail.append("M" + "XZ"[q in data_qubits_z], [q])
    tail.append("MX",  data_qubits_x)
    tail.append("M",  data_qubits_z)
    
    # Detectors
    for measure in chosen_basis_measure_index: ## Any need to sort this first? #sorted(chosen_basis_measure_coords, key=lambda c: (c.real, c.imag)):
        detectors: List[int] = []
        for delta in Map[measure]:
            # data = measure + delta
            # if data in p2q:
            detectors.append(-len(data_qubits)-len(measurement_qubits)+delta) ## This is the 4 data qubits associated with check qubit: 'measure'
        detectors.append(-len(data_qubits) - len(measurement_qubits) + measure) ## What is this? The check qubit added to the list of 4 data qubits? Parity check for all the 4 data qubits+ the corresponding check qubit?
        detectors.sort(reverse=True) ## Needed?
        tail.append("DETECTOR", [stim.target_rec(x) for x in detectors], [check_snake[measure].real, check_snake[measure].imag, 1.0]) ## -x? Maybe not since - has been added above already.

    tail.append("TICK", []) 

    # Logical observable
    obs_inc: List[int] = []
    for q in chosen_basis_observable:
        obs_inc.append(-len(data_qubits) + q) # the q here is 24+index of the data qubit. wait wait... It's not. chose_basis_observable, even though is data qubits, is indexed from 0 to 25
    obs_inc.sort(reverse=True) ## Why?
    tail.append("OBSERVABLE_INCLUDE", [stim.target_rec(x) for x in obs_inc], 0.0) ## The rec's are measurements, which need to be accounted for when indexing the measurments for the detectors.











    #####--FINAL--####################################################
    # Combine to form final circuit.
    return head + (body_1+body_2) * ((params.rounds - 2)//2) +(params.rounds%2)*body_1 + tail ## (rounds-2) factor requires rounds >= 2.