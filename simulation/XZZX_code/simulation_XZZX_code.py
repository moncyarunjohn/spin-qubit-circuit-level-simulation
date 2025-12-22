import stim
import sinter
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../spin-qubit-circuit-level-simulation'))
from circuits.ThreeNArray_XZZX_surface_code_architecture import create_XZZX_surface_code_architecture, CircuitGenParametersXZZX
 
"""
Rotated XZZX surface code simulation for the spin qubit architecture with 3N Array.
"""

# Generates surface code circuit tasks using Stim's circuit generation.
def generate_example_tasks(is_memory_H=False):
    etaXs = [0.5, 10, 1000]
    etas = [0.5,10,1000]#[0.5, 1, 10, 100, 1000, 10000]
    probabilities = [0.002,0.003,0.004]
    distances = [11,13,15,17]
    for etaX in etaXs:
        for eta in etas:   
            for p in probabilities:
                for d in distances:            
                    rounds = 3 * d     
                    params = CircuitGenParametersXZZX(
                                                        rounds=rounds,
                                                        distance=d,
                                                        after_clifford1_depolarization = p/10,
                                                        before_round_data_bias_probability= (p/10, eta),
                                                        before_measure_flip_probability = 2 * p,
                                                        after_reset_flip_probability =  2 * p,
                                                        after_clifford2_depolarization=p,                                    
                                                        pshuttle_biased = (0.1*p,etaX),
                                                        # nswaps=(3,2), # (Ny,Nx) in the main text, defines the swaps of checks and datas (per 2 qubit gate)
                                                    )
                    circuit = create_XZZX_surface_code_architecture(params, is_memory_H=is_memory_H)
                    # --- Add this line to print the STIM circuit ---
                    # print(circuit)
                    
                    yield sinter.Task(
                        circuit=circuit,
                        decoder=None,
                        # detector_error_model=decoder_dem,
                        json_metadata={
                            'p': p,
                            'd': d,
                            "eta": eta,
                            "etaX": etaX,
                            "params": params.__dict__,
                            "memory": ["V", "H"][is_memory_H]}       
                            )               

def main():
    filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "3NArray")#"FootprintSwap08_p_32")
    # Collect the samples (takes a few minutes).
    samples = sinter.collect(
        num_workers=multiprocessing.cpu_count()-1,
        max_shots=2_000_000,
        max_errors=2000,
        tasks=[task for task in generate_example_tasks(is_memory_H=False)] + [task for task in generate_example_tasks(is_memory_H=True)],
        decoders=["pymatching"],
        #count_detection_events=True,
        save_resume_filepath= os.path.join(filepath, "results.csv")
    )
    # Print samples as CSV data.
    print(sinter.CSV_HEADER)
    for sample in samples:
        print(sample.to_csv_line())


# NOTE: This is actually necessary! If the code inside 'main()' was at the
# module level, the multiprocessing children spawned by sinter.collect would
# also attempt to run that code.
if __name__ == '__main__':
    main()