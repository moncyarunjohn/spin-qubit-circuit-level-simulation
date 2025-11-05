The code is divided into different folders

1. **Circuits**: Design the architectures here.(EDIT)
    Files:
         rotated surface code layout (this file is used in the other two files. Can be used further as well)
         CSS Surface code architecture
         XZZX surface code architecture
   
   The different architectures are designed here.
   CSS, XZZX, GB code needs to be written (essentially, different QECC in different architectures - 2N, 3N Railways, instead of MEC).

2. **Simulation**:Code run with sinter and RESULTS saved here. (CHANGE PARAMATERS ONLY) 
   Files:
         CSS code (folder) - Contain plotting files and Results Folder
         XZZX code (folder) - Contain plotting files and Results Folder
   Run this to get the outputs (phy error rate vs logical error rate for different code distances)
   Change the probabilities, code distance, number of shots, etc here.
   No need to change the code much, except to alter the parameters.

3. **Visualize**: (DO NOT EDIT)
   The syndrome extraction + error correction (?) circuit is visualized. It's shown in the topological space (?). *How will this look for GB code?* *Or no need to write code to visualise it?* Can be done by Reza's concentric circle thing.



Sinter: This is the package which uses circuits from STIM and does the error rate plots, saving data and all other stuff. Hence it comes in simulation part. This need not be changed.
Knowing its working will simplify a lot of confusions and let us focus on STIM.

What are seaborn, sys, os?



Detectors in the code (Circuits folder).

**Message from Josu**: 
"dem = circuit.detector\_error\_model(decompose\_errors=True)

    bm = detector\_error\_model\_to\_check\_matrices(dem, allow\_undecomposed\_hyperedges=False)

from beliefmatching import BeliefMatching, detector\_error\_model\_to\_check\_matrices"

Detectors for the circuit need not be written explicitly. *How does it find the detector though?*



*Inside, simulation (to run) the file 'plot\_utils.py' does what? sounds like helps plot the phy error rate vs log error rate. doesn't sinter do that?*



Inside the circuits folder,
    rotated\_surface\_code is the default architecture. The other two (CSS and XZZX) are derived from it.

   This is because both are rotated surface codes (rotated from 2D Toric code). But the stabilizers are different.

   Hence the main work (stabilizers and measurement qubits) are defined inside the CSS and XZZX.

   XZZX comment says it's MEC based. *Not sure yet where.*




### Comments - After reading the code #######

The p2q and cnot_target lists need to be re-written as in another form. Like, integer to integer.
The snake pattern should give one order (data qubit).
The check qubits also has snake patterns.
But how to create a map?
- Follow Josu's indexing to fix the qubits in the 2D
