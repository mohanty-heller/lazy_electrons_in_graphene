# lazy_electrons_in_graphene
Public repository for "Lazy electrons in graphene" by Vaibhav Mohanty and Eric J. Heller, published in the Proceedings of the National Academy of Sciences.

Instructions for use:
'main_script.m' contains the same code as the Mathematica notebook 'main.nb', but suitable for use on a supercomputering cluster. This calculates matrix elements for the Hamiltonian for all discrete times specified. The Hamiltonian is then constructed by the first half of 'process_outputs.m' in MATLAB, which then constructs the Hamiltonian in matrix form and calculates the eigenenergies and eigenstates for all times. Specifying the intial ABO state, the TDSE is solved by 'ODE_solve.nb' in Mathematica. The TDSE wavefunction output at all times is exported and can be processed by the second half of 'process_outputs.m' to produce all relevant results, including overlap between TDSE and ABO states, autocorrelation functions, etc.
