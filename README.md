# CAFCA
Celular Automata for Fluid Charge Adjustments

Note to myself: 
1. Don't forget that the python script that calculates the number of electrons per radii is using a fixed dx for discretization, it must be generalized
2. Python calculates qpercube but it's rounded to the 8th decimal for clarity. This can be extended to more decimals, but output is visually nasty
3. When placing the electrons around the nuclei, if it is placed where a nucleus resides, the charge of the nucleus will dicrease. Bug.
