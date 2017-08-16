# CAFCA
Celular Automata for Fluid Charge Adjustments

Note to myself: 
1. Don't forget that the python script that calculates the number of electrons per radii is using a fixed dx for discretization, it must be generalized
2. Python calculates qpercube but it's rounded to the 8th decimal for clarity. This can be extended to more decimals, but output is visually nasty
3. When placing the electrons around the nuclei, if it is placed where a nucleus resides, the charge of the nucleus will dicrease. Bug.
4. Efficiency can be improved by calculating the acceleration only where there is charge in the cell.
5. If the maximum velocity is reached, the previous one is kept. Bug.
6. The molecule initialization puts the atoms in a grid where each cube corresponds to the coordinates of one vertex when it should correspond to the coordinates of the center of the cube
7. The input mol2 file must have atoms names with the atom symbol followed by at least a number
