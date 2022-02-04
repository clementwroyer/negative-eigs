# Results directory

Contains the output of our NES experiments

## Contents

Two type of files: the .mat files contain data about the runs (useful variables, eigenvalues, etc) while the ResNES files present the information in a table format.

## Naming convention

NFD: Indicates that matrices produced via Newton iterations and finite differences were used
Heuristics/AllOrders: Indicates which ordering were considered (4 heuristics or all possible ones)
RandPerm (optional): If present, indicates that a random permutation was applied to every matrix
dimsXtoY: Problems from dimensions X to Y were considered for that run
