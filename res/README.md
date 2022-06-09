# Results directory

Contains the output of our NES experiments

## Contents

Text files presenting the results of various runs.

## Naming convention

- *Heuristics/Comparison*: Indicates whether the file runs 4 heuristics or compares the best heuristic with all orderings.

- *Exact/FD/ExactFD*: Indicates whether actual Hessian matrices, finite-difference approximations thereof or both are considered.

- *All/1eMx* (Optional): When finite-difference matrices are used, indicates what tolerances are used (All of them or 10^{-x} with x a positive integer).

- *Determ/RandOrthog/RandPerm*: Indicates whether actual matrices are used (Exact) or whether these matrices are transformed via a random orthogonal matrix (RandOrthog) or a random perturbation of the indices (RandPerm).
