This directory contains all test cases used in our experiments associated with the paper titled:
"A Novel Robust and Efficient Method for Computing High-Quality AC OPF Solutions Based on Trust-Tech and Holomorphic Embedding".

The cases are collected from the following sources:

- [MATPOWER 8.1](https://matpower.org)
- [PGLIB-OPF](https://github.com/power-grid-lib/pglib-opf)
- [OptEnergyLocalOpt](https://webhomes.maths.ed.ac.uk/OptEnergy/LocalOpt/) – these cases are known to possess multiple local optimal solutions.

#### 1. `Case_mod/`

This folder contains modified cases adapted from [OptEnergyLocalOpt](https://webhomes.maths.ed.ac.uk/OptEnergy/LocalOpt/), along with one newly constructed benchmark:

- **LMBM3_50**: A 3‑bus system with ±10% voltage bounds, exhibiting two distinct local optima.
- **case9mod**: Modified 9‑bus system with reactive power generation lower bounds set to -3 MVars and demand scaled to 60%.
- **case118mod**: Modified 118‑bus system with generator real and reactive power bounds scaled by a factor of 7.

The above three cases are sourced from the [OptEnergyLocalOpt](https://webhomes.maths.ed.ac.uk/OptEnergy/LocalOpt/) repository.

- **case9mod_self**: A newly constructed 9‑bus benchmark with voltage bounds [0.98, 1.02], demand scaled to 60%, and reactive power lower bounds set to -12 MVars.


#### 2. `Cases(from MATPOWER 8.1 and PGLib)/` 

This folder contains standard benchmark cases from established repositories:

- **pglib_opf_case200_activ__api** – from PGLib OPF v23.07
- **pglib_opf_case240_pserc** – from PGLib OPF v23.07

The above PGLib cases correspond to standard congested and typical operating conditions.

- **case2383wp** – from MATPOWER 8.1
