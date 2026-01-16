This repository contains the code and results associated with the paper titled:
"A Novel Robust and Efficient Method for Computing High-Quality AC OPF Solutions Based on Trust-Tech and Holomorphic Embedding"



### 📌 Case Data (`../TableVI/casename/data/`)
Contains **ACOPF test case** collected from:
- [Matpower 8.1](https://matpower.org)
- [PGLIB-OPF](https://github.com/power-grid-lib/pglib-opf)
- [OptEnergyLocalOpt](https://webhomes.maths.ed.ac.uk/OptEnergy/LocalOpt/) : Cases satisfy multipue local optimal solutions (LOSs)
  
Contains **initial point** : 

The algorithm contains two types of initial points: random infeasible initial point and random feasible initial point.

**Random infeasible initial point** refer to initializations that are bound-feasible but infeasible with respect to the nonlinear AC OPF constraints, where all decision variables are randomly sampled within their prescribed lower and upper bounds.

**Random feasible initial point** are sampled directly from the feasible region of the AC OPF problem, i.e., all nonlinear power flow and operational constraints are satisfied at initialization.

---
## ⚙️ Environment Setup

### Required Toolboxes
- **MATPOWER 8.1**
- **MIPS** (version **1.5.1**)
- **IPOPT** (recommended version **v3.14+** or **MexIPOPT**) 
- **HSL** (**ma57 or ma97**)
- **KNITRO** (recommended version **15.0++**)

---



<!--

### ⚙️ IPOPT-Based Implementation (`/`)
- `ipoptopf_main.m`  
  Main function for **feasible region characterization** using IPOPT.

- `ipopt_convergence_main.m`  
  Automated **convergence testing** over all test cases.

### ⚙️ KNITRO-Based Implementation (`/`)
- `ipoptopf_main.m`  
  Main function for **feasible region characterization** using IPOPT.

- `ipopt_convergence_main.m`  
  Automated **convergence testing** over all test cases.

### ⚙️ MIPS-Based Implementation (`/`)
- `mips_opf_main.m`  
  Feasible region characterization using **MIPS**.

- `mips_convergence_main.m`  
  Automated convergence testing for all cases.

> **Note:** Due to the weaker convergence robustness of MIPS compared with IPOPT, **extensive random sampling** is required.

### 🧪 Reviewer Verification Script
- `reviewer #4_comment4.m`  
  Implements the **WB2 analytic feasible region formulation**, used to verify theoretical results in response to **Reviewer #1**.
-->
