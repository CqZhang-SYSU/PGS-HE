This repository contains the code and results associated with the paper titled:
"A Novel Robust and Efficient Method for Computing High-Quality AC OPF Solutions Based on Trust-Tech and Holomorphic Embedding"



### 📌 Case Data (`data/`)
Contains **ACOPF test cases** collected from:
- [Matpower 8.1](https://matpower.org)
- [PGLIB-OPF](https://github.com/power-grid-lib/pglib-opf)
- [OptEnergyLocalOpt](https://webhomes.maths.ed.ac.uk/OptEnergy/LocalOpt/) : Cases satify multipue local optimal solutions (LOSs)
  
Contains **initial point**  from:

---
## ⚙️ Environment Setup

### Required Toolboxes
- **MATPOWER 8.1**
- **MIPS**
- **IPOPT** (recommended version **v3.14+**)

---

### MATLAB Path Configuration
```matlab
addpath(genpath('/path/to/matpower8.1'));
addpath('/path/to/data');
```
---
### ⚙️ IPOPT-Based Implementation (`/`)
- `ipoptopf_main.m`  
  Main function for **feasible region characterization** using IPOPT.

- `ipopt_convergence_main.m`  
  Automated **convergence testing** over all test cases.

---

### ⚙️ KNITRO-Based Implementation (`/`)
- `ipoptopf_main.m`  
  Main function for **feasible region characterization** using IPOPT.

- `ipopt_convergence_main.m`  
  Automated **convergence testing** over all test cases.

---

### ⚙️ MIPS-Based Implementation (`/`)
- `mips_opf_main.m`  
  Feasible region characterization using **MIPS**.

- `mips_convergence_main.m`  
  Automated convergence testing for all cases.

> **Note:** Due to the weaker convergence robustness of MIPS compared with IPOPT, **extensive random sampling** is required.

---

### 🧪 Reviewer Verification Script
- `reviewer #4_comment4.m`  
  Implements the **WB2 analytic feasible region formulation**, used to verify theoretical results in response to **Reviewer #1**.








