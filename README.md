This repository contains the code and results associated with the paper titled:
"A Novel Robust and Efficient Method for Computing High-Quality AC OPF Solutions Based on Trust-Tech and Holomorphic Embedding"

### 📌 Case Data (`case/`)
Contains **22 state-of-the-art ACOPF test cases** collected from:
- **MATPOWER 7.1**
- **OptEnergy LocalOpt**
- **PGLIB-OPF** (for advanced stress testing; not included by default)

---

### ⚙️ IPOPT-Based Implementation (`ipopt_func/`)
- `ipoptopf_main.m`  
  Main function for **feasible region characterization** using IPOPT.

- `ipopt_convergence_main.m`  
  Automated **convergence testing** over all test cases.

---

### ⚙️ MIPS-Based Implementation (`mips_func/`)
- `mips_opf_main.m`  
  Feasible region characterization using **MIPS**.

- `mips_convergence_main.m`  
  Automated convergence testing for all cases.

> **Note:** Due to the weaker convergence robustness of MIPS compared with IPOPT, **extensive random sampling** is required.

---

### 🧪 Reviewer Verification Script
- `review1_comment4.m`  
  Implements the **WB2 analytic feasible region formulation**, used to verify theoretical results in response to **Reviewer #1**.

---

## 📊 Feasibility Data

Precomputed feasibility data used for validation and visualization:

- `WB2-feasible region.mat`  
  - Feasible voltage magnitudes  
  - Active/reactive power bounds  
  - Constraint satisfaction metrics  

- `case9mod-feasible region.mat`  
  - Feasible voltage magnitudes  
  - Active/reactive power bounds  
  - Constraint satisfaction metrics  

---

## ⚙️ Environment Setup

### Required Toolboxes
- **MATPOWER 4.1**
- **MIPS**
- **IPOPT** (recommended version **v3.12+**)

### MATLAB Path Configuration
```matlab
addpath(genpath('/path/to/matpower4.1'));
addpath('/path/to/ipopt_func');
addpath('/path/to/mips_func');
addpath('/path/to/case');
