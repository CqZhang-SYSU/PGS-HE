This repository contains the code and results associated with the paper titled:
"A Novel Robust and Efficient Method for Computing High-Quality AC OPF Solutions Based on Trust-Tech and Holomorphic Embedding"



### 📌 Case Data (`../Cases/`)

This directory contains all AC OPF test cases used in our experiments. The cases are collected from the following sources:

- [MATPOWER 8.1](https://matpower.org)
- [PGLIB-OPF](https://github.com/power-grid-lib/pglib-opf)
- [OptEnergyLocalOpt](https://webhomes.maths.ed.ac.uk/OptEnergy/LocalOpt/) – cases known to possess multiple local optimal solutions (LOSs).

That are:

- **LMBM3_50**: A 3‑bus system with ±10% voltage bounds, exhibiting two distinct local optima.
- **case9mod**: Modified 9‑bus system with reactive power generation lower bounds set to -3 MVars and demand scaled to 60%.
- **case118mod**: Modified 118‑bus system with generator real and reactive power bounds scaled by a factor of 7.
- **case9mod_self**: Newly constructed 9‑bus benchmark with voltage bounds [0.98, 1.02], demand scaled to 60%, and reactive power lower bounds set to -12 MVars.

The **PGLib cases** correspond to standard congested and typical operating conditions from PGLib OPF v23.07.


### 📌Table VI - AC OPF Solver Comparison (IPOPT vs. KNITRO)

This folder contains all the necessary data, source code, and benchmark results to reproduce **Table VI** in our paper. 

Table VI presents a comprehensive performance comparison between two state-of-the-art nonlinear programming solvers—**IPOPT** and **KNITRO**—applied to the Alternating Current Optimal Power Flow (AC OPF) problem. Crucially, both solvers are initialized from **identical initial points** to ensure a fair and rigorous evaluation of their convergence behavior and computational efficiency.

#### Initial Points Description

To investigate solver robustness and local convergence properties, we design two distinct categories of starting points:

#### 1. Random Infeasible Initial Points
- **Definition**: These points are **bound-feasible** but **violate** the nonlinear AC power flow equations and operational constraints (e.g., voltage magnitude, reactive power limits).
- **Generation**: All decision variables (voltage magnitudes `V`, angles `θ`, active/reactive power injections `P`/`Q`) are independently and uniformly sampled within their prescribed physical lower and upper bounds.
- **Purpose**: To evaluate the solver's ability to find a feasible stationary point from a strictly non-physical starting zone without prior constraint satisfaction.

#### 2. Random Feasible Initial Points
- **Definition**: These points strictly satisfy **all** nonlinear equality (power balance) and inequality constraints of the AC OPF model.
- **Generation**: Sampled directly from the feasible region using advanced constraint-satisfaction techniques (e.g., repeatedly solving a feasibility problem with perturbed targets).
- **Purpose**: To evaluate the local convergence rate and optimality gap when starting from a "warm" engineering guess.

> **Note**: For each category, we generate **N** distinct points (e.g., N=100). Both solvers are run from the **exact same set** of initial vectors to ensure a direct apples-to-apples comparison.

## ⚙️ Solvers Configuration

| Solver | Version | Key Options |
| :--- | :--- | :--- |
| **IPOPT** | 3.14+ | `linear_solver = mumps or ma57`,<br>`max_iter = 3000`,<br>`acceptable_compl_inf_tol = 0.001`,<br>`acceptable_constr_viol_tol = 5e-06`,<br>`acceptable_tol = 1e-08`,<br>`compl_inf_tol = 1e-05`,<br>`constr_viol_tol = 5e-06`,<br>`dual_inf_tol = 0.1`,<br>`tol = 1e-08`,<br>`mu_strategy = adaptive` |
| **KNITRO** | 15.0+ | `algorithm = 1` (Interior/Direct),<br>`maxit = 3000`,<br>`feastol = 1e-06`,<br>`feastol_abs = 5e-06`,<br>`ftol = 0.0001`,<br>`hessopt = 1`,<br>`opttol = 1e-06`,<br>`opttol_abs = 0.001`,<br>`xtol = 0.0001` |


## 🚀 How to Reproduce

### Prerequisites
- **For MATLAB**: Install [MATPOWER](https://matpower.org/) (or equivalent AC OPF framework) and the corresponding solver interfaces (IPOPT for MATLAB, KNITRO for MATLAB).


### Steps
1. **Clone the repository** and navigate to this folder:
   ```bash
   cd /path/to/Table_VI
---
## ⚙️ Environment Setup

### Required Toolboxes
- **MATPOWER 8.1** (https://matpower.org/)
- **MIPS** (version **1.5.1**)
- **MexIPOPT** (or **IPOPT** recommended version **v3.14+**) (https://github.com/ebertolazzi/mexIPOPT)
- **HSL** (**ma57 or ma97**)
- **KNITRO** (recommended version **15.0++**) (https://www.artelys.com/solvers/knitro/)

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
