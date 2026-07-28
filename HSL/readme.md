# HSL Math Library Dependency (IPOPT MATLAB Interface)

This directory contains precompiled Windows dynamic-link libraries (`.dll`) of the **HSL** (Harwell Subroutine Library) mathematical library. They are used to enable **IPOPT** to call the MA57 or MA97 sparse linear solver through the **MexIPOPT** interface in MATLAB.

## 📁 File List

| File | Description |
|------|-------------|
| `libhsl.dll` | Main HSL library, providing the MA57 sparse symmetric indefinite solver |
| `libmetis.dll` | METIS graph partitioning library for matrix reordering optimisation |
| `libgfortran-5.dll` | GCC Fortran runtime library (required dependency) |
| `libgomp-1.dll` | GCC OpenMP parallel support library (required dependency) |

> These are **Windows 64‑bit** precompiled binaries. `libhsl.dll` is specifically used when IPOPT is configured with the `ma57` (or `ma97`)linear solver.

## 🔧 Purpose

When solving large‑scale sparse nonlinear optimisation problems in MATLAB via IPOPT, you can specify the linear solver as `ma57` (or `ma97`):

```matlab
options.ipopt.linear_solver = 'ma57';
```

IPOPT will then dynamically load `libhsl.dll` to perform all sparse linear factorisations. The DLLs in this folder provide exactly that capability.

## 📦 Dependencies

- **MexIPOPT** – MATLAB interface for IPOPT  
  Recommended: [ebertolazzi/mexIPOPT](https://github.com/ebertolazzi/mexIPOPT)
- **MATLAB** R2018b or later (64‑bit)
- **Operating system**: Windows 64‑bit only

> For Linux or macOS, you must compile HSL yourself or switch to another linear solver (e.g., MUMPS).

## 🚀 Configuration

### 1. Add the DLL directory to the system path

In MATLAB, run the following command (assuming your current working directory is the project root):

```matlab
setenv('PATH', [getenv('PATH') ';' fullfile(pwd, 'HSL')]);
```

Alternatively, manually add the full path of the `HSL` folder to your system `PATH` environment variable.

### 2. Verify that the library is loaded correctly

Test with a simple IPOPT call:

```matlab
ipopt('example', struct('linear_solver', 'ma57'));
```

If IPOPT runs without reporting that `libhsl.dll` cannot be found, the configuration is successful.

## ⚠️ Important Notes

1. **License**: HSL is a commercial library and requires a valid licence. The DLLs provided here are intended **only** for reproducing the experiments of the associated paper. Please ensure you have the proper HSL licence before using them.
2. **Platform**: These `.dll` files are **Windows 64‑bit only**. For other operating systems, you need to obtain the appropriate HSL binaries from the HSL website or use an open‑source alternative like MUMPS.
3. **Runtime dependencies**: `libgfortran-5.dll` and `libgomp-1.dll` must be in the same directory as `libhsl.dll` or available in the system `PATH`, otherwise loading will fail.
4. **Version compatibility**: This HSL build has been verified with IPOPT 3.12+ and MexIPOPT. Other versions may cause symbol mismatches.

## 📖 Related Resources

- HSL official website: [http://www.hsl.rl.ac.uk/](http://www.hsl.rl.ac.uk/)
- MexIPOPT repository: [https://github.com/ebertolazzi/mexIPOPT](https://github.com/ebertolazzi/mexIPOPT)
- IPOPT documentation: [https://coin-or.github.io/Ipopt/](https://coin-or.github.io/Ipopt/)

## 📝 Additional Information

The files in this directory are part of the code accompanying the paper:

> *"A Novel Robust and Efficient Method for Computing High‑Quality AC OPF Solutions Based on Trust‑Tech and Holomorphic Embedding"*

They are provided for academic research purposes only. For any issues, please open an Issue in the main repository.
