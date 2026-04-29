# Data-driven discovery of bounded polynomial ODEs

Welcome to SILAS! This repository contains code and data used to produce results in the following paper:

[paper-here]

Most of the code is in MATLAB, with some Python scripts to generate data and post-process results. Details on the repository structure and a short description of each file are given below. To replicate the examples of the paper, you can simply run the scripts `ex01_discovery.m`, `ex02_main.m` and `ex03_main.m`, adjusting parameters as needed to match those given in the paper.

The code is provided "as-is" and without any guarantees. See the license file for more information.

**Authors:** Albert Alcalde, Giovanni Fantuzzi (FAU Erlangen-Nuernberg)

## Dependencies
To run the code, you must first install the following MATLAB packages:
* YALMIP
* MOSEK (or another SDP solver compatible with YALMIP)
* ChebFun

## Cite us
If you find this code useful in your work and publish your results, please cite the arXiv preprint above. An sample .bib code follows.

```
    @article{FantuzziAlcalde2026,
        AUTHOR = {Alcalde, Albert and Fantuzzi, Giovanni},
        TITLE = {Data-driven discovery of polynomial ODEs with provably bounded solutions},
        YEAR = {2026},
        HOWPUBLISHED = {arXiv:????.????},
        URL = {https://arxiv.org/abs/????.????}
    }
```

## Repository Structure

```text
silas-main/
├── LICENSE.txt
├── README.md
├── ex01-cubic/
├── ex02-dysts/
├── ex03-pde/
└── utils/
```

## Files

### Root

| File | Description |
|---|---|
| `LICENSE.txt` | Project license. |
| `README.md` | Original project documentation. |

### `ex01-cubic`

| File | Description |
|---|---|
| `ex01_simulate.m` | Simulates the cubic example and generates data. |
| `ex01_discovery.m` | Runs model discovery on the cubic example. |

### `ex02-dysts`

| File | Description |
|---|---|
| `ex02_main.m` | Main script for learning models from Dysts dynamical-system data. |
| `ex02_preprocess_data.m` | Prepares raw Dysts data for training and testing. |
| `ex02_attractor_plots.m` | Generates attractor plots for Dysts systems. |

### `ex02-dysts/private`

| File | Description |
|---|---|
| `exact_AtmosphericRegime.m` | Exact reference equations for the AtmosphericRegime system. |
| `exact_Bouali2.m` | Exact reference equations for the Bouali2 system. |
| `exact_GuckenheimerHolmes.m` | Exact reference equations for the Guckenheimer-Holmes system. |
| `exact_HyperRossler.m` | Exact reference equations for the HyperRossler system. |
| `exact_SprottA.m` | Exact reference equations for the SprottA system. |
| `exact_Thomas.m` | Exact reference equations for the Thomas system. |
| `exact_YuWang.m` | Exact reference equations for the YuWang system. |
| `plot_lyap.m` | Plots learned Lyapunov functions. |

### `ex02-dysts/python-scripts`

| File | Description |
|---|---|
| `data_generate.ipynb` | Notebook for generating Dysts datasets. |
| `process_results.ipynb` | Notebook for processing and summarizing experiment results. |

### `ex03-pde`

| File | Description |
|---|---|
| `ex03_main.m` | Main MATLAB script for the PDE reduced-order modeling example. |
| `ex03_helpers.py` | Python helper functions for PDE simulations and preprocessing. |
| `ex03_simulation.py` | Generates PDE simulation data. |

### `ex03-pde/private`

| File | Description |
|---|---|
| `lyap_plotProjDeg2.m` | Plots projected degree-2 Lyapunov functions. |
| `plot_attractors.m` | Plots attractors for the PDE reduced model. |
| `plot_pde_solution.m` | Visualizes PDE solution trajectories. |

### `utils`

| File | Description |
|---|---|
| `build_initial_lyap_constrained.m` | Builds an initial constrained Lyapunov candidate. |
| `build_initial_lyap_penalized.m` | Builds an initial penalized Lyapunov candidate. |
| `build_lyap_model.m` | Constructs a Lyapunov model. |
| `build_model_lyap.m` | Builds the combined learned model and Lyapunov structure. |
| `cheb_evaluate.m` | Evaluates Chebyshev polynomial expansions. |
| `cheb_gradient.m` | Computes gradients of Chebyshev expansions. |
| `cheb_mommat.m` | Builds Chebyshev moment matrices. |
| `cheb_multiply.m` | Multiplies Chebyshev polynomial terms. |
| `lyap_initialize.m` | Initializes Lyapunov-function learning. |
| `lyap_load.m` | Loads saved Lyapunov-function results. |
| `make_parameters.m` | Creates parameter settings for experiments. |
| `model_cheb2mon.m` | Converts models from Chebyshev to monomial form. |
| `model_evaluate.m` | Evaluates learned dynamical models. |
| `model_initialize.m` | Initializes model-learning variables. |
| `model_learn.m` | Learns a dynamical model from data. |
| `model_load.m` | Loads saved learned models. |
| `multi_indices.m` | Generates polynomial multi-index sets. |
| `save_results.m` | Saves experiment outputs. |
| `select_fekete.m` | Selects interpolation/sample points using Fekete-style criteria. |
| `viridis.m` | Provides the Viridis colormap for plots. |
