Redu3ECC v1.0
=====

Computing minimum edge clique covers (ECCs) using data **redu**ction for ECC, problem **redu**ction to VCC, and data **redu**ction for VCC.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


Citation
========

This is this code for experiments in the article 

**"Solving Edge Clique Cover Exactly via Synergistic Data Reduction"
by Anthony Hevia, Benjamin Kallus, Summer McClintic, Samantha Reisner, Darren Strash, and Johnathan Wilson.**

The article is to appear at the 31st Annual European Symposium on Algorithms (ESA 2023), which takes place September 4 to September 6, 2023. Until it is officially published, there is a preprint available at (https://arxiv.org/abs/2306.17804), which is citable as

```
@misc{hevia2023solving,
      title={Solving Edge Clique Cover Exactly via Synergistic Data Reduction},
      author={Anthony Hevia and Benjamin Kallus and Summer McClintic and Samantha Reisner and Darren Strash and Johnathan Wilson},
      year={2023},
      eprint={2306.17804},
      archivePrefix={arXiv},
      primaryClass={cs.DS}
}
```

Installation Notes
=====

Before you can start you need to install the following software packages:

- Scons (http://www.scons.org/)
- OpenMPI (https://www.open-mpi.org/) 
- Gurobi (https://www.gurobi.com/)

After installing the packages, run `scons program=redu3ecc variant=optimized` to build.

## Running

This package contains 3 different algorithms: **Redu3BnR**, **ReduILP**, **ReduIG**, which can be run as follows:

**Redu3BnR**
`./optimized/redu3ecc --preconfiguration=fsocial --k=2 --run_type="Redu3BnR" <input graph>`

**ReduILP**
`./optimized/redu3ecc --preconfiguration=fsocial --k=2 --run_type="ReduILP" <input graph>`

**ReduIG**
`./optimized/redu3ecc --preconfiguration=fsocial --k=2 [--mis=<independent set size>] --run_type="ReduIG" <input graph>`

Each algorithm can be provided with time limits, including:
- a time limit for the solver (in seconds) with `--solver_time_limit=<time in seconds>` and 
- a time limit for the overall running time with `--time_limit=<time in seconds>`

## Input Format

Redu3ECC uses **The unweighted METIS format**, which consists of

   `<# vertices> <# edges>`

   followed by `<# vertices>` lines of space-separated vertices,  where the `i`-th line consists of 
   all neighbors of `i`. All vertices range from `1` to `<# vertices>`

Loops and directed edges are not supported.

## Data Sets

You will find an example graph in the directory `examples`
