ReduVCC v1.0
=====

Computing minimum vertex clique covers using data reduction.

![build](https://github.com/darrenstrash/fast-vcc/actions/workflows/c-cpp.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Installation Notes
=====

Before you can start you need to install the following software packages:

- Scons (http://www.scons.org/)
- OpenMPI (https://www.open-mpi.org/) 

After installing the packages, run `scons program=vcc variant=optimized` to build.

## Running

This package contains 5 different algorithms: **Chalupa**, **Redu**, **ReduVCC**, **BnR**, and **EdgeBnR**, which can be run as follows:

**Chalupa**
`./optimized/vcc --preconfiguration=fsocial --k=2 --mis=<independent set size> --run_type="Chalupa" <input graph>`

**Redu**
`./optimized/vcc --preconfiguration=fsocial --k=2 --run_type="Redu" <input graph>`

**ReduVCC**
`./optimized/vcc --preconfiguration=fsocial --k=2 --mis=<independent set size> --run_type="ReduVCC" <input graph>`

**BnR**
`./optimized/vcc --preconfiguration=fsocial --k=2 --run_type="bnr" <input graph>`

**EdgeBnR**
`./optimized/vcc --preconfiguration=fsocial --k=2 --run_type="edge_bnr" <input graph>`

## Input Format

ReduVCC uses **The unweighted METIS format**, which consists of

   `<# vertices> <# edges>`

   followed by `<# vertices>` lines of space-separated vertices,  where the `i`-th line consists of 
   all neighbors of `i`. All vertices range from `1` to `<# vertices>`

Loops and directed edges are not supported.

## Data Sets

You will find an example graph in the directory `examples`
