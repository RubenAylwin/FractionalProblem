# FractionalProblem

[![Build & Test](https://github.com/RubenAylwin/FractionalProblem/actions/workflows/build.yml/badge.svg)](https://github.com/RubenAylwin/FractionalProblem/actions/workflows/build.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

C++ code for the resolution of problems involving the Riemann-Liouville fractional derivative and Reduced Basis for the same problem.

## Background

The Riemann-Liouville derivative generalizes ordinary differentiation to non-integer order and appears in models of anomalous diffusion, viscoelastic materials, and long-range interactions. Because it is a non-local operator, discretizing it directly produces dense (rather than sparse) stiffness matrices, and a new full finite-element solve is normally needed for every new choice of parameters. This project provides both a direct solver (`Problem 1`) and a Reduced Basis approach (`Problem 2`/`Problem 3`) that builds a low-dimensional approximation via a weak greedy algorithm, so repeated queries across a parameter range don't each require a full assembly.

## Dependencies

This project employs the following open-source libaries:
- [Boost C++ Libraries](https://www.boost.org/)
- [Eigen](https://eigen.tuxfamily.org/)
- [Sobol](https://people.sc.fsu.edu/~jburkardt/cpp_src/sobol/sobol.html) - included in src/QMC/
- [Halton](https://people.sc.fsu.edu/~jburkardt/cpp_src/halton/halton.html) - included in src/QMC/
- [ProgressBar](https://github.com/gipert/progressbar) - included in src/progressbar

## Usage
### Installation
Do:
```
mkdir build
cd build
cmake ../
make
```
and a binary named "VP" should appear.

### To solve a RL problem
The following options need to be given to the binary with double dash (--):
- "Problem" -> 1 for solving the RL problem..
- "dt" -> "poly", "trig" or "pw" indicates whether the coming expansion for the diffusion coefficient will be a taylor series, trigonometric expansion or piecewise partition of the coefficient.
- "qt" -> "poly", "trig" or "pw" indicates whether the coming expansion for the reaction coefficient will be a taylor series, trigonometric expansion or piecewise partition of the coefficient.
- "d" -> A space separated sequence of real numbers indicating the coefficients for an expansion of the diffusion coefficient.
- "q" -> A space separated sequence of real numbers indicating the coefficients for an expansion of the reaction coefficient.
- "rhs"-> A space separated sequence of real numbers indicating the coefficients for polynomial expansion of the right hand side.
- "fo"-> A number between 50 and 100 (non inclusive) indicating the half-order of the fractional derivative $\times 100$ (if we wish to consider, e.g., s = 1.4 we take --fo 70).
- "mn" -> Number of mesh elements (uniform mesh of the unit interval)

E.g. To solve the problem with the coefficients $d(x) = 4+\sin(2\pi x)$ and $q(x) = \cos(2\pi x)$ and the right hand side $f(x)=x(1-x)$, for $s=1.5$ and $200$ elements in the mesh, we take
```
./VP --Problem 1 --dt trig --qt trig --d 4. 2. --q 0. 0. 1. --rhs 0. 1. -1. --fo 75 --mn 200
```
This produces several files from which plots of the parameters and solution can be constructed, in the format (x, f(x)) for $x\in (0,1)$.

### To evaluate a RB for the RL problem
All previous parameters must be present. Moreover,
- "Problem" -> set to 2 to evaluate a RB.
- "dv" -> For each parameter in the expansion indicated in "d", indicate a $\pm$ variation for it.
- "qv" -> For each parameter in the expansion indicated in "q", indicate a $\pm$ variation for it.
- "rbp" -> Number of quadrature points for the greedy algorithm in each dimension.

E.g. To evaluate a RB for the problem with the coefficients $d(x) = \alpha+\beta\sin(2\pi x)$ and $q(x) = \gamma\cos(2\pi x)$, where $\alpha\in (3, 5),\ \beta\in (0, 1)$ and $\gamma\in(0, 2)$, and the right hand side $f(x)=x(1-x)$, for $s=1.5$, $200$ elements in the mesh and $10$ quadrature points in each dimension, we take
```
./VP --Problem 2 --dt trig --qt trig --d 4. .5 --q 0. 0. 1. --dv 1. .5 --qv 0. 0. 1. --rhs 0. 1. -1. --fo 75 --mn 200 --rbp 10
```
This will print out messages with the approximation information through the weak greedy algorithm. The environment variable "RIE_PROB" must be set to a value equal or higher than $1$.

### To evaluate RBs for the RL problem for a series of orders
All previous parameters must be present, except "fo", since now we add a new one for a list of orders. Moreover,
- "Problem" -> set to 3 to evaluate a sequence of RBs.
- "slist" -> List of orders to evaluate RBs on. Same format as "fo" before with integers between 50 and 100 (non inclusive) representing the half order $\times 100$.

E.g. To evaluate RBa for the problem with the coefficients $d(x) = \alpha+\beta\sin(2\pi x)$ and $q(x) = \gamma\cos(2\pi x)$, where $\alpha\in (3, 5),\ \beta\in (0, 1)$ and $\gamma\in(0, 2)$, and the right hand side $f(x)=x(1-x)$, for $s=1.2$, $1.5$ and $1.8$, $200$ elements in the mesh and $10$ quadrature points in each dimension, we take
```
./VP --Problem 3 --dt trig --qt trig --d 4. .5 --q 0. 0. 1. --dv 1. .5 --qv 0. 0. 1. --rhs 0. 1. -1. --slist 60 75 90 --mn 200 --rbp 10
```
This will print out messages with the approximation information for each order.
We recommend the following environment variable to be set to "true": PARALLELIZE_MATRIX_CONSTRUCTION. This allows parallel computation of stiffness matrices in the FEM.
