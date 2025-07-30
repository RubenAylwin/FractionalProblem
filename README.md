# FractionalProblem
C++ code for the resolution of problems involving the Riemann-Liouville fractional derivative and Reduced Basis for the same problem.

## Dependencies

This project employs the following open-source libaries:
- [Boost C++ Libraries](https://www.boost.org/)
- [Eigen](https://eigen.tuxfamily.org/) - Location must be indicated in CMakeLists.txt
- [Sobol](https://people.sc.fsu.edu/~jburkardt/cpp_src/sobol/sobol.html) - included in src/QMC/
- [Halton](https://people.sc.fsu.edu/~jburkardt/cpp_src/halton/halton.html) - included in src/QMC/
- [ProgressBar](https://github.com/gipert/progressbar) - included in src/progressbar

## Notes
- CMakeLists.txt should indicate where to find the files for the Eigen library.

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
- "Problem" -> 1 for solving the RL problem, 2 for constructing a Reduced Basis.
- "dt" -> "poly", "trig" or "pw" indicates whether the coming expansion for the diffusion coefficient will be a taylor series, trigonometric expansion or piecewise partition of the coefficient.
- "qt" -> "poly", "trig" or "pw" indicates whether the coming expansion for the reaction coefficient will be a taylor series, trigonometric expansion or piecewise partition of the coefficient.
- "d" -> A space separated sequence of real numbers indicating the coefficients for an expansion of the diffusion coefficient.
- "q" -> A space separated sequence of real numbers indicating the coefficients for an expansion of the reaction coefficient.
- "rhs"-> A space separated sequence of real numbers indicating the coefficients for polynomial expansion of the right hand side.
- "fo"-> A number between 50 and 100 (non inclusive) indicating the half-order of the fractional derivative x100 (if we wish to consider, e.g., s = 1.4 we take --fo 70).
- "mn" -> Number of mesh elements (uniform mesh of the unit interval).

E.g. To solve the problem with the coefficients $d(x) = 4+\sin(2\pi x)$ and $q(x) = \cos(2\pi x)$ and the right hand side $f(x)=x(1-x)$, for $s=1.5$ and $200$ elements in the mesh, we take
```
./VP --Problem 1 --dt trig --qt trig --d 4. 2. --q 0. 0. 1. --rhs 0. 1. -1. --fo 75 --mn 200
```
This produces several files from which plots of the parameters and solution can be constructed, in the format (x, f(x)) for $x\in (0,1)$.