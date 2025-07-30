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
