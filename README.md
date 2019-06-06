# Burgers-Equation
Implementation of the split form of Burgers' Equation in C++ using the Discontinuous Galerkin Scheme.
Currently the scheme is in split, weak form, with collocated Gauss-Lobatto nodes for quadrature and interpolation.
The code is templated so it can be run in any dimension. To build the code, run the following:

```
rm -rf ./build
cd build
cmake -DDEAL_II_DIR=/path/to/dealii/installation ..
make -j<N>
./main
```
