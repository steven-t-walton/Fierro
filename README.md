# Fierro

**Fierro** (LANL code number C21030) is a modern C++ code designed to simulate quasi-static solid mechanics problems and transient, compressible material dynamic problems with Lagrangian methods, which have meshes with constant mass elements that move with the material, or with Eulerian methods, which have stationary meshes.  **Fierro** is designed to aid a) material model research that has historically been done using commercial implicit and explicit finite element codes, b) numerical methods research, and c) computer science research.  **Fierro** will soon support user developed material models that adhere to several industry standard formats by using a C++ to Fortran interface to couple the model to the numerical solvers.  **Fierro** is built on the **ELEMENTS** library that supports a diverse suite of element types, including high-order elements, and quadrature rules. The mesh class within the **ELEMENTS** library is designed for efficient calculations on unstructured meshes and to minimize memory usage.  **Fierro** is designed to readily accommodate a range of numerical methods including continuous finite element, finite volume, and discontinuous Galerkin methods.  **Fierro** is designed to support explicit and implicit time integration methods as well as implicit optimization methods.  


## Computer implementation
**Fierro** is implemented in C++ following modern programming practices.  **Fierro** leverages the unique features of the **ELEMENTS** library, as such, this code serves as an example for solving a system of partial differential equations using the mesh class and geometric functions within the **ELEMENTS** library.  **Fierro** registers state at material points within the element, registers polynomial fields in the element, and registers kinematic variables at element vertices.  The routines for the state are designed for highly efficient computations and to minimize memory usage.  The loops are written to aid fine-grained parallelization and to allow vectorization. **Fierro** is a light-weight software application, and cleanly written following modern programming practices, so it useful for researching computer science based technologies for software performances portability.  

## Spatial discretization methods 
**Fierro** has an established conservative low-order Lagrangian finite element method, a low-order Lagrangian cell-centered finite volume method, and an arbitrary-order Lagrangian Discontinuous Galerkin method for solving the governing equations (e.g., mass, momentum, and energy evolution equations) for compressible material dynamics using unstructured hexahedral meshes.  These methods are combined with a multidirectional approximate Riemann solver (MARS) for improved accuracy on smooth flows and stable solutions near velocity discontinuities and large gradients in a flow. **Fierro** is designed for both low and high-order Lagrangian methods research and development but other types of numerical methods can be readily added to the code.  Likewise, other high-order methods can be studied within the code because it is built upon the **ELEMENTS** library that supports high-order elements and high-order quadrature rules.  Numerical methods are being added to **Fierro** to simulate quasi-static solid mechanics problems.  Likewise, direct Eulerian hydrodynamic methods can be investigated using **Fierro**.

## Temporal discretization methods 
**Fierro** supports a range of published multi-step time integration methods.  The code has an explicit multi-step Runge Kutta time integration method.  Implicit time integration methods can be implemented in **Fierro**.

## Material models  
The classical ideal gas model is the only material model implemented in the code, and it is useful for verification tests of the software and simple gas dynamic simulations.  A forthcoming version of the code will have C++ to Fortran interfaces to enable code users the ability to use their own material models and test them on quasi static problems or material dynamic applications.  The interfaces follow an industry standard format so that **Fierro** can be used for model research and development that has historically been done with commercial implicit or explicit finite element codes. 

## Cloning the code
If the user has set up ssh keys with GitHub, type
```
git clone --recursive ssh://git@github.com/lanl/Fierro.git
```
The code can also be cloned using
```
git clone --recursive https://github.com/lanl/Fierro.git
```

## Building the code
The user should create a new directory where the compiled code will reside.  
```
mkdir bin
```
Next, go to the folder and, for a default build, type
```
cmake ..
```
To compile the code type
```
make -j
```
The fierro executable will be in the bin/test folder. The following are the possible cmake build variables with their default values shown
```
BUILD_ELEMENTS=ON (Tells cmake whether to build the Elements libraries. Otherwise the user must compile them in Elements/build as instructed by the Elements readme)

BUILD_EXPLICIT_SOLVER=ON (Tells cmake whether to build the explicit solver components of Fierro)

BUILD_IMPLICIT_SOLVER=OFF (Tells cmake whether to build the implicit solver components of Fierro. This requires the user to build Trilinos in the folder Fierro/Trilinos/build)
```

## Building the explicit Lagrangian methods with Kokkos
Explicit Lagrangian codes are being added to the repository that are written using MATAR+Kokkos and run with fine-grained parallellism on multi-core CPUs and GPUs.  Build scripts are provided for each Lagrangian code, and those scripts follow those used in the [MATAR](https://github.com/lanl/MATAR/) GitHub repository. The scripts to build the Lagrangian codes (that use MATAR+Kokkos) are in the scripts folder.  The user must update the modules loaded by the build scripts (for the compiler etc.), and then type
```
source build-it.sh
```
The build-it.sh script sources the other scripts in the folder.  The compiled code will be in a folder (named after the explicit Lagrangian method) in the Fierro directory.  A range of scripts are provided for many architectures; however, they might not be correctly configured for the user's hardware.  The CPU architecture information needs to be listed if running with the Kokkos serial, OpenMP, and pthreads backends; GPU architecture information must be listed if using a Kokkos GPU backend. We refer the user to Kokkos compiling page to see the large list of compilation options,
```
https://github.com/kokkos/kokkos/wiki/Compiling
```
If the scripts fail to build a Lagrangian code, then carefully review the modules used and the computer architecture settings.  A more lenghtly discussion of the build scripts is provided in the MATAR GitHub repository.  



## Trilinos Dependencies to Install
OpenMPI
g++
gfortran
BLAS
LAPACK


## Running the Fierro code using explicit Lagrangian methods 
To run the fierro exectuable (see subsection above here to make the executable) go to bin/test and type
```
./fierro my_mesh.geo
```
The user must supply a mesh when executing the code, a range of meshes are provided in the meshes/ folder in the repository.


## Updating submodules
The ELEMENTS library and MATAR library can be updated to the newest release using
```
git submodule update --remote --merge
```







