# GraMPM
A C++ MPM simulation library for granular materials parallelised for CPUs with OpenMP. The library is a thin header-only library which gives you accesses to classes and methods which you can customize for your own simulations.

I wrote this mostly because I was curious about MPM and comparing it to SPH, as well as getting some practice implementing OOP practices.

The algorithms coded here are based on the description provided in the book: The Material Point Method by Zhang et al. (2017).
If you're looking to get started learning MPM, this is an ok start. But there are probably better books out there.

## Requirements

You will need a working C++ compiler (compilers known to work are g++, and nvc++).

## Getting started

After cloning this repo, build the executable with

```bash
make
```

which builds the `mpm.x` binary.

Run the example case with

```bash
./mpm.x
```

using the `OMP_NUM_THREADS` environment variable to control the number of parallel threads to use.

## API

### Controlling simulation geometry

The MPM grid extents are defined manually, and passed to the `GraMPM::MPM_system` class constructor, along with the grid cell size. These grid dimensions cannot be changed for the duration of the program.
When defining the grid extents, you must manually account for any buffer grid nodes that are needed for the mapping of data between the particles and grid. This is relevant moreso when using kernels that have radius greater than 1 grid cell (like the cubic B-spline).

The particles can be generated using the vector-like `p_push_back` method, which accepts an instance of the `particle` struct. As a minimum, the particle position, density, and mass need to be defined.

This is shown in the `src/main.cpp` example code.

### Defining boundary conditions

MPM boundary conditions are defined by the user as a function. The function signature is (using the momentum_boundary as an example):

```cpp
void momentum_boundary(GraMPM::MPM_system<double> &self, const size_t &timestep, const double &dt) {
```

where `self` is the instance of the MPM_system that the boundary conditions will be applied to, `timestep` is the current timestep, and `dt` is the timestep size (constant throughout the simulation).
The latter two can be used when, for example, wanting to define a moving boundary.

The definition of `momentum_boundary` and `force_boundary` in `src/main.cpp` shows how fully-fixed and free-slip boundaries can be defined in these functions.
These functions are then passed to an instance of the `MPM_system` on line 96 and 97 of `src/main.cpp`.

Note that OpenMP parallelism must be implemented within the user-defined function for these functions to be parallelised.

### Defining a stress update function

Like the force and momentum boundaries, the stress update of all the MPM particles is performed through a user-defined function. These functions need to have a signature:

```cpp
void stress_update(GraMPM::MPM_system<double> &self, const double dt)
```

Both Generalized Hooke's Law and Drucker Prager elasto-plasticity models are defined in `include/grampm_stress_update_functions.hpp` as template functions which shows how the strainrate and spinrate information can be used to update the stress state of each particle.
`src/main.cpp` is currently setup to use the Drucker Prager elasto-plasticity model to simulate the collapse of a granular column. To change the model used to Hooke's Law, change the function passed to `set_stress_update_function` on line 95 of `src/main.cpp` to `GraMPM::stress_update::hookes_law<double>`.

Parameters can be stored withing the instance of the `MPM_system` class by using the `set_stress_update_param`, which is a thin setter function which updates a `std::unordered_map`. 
The corresponding `get_stress_update_param` getter function can be used in the user-defined stress update functions.

Note that OpenMP parallelism must be implemented within the user-defined function for these functions to be parallelised.
