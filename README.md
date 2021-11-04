# Simple Lattice Boltzmann

A toy implementation of the D2Q9 Lattice Boltzmann Method, generalised across
the input datafile. This project is originally intended to help me learn the
basics of LBM, along with performance gains from parallel implementations.

## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes.

### Prerequisites

- HDF5 Library. *Note: On the UQ `goliath` and `getafix` clusters, this
   library can be loaded in with*
```bash
module load hdf5-serial
module load cuda
```

### Building

The location of the HDF5 library must provided in the variable `HDF5HOME`, or
it must be manually included in the [makefile](./Makefile).

```makefile
# if not already defined, set the HDF5 Home
HDF5HOME ?= /home/alex/hdf5
```

Once the library is specified, the code can be built with `make cpu` or `make gpu`.

### Example

Once built, the file can be called with
```bash
./bin/main "taylor_green.h5"
./bin/main_gpu "taylor_green.h5"
```

## Authors

  - **Alex Muirhead** - [Alex-Muirhead](https://github.com/Alex-Muirhead/)

## Acknowledgments

  - My supervisor Chris Leonardi who put me on the path of LBM
