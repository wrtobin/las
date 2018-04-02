## Overview

LAS provides a zero-overhead API to use the general operations of a linear algebra library operating on large (parallel) matrices and vectors. Currently the library only allows double-values as scalars.

The solution process and multiplication process are not zero-overhead calls since the cost of these function calls is significantly less than the cost of the solve or multiplication procedure.

The operations currently exposed through the (simple) API are:
 - zero a vector or matrix
 - assemble (add) values to multiple locations in a vector or matrix
 - set values at multiple locations in a vector or matrix
 - calculate the norm of a vector
 - calculate the dot product of two vectors
 - a general axpy procedure (y = a*x + y) for vectors x and y and scalar a
 - retrieval of the (local) underlying array of a vector
 - restore the (local) underlying array of a vector
 - matrix-vector multiplication (Ax = y) (not zero-overhead)
 - solving a system (Ax = b) (not zero-overhead)

## Development Guidelines

In order to implement a new backend for the API, a developer must create (minimally) two files:
lasXXX.h and lasXXX_impl.h. The first of which declares a class as follows:

class xxxOps : public LasOps<xxxOps>
{ ... };

**The implementation of the class interface functions may occur inline in the function declaration in lasXXX.h for ONLY FOR simple functions which DO NOT DIRECTLY CALL ANY BACKEND API functions.**

**Anything which explicitly makes calls to the backend API is implemented in the lasXXX_impl.h file, and is defined as an inline member function.**

Following these guidelines (see lasSparskit.h and lasSparskit_impl.h for an example) allows the use of these functions to be zero-overhead when building in release mode using -O3.

## Style Guide

Until I get around to an actual style guide refer to the SCOREC/core [style guide](https://github.com/SCOREC/core/blob/develop/STYLE.md) and mimic the rest of the codebase. Specifically: spaces not tabs, 2-space indentation, no trailing whitespace, no empty lines.

Always have a line-break for scope-based brackets, since the code style doesn't allow empty lines (vertical whitespace), this causes code blocks to stand out well.

```
for (i = 0; i < 3; ++i)
{
  b[i] = a[i];
}
```

## Dependencies
Currently the library has only optional dependencies, but isn't particularly useful without at least one linear algebra backend, though sparskit is included and built by default.
 * [MPI](https://www.open-mpi.org/) when using the linear algebra backends in parallel on HPC architectures. Basically always build with this if you have it available by setting CC and CXX for the cmake scripts.
 * [SCOREC/core](https://github.com/SCOREC/core) is used to interface with unstructured meshes / tensor fields / and dof numberings of those tensor fields to develop algorithms to define the sparsity patterns of matrices.
 * [PETSc](https://www.mcs.anl.gov/petsc/) provides a robust parallel linear algebra library for use in MPI-parallel applications
 * [cuSparse](http://docs.nvidia.com/cuda/cusparse/index.html) is a sparse matrix linear algebra library for use on CUDA-enable devices
