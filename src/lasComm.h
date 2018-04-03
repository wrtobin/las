#ifndef LAS_COMM_H_
#define LAS_COMM_H_
#ifdef MPI_VERSION
#include <mpi.h>
extern MPI_Comm LAS_COMM_WORLD;
#else
#define MPI_Comm int
#endif
#endif
