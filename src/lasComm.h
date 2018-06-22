#ifndef LAS_COMM_H_
#define LAS_COMM_H_
#define USE_MPI 1
#if USE_MPI == 1
#include <mpi.h>
#else
typdef int MPI_Comm
#endif
extern MPI_Comm LAS_COMM_WORLD;
#endif
