#include "lasComm.h"
#ifdef MPI_VERSION
MPI_Comm LAS_COMM_WORLD = MPI_COMM_WORLD;
#endif
