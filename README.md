LAS provides a zero-overhead API to use the general operations of a linear algebra library operating on large (parallel) matrices and vectors. Currently the library only allows double-values as scalars.

Yhe solution process and multiplication process are not zero-overhead calls since the cost of these function calls is significantly less than the cost of the solve or multiplication proceedure.

The operations currently exposed through the (simple) API are:
 - zero a vector or matrix
 - assemble (add) values to multiple locations in a vector or matrix
 - set values at multiple locations in a vector or matrix
 - calculate the norm of a vector
 - calculate the dot product of two vectors
 - a general axpy proceedure (y = a*x + y) for vectors x and y and scalar a
 - retrieval of the (local) underlying array of a vector
 - restore the (local) underlying array of a vector

In order to implement a new backend for the API, a developer must create (minimally) two files:
lasXXX.h and lasXXX_impl.h. The first of which declares a class as follows:

class xxxOps : public LasOps<xxxOps>
{ ... };

The implementation of the class interface functions happens either inline in the function declaration in lasXXX.h for simple functions which DO NOT DIRECTLY CALL ANY BACKEND API functions.

Anything which explicitly makes calls to the backend API is implemented in the lasXXX_impl.h file, and is implemented as an inline function.

Following these guidelines (see lasSparskit.h and lasSparskit_impl.h for an example) allows the use of these functions to be zero-overhead when building in release mode using -O3.
