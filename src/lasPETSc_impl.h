#ifndef LAS_PETSC_IMPL_H_
#define LAS_PETSC_IMPL_H_
#include "lasComm.h"
#include "lasDebug.h"
#include <petsc.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <iostream>
namespace las
{
  inline void initPETScLAS(int * argc, char ** argv[], MPI_Comm cm = LAS_COMM_WORLD)
  {
    LAS_COMM_WORLD = cm;
    PETSC_COMM_WORLD = cm;
    PetscInitialize(argc,argv,PETSC_NULL,PETSC_NULL);
  }
  inline void finalizePETScLAS()
  {
    PetscFinalize();
  }
  inline ::Mat * getPetscMat(las::Mat * m)
  {
    return reinterpret_cast<::Mat*>(m);
  }
  inline ::Vec * getPetscVec(las::Vec * v)
  {
    return reinterpret_cast<::Vec*>(v);
  }
  inline void PetscOps::_zero(las::Mat * m)
  {
    MatZeroEntries(*getPetscMat(m));
  }
  inline void PetscOps::_zero(las::Vec * v)
  {
    VecZeroEntries(*getPetscVec(v));
  }
  inline void PetscOps::_assemble(las::Vec * v, int cnt, int * rws, double * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,ADD_VALUES);
  }
  inline void PetscOps::_assemble(las::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,ADD_VALUES);
  }
  inline void PetscOps::_set(las::Vec * v, int cnt, int * rws, double * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,INSERT_VALUES);
  }
  inline void PetscOps::_set(las::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,INSERT_VALUES);
  }
  inline double PetscOps::_norm(las::Vec * v)
  {
    double n = 0.0;
    VecNorm(*getPetscVec(v),NORM_2,&n);
    return n;
  }
  inline double PetscOps::_dot(las::Vec * v0, las::Vec * v1)
  {
    double d = 0.0;
    VecDot(*getPetscVec(v0),*getPetscVec(v1),&d);
    return d;
  }
  inline void PetscOps::_axpy(double a, Vec * x, Vec * y)
  {
    VecAXPY(*getPetscVec(y),a,*getPetscVec(x));
  }
  inline void PetscOps::_get(las::Vec * v, double *& vls)
  {
    VecGetArray(*getPetscVec(v),&vls);
  }
  inline void PetscOps::_restore(las::Vec * v, double *& vls)
  {
    VecRestoreArray(*getPetscVec(v),&vls);
  }
  // todo : create variant using nnz structure
  inline las::Mat * createPetscMatrix(int g, int l, MPI_Comm cm = LAS_COMM_WORLD)
  {
    ::Mat * m = new ::Mat;
    MatCreateAIJ(cm,
                 l,l,g,g,
                 sqrt(g), // dnz approximation
                 PETSC_NULL,
                 sqrt(g), // onz approximation
                 PETSC_NULL,
                 m);
    MatSetOption(*m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    return reinterpret_cast<las::Mat*>(m);
  }
  inline las::Vec * createPetscVector(int g, int l, MPI_Comm cm = LAS_COMM_WORLD)
  {
    ::Vec * v = new ::Vec;
    VecCreateMPI(cm,l,g,v);
    VecSetOption(*v,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    return reinterpret_cast<las::Vec*>(v);
  }
  inline LasOps<PetscOps> * getPetscOps()
  {
    static PetscOps * ops = NULL;
    if(ops == NULL)
      ops = new PetscOps;
    return ops;
  }
  class PetscLUSolve : public LasSolve
  {
    MPI_Comm cm;
    ::KSP ksp;
  public:
    PetscLUSolve(MPI_Comm c = LAS_COMM_WORLD)
      : cm(c)
    {
      KSPCreate(cm,&ksp);
    }
    ~PetscLUSolve()
    {
      KSPDestroy(&ksp);
    }
    virtual void solve(las::Mat * k, las::Vec * u, las::Vec * f)
    {
      ::Mat * pk = getPetscMat(k);
      ::Vec * pu = getPetscVec(u);
      ::Vec * pf = getPetscVec(f);
      DBG(
        MPI_Comm kspCm = MPI_COMM_NULL;
        MPI_Comm matCm = MPI_COMM_NULL;
        PetscObjectGetComm((PetscObject)ksp,&kspCm);
        PetscObjectGetComm((PetscObject)pk,&matCm);
        assert(kspCm == matCm);
        )
        VecAssemblyBegin(*pf);
      VecAssemblyEnd(*pf);
      MatAssemblyBegin(*pk,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(*pk,MAT_FINAL_ASSEMBLY);
      KSPSetOperators(ksp,*pk,*pk);
      KSPSetFromOptions(ksp);
      KSPSolve(ksp,*pf,*pu);
    }
  };
  inline LasSolve * createPetscLUSolve(MPI_Comm cm = LAS_COMM_WORLD)
  {
    return new PetscLUSolve(cm);
  }
  /*
    PetscErrorCode PetscIterate(SNES ,::Vec x,::Vec f,void * i)
    {
    LasOps * ops = getPetscOps();
    // x vector needs to be applied to simulation prior to assembling force vector...
    // typically we update our fields just after a solve, which gets us into position
    // for the next iteration
    void * hd = i;
    apf::Numbering * fld = NULL;
    hd = byte_read(hd,fld);
    double * sol = NULL;
    las::Vec * lx = reinterpret_cast<las::Vec*>(&x);
    ops->get(lx,sol);
    amsi::AccumOp op;
    amsi::ApplyVector appsol(fld,apf::getField(fld),&sol[0],3,&op);
    appsol.run();
    // get the iteration object to run to assemble the force vector
    amsi::Iteration * itr = NULL;
    hd = byte_read(hd,itr);
    // also we need to assemble the values INTO f, which is provided by amsi, we can't assume
    // that it is the same f we passed into the PetscQNSolve::solve() function,
    // so we pass a pointer to the internal iteration pointer to the vector to assemble into,
    // and then we set that to point to the f in this function
    las::Vec ** lf = NULL;
    hd = byte_read(hd,lf); // lf points to f inside of itr
    *lf = reinterpret_cast<las::Vec*>(&f);
    itr->iterate();
    // set *lf to NULL since the local f might no longer be valid after this function ends
    *lf = NULL;
    return(0);
    }
  */
  class PetscQNSolve : public LasSolve
  {
  private:
    void * args;
  public:
    PetscQNSolve(void * a) : args(a)
    {}
    virtual void solve(las::Mat * k, las::Vec * u, las::Vec * f)
    {
      ::Mat * pk = getPetscMat(k);
      ::Vec * pu = getPetscVec(u);
      ::Vec * pf = getPetscVec(f);
      ::SNES snes;
      SNESCreate(PETSC_COMM_WORLD,&snes);
      // set snes options
      SNESSetType(snes,SNESQN);
      SNESSetJacobian(snes,*pk,*pk,NULL,NULL);
      //SNESSetFunction(snes,*pf,&PetscIterate,args);
      SNESSolve(snes,*pf,*pu);
      int it = -1;
      SNESGetIterationNumber(snes,&it);
      std::cout << "BFGS QN Solve complete in " << it << " iterations." << std::endl;
      SNESDestroy(&snes);
    }
  };
  inline LasSolve * createPetscQNSolve(void * a)
  {
    return new PetscQNSolve(a);
  }
  class PetscMultiply : public LasMultiply
  {
    void exec(Mat * x, Vec * a, Vec * b)
    {
      ::Mat * px = getPetscMat(x);
      ::Vec * pa = getPetscVec(a);
      ::Vec * pb = getPetscVec(b);
      MatMult(*px,*pa,*pb);
    }
  };
  inline LasMultiply * createPetscMultiply()
  {
    return new PetscMultiply;
  }
}
#endif
