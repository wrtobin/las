#ifndef LAS_PETSC_IMPL_H_
#define LAS_PETSC_IMPL_H_
#include "lasComm.h"
#include "lasDebug.h"
#include "lasNNZ.h"
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
  inline void PetscOps::_assemble(las::Vec * v, int cnt, int * rws, scalar * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,ADD_VALUES);
  }
  inline void PetscOps::_assemble(las::Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,ADD_VALUES);
  }
  inline void PetscOps::_set(las::Vec * v, int cnt, int * rws, scalar * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,INSERT_VALUES);
  }
  inline void PetscOps::_set(las::Mat * m, int cntr, int * rws, int cntc, int * cls, scalar * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,INSERT_VALUES);
  }
  inline scalar PetscOps::_norm(las::Vec * v)
  {
    scalar n = 0.0;
    VecNorm(*getPetscVec(v),NORM_2,&n);
    return n;
  }
  inline scalar PetscOps::_dot(las::Vec * v0, las::Vec * v1)
  {
    scalar d = 0.0;
    VecDot(*getPetscVec(v0),*getPetscVec(v1),&d);
    return d;
  }
  inline void PetscOps::_axpy(scalar a, Vec * x, Vec * y)
  {
    VecAXPY(*getPetscVec(y),a,*getPetscVec(x));
  }
  inline void PetscOps::_get(las::Vec * v, scalar *& vls)
  {
    VecGetArray(*getPetscVec(v),&vls);
  }
  inline void PetscOps::_restore(las::Vec * v, scalar *& vls)
  {
    VecRestoreArray(*getPetscVec(v),&vls);
  }
  inline las::Mat * createPetscMatrix(unsigned l, unsigned g, unsigned bs = 1, Sparsity * sprs = nullptr, MPI_Comm cm = LAS_COMM_WORLD)
  {
    bool have_sparsity = sprs != nullptr;
    NNZ * nnz = have_sparsity ? reinterpret_cast<NNZ*>(sprs) : nullptr;
    const char * mat_tps[][2] = { {MATSEQAIJ, MATSEQBAIJ}, {MATMPIAIJ, MATMPIBAIJ} };
    bool is_par = cm != MPI_COMM_SELF;
    bool blk = bs > 1;
    ::Mat * m = new ::Mat;
    MatCreate(cm,m);
    MatSetType(*m, mat_tps[is_par][blk]);
    MatSetSizes(*m, l, l, g, g);
    MatSetBlockSize(*m, bs);
    if(have_sparsity)
    {
      if(!is_par)
      {
        if(nnz->dnnz.size() == l)
        {
          if(nnz->onnz.size() == l)
          {
            for(unsigned ii = 0; ii < l; ++ii)
              nnz->dnnz[ii] += nnz->onnz[ii];
          }
          if(!blk)
            MatSeqAIJSetPreallocation(*m,0,&nnz->dnnz[0]);
          else if(blk)
            MatSeqBAIJSetPreallocation(*m,bs,0,&nnz->dnnz[0]);
        }
      }
      else
      {
        if(nnz->dnnz.size() == l && nnz->onnz.size() == l)
        {
          if(!blk)
            MatMPIAIJSetPreallocation(*m,0,&nnz->dnnz[0],0,&nnz->onnz[0]);
          else if(blk)
            MatMPIBAIJSetPreallocation(*m,bs,0,&nnz->dnnz[0],0,&nnz->onnz[0]);
        }
      }
      MatSetOption(*m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
    }
    if(!blk)
      MatSetOption(*m,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);
    return reinterpret_cast<las::Mat*>(m);
  }
  inline las::Vec * createLHSVec(las::Mat * m)
  {
    ::Vec * v = new ::Vec;
    MatCreateVecs(*getPetscMat(m),v,nullptr);
    return reinterpret_cast<las::Vec*>(v);
  }
  inline las::Vec * createRHSVec(las::Mat * m)
  {
    ::Vec * v = new ::Vec;
    MatCreateVecs(*getPetscMat(m),nullptr,v);
    return reinterpret_cast<las::Vec*>(v);
  }
  inline las::Vec * createPetscVector(unsigned l, unsigned g, unsigned bs, MPI_Comm cm = LAS_COMM_WORLD)
  {
    ::Vec * v = new ::Vec;
    VecCreateMPI(cm,l,g,v);
    VecSetBlockSize(*v,bs);
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
    scalar * sol = NULL;
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
    //void * args;
  public:
    PetscQNSolve(void *)// : args(a)
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
  inline mat_builder getPetscMatBuilder()
  {
    return
      [](unsigned lcl,
         unsigned gbl,
         unsigned bs,
         Sparsity * sprs,
         MPI_Comm cm)
   {
     return createPetscMatrix(lcl,gbl,bs,sprs,cm);
   };
  }
  inline vec_builder getPetscVecBuilder()
  {
    return
      [](unsigned lcl,
         unsigned gbl,
         unsigned bs,
         MPI_Comm cm)
    {
      return createPetscVector(lcl,gbl,bs,cm);
    };
  }
}
#endif
