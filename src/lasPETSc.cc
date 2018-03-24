#include "lasPETSc.h"
#include <petsc.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <iostream>
namespace las
{
  class PetscOps : public LasOps
  {
  public:
    virtual void zero(Mat * m);
    virtual void zero(Vec * v);
    virtual void assemble(Vec * v, int cnt, int * rws, double * vls);
    virtual void assemble(Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls);
    virtual void set(Vec * v, int cnt, int * rws, double * vls);
    virtual void set(Mat * v, int cntr, int * rws, int cntc, int * cls, double * vls);
    virtual double norm(Vec * v);
    virtual double dot(Vec * v0, Vec * v1);
    virtual void axpy(double a, Vec * x, Vec * y);
    virtual void get(Vec * v, double *& vls);
    virtual void restore(Vec * v, double *& vls);
  };
  // todo : create variant using nnz structure
  las::Mat * createPetscMatrix(int g, int l)
  {
    ::Mat * m = new ::Mat;
    MatCreateAIJ(PETSC_COMM_WORLD,
                 l,l,g,g,
                 sqrt(g), // dnz approximation
                 PETSC_NULL,
                 sqrt(g), // onz approximation
                 PETSC_NULL,
                 m);
    MatSetOption(*m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    return reinterpret_cast<las::Mat*>(m);
  }
  las::Vec * createPetscVector(int g, int l)
  {
    ::Vec * v = new ::Vec;
    VecCreateMPI(PETSC_COMM_WORLD,l,g,v);
    VecSetOption(*v,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    return reinterpret_cast<las::Vec*>(v);
  }
  inline ::Mat * getPetscMat(las::Mat * m)
  {
    return reinterpret_cast<::Mat*>(m);
  }
  inline ::Vec * getPetscVec(las::Vec * v)
  {
    return reinterpret_cast<::Vec*>(v);
  }
  LasOps * getPetscOps()
  {
    static PetscOps * ops = NULL;
    if(ops == NULL)
      ops = new PetscOps;
    return ops;
  }
  void PetscOps::zero(las::Mat * m)
  {
    MatZeroEntries(*getPetscMat(m));
  }
  void PetscOps::zero(las::Vec * v)
  {
    VecZeroEntries(*getPetscVec(v));
  }
  void PetscOps::assemble(las::Vec * v, int cnt, int * rws, double * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,ADD_VALUES);
  }
  void PetscOps::assemble(las::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,ADD_VALUES);
  }
  void PetscOps::set(las::Vec * v, int cnt, int * rws, double * vls)
  {
    VecSetValues(*getPetscVec(v),cnt,rws,vls,INSERT_VALUES);
  }
  void PetscOps::set(las::Mat * m, int cntr, int * rws, int cntc, int * cls, double * vls)
  {
    MatSetValues(*getPetscMat(m),cntr,rws,cntc,cls,vls,INSERT_VALUES);
  }
  double PetscOps::norm(las::Vec * v)
  {
    double n = 0.0;
    VecNorm(*getPetscVec(v),NORM_2,&n);
    return n;
  }
  double PetscOps::dot(las::Vec * v0, las::Vec * v1)
  {
    double d = 0.0;
    VecDot(*getPetscVec(v0),*getPetscVec(v1),&d);
    return d;
  }
  void PetscOps::axpy(double a, Vec * x, Vec * y)
  {
    VecAXPY(*getPetscVec(y),a,*getPetscVec(x));
  }
  void PetscOps::get(las::Vec * v, double *& vls)
  {
    VecGetArray(*getPetscVec(v),&vls);
  }
  void PetscOps::restore(las::Vec * v, double *& vls)
  {
    VecRestoreArray(*getPetscVec(v),&vls);
  }
  class PetscLUSolve : public LasSolve
  {
  public:
    virtual void solve(las::Mat * k, las::Vec * u, las::Vec * f)
    {
      ::Mat * pk = getPetscMat(k);
      ::Vec * pu = getPetscVec(u);
      ::Vec * pf = getPetscVec(f);
      ::KSP s;
      KSPCreate(PETSC_COMM_WORLD,&s);
      VecAssemblyBegin(*pf);
      VecAssemblyEnd(*pf);
      MatAssemblyBegin(*pk,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(*pk,MAT_FINAL_ASSEMBLY);
      KSPSetOperators(s,*pk,*pk);
      KSPSetFromOptions(s);
      KSPSolve(s,*pf,*pu);
      KSPDestroy(&s);
    }
  };
  LasSolve * createPetscLUSolve()
  {
    return new PetscLUSolve;
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
  LasSolve * createPetscQNSolve(void * a)
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
  LasMultiply * createPetscMultiply()
  {
    return new PetscMultiply;
  }
}
