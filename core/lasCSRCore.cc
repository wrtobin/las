#include "lasCSRCore.h"
#include "lasCSRBuilder.h"
#include <cassert>
namespace las
{
  class CSRBuilder
  {
  protected:
    apf::Numbering * nm;
    apf::MeshEntity * ment;
    int nedofs;
    int ndofs;
    int nnz;
    std::vector<int> rws;
    std::vector<int> cls;
  public:
    CSRBuilder(apf::Numbering * n, int nd)
      : nm(n)
      , ment(NULL)
      , nedofs(0)
      , ndofs(nd)
      , nnz(0)
      , rws(ndofs+1)
      , cls(ndofs*ndofs)
    {
      assert(nm);
      rws.assign(ndofs+1,1);
      cls.assign(ndofs*ndofs,0);
    }
    void apply()
    {
      apf::Mesh * msh = apf::getMesh(apf::getField(nm));
      // find the highest dimension with mesh entities
      int dim = -1;
      for(dim = 3; dim >= 0; --dim)
        if(msh->count(dim) != 0)
          break;
      apf::MeshEntity * ent = NULL;
      apf::MeshIterator * it = NULL;
      for(it = msh->begin(dim); (ent = msh->iterate(it));)
      {
        apf::NewArray<int> dofs;
        int nedofs = apf::getElementNumbers(nm, ent, dofs);
        for(int ii = 0; ii < nedofs; ii++)
          for(int jj = 0; jj < nedofs; jj++)
            if(addNonzero(&rws[0],dofs[ii],&cls[0],dofs[jj],ndofs))
              nnz++;
      }
    }
    CSR * finalize()
    {
      return new CSR(ndofs,nnz,&rws[0],&cls[0]);
    }
  };
  CSR * createCSR(apf::Numbering * num, int ndofs)
  {
    CSRBuilder bldr(num,ndofs);
    bldr.apply();
    return bldr.finalize();
  }
}
