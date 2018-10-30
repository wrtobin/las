#include "lasCSRCore.h"
#include "lasCSRBuilder.h"
#include <cassert>
namespace las
{
  class CSRFromNumbering : public CSRBuilder
  {
  protected:
    apf::Numbering * nm;
    apf::MeshEntity * ment;
    int nedofs;
    int ndofs;
  public:
    CSRFromNumbering(apf::Numbering * n, int nd)
      : CSRBuilder(nd,nd)
      , nm(n)
      , ment(NULL)
      , nedofs(0)
      , ndofs(nd)
    {
      assert(nm);
    }
    void run()
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
            add(dofs[ii],dofs[jj]);
      }
    }
  };
  class IdentityCSR : public CSRBuilder
  {
  protected:
    int ndofs;
  public:
    IdentityCSR(int ndofs)
      : CSRBuilder(ndofs,ndofs)
      , ndofs(ndofs)
    {
    }
    void run()
    {
      for(int i=0; i<ndofs;++i)
      {
        add(i,i);
      }
    }
  };
  Sparsity * createCSR(apf::Numbering * num, int ndofs)
  {
    CSRFromNumbering bldr(num,ndofs);
    bldr.run();
    return bldr.finalize();
  }
  Sparsity * createIdentityCSR(int ndofs)
  {
    IdentityCSR bldr(ndofs);
    bldr.run();
    return bldr.finalize();
  }
}
