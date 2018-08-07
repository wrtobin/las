#include "lasCore.h"
#include "lasConfig.h"
#include <apf.h>
#include <apfAggregateNumbering.h>
#include <apfMesh.h>
namespace las
{
  template <class B>
  Mat * createMatrix(apf::AggNumbering * num, MPI_Comm cm, int bld_id)
  {
    apf::Field * fld = apf::getField(num);
    apf::Mesh * msh = apf::getMesh(fld);
    int dim = msh->getDimension();
    bool is_par = cm != MPI_COMM_SELF;
    int own_nds = 0;
    for(int dd = 0; dd < dim; ++dd)
      own_nds += apf::countOwnedNodes(msh,dd,apf::getShape(fld));
    int cmps = apf::countComponents(fld);
    int bs = apf::countDOFsPerBlock(num);
    Sparsity * sps = createSparsity<B>(num,cm,bld_id);
    // get num local dofs from the field
    // get block size from field as well... really need a numbering for this
    // get the sparsity pattern from the field also...
    LasCreateMat * mb = getMatBuilder<B>(bld_id);
    return mb->create(lcl,bs,sps,cm);
  }
  void destroyMatrix(Mat * m)
  {

  }
}
