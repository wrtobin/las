#include "las.h"
#include "lasDebug.h"
#include <iostream>
namespace las
{
  Vec * LasCreateVec::createRHS(Mat*)
  {
    DBG(std::cerr << "[LAS] : createRHS(Mat*) not implemented for configured backend" << std::endl);
    return nullptr;
  }
  Vec * LasCreateVec::createLHS(Mat*)
  {
    DBG(std::cerr << "[LAS] : createLHS(Mat*) not implemented for configured backend" << std::endl);
    return nullptr;
  }
}
