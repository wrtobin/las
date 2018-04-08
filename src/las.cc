#include "las.h"
#include <iostream>
namespace las
{
  Vec * createRHS(Mat * m)
  {
    DBG(std::cerr << "[LAS] : createRHS(Mat*) not implemented for configured backend" << std::endl);
  }
  Vec * createLHS(Mat * m)
  {
    DBG(std::cerr << "[LAS] : createLHS(Mat*) not implemented for configured backend" << std::endl);
  }
}
