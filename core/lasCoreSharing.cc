#include "lasCoreSharing.h"
namespace apf
{
  AllSharing::AllSharing(Mesh * m)
    : msh(m)
  { }
  bool AllSharing::isShared(MeshEntity * e)
  {
    return msh->isShared(e);
  }
}
