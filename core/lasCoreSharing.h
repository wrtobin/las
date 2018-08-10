#include <apfMesh.h>
#include <PCU.h>
namespace apf
{
  class AllSharing : public Sharing
  {
  public:
    AllSharing(Mesh * m);
    virtual bool isOwned(MeshEntity *) { return true; }
    virtual int getOwner(MeshEntity *) { return PCU_Comm_Self(); }
    virtual void getCopies(MeshEntity *, CopyArray &) {};
    virtual bool isShared(MeshEntity *);
  protected:
    Mesh * msh;
  };
}
