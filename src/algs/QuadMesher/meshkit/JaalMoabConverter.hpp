#ifndef MOAB_QUADMESH_H
#define MOAB_QUADMESH_H

#include <meshkit/Mesh.hpp>

#ifdef HAVE_IMESH
#include <SimpleArray.hpp>
#include <iMesh.h>
#include <MBInterface.hpp>
#endif

using namespace Jaal;

class JaalMoabConverter
{
   public:

    void clear()
    {
        moabnode.clear();
        moabface.clear();
        jaalnode.clear();
        jaalface.clear();
    }
    //  Converts the mesh into MOAB data structures.
    int toMOAB(Mesh *mesh, iMesh_Instance &imesh, iBase_EntitySetHandle eset = 0);

    //  Fill the mesh from MOAB..
    Mesh *fromMOAB(iMesh_Instance imesh, iBase_EntitySetHandle eset = 0);


private:
    Mesh *jmesh;

    iBase_EntityHandle new_MOAB_Handle(iMesh_Instance imesh, Vertex *vertex);
    iBase_EntityHandle new_MOAB_Handle(iMesh_Instance imesh, Face *f);

    std::map<PNode, iBase_EntityHandle> moabnode;
    std::map<PFace, iBase_EntityHandle> moabface;
    std::map<iBase_EntityHandle, PNode> jaalnode;
    std::map<iBase_EntityHandle, PFace> jaalface;
};

#endif
