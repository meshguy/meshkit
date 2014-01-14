#ifndef MOAB_QUADMESH_H
#define MOAB_QUADMESH_H

#include "Mesh.hpp"

#include <meshkit/SimpleArray.hpp>
#include <meshkit/iMesh.hpp>

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
    int toMOAB(Jaal::Mesh *mesh, iMesh_Instance &imesh, iBase_EntitySetHandle eset = 0);

    //  Fill the mesh from MOAB..
    Jaal::Mesh *fromMOAB(iMesh_Instance imesh, Jaal::Mesh *m = NULL, iBase_EntitySetHandle eset = 0);

private:
    Jaal::Mesh *jmesh;

    iBase_EntityHandle new_MOAB_Handle(iMesh_Instance imesh, Vertex *vertex);
    iBase_EntityHandle new_MOAB_Handle(iMesh_Instance imesh, Face *f);

    std::map<PNode, iBase_EntityHandle> moabnode;
    std::map<PFace, iBase_EntityHandle> moabface;
    std::map<iBase_EntityHandle, PNode> jaalnode;
    std::map<iBase_EntityHandle, PFace> jaalface;
};

#endif
