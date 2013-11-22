#ifndef NG_EDGMESH_H
#define NG_EDGMESH_H

#include <meshing.hpp>

#include <iGeom.h>

#include <vector>

using namespace std;
using namespace netgen;

#include "SimpleArray.hpp"
#include "SearchUV.h"

#define TCL_OK 0
#define TCL_ERROR 1

#define DIVIDE_EDGE_SECTIONS 1000
#define IGNORE_CURVE_LENGTH  1e-4

class ITAP_NetGen_EdgeMesh
{
   public:
      ITAP_NetGen_EdgeMesh(iGeom_Instance &g, Mesh &m)
       {   geometry = g; mesh = m; global_meshedge_id = 1;}

     void execute();
   private:
     Mesh  mesh;
     iGeom_Instance geometry;

     iBase_TagHandle  id_tag;

     int global_meshedge_id;

     void discretize_close_edge(iBase_EntityHandle edgeHandle, int numEdges);
     void discretize_open_edge( iBase_EntityHandle edgeHandle, int numEdges);
     void discretize( iBase_EntityHandle edgeHandle, int numEdges);
};

#endif


