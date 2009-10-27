#ifndef ITAPS_H
#define ITAPS_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>
#include <math.h>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include "SimpleArray.h"

using namespace std;


///////////////////////////////////////////////////////////////////////////////

class EdgeMesher;
class SurfaceMesher;
class VolumeMesher;

struct GeomMesh
{
    GeomMesh()
    {
        char *options = NULL;
        int err, optlen = 0;
        iMesh_newMesh(options, &mesh, &err, optlen);
        iMesh_getRootSet(mesh, &meshRootSet, &err);
        iRel_newAssoc(0, &assoc, &err, 0);
        surfMesher  = NULL;
        volMesher   = NULL;
        edgeMesher  = NULL;
    }

    void setGeometry(iGeom_Instance & g)
    {
        int err;
        geom = g;
        iRel_createAssociation(assoc,
                               geom, 0, iRel_IGEOM_IFACE,
                               mesh, 0, iRel_IMESH_IFACE, &relation00, &err);
        iRel_createAssociation(assoc,
                               geom, 0, iRel_IGEOM_IFACE,
                               mesh, 2, iRel_IMESH_IFACE, &relation02, &err);
        iGeom_getRootSet(geom, &geomRootSet, &err);
    }

    iMesh_Instance mesh;
    iGeom_Instance geom;
    iRel_Instance assoc;
    iRel_RelationHandle relation00, relation02;
    iBase_EntitySetHandle meshRootSet, geomRootSet;

    EdgeMesher    *edgeMesher;
    SurfaceMesher *surfMesher;
    VolumeMesher  *volMesher;
};

#include "Mesh.h"

#include "EdgeMesher.h"
#include "SurfMesher.h"
#include "VolMesher.h"

void generate_spectral_elements(GeomMesh &geomesh, int numNodes);

#endif

