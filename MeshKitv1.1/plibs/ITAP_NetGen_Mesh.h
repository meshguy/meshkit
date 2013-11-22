#ifndef ITAP_NETGEN_H
#define ITAP_NETGEN_H

#include <meshing.hpp>
#include <iGeom.h>

#include "ITAP_NetGen_EdgeMesh.h"

#include "ITAP_NetGen_SurfMesh.h"

int ITAP_NetGen_GenerateMesh ( iGeom_Instance & geom, Mesh &mesh, 
                               int perfstepsstart, int perfstepsend);

#endif
