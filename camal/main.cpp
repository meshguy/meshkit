/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
/**
 * \file main.cpp
 *
 * \brief main for camel executable
 *
 */

#include "iRel.h"
#include "iMesh.h"
#include "iGeom.h"

#include "camel.hpp"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <math.h>
#include <assert.h>
#include <string.h>

bool debug = false;

int main( int argc, char *argv[] )
{
    // Check command line arg
  std::string geom_filename, mesh_filename;
  double mesh_size = -1.0;
  int mesh_intervals = -1;
  bool force_intervals = false;

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <geom_filename> <mesh_filename> [-s <uniform_size>] [-i <uniform_int>] [-f] "
              << std::endl
              << "  -s <uniform_size> = mesh with this size" << std::endl
              << "  -i <uniform_int> = mesh curves with this # intervals" << std::endl
              << "  -f = force these size/interval settings even if geometry has interval settings" << std::endl
              << std::endl;
    
    return 0;
  }
  else {
    geom_filename = argv[1];
    mesh_filename = argv[2];
    int argno = 3;
    while (argno < argc) {
      if (!strcmp(argv[argno], "-s")) {
        argno++;
        sscanf(argv[argno], "%lf", &mesh_size);
        argno++;
      }
      else if (!strcmp(argv[argno], "-i")) {
        argno++;
        sscanf(argv[argno], "%d", &mesh_intervals);
        argno++;
      }
      else if (!strcmp(argv[argno], "-f")) {
        argno++;
        force_intervals = true;
      }
      else {
        std::cerr << "Unrecognized option: " << argv[argno] << std::endl;
        return 1;
      }
    }
  }
  
    // initialize the geometry and mesh interface instances
  int err;
  
  iGeom_Instance geom;
  iGeom_newGeom(0, &geom, &err, 0);
  assert(iBase_SUCCESS == err);
  
  iMesh_Instance mesh;
  iMesh_newMesh(0, &mesh, &err, 0);
  assert(iBase_SUCCESS == err);
  
  iRel_Instance relate;
  iRel_newAssoc(0, &relate, &err, 0);
  assert(iBase_SUCCESS == err);

    // create an association pair
  iRel_RelationHandle classification;
  iRel_createAssociation( relate, geom, 0, iRel_IGEOM_IFACE,
                                  mesh, 1, iRel_IMESH_IFACE,
                                  &classification, &err );
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR : can not create an assoc pair." << std::endl;
    return 1;
  }
  

    // read geometry
  iGeom_load( geom, geom_filename.c_str(), 0, &err, geom_filename.size(), 0 );
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR : can not load a geometry from " << geom_filename << std::endl;
    return 1;
  }

    // mesh the geometry
  CMEL cmel(geom, mesh, relate, classification);
  bool success = cmel.mesh_geometry(mesh_size, mesh_intervals, force_intervals);
  if (!success)
    std::cerr << "Problems meshing." << std::endl;
  
    // write the mesh
  iBase_EntitySetHandle root;
  iMesh_getRootSet( mesh, &root, &err );
  assert(iBase_SUCCESS == err);
  iMesh_save(mesh, root, mesh_filename.c_str(), 0, &err, mesh_filename.length(), 0);
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR saving mesh to " << mesh_filename << std::endl;
    return 1;
  }

  return !success;
}
