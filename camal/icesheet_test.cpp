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
#include "MBRange.hpp"
#include "MBSkinner.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBiMesh.hpp"

#include <iostream>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <map>
#include <assert.h>
#include <string.h>

extern bool debug;

#define ERRORR(a) {if (iBase_SUCCESS != err) {printf(a); return 1;}}
#define MBERRORR(a) {if (MB_SUCCESS != rval) {std::cerr << a << std::endl; return 1;}}

int main( int argc, char *argv[] )
{
  // Check command line arg
  std::string input_filename, output_filename;
  double mesh_size = -1.0;
  int mesh_intervals = -1;
  bool force_intervals = false;

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <input_mesh_filename> <output_mesh_filename> [-s <uniform_size>] [-i <uniform_int>] [-f] "
              << std::endl
              << "  -s <uniform_size> = mesh with this size" << std::endl
              << "  -i <uniform_int> = mesh curves with this # intervals" << std::endl
              << "  -f = force these size/interval settings even if geometry has interval settings" << std::endl
              << std::endl;
    
    return 0;
  }
  else {
    input_filename = argv[1];
    output_filename = argv[2];
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
  
  // initialize the mesh interface instances
  int i, err;
  iGeom_Instance geom = NULL;
  iGeom_newGeom(0, &geom, &err, 0);
  assert(iBase_SUCCESS == err);
  
  iMesh_Instance mesh = NULL;
  iMesh_newMesh(0, &mesh, &err, 0);
  assert(iBase_SUCCESS == err);

  iMesh_Instance out_mesh = NULL;
  iMesh_newMesh(0, &out_mesh, &err, 0);
  assert(iBase_SUCCESS == err);
  
  iRel_Instance relate;
  iRel_newRel(0, &relate, &err, 0);
  assert(iBase_SUCCESS == err);

  iBase_EntitySetHandle mesh_root_set;
  iMesh_getRootSet(mesh, &mesh_root_set, &err);
  if (iBase_SUCCESS != err) {
    printf("Failed to return a root set.\n");
    return 1;
  }

  // create an association pair
  iRel_PairHandle classification;
  iRel_createPair(relate, geom, 0, iRel_IGEOM_IFACE,
			 out_mesh, 1, iRel_IMESH_IFACE,
			 &classification, &err );
  ERRORR("ERROR : can not create an assoc pair.");

  // read surface mesh
  iMesh_load(mesh, mesh_root_set, input_filename.c_str(), NULL, &err, input_filename.size(), 0);
  ERRORR("ERROR : can not load input mesh from file");

  // get surfaces, edges and verticies
  iBase_EntityHandle *triangles = NULL;
  int triangles_size, triangles_alloc = 0;
  iMesh_getEntities(mesh, mesh_root_set, iBase_ALL_TYPES,
		    iMesh_TRIANGLE, &triangles, &triangles_alloc, 
		    &triangles_size, &err);
  ERRORR("Failed to get triangles.\n");

  MBRange surface_ents, edge_ents, loop_range;
  for (i = 0; i < triangles_size; i++) {
    surface_ents.insert(reinterpret_cast<MBEntityHandle> (triangles[i]));
  }

  MBSkinner tool(reinterpret_cast<MBiMesh*> (mesh)->mbImpl);
  MBErrorCode rval = tool.find_skin(surface_ents, 1, edge_ents);
  MBERRORR("Skinner failure.");

  iBase_EntityHandle* edge_entities = (iBase_EntityHandle*)malloc(edge_ents.size()*sizeof(iBase_EntityHandle));
  MBEntityHandle *tmp_ptr = reinterpret_cast<MBEntityHandle*>(edge_entities);
  std::copy(edge_ents.begin(), edge_ents.end(), tmp_ptr);

  iBase_EntityHandle* verts = NULL;
  int verts_alloc = 0, verts_size;
  iMesh_getEntAdj(mesh, edge_entities[0], iBase_VERTEX,
		  &verts, &verts_alloc, &verts_size, &err);
  ERRORR("Problem getting edge vertex.\n");

  if (debug) {
    double x, y, z;
    iMesh_getVtxCoord(mesh, verts[0], &x, &y, &z, &err);
    std::cout << "vertex coords=" << x << ", " << y << ", " << z << std::endl;
  }
  
  // make all geometry mesh sets and add entities
  iBase_EntitySetHandle sets[3];
  for (i = 0; i < 3; i++) {
    iMesh_createEntSet(mesh, false, &sets[i], &err);
    ERRORR("Problem creating geometry entityset.\n");
  }
  iMesh_addEntArrToSet(mesh, verts, 1, sets[0], &err);
  ERRORR("Failed to add vertex in entity set\n");
  iMesh_addEntArrToSet(mesh, edge_entities, edge_ents.size(), sets[1], &err);
  ERRORR("Failed to add edges in entity set.\n");
  iMesh_addEntArrToSet(mesh, triangles, triangles_size, sets[2], &err);
  ERRORR("Failed to add triangles in entity set.\n");
  iMesh_addPrntChld(mesh, sets[2], sets[1], &err);
  ERRORR("Problem add parent in entity_sets_test.\n");
  iMesh_addPrntChld(mesh, sets[1], sets[0], &err);
  ERRORR("Problem add parent in entity_sets_test.\n");

  if (debug) {
    iMesh_save(mesh, sets[1], "boundary_edges.vtk", 0, &err, 18, 0);
    if (iBase_SUCCESS != err) {
      std::cerr << "ERROR saving mesh to boundary_edges.vtk" << std::endl;
      return 1;
    }
  }

  // set geometry tag to geometry sets
  iBase_TagHandle geom_tag_handle;
  const char* geom_tag_name = "GEOM_DIMENSION";
  iBase_TagHandle id_tag_handle;
  const char* id_tag_name = "GLOBAL_ID";
  iMesh_getTagHandle(mesh, geom_tag_name, &geom_tag_handle, &err, 14);
  ERRORR("Couldn't get geometry tag handle.\n");
  iMesh_getTagHandle(mesh, id_tag_name, &id_tag_handle, &err, 9);
  ERRORR("Couldn't get geometry tag handle.\n");
  for (i = 0; i < 3; i++) {
    iMesh_setEntSetIntData(mesh, sets[i], geom_tag_handle, i, &err);
    ERRORR("Failed to set geom tag values.\n");
    iMesh_setEntSetIntData(mesh, sets[i], id_tag_handle, 1, &err);
    ERRORR("Failed to set geom tag values.\n");
  }
  // set some senses

  if (debug) {
    iBase_EntityHandle *entities;
    int entities_alloc, entities_size;
    for (i = 0; i < iMesh_ALL_TOPOLOGIES; i++) {
      entities = NULL;
      entities_alloc = 0;
      iMesh_getEntities(mesh, mesh_root_set, iBase_ALL_TYPES,
			i, &entities, &entities_alloc, 
			&entities_size, &err);
      ERRORR("Failed to get entities\n");
      std::cout << "type=" << i << ", number=" << entities_size << std::endl;
    }
    
    iBase_EntitySetHandle *esets = NULL;
    int esets_alloc = 0, esets_size;
    iMesh_getEntSets(mesh, mesh_root_set, 1, &esets, &esets_alloc, &esets_size, &err);
    ERRORR("Failed to get entity sets\n");
    std::cout << "entity set number=" << esets_size << std::endl;
  }
  moab::GeomTopoTool geomTool(reinterpret_cast<MBiMesh*> (mesh)->mbImpl);
  int sense = 1;//
  std::vector<int> senses;
  senses.push_back(sense);
  std::vector<MBEntityHandle> setFaces;
  //MBEntityHandle setFace((MBEntityHandle) sets[2]);
  setFaces.push_back( (MBEntityHandle)sets[2] );
  //MBEntityHandle setCurve((MBEntityHandle) sets[1]);
  // GeomTopoTool::set_senses(EntityHandle edge,  std::vector<EntityHandle, std::vector<int> &senses)
  geomTool.set_senses((MBEntityHandle)sets[1], setFaces, senses);// only one set of dimension 1

  // write the mesh
  iMesh_save(mesh, mesh_root_set, "surface_geom.h5m", 0, &err, 16, 0);
  ERRORR("Couldn't get geometry tag handle.\n");

  // read mesh geometry
  iGeom_load(geom, "surface_geom.h5m", 0, &err, 16, 0);
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR : can not load a input mesh geometry from "
	      << input_filename << std::endl;
    return 1;
  }
  
  // mesh the geometry
  CMEL cmel(geom, out_mesh, relate, classification);
  bool success = cmel.mesh_geometry(mesh_size, mesh_intervals, force_intervals);
  if (!success) {
    std::cerr << "Problems meshing." << std::endl;
  }

  // write the mesh
  iBase_EntitySetHandle out_root_set;
  iMesh_getRootSet(out_mesh, &out_root_set, &err);
  assert(iBase_SUCCESS == err);
  iMesh_save(out_mesh, out_root_set, output_filename.c_str(), 0, &err, output_filename.length(), 0);
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR saving mesh to " << input_filename << std::endl;
    return 1;
  }

  return !success;
}
