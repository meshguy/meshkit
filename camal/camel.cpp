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
//
// mesh something with camel
//
#include "camel.hpp"

#include "camal_interface.hpp"

/*
#include "iRel.h"
#include "iMesh.h"
#include "iGeom.h"
#include <math.h>
*/
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <assert.h>
#include <string.h>

extern bool debug ;

CMEL::CMEL( iGeom_Instance geom, 
            iMesh_Instance mesh, 
            iRel_Instance relate,
            iRel_PairHandle rel ) 
{
  int err;
  geomIface = geom;
  if (NULL == mesh) iMesh_newMesh(NULL, &meshIface, &err, 0);
  else meshIface = mesh;
  
  if (NULL == relate) iRel_newRel(NULL, &relateIface, &err, 0);
  else relateIface =  relate; 
  
  if (NULL == relate || NULL == rel)
    iRel_createPair( relateIface, geom, 1, iRel_IGEOM_IFACE,
                                         mesh, 1, iRel_IMESH_IFACE,
                                         &pairHandle, &err );
  else
    pairHandle = rel;
}

bool CMEL::mesh_geometry(double mesh_size, int mesh_intervals, const bool force_intervals,
      const bool quadMesh)
{
    // get the maximal-dimension entities
  std::vector<iBase_EntityHandle> new_entities;
  int result = iBase_SUCCESS;
  
  for (int dim = 3; dim >= 0; dim--) {
    iBase_EntityHandle *these_gents = NULL;
    int these_gents_size = 0, these_gents_alloc = 0;
        // get all entities of this dimension
    iGeom_getEntities(geomIface, 0, dim, 
                      &these_gents, &these_gents_alloc, &these_gents_size, &result);
    if (iBase_SUCCESS != result) {
      std::cerr << "Trouble getting gentities of dimension" << dim << std::endl;
      return false;
    }
  
    if (0 != these_gents_size)
      return mesh_entities(these_gents, these_gents_size, mesh_size, 
                           mesh_intervals, force_intervals, quadMesh);
  }

    // if we got here, we didn't mesh anything
  return false;
}

bool CMEL::mesh_entities(iBase_EntityHandle *gentities, 
                         const int num_geom_entities, 
                         double mesh_size, int mesh_intervals,
                         const bool force_intervals,
                         const bool quadMesh)
{
  bool result = true;
  std::vector<iBase_EntityHandle> new_entities;
  for (int i = 0; i < num_geom_entities; i++) {
    new_entities.clear();
    bool tmp_result = mesh_entity(gentities[i], mesh_size, mesh_intervals,
                                  force_intervals, new_entities, quadMesh);
    if (!tmp_result) result = tmp_result;
  }
  
  return result;
}

bool CMEL::mesh_boundary(iBase_EntityHandle gentity, 
                         double mesh_size, int mesh_intervals, const bool force_intervals,
                         iBase_EntityHandle **bounding_ents, int *bounding_ent_size,
                         const bool quadMesh)
{
    // get dimension
  int dim = -1;
  int gtype, result;
  iGeom_getEntType(geomIface, gentity, &gtype, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting dimension of gentity." << std::endl;
    return false;
  }
  else 
    dim = gtype;

  if (iBase_VERTEX == gtype) return true;
  
    // get bounding entities; optionally pass back these entities
  iBase_EntityHandle *tmp_gents = NULL;
  iBase_EntityHandle **adj_gents = (bounding_ents ? bounding_ents : &tmp_gents);
  int tmp_size = 0;
  int *adj_gents_size = (bounding_ent_size ? bounding_ent_size : &tmp_size);
  iGeom_getEntAdj(geomIface, gentity, gtype-1, 
                  adj_gents, adj_gents_size, adj_gents_size, &result);
  if (iBase_SUCCESS != result) return false;
  
  bool success = true;
  
  for (int g = 0; g < *adj_gents_size; g++) {
    
      // mesh this entity if it isn't already meshed
    std::vector<iBase_EntityHandle> ments_vec;

    if (!is_meshed((*adj_gents)[g])) {
      bool tmp_success = mesh_entity((*adj_gents)[g], mesh_size, 
                   mesh_intervals, force_intervals, ments_vec, quadMesh);
      if (!tmp_success) {
        success = tmp_success;
        continue;
      }
    }
  }

  if (NULL == bounding_ents) free(*adj_gents);

  return success;
}
  
bool CMEL::mesh_entity(iBase_EntityHandle gentity, 
                       double mesh_size, int mesh_intervals, 
                       const bool force_intervals,
                       std::vector<iBase_EntityHandle> &new_entities,
                       const bool quadMesh)
{
    // check to make sure this entity isn't meshed
  if (is_meshed(gentity)) return true;

  // if gentity is surface (one surface) and we have some new boundary established, mesh
  // using a new method (rather quirky design)

   int gtype, result;
   iGeom_getEntType(geomIface, gentity, &gtype, &result);
   if (iBase_SUCCESS != result) {
      std::cerr << "Trouble getting dimension of gentity." << std::endl;
      return false;
   }

   if (iBase_FACE == gtype && this->newBoundary.size() > 0)
   {
      return CAMAL_mesh_trimmed_surface(this, gentity,
             mesh_size,  new_entities,
             newBoundary,  quadMesh  );
   }
  //
  bool success = mesh_boundary(gentity, mesh_size, mesh_intervals, force_intervals, NULL, NULL, quadMesh);

    // get the size or interval #
  if (!force_intervals) {
    int this_int;
    double this_size;
    bool success = get_attrib_intervals(gentity, this_size, this_int);
    if (success) {
      if (0 < this_size) mesh_size = this_size;
      else if (0 < this_int) mesh_intervals = this_int;
    }
  }

  success = CAMAL_mesh_entity(this, gentity, mesh_size, mesh_intervals, 
                              force_intervals, new_entities, quadMesh);
  assert(!success || is_meshed(gentity));
  
  return success;
}

bool CMEL::is_meshed(iBase_EntityHandle gentity) 
{
    // get dimension of this gentity
  int dimension = -1;

  int result, this_type;
  iGeom_getEntType(geomIface, gentity, &this_type, &result);
  if (iBase_SUCCESS != result) return false;
  
  dimension = this_type - iBase_VERTEX;

  iBase_EntityHandle *ments = NULL;
  int ments_size = 0, ments_alloc = 0;
  iRel_getEntEntArrRelation(relateIface, pairHandle, 
                               gentity, false, 
                               &ments, &ments_alloc, &ments_size, &result);
  if (result != iBase_TAG_NOT_FOUND && iBase_SUCCESS != result) {
    char descr[100];
    iRel_getDescription ( relateIface, descr, &result, 99);
    std::cerr << "is_meshed: returned error, message: "
              << descr << std::endl;
    assert(false);
  }
  
  return (0 != ments_size ? true : false);
}

bool CMEL::get_attrib_intervals(iBase_EntityHandle gentity,
                                double &mesh_size,
                                int &mesh_interval) 
{
    // look for mesh size or interval settings
  mesh_size = -1.0;
  mesh_interval = -1;
  int result;
  
  bool have_interval = false;
  iBase_TagHandle int_tag;
  iGeom_getTagHandle(geomIface, "MESH_INTERVAL", &int_tag, &result, strlen("MESH_INTERVAL"));
  if (iBase_SUCCESS == result) {
    int interval;
    iGeom_getIntData( geomIface, gentity, int_tag, &interval, &result );
    if (iBase_SUCCESS == result) {
      mesh_interval = interval;
      have_interval = true;
    }
  }
  
  bool have_size = false;
  iBase_TagHandle size_tag;
  iGeom_getTagHandle(geomIface, "MESH_SIZE", &size_tag, &result, strlen("MESH_SIZE"));
  if (iBase_SUCCESS == result) {
    double size;
    iGeom_getDblData( geomIface, gentity, size_tag, &size, &result );
    if (iBase_SUCCESS == result) {
      mesh_size = size;
      have_size = true;
    }
  }
  
  return have_interval || have_size;
}

/*
bool CMEL::get_bounding_mesh_elements(iBase_EntityHandle gentity, 
                                      double mesh_size,
                                      int mesh_intervals, const bool force_intervals,
                                      std::vector<iBase_EntityHandle> &bdy_mesh, 
                                      std::vector<int> &bdy_orientations) 
{
    // get dimension of this gentity
  int dimension = -1;
  int result, ent_type;
  iGeom_getEntType(geomIface, gentity, &ent_type, &result);
  if (iBase_SUCCESS != result) return false;
  
  dimension = ent_type - iBase_VERTEX;

    // mesh boundary; get back bounding entities too
  iBase_EntityHandle *adj_gents = NULL;
  int adj_gents_size = 0;
  bool success = mesh_boundary(gentity, mesh_size, mesh_intervals, force_intervals,
                               &adj_gents, &adj_gents_size);
  if (!success) return false;
  
  if (0 == dimension) return true;
  else if (-1 == dimension) return false;

  std::vector<iBase_EntityHandle> new_ents;
  
  for (int g = 0; g < adj_gents_size; g++) {
    
      // mesh this entity if it isn't already meshed
    iBase_EntityHandle *ments = NULL;
    int ments_size = 0;

      // get the mesh on that entity
    iRel_getGeomOwnedMesh(relateIface, relationHandle, 
                           adj_gents[g], dimension-1,
                           &ments, &ments_size, &ments_size, &result);
    if (iBase_SUCCESS != result) return false;
    
      // decide whether to store forward or reverse flags
    int this_sense = 1;
    if (dimension > 1) {
      if (dimension == 2) this_sense = iGeom_getEgFcSense(geomIface, adj_gents[g], gentity);
      else this_sense = iGeom_getEntNrmlSense(geomIface, adj_gents[g], gentity);
    }
    int old_size = bdy_mesh.size();

      // might need to reverse list if it's a curve and reverse flag is set
    if (-1 == this_sense && dimension == 2)
      std::reverse(ments, ments+ments_size);
      
    std::copy(ments, ments+ments_size, std::back_inserter(bdy_mesh));

    free(ments);
    
    bdy_orientations.resize(bdy_mesh.size());
    std::fill(&bdy_orientations[old_size], &bdy_orientations[bdy_mesh.size()], this_sense);
  }

  return success;
}

bool CMEL::get_bounding_mesh_connect(iBase_EntityHandle gentity,
                                     double mesh_size,
                                     int mesh_intervals, const bool force_intervals,
                                     std::vector<iBase_EntityHandle> &bdy_verts,
                                     std::vector<double> &coords, 
                                     std::vector<int> &connect) 
{
  std::vector<iBase_EntityHandle> bdy_mesh;
  std::vector<int> bdy_orientations;
  bool success = get_bounding_mesh_elements(gentity, mesh_size, mesh_intervals, 
                                            force_intervals, bdy_mesh, bdy_orientations);
  if (!success) {
    std::cerr << "Couldn't get unordered bounded mesh." << std::endl;
    return success;
  }
  
    // arrange vertex coordinates in a sequential array, and connectivity of bdy mesh
    // in integer-indexed array
  iBase_EntityHandle *bdy_verts_vec = NULL;
  int bdy_verts_size = 0;
  int *offsets = NULL;
  int offsets_size = 0;
  
  int result;
  iMesh_getEntArrAdj(meshIface, &bdy_mesh[0], bdy_mesh.size(),
                     iMesh_VERTEX,
                     &bdy_verts_vec, &bdy_verts_size, &bdy_verts_size,
                     &offsets, &offsets_size, &offsets_size, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't get bounding vertices." << std::endl;
    return false;
  }

    // make unique; need to create copy because we need the original for
    // assigning indices in connect
  bdy_verts.resize(bdy_verts_size);
  std::copy(bdy_verts_vec, bdy_verts_vec+bdy_verts_size, &bdy_verts[0]);
  std::sort(bdy_verts.begin(), bdy_verts.end());
  std::vector<iBase_EntityHandle>::iterator new_end = 
    std::unique(bdy_verts.begin(), bdy_verts.end());
  bdy_verts.erase(new_end, bdy_verts.end());

    // assign id sequence & get ids in connectivity array
  std::vector<int> indices(bdy_verts.size());
  for (unsigned int i = 0; i < bdy_verts.size(); i++) indices[i] = i;
  iBase_TagHandle index_tag;
  connect.resize(bdy_verts_size);
  int *connect_ptr = &connect[0];
  int connect_size = bdy_verts_size;
  iMesh_createTag(meshIface, "CMEL_index", 1, iBase_INTEGER,
                  &index_tag, &result, strlen("CMEL_index"));
  if (iBase_SUCCESS != result && iBase_TAG_ALREADY_EXISTS == result) return false;
  
  iMesh_setIntArrData(meshIface, &bdy_verts[0], bdy_verts.size(), index_tag, &indices[0], bdy_verts.size(), &result);
  if (iBase_SUCCESS != result) return false;
  iMesh_getIntArrData(meshIface, bdy_verts_vec, bdy_verts_size, index_tag, 
                      &connect_ptr, &connect_size, &connect_size, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't set indices of bounding vertices or get connect indices." << std::endl;
    return false;
  }

    // reverse the indices for reverse-sense bdy entities
  for (unsigned int i = 0; i < bdy_mesh.size(); i++) {
    if (bdy_orientations[i] < 0) 
      std::reverse(&connect[offsets[i]], &connect[offsets[i+1]]);
  }
  
    // get vertex coordinates
  int coords_size = 3*bdy_verts.size();
  coords.resize(coords_size);
  double *coords_ptr = &coords[0];
  iBase_StorageOrder this_order = iBase_INTERLEAVED;
  iMesh_getVtxArrCoords(meshIface, &bdy_verts[0], bdy_verts.size(), 
                        &this_order, &coords_ptr, &coords_size, &coords_size,
                        &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't get coords of bounding vertices." << std::endl;
    return false;
  }
  
    // finally, destroy the indices tag
  iMesh_destroyTag(meshIface, index_tag, true, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't destroy index tag." << std::endl;
    return false;
  }
  
  return true;
}

bool CMEL::get_bounding_curve_mesh(iBase_EntityHandle gentity, double mesh_size,
                                   int mesh_intervals, const bool force_intervals,
                                   std::vector<iBase_EntityHandle> &bdy_verts,
                                   std::vector<double> &bdy_coords,
                                   std::vector<int> &loop_sizes, std::vector<int> &loops) 
{
  std::vector<int> connect;
  
  bool success = get_bounding_mesh_connect(gentity, mesh_size, mesh_intervals,
                                           force_intervals,
                                           bdy_verts, bdy_coords, connect);
  if (!success) return success;
  
  bool done = false;
  std::vector<int>::iterator cit = connect.begin();
  std::vector<iBase_EntityHandle>::iterator eit = bdy_verts.begin();
  int last_size = 0;
  while (!done && cit != connect.end()) {
    int start_v = *cit;
    loops.push_back(start_v);
    int next_v;
    do {
      next_v = *(cit+1);
      if (next_v != start_v) loops.push_back(next_v);
      cit += 2;
    }
    while (next_v != start_v && cit != connect.end());
    
      // have a complete loop; compute the size
    loop_sizes.push_back(loops.size()-last_size);
    last_size = loops.size();
  }
  
    // sanity check on # loops
  if (loops.size() != bdy_verts.size()) {
    std::cerr << "Wrong # edges bounding this surface." << std::endl;
    return false;
  }

  return true;
}

*/

bool CMEL::mesh_vertex(iBase_EntityHandle gentity, 
                       std::vector<iBase_EntityHandle> &new_entities) 
{
    // just get coordinates of the vertex, create a mesh, assign, return
  std::vector<double> coords(3);
  
  int result;
  iGeom_getVtxCoord(geomIface, gentity, &coords[0], &coords[1], &coords[2], &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting coords of gentity." << std::endl;
    return false;
  }

  new_entities.clear();
  std::vector<int> connect;
  return create_vertices_elements(gentity, new_entities, coords, connect,
                                  iMesh_ALL_TOPOLOGIES,
                                  new_entities);
}

bool CMEL::mesh_curve(iBase_EntityHandle gentity, 
                      double mesh_size,
                      int mesh_intervals,
                      std::vector<iBase_EntityHandle> &new_entities) 
{
  std::vector<iBase_EntityHandle> bdy_mesh;

  bool success = bdy_nodes(gentity, bdy_mesh);
  if (!success) 
    return success;
  
    // get the parametric range of the edge
  double uv_min[2], uv_max[2];
  int result;
  iGeom_getEntURange(geomIface, gentity, uv_min, uv_max, &result);
  if (iBase_SUCCESS != result) return false;

  double length;
  double *lengths = &length;
  int lengths_size = 1, lengths_alloc = 1;
  iGeom_measure(geomIface, &gentity, 1, &lengths, &lengths_alloc, &lengths_size, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting uv range for curve." << std::endl;
    return false;
  }

  int nint;
  if (0 < mesh_size) nint = (int) ceil(length/mesh_size);
  else if (0 < mesh_intervals) nint = mesh_intervals;
  else return false;

  if (0 == nint) nint = 1;

    // compute parametric interval length and parametric positions
  double dp = (uv_max[0] - uv_min[0]) / (double) nint;
  std::vector<double> vert_uvs(2*(nint-1));
  for (int i = 1; i < nint; i++)
    vert_uvs[2*(i-1)] = uv_min[0] + i*dp;

    // get the coordinates for the new vertices
  std::vector<double> vert_xyzs(3*(nint-1));
  for (int i = 0; i < nint-1; i++) {
    int tmp_result; 
    iGeom_getEntUtoXYZ(geomIface, gentity, vert_uvs[2*i], 
                        &vert_xyzs[3*i], &vert_xyzs[3*i+1], &vert_xyzs[3*i+2],
                        &tmp_result);
    if (iBase_SUCCESS != tmp_result) result = tmp_result;
  }
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting xyzs for curve vertices." << std::endl;
    return false;
  }

    // pack an integer array with proper vertex indices
  std::vector<int> edge_verts(2*nint);
  
  edge_verts[0] = 0;
  edge_verts[1] = (bdy_mesh.size() == 2 ? 2 : 1);
  for (int i = 1; i < nint; i++) {
    edge_verts[2*i] = edge_verts[2*i-1];
    edge_verts[2*i+1] = edge_verts[2*i-1]+1;
  }
  edge_verts[2*nint-1] = (bdy_mesh.size() == 2 ? 1 : 0);
  
    // pass the whole thing in to a function to create the mesh
  success = create_vertices_elements(gentity, bdy_mesh, vert_xyzs, edge_verts,
                                     iMesh_LINE_SEGMENT,
                                     new_entities);

  if (debug) {
    iBase_EntityHandle *adj_gents = NULL;
    int adj_gents_size = 0, adj_gents_alloc = 0;
    iGeom_getEntAdj(geomIface, gentity, iBase_VERTEX,
                  &adj_gents, &adj_gents_alloc, &adj_gents_size, &result);
  
    std::cout << "Meshed curve " << get_gentity_id(gentity) << ", bounding vertices "
              << get_gentity_id(adj_gents[0]) << " and " 
              << get_gentity_id(adj_gents[(adj_gents_size == 2 ? 1 : 0)])
              << " with " << new_entities.size()
              << " edges." << std::endl
              << "New vertex coords: " << std::endl;
    for (std::vector<double>::iterator vit = vert_xyzs.begin(); vit != vert_xyzs.end(); vit+=3)
      std::cout << *vit << ", " << *(vit+1) << ", " << *(vit+2) << std::endl;
    std::cout << std::endl;
    iBase_EntityHandle* adjacentEntityHandles = NULL;
    int adj_entity_handles_size = 0;
    int  adjacentEntityHandles_allocated = 0;
    int* offset = NULL;
    int offset_allocated = 0, offset_size=0;

    iMesh_getEntArrAdj(meshIface, &new_entities[0], new_entities.size(),
          iMesh_POINT, &adjacentEntityHandles,
                              /*inout*/ &adjacentEntityHandles_allocated,
                              /*out*/ &adj_entity_handles_size,
                              /*inout*/ &offset,
                              /*inout*/ &offset_allocated,
                              /*out*/ &offset_size,
                              /*out*/ &result);
    for (unsigned int i=0; i<new_entities.size(); i++)
       std::cout<<" edge:" << std::hex<< new_entities[i] <<" v:" << std::dec<<
          adjacentEntityHandles[2*i] << " " <<  adjacentEntityHandles[2*i+1] << std::endl;

  }
  
  return success;
}

bool CMEL::assign_mesh(iBase_EntityHandle gentity, std::vector<iBase_EntityHandle> &mesh) 
{
  int result;
  iRel_setEntEntArrRelation(relateIface, pairHandle, gentity, false,
                               &mesh[0], mesh.size(), &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't set classification." << std::endl;
    return false;
  }

  return true;
}

bool CMEL::create_vertices_elements(iBase_EntityHandle gentity, 
                                    std::vector<iBase_EntityHandle> &bdy_verts, 
                                    std::vector<double> &coords, std::vector<int> &connect, 
                                    int ent_type,
                                    std::vector<iBase_EntityHandle> &new_entities) 
{
  return create_vertices_elements(gentity, bdy_verts, &coords[0], coords.size()/3,
                                  connect, ent_type, new_entities);
}

bool CMEL::get_mesh_set(iBase_EntityHandle gentity, iBase_EntitySetHandle &mesh_set,
                        bool create_if_missing) 
{
  int result;
  iRel_getEntSetRelation(relateIface,
                            pairHandle,
                            gentity, false, 
                            &mesh_set, &result);
  if (iBase_SUCCESS != result && create_if_missing) {
      // get the dimension of the geometry entity
    int this_type;
    iGeom_getEntType(geomIface, gentity, &this_type, &result);
    if (iBase_SUCCESS != result) return false;
    assert(iBase_VERTEX <= this_type && iBase_REGION >= this_type);
    
      // create the set
    bool isList = false;
    if (this_type == iBase_EDGE)
       isList = true;
    iMesh_createEntSet(meshIface, isList, &mesh_set, &result);
    if (iBase_SUCCESS != result) return false;

      // set tags for global id and telling that it's a geom topo set
    int this_id = get_gentity_id(gentity);
    iBase_TagHandle gid_tag;
    iMesh_getTagHandle(meshIface, "GLOBAL_ID", &gid_tag, &result, 10);
    assert (0 != gid_tag);
    iMesh_setEntSetIntData(meshIface, mesh_set, gid_tag, this_id, &result);
    if (iBase_SUCCESS != result) return false;
    
    iBase_TagHandle geom_tag;
    iMesh_getTagHandle(meshIface, "GEOM_DIMENSION", &geom_tag, &result, 15);
    if (0 == geom_tag) {
      iMesh_createTag(meshIface, "GEOM_DIMENSION", 1, iBase_INTEGER,
                      &geom_tag, &result, 15);
      if (iBase_SUCCESS != result) return false;
    }
    assert (0 != geom_tag);
    iMesh_setEntSetIntData(meshIface, mesh_set, geom_tag, 
                     (int)(this_type - iBase_VERTEX), &result);
    if (iBase_SUCCESS != result) return false;

      // now associate it with gentity
    iRel_setEntSetRelation(relateIface,
                               pairHandle,
                               gentity, 
                               mesh_set,
                               &result);
    if (iBase_SUCCESS != result) return false;
  }
  
  return (result == iBase_SUCCESS);
}

void CMEL::print_meshed_entity( iBase_EntityHandle gentity,
                                int elem_type,
                                int num_node,
                                int num_conn )
{
  static const char* const typenames[] = 
    { "vertex", "curve", "surface", "volume", "<invalid>" };
  int err, type, id = get_gentity_id( gentity );
  iGeom_getEntType( geomIface, gentity, &type, &err );
  if (iBase_SUCCESS != err || type < 0 || type > 3)
    type = 4;
    
  int vert_per_elem;
  const char* elem_type_name = 0;
  switch (elem_type) {
    case iMesh_LINE_SEGMENT:  vert_per_elem = 2; elem_type_name = "edge"; break;
    case iMesh_TRIANGLE:      vert_per_elem = 3; elem_type_name = "tri";  break;
    case iMesh_QUADRILATERAL: vert_per_elem = 4; elem_type_name = "quad"; break;
    case iMesh_TETRAHEDRON:   vert_per_elem = 4; elem_type_name = "tet";  break;
    case iMesh_HEXAHEDRON:    vert_per_elem = 8; elem_type_name = "hex";  break;
    default:                  vert_per_elem = 0; elem_type_name = 0;      break;
  }
  
  std::cout << "Meshed " << typenames[type] << " " << id << " with " 
            << num_node << " node(s)";
  if (num_conn) {
    if (elem_type_name) 
      std::cout << " and " << num_conn/vert_per_elem << " " << elem_type_name << "(s)";
    else
      std::cout << " and elements of an unknown type (type = " << elem_type
                << ", num_conn = " << num_conn << ")";
  }
  std::cout << std::endl;
}

bool CMEL::create_vertices_elements(iBase_EntityHandle gentity, 
                                    std::vector<iBase_EntityHandle> &bdy_verts, 
                                    double *coords, const int num_points,
                                    std::vector<int> &connect, 
                                    int ent_type,
                                    std::vector<iBase_EntityHandle> &new_entities) 
{
  iBase_EntitySetHandle mesh_set;
  bool success = get_mesh_set(gentity, mesh_set, true);
  if (!success) {
    std::cerr << "Couldn't get mesh set for a gentity." << std::endl;
    return false;
  }

  if (debug)
    print_meshed_entity( gentity, ent_type, num_points, connect.size() );
  
    // create new vertices
  std::vector<iBase_EntityHandle> new_verts(num_points);
  iBase_EntityHandle *new_verts_ptr = &new_verts[0];
  int new_verts_size = num_points, new_verts_alloc = num_points;
  int result;
  iMesh_createVtxArr(meshIface, num_points, iBase_INTERLEAVED,
                     coords, 3*num_points, 
                     &new_verts_ptr, &new_verts_alloc, &new_verts_size,
                     &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't create new vertices." << std::endl;
    return false;
  }
  
    // add vertices to the set
  iMesh_addEntArrToSet(meshIface, new_verts_ptr, num_points,
                       mesh_set, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't assign mesh from create_vertices_elements." << std::endl;
    return false;
  }

  if (connect.empty()) {
      // if we're just doing vertices, return new verts in new_entities
    std::copy(new_verts.begin(), new_verts.end(), std::back_inserter(new_entities));
    return true;
  }
  
    // assemble connectivity of new elements
  std::vector<iBase_EntityHandle> new_connect(connect.size());
  int old_verts = bdy_verts.size();
  for (unsigned int i = 0; i < connect.size(); i++) {
    if (connect[i] >= old_verts)
      new_connect[i] = new_verts[connect[i]-old_verts];
    else
      new_connect[i] = bdy_verts[connect[i]];
  }
  
    // create them
  int *status = NULL;
  int status_size = 0, status_alloc = 0;
  iBase_EntityHandle *new_ments = NULL;
  int new_ments_size = 0, new_ments_alloc = 0;
  iMesh_createEntArr(meshIface, ent_type, &new_connect[0], connect.size(),
                     &new_ments, &new_ments_alloc, &new_ments_size,
                     &status, &status_alloc, &status_size,
                     &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't create new elements." << std::endl;
    return false;
  }
  
    // put them into the new_entities vector
  new_entities.resize(new_ments_size);
  std::copy(new_ments, new_ments+new_ments_size, &new_entities[0]);

    // add elements to the geometric owner's set
  iMesh_addEntArrToSet(meshIface, new_ments, new_ments_size,
                       mesh_set, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Couldn't assign mesh from create_vertices_elements." 
              << std::endl;
    return false;
  }
  
  return true;
}

int CMEL::get_gentity_id(iBase_EntityHandle gentity) 
{
  iBase_TagHandle gid_tag;
  int result, err;
  iGeom_getTagHandle(geomIface, "GLOBAL_ID", &gid_tag, &result, strlen("GLOBAL_ID"));
  if (iBase_SUCCESS != result) return -1;
  
  iGeom_getIntData(geomIface, gentity, gid_tag, &result, &err);
  return (iBase_SUCCESS == err) ? result : -1;
}

bool CMEL::bdy_geom_grouped(iBase_EntityHandle gentity, 
                            std::vector<iBase_EntityHandle> &groups,
                            std::vector<int> &group_sizes) 
{
    // for a given surface, return the bounding edges in the form of edge groups,
    // oriented ccw around surface
  int result, ent_type;
  
  iGeom_getEntType(geomIface, gentity, &ent_type, &result);
  if (iBase_SUCCESS != result) return false;

    // shouldn't be calling this function if we're a vertex
  if (ent_type == iBase_VERTEX) return false;
  
    // get adj ents & senses
  iBase_EntityHandle *adj_gents = NULL;
  int adj_gents_size = 0, adj_gents_alloc = 0;
  iGeom_getEntAdj(geomIface, gentity, ent_type-1,
                  &adj_gents, &adj_gents_alloc, &adj_gents_size, &result);
  if (iBase_SUCCESS != result) return false;

    // special handling for edges and entities with single bounding entity - 
    // will just have 1 or 2 adj entities, set those directly and return
  if (ent_type == iBase_EDGE || adj_gents_size == 1) {
    for (int i = 0; i < adj_gents_size; i++) {
      groups.push_back(adj_gents[i]);
      group_sizes.push_back(1);
    }
    return true;
  }
  
  std::vector<iBase_EntityHandle> intersected_ents(adj_gents_size);

    // get adjacent entities into a sorted, mutable list
  std::vector<iBase_EntityHandle> b_ents;
  std::set<iBase_EntityHandle> dbl_curves, sgl_curves;
  std::copy(adj_gents, adj_gents+adj_gents_size, std::back_inserter(b_ents));
  std::sort(b_ents.begin(), b_ents.end());

  while (!b_ents.empty()) {
      // get 1st in group
    
    std::vector<iBase_EntityHandle> group_stack, this_group;
    group_stack.push_back(b_ents.front());
    
      // while there are still entities on the stack
    while (!group_stack.empty()) {
      
        // pop one off & get its sense
      iBase_EntityHandle this_entity = group_stack.back(); group_stack.pop_back();
      
      int this_sense;
      if (ent_type == iBase_FACE)
        iGeom_getEgFcSense(geomIface, this_entity, gentity, &this_sense, &result);
      else if (ent_type == iBase_REGION)
        iGeom_getEntNrmlSense(geomIface, this_entity, gentity, &this_sense, &result);
      else return false;

        // if we already have this one, continue
      if ((0 != this_sense && 
           std::find(this_group.begin(), this_group.end(), this_entity) != this_group.end()) ||
          std::find(dbl_curves.begin(), dbl_curves.end(), this_entity) != dbl_curves.end())
        continue;

        // either way we need the d-2 entities
      iBase_EntityHandle *bridges = NULL;
      int bridges_size = 0, bridges_alloc = 0;
      iGeom_getEntAdj(geomIface, this_entity, ent_type-2, 
                      &bridges, &bridges_alloc, &bridges_size, &result);
      if (iBase_SUCCESS != result) return false;

        // only remove from the list of candidates if it's not dual-sensed
      if (0 != this_sense) 
        b_ents.erase(std::remove(b_ents.begin(), b_ents.end(), this_entity), b_ents.end());

        // if it's double-sensed and this is the first time we're seeing it, find the right sense
      if (0 == this_sense) {
        if (std::find(sgl_curves.begin(), sgl_curves.end(), this_entity) == sgl_curves.end()) {
          sgl_curves.insert(this_entity);
          if (!this_group.empty()) {
            iBase_EntityHandle common_v = shared_entity(this_entity, this_group.back(), iBase_VERTEX);
            if (common_v == bridges[0]) this_sense = 1;
            else if (common_v == bridges[1]) this_sense = -1;
            else return false;
          }
            // else, if this is the first one in the loop, just choose a sense
          else {
            this_sense = 1;
          }
        }
        else {
            // else this is the second time we're seeing it, move it to the dbl_curves list and remove
            // from the candidates list
          dbl_curves.insert(this_entity);
          sgl_curves.erase(this_entity);
          b_ents.erase(std::remove(b_ents.begin(), b_ents.end(), this_entity), b_ents.end());
        }
      }
          
        // it's in the group; put on group & remove from untreated ones
      this_group.push_back(this_entity);
      std::vector<iBase_EntityHandle> tmp_from, tmp_adjs;

        // if we're on a face and we're the first in a group, check sense of this first
        // edge; make sure "next" in loop sense is last on list
      if (ent_type == iBase_FACE && this_group.size() == 1) {

          // get vertex which we know is shared by the "right" next edge; first get the vertices
          // if sense of current edge is forward, it's the 2nd vertex we want,
          // otherwise the first
        if (1 == this_sense && bridges_size > 1) tmp_from.push_back(bridges[1]);
        else tmp_from.push_back(bridges[0]);
        tmp_adjs = b_ents;
        get_adjs_bool(tmp_from, ent_type-1, tmp_adjs, CMEL::INTERSECT);
      }
      else {
        std::copy(bridges, bridges+bridges_size, std::back_inserter(tmp_from));
        std::vector<iBase_EntityHandle> tmp_adjs2;
        get_adjs_bool(tmp_from, ent_type-1, tmp_adjs2, CMEL::UNION);
        std::sort(tmp_adjs2.begin(), tmp_adjs2.end());
        tmp_adjs.resize(tmp_adjs2.size());
        tmp_adjs.erase(std::set_intersection(b_ents.begin(), b_ents.end(),
                                             tmp_adjs2.begin(), tmp_adjs2.end(),
                                             tmp_adjs.begin()), tmp_adjs.end());
      }

      if (ent_type == iBase_FACE && tmp_adjs.size() > 1) {
          // more than one adjacent edge - need to evaluate winding angle to find right one
        iBase_EntityHandle next_ent = next_winding(this_entity, gentity, this_sense, tmp_adjs);
        if (0 == next_ent) return false;
        group_stack.push_back(next_ent);
      }
      else if (!tmp_adjs.empty()) 
        std::copy(tmp_adjs.begin(), tmp_adjs.end(), std::back_inserter(group_stack));
      
      free(bridges);
    }
    
    // put group in group list
    std::copy(this_group.begin(), this_group.end(), std::back_inserter(groups));
    group_sizes.push_back(this_group.size());
  }
  
  return true;
}

iBase_EntityHandle CMEL::shared_entity(iBase_EntityHandle ent1,
                                 iBase_EntityHandle ent2,
                                 int to_type) 
{
    // find the shared entity between the two entities of the prescribed dimension
  std::vector<iBase_EntityHandle> from_ents(2), to_ents;
  from_ents[0] = ent1;
  from_ents[1] = ent2;
  bool success = get_adjs_bool(from_ents, to_type, to_ents, CMEL::INTERSECT);
  if (!success || to_ents.empty()) return 0;
  else return to_ents[0];
}

bool CMEL::get_adjs_bool(std::vector<iBase_EntityHandle> &from_ents,
                         int to_type,
                         std::vector<iBase_EntityHandle> &to_ents,
                         CMEL::BooleanType op_type) 
{
  if (from_ents.empty()) {
    to_ents.clear();
    return true;
  }

  int result;
  iBase_EntityHandle *bridges = NULL;
  int bridges_size = 0, bridges_alloc = 0;

  std::vector<iBase_EntityHandle>::iterator from_it = from_ents.begin();
  if (to_ents.empty() && op_type == CMEL::INTERSECT) {
    iGeom_getEntAdj(geomIface, from_ents.front(), to_type,
                    &bridges, &bridges_alloc, &bridges_size, &result);
    if (iBase_SUCCESS != result) return false;
    std::copy(bridges, bridges+bridges_size, std::back_inserter(to_ents));
    from_it++;
    free(bridges);
    bridges = NULL;
    bridges_alloc = 0;
  }

  std::sort(to_ents.begin(), to_ents.end());
  std::vector<iBase_EntityHandle> result_ents(to_ents.size());
  for (; from_it != from_ents.end(); from_it++) {
    iGeom_getEntAdj(geomIface, *from_it, to_type,
                    &bridges, &bridges_alloc, &bridges_size, &result);
    if (iBase_SUCCESS != result) return false;
    if (op_type == CMEL::INTERSECT) {
      std::sort(bridges, bridges+bridges_size);
      result_ents.erase(std::set_intersection(to_ents.begin(), to_ents.end(),
                                              bridges, bridges+bridges_size,
                                              result_ents.begin()), result_ents.end());
    }
    else {
      std::copy(bridges, bridges+bridges_size, std::back_inserter(result_ents));
    }
    
    to_ents = result_ents;
    free(bridges);
    bridges = NULL;
    bridges_alloc = 0;
  }
  
  if (op_type == CMEL::UNION)
    std::sort(to_ents.begin(), to_ents.end());
  
  to_ents.erase(std::unique(to_ents.begin(), to_ents.end()), to_ents.end());
  
  return true;
}

iBase_EntityHandle CMEL::next_winding(iBase_EntityHandle this_edge, 
                                iBase_EntityHandle gface, 
                                int this_sense, 
                                std::vector<iBase_EntityHandle> &tmp_adjs) 
{
    // given this_entity, a d-1 entity bounding gentity, and optional "next"
    // entities in tmp_adjs, find the next one based on windings around the shared
    // vertex
  
    // first, get the shared vertex 
  iBase_EntityHandle *verts = NULL;
  int result, verts_size = 0, verts_alloc = 0;
  iGeom_getEntAdj(geomIface, this_edge, iBase_VERTEX,
                  &verts, &verts_alloc, &verts_size, &result);
  if (iBase_SUCCESS != result) return 0;

  iBase_EntityHandle shared_vert = verts[0];
  if (this_sense == 1 && verts_size > 1) shared_vert = verts[1];

    // get locations just before the vertex, at the vertex, and just after the vertex
  double v1[3], v2[3], v3[3];
  double umin, umax;
  iGeom_getEntURange(geomIface, this_edge, &umin, &umax, &result);
  if (iBase_SUCCESS != result) return 0;
  double utgt;
  if (1 == this_sense) utgt = umin + 0.9 * (umax - umin);
  else utgt = umin + 0.1 * (umax - umin);
  iGeom_getEntUtoXYZ(geomIface, this_edge, utgt,
                     v1, v1+1, v1+2, &result);
  if (iBase_SUCCESS != result) return 0;
  
  iGeom_getVtxCoord(geomIface, shared_vert, v2, v2+1, v2+2, &result);
  if (iBase_SUCCESS != result) return 0;

  double v21[] = {v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]};
/*#define DOT(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define LENGTH_SQ(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define NORMALIZE(a) {double lsq = LENGTH_SQ(a); \
                      if (lsq == 0.0) return 0; \
                      lsq = 1.0 / sqrt(lsq); \
                      for (int i = 0; i < 3; i++) a[i] *= lsq;}
#define CROSS(a, b, c) {c[0] = a[1]*b[2] - a[2]*b[1]; \
                        c[1] = a[2]*b[0] - a[0]*b[2]; \
                        c[2] = a[0]*b[1] - a[1]*b[0];}
  const double PI = acos(-1.0);
  const double TWO_PI = 2.0 * PI*/;
  
  NORMALIZE(v21);

    // get the normal vector at the vertex
  double normal[3];
  iGeom_getEntNrmlXYZ(geomIface, gface, v2[0], v2[1], v2[2],
                      normal, normal+1, normal+2, &result);
  if (iBase_SUCCESS != result) return 0;
  
    // now loop over candidates, finding magnitude of swept angle
  iBase_EntityHandle other = 0;
  double angle = TWO_PI;
  int is_adj;
  for (std::vector<iBase_EntityHandle>::iterator vit = tmp_adjs.begin(); vit != tmp_adjs.end(); vit++) {
      // if we're here, we have multiple candidates, therefore don't choose the same one
    if (*vit == this_edge)
      continue;
    iGeom_isEntAdj(geomIface, *vit, shared_vert, &is_adj, &result);
    if (iBase_SUCCESS != result)
      return 0;
    if (!is_adj)
      continue;

      // get param range
    iGeom_getEntURange(geomIface, *vit, &umin, &umax, &result);
    if (iBase_SUCCESS != result) return 0;
      // get sense
    int tmp_sense;
    iGeom_getEgFcSense(geomIface, *vit, gface, &tmp_sense, &result);
    if (iBase_SUCCESS != result)
      return 0;
    if (1 == tmp_sense) utgt = umin + 0.1 * (umax - umin);
    else utgt = umin + 0.9 * (umax - umin);
    iGeom_getEntUtoXYZ(geomIface, *vit, utgt,
                       v3, v3+1, v3+2, &result);
    if (iBase_SUCCESS != result) return 0;
    double v23[] = {v3[0]-v2[0], v3[1]-v2[1], v3[2]-v2[2]};
    NORMALIZE(v23);

    CROSS(normal, v23, v1);
    
    CROSS(v1, normal, v3);

    double x = DOT(v21, v3);
    double y = DOT(v21, v1);

    assert(x != 0.0 || y != 0.0);
    double this_angle = atan2(y, x);
    if (this_angle < 0.0)
    {
      this_angle += TWO_PI;
    }
    if (this_angle < angle) {
      other = *vit;
      angle = this_angle;
    }
  }
  
  return other;
}
  
// this will return the nodes that are at the ends of the curve
// one node means periodic edge
bool CMEL::bdy_nodes(iBase_EntityHandle gentity,
                               std::vector<iBase_EntityHandle> &end_nodes)
{

  std::vector<int> group_sizes;
 
  // get bounding geometry groups; each size is 1
  // these are vertex sets
  std::vector<iBase_EntityHandle> gent_groups;
  std::vector<int> gent_group_sizes;
  bool success = bdy_geom_grouped(gentity, gent_groups, gent_group_sizes);
  if (!success) return success;
  // if we have only one group, it is fine, just get the mesh set , and a vertex
  std::vector<iBase_EntitySetHandle> mesh_sets;
  for (unsigned int i=0; i<gent_groups.size(); i++)
  {
     iBase_EntitySetHandle mesh_set;
     get_mesh_set(gent_groups[i], mesh_set, false);// do not create if missing
     mesh_sets.push_back(mesh_set);
  }
  // we should have only one node at each vertex, so in each set
  // get the nodes from the first, and possibly, unique mesh set
  iBase_EntityHandle *ments = NULL;
  int result;
  int ments_size = 0, ments_alloc = 0;

    // get just the 0-dimensional entities from the first mesh set
  iMesh_getEntities(meshIface, mesh_sets[0],
                    0,
                    iMesh_ALL_TOPOLOGIES,
                    &ments, &ments_alloc, &ments_size,
                    &result);
  if (iBase_SUCCESS != result || ments_size!=1) return false;
  if (gent_groups.size()==1)
  {
     end_nodes.push_back(ments[0]); // really, just one node in the boundary
     free(ments);
     return true;
  }
  // else, we have to see which one is the first vertex, and which one the second
  assert(gent_groups.size()==2);
  //
  // get the second set of nodes from vertices
  iBase_EntityHandle *ments2 = NULL;
  int ments_size2 = 0, ments_alloc2 = 0;

      // get just the 0-dimensional entities from the first mesh set
  iMesh_getEntities(meshIface, mesh_sets[1],
                      0,
                      iMesh_ALL_TOPOLOGIES,
                      &ments2, &ments_alloc2, &ments_size2,
                      &result);
  if (iBase_SUCCESS != result || ments_size2!=1) return false;

  int sense_out=1;
  iGeom_getEgVtxSense( geomIface, gentity, gent_groups[0],gent_groups[1],
                              &sense_out, &result);
  if (sense_out ==-1)
  {
     // start with second set
     end_nodes.push_back(ments2[0]);
     end_nodes.push_back(ments[0]); // exactly 2 nodes
  }
  else // assert(sense_out == 1);
  {
     end_nodes.push_back(ments[0]);
     end_nodes.push_back(ments2[0]); // exactly 2 nodes
  }

  //
  free(ments);
  free(ments2);

  return true;

}

bool CMEL::bdy_elements_senses(iBase_EntityHandle gentity,
                               std::vector<iBase_EntityHandle> &elements,
                              std::vector<int> &element_senses)
{
   std::vector<int> group_sizes;
   return bdy_elements_senses_grouped(gentity, elements, element_senses, group_sizes);
}

bool CMEL::bdy_elements_senses_grouped(iBase_EntityHandle gentity,
                                       std::vector<iBase_EntityHandle> &elements,
                                       std::vector<int> &element_senses,
                                       std::vector<int> &group_sizes) 
{
    // get bounding geometry groups
  std::vector<iBase_EntityHandle> gent_groups;
  std::vector<int> gent_group_sizes;
  bool success = bdy_geom_grouped(gentity, gent_groups, gent_group_sizes);
  if (!success) return success;

    // get dimension of gentity
  int result, gentity_type;
  iGeom_getEntType(geomIface, gentity, &gentity_type, &result);
  if (iBase_SUCCESS != result) return false;
  int gentity_dim = gentity_type - iBase_VERTEX;

    // now get elements & senses for each of the groups
  std::vector<iBase_EntityHandle>::iterator gent_it = gent_groups.begin();
  
    // loop over groups:
  for (unsigned int grp = 0; grp < gent_group_sizes.size(); grp++) {

    int grp_size = elements.size();
    
      // loop over entities in this group:
    for (int i = 0; i < gent_group_sizes[grp]; i++) {
      iBase_EntityHandle this_gent = *gent_it++;

      iBase_EntityHandle *ments = NULL;
      int result;
      int ments_size = 0, ments_alloc = 0;
      iBase_EntitySetHandle this_set = 0;
      iRel_getEntSetRelation(relateIface, pairHandle, this_gent,
                                false, &this_set, &result);
      if (iBase_SUCCESS != result) return false;

        // get just the (d-1)-dimensional entities from that set
      iMesh_getEntities(meshIface, this_set, 
                        (gentity_dim-1),
                        iMesh_ALL_TOPOLOGIES,
                        &ments, &ments_alloc, &ments_size,
                        &result);
      if (iBase_SUCCESS != result) return false;
      
        // get sense of this gent wrt gentity
      int this_sense = 1;
      result = iBase_SUCCESS;
      if (gentity_dim == 3) 
        iGeom_getEntNrmlSense(geomIface, this_gent, gentity, &this_sense, &result);
      else if (gentity_dim == 2) 
        iGeom_getEgFcSense(geomIface, this_gent, gentity, &this_sense, &result);
      if (iBase_SUCCESS != result)
        return false;

      if (0 == this_sense) get_both_sense(this_gent, 
                                          (grp == 0 ? *gent_it : *(gent_it-2)),
                                          grp == 0, this_sense);
      
        // if necessary, reverse the list itself
      if (-1 == this_sense && gentity_dim == 2)
        std::reverse(ments, ments+ments_size);
        
        // append the elements to the group
      std::copy(ments, ments+ments_size, std::back_inserter(elements));

        // fill in senses for these elements
      element_senses.resize(elements.size(), this_sense);
    }

    grp_size = elements.size() - grp_size;
    group_sizes.push_back(grp_size);
  }

  return true;
}

bool CMEL::get_both_sense(iBase_EntityHandle this_gent, 
                          iBase_EntityHandle other_gent,
                          bool next_ent, int &this_sense) 
{
    // get the proper sense at this point for this_gent, assuming
    // other_gent is last (next_ent = false) or next (next_ent = true)
    // entity along the loop

    // get the adjacent entities to each
  iBase_EntityHandle *this_adjs = NULL, *other_adjs = NULL;
  
  int result, this_alloc = 0, other_alloc = 0, this_size, other_size;
  
  iGeom_getEntAdj(geomIface, this_gent, 0,
                  &this_adjs, &this_alloc, &this_size,
                  &result);
  if (iBase_SUCCESS != result) return false;
  
  iGeom_getEntAdj(geomIface, other_gent, 0,
                  &other_adjs, &other_alloc, &other_size, &result);
  if (iBase_SUCCESS != result) return false;
  
    // should never have single-vertex curve with 'both'-type sense
  if (this_size == 1) return false;
  
    // if other curve is at start of this, assume forward
  if (this_adjs[0] == other_adjs[0] || 
      (other_size > 1 && this_adjs[0] == other_adjs[1]))
    this_sense = 1;
  else this_sense = -1;
  
    // now switch depending on whether we used preceding or next
    // curve for other
  if (next_ent) this_sense *= -1;
  
  return true;
}

bool CMEL::assign_tmp_indices(std::vector<iBase_EntityHandle> &bdy_verts,
                              iBase_TagHandle &index_tag, int offset) 
{
  std::vector<int> indices(bdy_verts.size());
  int idx = offset;
  for (unsigned int i = 0; i < bdy_verts.size(); i++) 
    indices[i] = idx++;
  if (debug)
  {
     std::cout<<"assign tmp indices:" << std::endl;
     for (unsigned int i =0; i<bdy_verts.size(); i++)
        std::cout<< "bvert: "<< (unsigned long int)bdy_verts[i]<< " indx:" << indices[i]<< std::endl ;
  }

  int result;
  iMesh_createTag(meshIface, "CMEL_index", 1, iBase_INTEGER,
                  &index_tag, &result, 11);
  if (iBase_SUCCESS != result && iBase_TAG_ALREADY_EXISTS != result) return false;
  
  iMesh_setIntArrData(meshIface, &bdy_verts[0], bdy_verts.size(), 
                      index_tag, &indices[0], bdy_verts.size(), &result);
  if (iBase_SUCCESS != result) return false;

  return true;
}
bool CMEL::trimSurface(const char * polygon_filename, int len)
{
   
   // read the file with the polygon user data
   std::ifstream datafile(polygon_filename, std::ifstream::in);
   if (!datafile) {
      std::cout << "can't read file\n";
      return 1;
   }
   //
   char temp[100];
   double direction[3];// normalized
   double gridSize;
   datafile.getline(temp, 100);// first line 

   // get direction and mesh size along polygon segments, from file
   sscanf(temp, " %lf %lf %lf %lf ", direction, direction+1, direction+2, &gridSize);
   NORMALIZE(direction);// just to be sure

   std::vector<double> xs, ys, zs;
   while (!datafile.eof()) {
      datafile.getline(temp, 100);
      //int id = 0;
      double x, y, z;
      int nr = sscanf(temp, "%lf %lf %lf", &x, &y, &z);
      if (nr==3)
      {
         xs.push_back(x);
         ys.push_back(y);
         zs.push_back(z);
      }
   }
   int sizePolygon = xs.size() ;
   if (sizePolygon < 3)
   {
      std::cerr<< " Not enough points in the polygon" << std::endl;
      return false;
   }
   // now advance in the direction of the polygon edge, with gridSize, and compute points
   // on the surface (it should be first surface); those points will form a new boundary on the
   // (first) surface, which will be used for camal;
   //int startPolygon = 0;
   double currentPosition[3] = {xs[0], ys[0], zs[0]};
   double currentDirection[3] = {xs[1] - xs[0], ys[1] - ys[0], zs[1] - zs[0]};
   int nextPolygonPt = 1;
   double nextPolyCoord[3] = {xs[1], ys[1], zs[1]};
   double currentLengthLeft = sqrt(LENGTH_SQ(currentDirection));
   NORMALIZE(currentDirection);

   int reachedEnd = 0;
   int err;
   while (!reachedEnd)
   {
      // march along given direction  ; it should be at most size / 3?

      iBase_EntityHandle * intersect_entity_handles = NULL;
      int intersect_entity_handles_allocated = 0, intersect_entity_handles_size = 0;
      double * intersect_coords = NULL;
      int intersect_coords_allocated =0 ,  intersect_coords_size = 0;
       double * param_coords = NULL;
      int param_coords_allocated = 0, param_coords_size =0;
      iGeom_getPntRayIntsct( geomIface,
            currentPosition[0], currentPosition[1], currentPosition[2],
            direction[0], direction[1], direction[2],
            &intersect_entity_handles, &intersect_entity_handles_allocated,
            &intersect_entity_handles_size, iBase_INTERLEAVED,
            &intersect_coords, &intersect_coords_allocated, &intersect_coords_size,
            &param_coords, &param_coords_allocated, &param_coords_size,
            &err );
      // get the first coordinate
      if (err != 0 || intersect_entity_handles_size ==0)
         return false;
      // consider only the first intersection point
      for (int j=0; j<3; j++)
         newBoundary.push_back( intersect_coords[j]);
      free(intersect_entity_handles);
      free(intersect_coords);
      free(param_coords);
      // decide a new starting point for the ray
      currentLengthLeft -= gridSize;
      if (currentLengthLeft > 0)
      {
         for (int i=0; i<3; i++)
            currentPosition[i] += gridSize*currentDirection[i];
      }
      else
      {
         double point[3];
         while (currentLengthLeft < 0)
         {
            nextPolygonPt++;
            nextPolygonPt = nextPolygonPt%sizePolygon;
            if (nextPolygonPt==1)
            {
               reachedEnd = 1;
               break;
            }
            point[0] = xs[nextPolygonPt];
            point[1] = ys[nextPolygonPt];
            point[2] = zs[nextPolygonPt];
            for (int i=0; i<3; i++)
               currentDirection[i] = point[i] - nextPolyCoord[i];
            double length = sqrt(LENGTH_SQ(currentDirection));
            currentLengthLeft = currentLengthLeft+length;
            NORMALIZE(currentDirection);
            for (int j=0; j<3; j++)
               nextPolyCoord[j] = point[j];
         }
         // compute current position from nextPolygonPt
         for (int i=0; i<3; i++)
            currentPosition[i] = nextPolyCoord[i] - currentLengthLeft*currentDirection[i];

      }

   }// end while

   // newBoundary has the double coordinates; create a set in the mesh file with these
   // coordinates, and add those edge segments
   int newNodes = newBoundary.size()/3; // interleaved for vertex creation

   iBase_EntityHandle * new_verts_ptr = NULL;
   int new_verts_alloc = 0, new_verts_size = 0;
   iMesh_createVtxArr(meshIface, newNodes, iBase_INTERLEAVED,
                        &newBoundary[0], 3*newNodes,
                        &new_verts_ptr, &new_verts_alloc, &new_verts_size,
                        &err);
   if (iBase_SUCCESS != err) {
      std::cerr << "Couldn't create new vertices." << std::endl;
      return false;
   }

   iBase_EntitySetHandle entity_set;
   iMesh_createEntSet(meshIface,  /*in const int isList */ 1,
                                  /*out*/ &entity_set,
                                  &err);
   // add vertices to the set
   iMesh_addEntArrToSet(meshIface, new_verts_ptr, newNodes, entity_set, &err);
   if (iBase_SUCCESS != err) {
      std::cerr << "Couldn't add new vertices to set."
            << std::endl;
      return false;
   }

   // create connectivity for edges
   std::vector<iBase_EntityHandle> connect(2*newNodes);
   for (int i=0; i<newNodes; i++)
   {
      int nextNodeIndex = (i+1)%newNodes;
      connect[2*i] = new_verts_ptr[i];
      connect[2*i+1] = new_verts_ptr[nextNodeIndex];
   }
   iBase_EntityHandle * new_elem_ptr = NULL;
   int new_entity_handles_allocated = 0, new_entity_handles_size =0;
   int * status = NULL;
   int status_allocated = 0, status_size =0;
   iMesh_createEntArr(meshIface, iMesh_LINE_SEGMENT,
                             &connect[0], 2*newNodes,
                             &new_elem_ptr,
                             &new_entity_handles_allocated,
                             &new_entity_handles_size,
                             &status,
                             &status_allocated,
                             &status_size,
                             &err);
   // add segments to the set
   iMesh_addEntArrToSet(meshIface, new_elem_ptr, new_entity_handles_size, entity_set, &err);
   if (iBase_SUCCESS != err) {
      std::cerr << "Couldn't add new segments to set."
            << std::endl;
      return false;
      }
   // temp file
   const char * name = "newBound.vtk";
   iMesh_save(meshIface,
         /*in*/ entity_set,
         /*in*/  name,
         /*in*/ NULL,
         /*out*/ &err,
         /*in*/ 12,
         /*in*/ 0);

   free(new_elem_ptr);
   free(new_verts_ptr);
   // maybe delete some entity arrays
   return true;
}

#ifdef TEST_CAMEL

bool test_bdy_geom_grouped(CMEL *cmel) 
{
    // given an entity, return groups of entities in loops/shells

    // get all faces, find ones with multiple loops
  iBase_EntityHandle *gents = NULL;
  int result, gents_alloc = 0, gents_size = 0;
  iGeom_getEntities(cmel->geomIface, 0, iBase_FACE, 
                    &gents, &gents_alloc, &gents_size, &result);
  if (iBase_SUCCESS != result) {
    std::cerr << "Trouble getting gfaces" << std::endl;
    return false;
  }

  std::vector<int> nsconnect_faces;
  for (int i = 0; i < gents_size; i++) {
    std::vector<iBase_EntityHandle> bdy_groups;
    std::vector<int> bdy_groups_sizes;
    bool success = cmel->bdy_geom_grouped(gents[i], bdy_groups, bdy_groups_sizes);
    if (!success) {
      std::cerr << "Trouble bdy_geom_grouped for face " << cmel->get_gentity_id(gents[i])
                << std::endl;
    }
    if (bdy_groups_sizes.size() > 1) 
      nsconnect_faces.push_back(cmel->get_gentity_id(gents[i]));

    if (debug) {
      std::cerr << "Face " << cmel->get_gentity_id(gents[i]) << " edge(s) ";
      int num_edges = 0, num_loops = 1, *loop_sizes = &bdy_groups_sizes[0];
      for (int j = 0; j < bdy_groups.size(); j++) {
        
        std::cerr << cmel->get_gentity_id(bdy_groups[j]) << "("
                  << iGeom_getEgFcSense(cmel->geomIface, bdy_groups[j], gents[i]) << ") ";
        num_edges++;
        if (num_edges == *loop_sizes && num_loops < bdy_groups_sizes.size()) {
          loop_sizes++;
          std::cerr << "/ ";
          num_edges = 0;
          num_loops++;
        }
      }
      std::cerr << std::endl;
    }
  }
  
  std::cerr << "Found " << nsconnect_faces.size() << " non-simply-connected face(s)";
  if (!nsconnect_faces.empty()) {
    std::cerr << " with id(s): ";
    for (int i = 0; i < nsconnect_faces.size(); i++) 
      std::cerr << nsconnect_faces[i] << " ";
  }
  std::cerr << std::endl;

  return true;
}

/*

  bool bdy_coords_connect(std::vector<iBase_EntityHandle> &bdy_mesh, 
                          std::vector<int> &bdy_orientations, 
                          std::vector<iBase_EntityHandle> &bdy_verts,
                          std::vector<double> &coords, 
                          std::vector<int> &connect);

  bool bdy_elements_senses_grouped(iBase_EntityHandle gentity,
                                   std::vector<iBase_EntityHandle> &elements,
                                   std::vector<int> &element_senses,
                                   std::vector<int> &group_sizes);
  
  bool bdy_elements_senses(iBase_EntityHandle gentity,
                           std::vector<iBase_EntityHandle> &elements,
                           std::vector<int> &element_senses);
  
*/

int main(int argc, char **argv) 
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << "<geom_file> [mesh_size]" << std::endl;
    return 1;
  }

  iGeom_Instance geom = iGeom_newGeom(NULL, 0);
  iMesh_Instance mesh;
  int result;
  iMesh_newMesh(NULL, 0, &mesh, &result);
  iRel_Instance relate = iRel_newAssoc(NULL, 0);

    // read geometry
  iBase_ErrorType result = iGeom_load(geom, argv[1], NULL, 0);
  if (iBase_SUCCESS != result) {
    std::cerr << "ERROR : can not load a geometry from " << argv[1] << std::endl;
    return 1;
  }

  CMEL cmel(geom, mesh, relate);

  bool success = test_bdy_geom_grouped(&cmel);
  if (!success) {
    std::cerr << "CMEL::bdy_geom_grouped test failed." << std::endl;
  }
  
  return 0;
}

#endif
