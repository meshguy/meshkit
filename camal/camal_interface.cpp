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
#include "camal_interface.hpp"
#include "camel.hpp"
#include "CAMALGeomEval.hpp"
#include "CAMALSizeEval.hpp"

#include <iostream>
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <math.h>

#if CAMAL_VERSION > 500
  #include "CMLDefines.h"
  #include "CMLTriAdvance.hpp"
  #include "CMLTriDelaunay.hpp"
  #include "CMLTetMesher.hpp"
#else
  #include "CMLTriMesher.hpp"
  #include "CMLTetMesher.hpp"
#endif

const bool debug = false;
iBase_TagHandle index_tag;

bool CAMAL_bdy_loops_coords(CMEL *cmel,
                            iBase_EntityHandle gentity,
                            std::vector<int> &loops,
                            std::vector<int> &loop_sizes,
                            std::vector<iBase_EntityHandle> &bdy_verts,
                            std::vector<double> &coords) 
{
    // get the groups of mesh entities and their senses
  std::vector<int> loop_senses;
  
  std::vector<iBase_EntityHandle> loop_ents;
  bool success = cmel->bdy_elements_senses_grouped(gentity, loop_ents, 
                                                   loop_senses, loop_sizes);
  if (!success) return success;
  assert(loop_ents.size() == loop_senses.size());
  
    // get all vertices used in all elements
  int result;
  
  iBase_EntityHandle *tmp_verts = NULL;
  int tmp_verts_size = 0, tmp_verts_alloc = 0;
  int *offsets = NULL, offsets_size = 0, offsets_alloc = 0;
  iMesh_getEntArrAdj(cmel->meshIface, &loop_ents[0], loop_ents.size(),
                     iBase_VERTEX,
                     &tmp_verts, &tmp_verts_alloc, &tmp_verts_size,
                     &offsets, &offsets_alloc, &offsets_size,
                     &result);
  if (iBase_SUCCESS != result) return result;
    
  std::copy(tmp_verts, tmp_verts+tmp_verts_size, std::back_inserter(bdy_verts));
  free(tmp_verts);
  free(offsets);

    // make sorted & unique
  std::sort(bdy_verts.begin(), bdy_verts.end());
  bdy_verts.erase(std::unique(bdy_verts.begin(), bdy_verts.end()), bdy_verts.end());

    // assign ids to vertices
  success = cmel->assign_tmp_indices(bdy_verts, index_tag);
  if (!success) return success;

    // get connect list for this group
  tmp_verts = NULL;
  tmp_verts_size = 0;
  tmp_verts_alloc = 0;
  offsets = NULL; 
  offsets_size = 0;
  offsets_alloc = 0;
  iMesh_getEntArrAdj(cmel->meshIface, &loop_ents[0], loop_ents.size(),
                     iBase_VERTEX,
                     &tmp_verts, &tmp_verts_alloc, &tmp_verts_size,
                     &offsets, &offsets_alloc, &offsets_size,
                     &result);
  if (iBase_SUCCESS != result) return result;
    
    // reverse if necessary
  for (unsigned int i = 0; i < loop_ents.size(); i++) {
    if (loop_senses[i] == -1) {
      std::reverse(tmp_verts+offsets[i], tmp_verts+offsets[i+1]);
    }
  }

    // read indices into that list in place of vertex handles
  int *connect_verts_ptr = (int*) tmp_verts;
  int connect_verts_size = tmp_verts_size, connect_verts_alloc = connect_verts_size;
  iMesh_getIntArrData(cmel->meshIface, tmp_verts, tmp_verts_size, index_tag, 
                      &connect_verts_ptr, &connect_verts_alloc, 
                      &connect_verts_size, &result);
  if (iBase_SUCCESS != result) 
    return false;

    // now, stuff into loop list
  for (int *i = connect_verts_ptr; i < connect_verts_ptr+connect_verts_size; i+=2) {
    loops.push_back(*i);
  }
  
    // lastly, get the vertex coords
  coords.resize(3*bdy_verts.size());
  double *coords_ptr = &coords[0];
  int coords_size = 3*bdy_verts.size(), coords_alloc = coords_size;
  iMesh_getVtxArrCoords(cmel->meshIface, &bdy_verts[0], bdy_verts.size(), 
                        iBase_INTERLEAVED, &coords_ptr, &coords_alloc, &coords_size,
                        &result);
  if (iBase_SUCCESS != result) 
    return false;

  return true;
}

bool CAMAL_bdy_connect_coords(CMEL *cmel,
                              std::vector<iBase_EntityHandle> &bdy_elems, 
                              std::vector<int> &bdy_elem_senses,
                              std::vector<int> &bdy_connect,
                              std::vector<iBase_EntityHandle> &bdy_verts, 
                              std::vector<double> &bdy_coords) 
{
    // get connect list
  iBase_EntityHandle *tmp_verts = NULL;
  int tmp_verts_size = 0, tmp_verts_alloc = 0;
  int *offsets = NULL, offsets_size = 0, offsets_alloc = 0;
  int result;
  iMesh_getEntArrAdj(cmel->meshIface, &bdy_elems[0], bdy_elems.size(),
                     iBase_VERTEX,
                     &tmp_verts, &tmp_verts_alloc, &tmp_verts_size,
                     &offsets, &offsets_alloc, &offsets_size,
                     &result);
  if (iBase_SUCCESS != result) return result;
    
    // reverse if necessary
  std::vector<int>::iterator sense_it = bdy_elem_senses.begin();
  unsigned int i = 0;
  for (; i < bdy_elem_senses.size(); i++, sense_it++) {
    if (*sense_it == -1) {
      std::reverse(tmp_verts+offsets[i], tmp_verts+offsets[i+1]);
    }
  }

    // copy, then make sorted & unique
  std::copy(tmp_verts, tmp_verts+tmp_verts_size, std::back_inserter(bdy_verts));
  std::sort(bdy_verts.begin(), bdy_verts.end());
  bdy_verts.erase(std::unique(bdy_verts.begin(), bdy_verts.end()), bdy_verts.end()); 

    // assign ids to vertices
  bool success = cmel->assign_tmp_indices(bdy_verts, index_tag);
  if (!success) return false;

    // read indices into that list in place of vertex handles
  bdy_connect.resize(tmp_verts_size);
  int *connect_verts_ptr = &bdy_connect[0];
  int connect_verts_size = tmp_verts_size, connect_verts_alloc = connect_verts_size;
  iMesh_getIntArrData(cmel->meshIface, tmp_verts, tmp_verts_size, index_tag, 
                      &connect_verts_ptr, &connect_verts_alloc, 
                      &connect_verts_size, &result);
  if (iBase_SUCCESS != result) 
    return false;
  
    // lastly, get the vertex coords
  bdy_coords.resize(3*bdy_verts.size());
  double *coords_ptr = &bdy_coords[0];
  int coords_size = 3*bdy_verts.size(), coords_alloc = coords_size;
  iMesh_getVtxArrCoords(cmel->meshIface, &bdy_verts[0], bdy_verts.size(), 
                        iBase_INTERLEAVED, &coords_ptr, &coords_alloc, &coords_size,
                        &result);
  if (iBase_SUCCESS != result) 
    return false;

  return true;
}
  
bool CAMAL_mesh_entity(CMEL *cmel,
                       iBase_EntityHandle gentity, 
                       double mesh_size, int mesh_intervals,
                       const bool force_intervals,
                       std::vector<iBase_EntityHandle> &new_entities) 
{
    // assemble the bounding edges into loops of nodes
    // first get all the nodes
  std::vector<iBase_EntityHandle> bdy_verts;
  std::vector<double> bdy_coords;
  std::vector<int> connect;
  bool success;
  int new_points;

  CAMALGeomEval geom_eval(cmel->geomIface, gentity);
  geom_eval.set_mesh_size(mesh_size);
  
  iMesh_EntityTopology etop = iMesh_ALL_TOPOLOGIES;

  if (geom_eval.get_dimension() == 0) {
    return cmel->mesh_vertex(gentity, new_entities);
  }
  
  else if (geom_eval.get_dimension() == 1) {
    return cmel->mesh_curve(gentity, mesh_size, mesh_intervals, 
                            new_entities);
  }
  
  else if (geom_eval.get_dimension() == 2) {
      // connect was returned as edge connectivity represented as vertex indices;
      // revise to make loops out of them
    std::vector<int> loops, loop_sizes;
    success = CAMAL_bdy_loops_coords(cmel, gentity, loops, loop_sizes,
                                     bdy_verts, bdy_coords);
    if (!success) {
      std::cerr << "Couldn't get bounding mesh for surface." << std::endl;
      return success;
    }

    // set boundary mesh size
    int n_bdy_verts = bdy_verts.size();
    if (n_bdy_verts != bdy_coords.size()/3) {
      std::cerr << "# of boundary vertices are not matched." << std::endl;
      return false;
    }


    if (mesh_size < 0.) {
      double tot_length = 0.;
      for (unsigned int i = 0; i < n_bdy_verts; i++) {
	tot_length += sqrt((bdy_coords[3*i] - bdy_coords[3*((i + 1)%n_bdy_verts)])*
			   (bdy_coords[3*i] - bdy_coords[3*((i + 1)%n_bdy_verts)]) +
			   (bdy_coords[3*i + 1] - bdy_coords[3*((i + 1)%n_bdy_verts) + 1])*
			   (bdy_coords[3*i + 1] - bdy_coords[3*((i + 1)%n_bdy_verts) + 1]) +
			   (bdy_coords[3*i + 2] - bdy_coords[3*((i + 1)%n_bdy_verts) + 2])*
			   (bdy_coords[3*i + 2] - bdy_coords[3*((i + 1)%n_bdy_verts) + 2]));
      }
      mesh_size = tot_length/n_bdy_verts;
      geom_eval.set_mesh_size(mesh_size);
    }

#if CAMAL_VERSION > 500
    CAMALSizeEval size_eval(mesh_size);
#endif

      // pass to CAMAL
    if (debug) {
      std::cout << "Surface " << cmel->get_gentity_id(gentity)
		<< ", mesh_size = " << mesh_size
                << ", boundary mesh: " << std::endl;
      std::cout << bdy_coords.size()/3 << "  " << loop_sizes.size() << std::endl;
      for (unsigned int i = 0; i < bdy_coords.size()/3; i++)
        std::cout << bdy_coords[3*i] << "  " 
                  << bdy_coords[3*i+1] << "  " 
                  << bdy_coords[3*i+2] << std::endl;
      
      for (std::vector<int>::iterator vit = loop_sizes.begin();
           vit != loop_sizes.end(); vit++)
        std::cout << *vit << "  ";
      
      std::cout << std::endl;
      for (std::vector<int>::iterator vit = loops.begin();
           vit != loops.end(); vit++)
        std::cout << *vit << std::endl;
      
      int err;
      iBase_EntitySetHandle outset;
      std::string outfile;
      std::stringstream ss;
      ss << "boundary";
      ss << mesh_intervals;
      ss >> outfile;
      outfile += ".vtk";
      iMesh_createEntSet(cmel->meshIface, false, &outset, &err);
      iMesh_addEntArrToSet(cmel->meshIface, &bdy_verts[0], bdy_verts.size(), outset, &err);
      iMesh_save(cmel->meshIface, outset, outfile.c_str(), 0,
		 &err, outfile.length(), 0);
    }
      
#if CAMAL_VERSION > 500
    CMLTriAdvance tri_mesher(&geom_eval);
#else
    CMLTriMesher tri_mesher(&geom_eval);
#endif

    success = tri_mesher.set_boundary_mesh(bdy_coords.size()/3, &bdy_coords[0],
                                           (int)loop_sizes.size(), 
                                           &loop_sizes[0], &loops[0]);
    if (!success) {
      std::cerr << "Failed setting boundary mesh" << std::endl;
      return success;
    }

#if CAMAL_VERSION > 500 
    tri_mesher.set_sizing_function(CML::LINEAR_SIZING);
#endif

      // generate the mesh
    int num_tris;
    success = tri_mesher.generate_mesh(new_points, num_tris);
    if (!success) {
      std::cerr << "Failed generating mesh" << std::endl;
      return success;
    }
    
    std::cout << "Meshed surface " << cmel->get_gentity_id(gentity) << " with "
	      << new_points << " new vertices and " << num_tris << " triangles." << std::endl;

      // get the generated mesh
    //bdy_coords.resize(3*(bdy_verts.size() + new_points));
    bdy_coords.resize(3* new_points);
    connect.resize(3*num_tris);
    success = tri_mesher.get_mesh(new_points, &bdy_coords[0], num_tris, &connect[0]);
    if (!success) {
      std::cerr << "Failed get generated mesh" << std::endl;
      return success;
    }

    etop = iMesh_TRIANGLE;
  }
  else if (geom_eval.get_dimension() == 3) {

      // get bdy elem handles
    std::vector<iBase_EntityHandle> bdy_elems;
    std::vector<int> bdy_elem_senses;
    success = cmel->bdy_elements_senses(gentity, bdy_elems, bdy_elem_senses);
    if (!success) {
      std::cerr << "Trouble getting bdy elements & senses." << std::endl;
      return success;
    }
    
      // convert to indexed list
    std::vector<int> bdy_conn;
    success = CAMAL_bdy_connect_coords(cmel, bdy_elems, bdy_elem_senses,
                                       bdy_conn, bdy_verts, bdy_coords);
    if (!success) {
      std::cerr << "Trouble getting bdy connect & coords." << std::endl;
      return success;
    }

    if (debug) {
      int result;
      iMesh_save(cmel->meshIface, 0, "test.h5m", 0, &result, 8, 0);
      if (iBase_SUCCESS != result) return false;
    }
    
    CMLTetMesher tet_mesher;
    success = tet_mesher.set_boundary_mesh(bdy_coords.size()/3, &bdy_coords[0],
					   bdy_conn.size()/3, &bdy_conn[0]);
    if (!success) {
      std::cerr << "Failed setting boundary mesh" << std::endl;
      return success;
    }
  
      // generate the mesh
    int num_tets;
    success = tet_mesher.generate_mesh(new_points, num_tets);
    if (!success) {
      std::cerr << "Failed generating mesh" << std::endl;
      return success;
    }
    
     // get the generated mesh
    //bdy_coords.resize(3*(bdy_verts.size() + new_points));
    bdy_coords.resize(3*new_points);
    connect.resize(4*num_tets);
    success = tet_mesher.get_mesh(new_points, &bdy_coords[0], num_tets, &connect[0]);
    if (!success) {
      std::cerr << "Failed get generated mesh" << std::endl;
      return success;
    }

    etop = iMesh_TETRAHEDRON;
  }

    // put new mesh back into interface
  success = cmel->create_vertices_elements(gentity, bdy_verts, 
                                           &bdy_coords[3*bdy_verts.size()], new_points-bdy_verts.size(),
                                           connect, etop, new_entities);
  if (!success)
    std::cerr << "Failed create generated mesh" << std::endl;

  return success;
}
