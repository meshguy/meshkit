/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at boundary
 *  option) any later version.
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
// the paver is here:
#include "CMLPaver.hpp"
//  #include "CMLTetMesher.hpp"
#else
#include "CMLTriMesher.hpp"
#include "CMLTetMesher.hpp"
#endif

bool debug = false;
iBase_TagHandle index_tag;

bool CAMAL_bdy_loops_coords(CMEL *cmel, iBase_EntityHandle gentity,
      std::vector<int> &loops, std::vector<int> &loop_sizes, std::vector<
            iBase_EntityHandle> &bdy_verts, std::vector<double> &coords) {
   // get the groups of mesh entities and their senses
   std::vector<int> loop_senses;

   std::vector<iBase_EntityHandle> loop_ents;
   bool success = cmel->bdy_elements_senses_grouped(gentity, loop_ents,
         loop_senses, loop_sizes);
   if (!success)
      return success;
   assert(loop_ents.size() == loop_senses.size());
   if (debug) {
      std::cout << " Surface: " << cmel->get_gentity_id(gentity)
            << " loop_ents size: " << loop_ents.size() << std::endl;
      for (unsigned int i = 0; i < loop_ents.size(); i++)
         std::cout << std::hex << loop_ents[i] << " " << std::dec << " sense: "
               << loop_senses[i] << std::endl;
   }

   // get all vertices used in all elements
   int result;

   iBase_EntityHandle *tmp_verts = NULL;
   int tmp_verts_size = 0, tmp_verts_alloc = 0;
   int *offsets = NULL, offsets_size = 0, offsets_alloc = 0;
   iMesh_getEntArrAdj(cmel->meshIface, &loop_ents[0], loop_ents.size(),
         iBase_VERTEX, &tmp_verts, &tmp_verts_alloc, &tmp_verts_size, &offsets,
         &offsets_alloc, &offsets_size, &result);
   if (iBase_SUCCESS != result)
      return result;

   std::copy(tmp_verts, tmp_verts + tmp_verts_size, std::back_inserter(
         bdy_verts));
   free(tmp_verts);
   free(offsets);

   // make sorted & unique
   std::sort(bdy_verts.begin(), bdy_verts.end());
   bdy_verts.erase(std::unique(bdy_verts.begin(), bdy_verts.end()),
         bdy_verts.end());

   // assign ids to vertices
   success = cmel->assign_tmp_indices(bdy_verts, index_tag);
   if (!success)
      return success;

   // get connect list for this group
   tmp_verts = NULL;
   tmp_verts_size = 0;
   tmp_verts_alloc = 0;
   offsets = NULL;
   offsets_size = 0;
   offsets_alloc = 0;
   iMesh_getEntArrAdj(cmel->meshIface, &loop_ents[0], loop_ents.size(),
         iBase_VERTEX, &tmp_verts, &tmp_verts_alloc, &tmp_verts_size, &offsets,
         &offsets_alloc, &offsets_size, &result);
   if (iBase_SUCCESS != result)
      return result;

   // reverse if necessary
   for (unsigned int i = 0; i < loop_ents.size(); i++) {
      if (loop_senses[i] == -1) {
         std::reverse(tmp_verts + offsets[i], tmp_verts + offsets[i + 1]);
      }
   }

   // read indices into that list in place of vertex handles
   int *connect_verts_ptr = (int*) tmp_verts;
   int connect_verts_size = tmp_verts_size, connect_verts_alloc =
         connect_verts_size;
   iMesh_getIntArrData(cmel->meshIface, tmp_verts, tmp_verts_size, index_tag,
         &connect_verts_ptr, &connect_verts_alloc, &connect_verts_size, &result);
   if (iBase_SUCCESS != result)
      return false;

   // now, stuff into loop list
   for (int *i = connect_verts_ptr; i < connect_verts_ptr + connect_verts_size; i
         += 2) {
      loops.push_back(*i);
   }

   // lastly, get the vertex coords
   coords.resize(3 * bdy_verts.size());
   double *coords_ptr = &coords[0];
   int coords_size = 3 * bdy_verts.size(), coords_alloc = coords_size;
   iMesh_getVtxArrCoords(cmel->meshIface, &bdy_verts[0], bdy_verts.size(),
         iBase_INTERLEAVED, &coords_ptr, &coords_alloc, &coords_size, &result);
   if (iBase_SUCCESS != result)
      return false;

   return true;
}

bool CAMAL_bdy_connect_coords(CMEL *cmel,
      std::vector<iBase_EntityHandle> &bdy_elems,
      std::vector<int> &bdy_elem_senses, std::vector<int> &bdy_connect,
      std::vector<iBase_EntityHandle> &bdy_verts,
      std::vector<double> &bdy_coords) {
   // get connect list
   iBase_EntityHandle *tmp_verts = NULL;
   int tmp_verts_size = 0, tmp_verts_alloc = 0;
   int *offsets = NULL, offsets_size = 0, offsets_alloc = 0;
   int result;
   iMesh_getEntArrAdj(cmel->meshIface, &bdy_elems[0], bdy_elems.size(),
         iBase_VERTEX, &tmp_verts, &tmp_verts_alloc, &tmp_verts_size, &offsets,
         &offsets_alloc, &offsets_size, &result);
   if (iBase_SUCCESS != result)
      return result;

   // reverse if necessary
   std::vector<int>::iterator sense_it = bdy_elem_senses.begin();
   unsigned int i = 0;
   for (; i < bdy_elem_senses.size(); i++, sense_it++) {
      if (*sense_it == -1) {
         std::reverse(tmp_verts + offsets[i], tmp_verts + offsets[i + 1]);
      }
   }

   // copy, then make sorted & unique
   std::copy(tmp_verts, tmp_verts + tmp_verts_size, std::back_inserter(
         bdy_verts));
   std::sort(bdy_verts.begin(), bdy_verts.end());
   bdy_verts.erase(std::unique(bdy_verts.begin(), bdy_verts.end()),
         bdy_verts.end());

   // assign ids to vertices
   bool success = cmel->assign_tmp_indices(bdy_verts, index_tag);
   if (!success)
      return false;

   // read indices into that list in place of vertex handles
   bdy_connect.resize(tmp_verts_size);
   int *connect_verts_ptr = &bdy_connect[0];
   int connect_verts_size = tmp_verts_size, connect_verts_alloc =
         connect_verts_size;
   iMesh_getIntArrData(cmel->meshIface, tmp_verts, tmp_verts_size, index_tag,
         &connect_verts_ptr, &connect_verts_alloc, &connect_verts_size, &result);
   if (iBase_SUCCESS != result)
      return false;

   // lastly, get the vertex coords
   bdy_coords.resize(3 * bdy_verts.size());
   double *coords_ptr = &bdy_coords[0];
   int coords_size = 3 * bdy_verts.size(), coords_alloc = coords_size;
   iMesh_getVtxArrCoords(cmel->meshIface, &bdy_verts[0], bdy_verts.size(),
         iBase_INTERLEAVED, &coords_ptr, &coords_alloc, &coords_size, &result);
   if (iBase_SUCCESS != result)
      return false;

   return true;
}

bool CAMAL_mesh_entity(CMEL *cmel, iBase_EntityHandle gentity,
      double mesh_size, int mesh_intervals, const bool force_intervals,
      std::vector<iBase_EntityHandle> &new_entities, const bool quadMesh) {
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
      return cmel->mesh_curve(gentity, mesh_size, mesh_intervals, new_entities);
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
      if (n_bdy_verts != bdy_coords.size() / 3) {
         std::cerr << "# of boundary vertices are not matched." << std::endl;
         return false;
      }

      if (mesh_size < 0.) {
         double tot_length = 0.;
         unsigned int start_current_loop = 0;
         for (unsigned int k = 0; k < loop_sizes.size(); k++) {
            // for each loop, compute the edge lengths individually
            int current_loop_size = loop_sizes[k];
            for (unsigned int i = 0; i < current_loop_size; i++) {
               unsigned int i1 = loops[start_current_loop + i];
               unsigned int i2 = loops[start_current_loop + (i + 1)
                     % current_loop_size];
               tot_length += sqrt((bdy_coords[3 * i1] - bdy_coords[3 * i2])
                     * (bdy_coords[3 * i1] - bdy_coords[3 * i2])
                     + (bdy_coords[3 * i1 + 1] - bdy_coords[3 * i2 + 1])
                           * (bdy_coords[3 * i1 + 1] - bdy_coords[3 * i2 + 1])
                     + (bdy_coords[3 * i1 + 2] - bdy_coords[3 * i2 + 2])
                           * (bdy_coords[3 * i1 + 2] - bdy_coords[3 * i2 + 2]));
            }
            start_current_loop += current_loop_size;
         }
         mesh_size = tot_length / n_bdy_verts;
         geom_eval.set_mesh_size(mesh_size);
      }

#if CAMAL_VERSION > 500
      CAMALSizeEval size_eval(mesh_size);
#endif

      // pass to CAMAL
      if (debug) {
         std::cout << "Surface " << cmel->get_gentity_id(gentity)
               << ", mesh_size = " << mesh_size << ", boundary mesh: "
               << std::endl;
         std::cout << bdy_coords.size() / 3 << "  " << loop_sizes.size()
               << std::endl;
         for (unsigned int i = 0; i < bdy_coords.size() / 3; i++)
            std::cout << bdy_coords[3 * i] << "  " << bdy_coords[3 * i + 1]
                  << "  " << bdy_coords[3 * i + 2] << std::endl;

         for (std::vector<int>::iterator vit = loop_sizes.begin(); vit
               != loop_sizes.end(); vit++)
            std::cout << *vit << "  ";

         std::cout << std::endl;
         for (std::vector<int>::iterator vit = loops.begin(); vit
               != loops.end(); vit++)
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
         iMesh_addEntArrToSet(cmel->meshIface, &bdy_verts[0], bdy_verts.size(),
               outset, &err);
         iMesh_save(cmel->meshIface, outset, outfile.c_str(), 0, &err,
               outfile.length(), 0);
      }

      if (quadMesh) {
         // start copy
         //CAMALSizeEval size_eval(mesh_size);
         CMLPaver pave_mesher(&geom_eval, &size_eval);
         // the loops are established, get the mesh sets from each edge, and create the boundary loops arrays
         // mesh faces;

         // set only num_points_out -1 , because the last one is repeated
         success = pave_mesher.set_boundary_mesh(bdy_coords.size() / 3,
               &bdy_coords[0], (int) loop_sizes.size(), &loop_sizes[0],
               &loops[0]);
         if (!success) {
            std::cerr << "Failed setting boundary mesh" << std::endl;
            return false;
         }

         pave_mesher.set_sizing_function(CML::LINEAR_SIZING);

         // generate the mesh
         int num_quads;
         success = pave_mesher.generate_mesh(new_points, num_quads);
         if (!success) {
            std::cerr << "Failed generating mesh" << std::endl;
            //return success;
         }

         std::cout << "Meshed surface with " << new_points
               << " new vertices and " << num_quads << " quadrilaterals."
               << std::endl;

         // get the generated mesh
         bdy_coords.resize(3 * new_points);
         //std::vector<double> new_coords;
         //new_coords.resize(3 * new_points);
         connect.resize(4 * num_quads);
         success = pave_mesher.get_mesh(new_points, &bdy_coords[0], num_quads,
               &connect[0]);
         if (!success) {
            std::cerr << "Failed to get generated mesh" << std::endl;
            //return success;
         }
         // end copy
         etop = iMesh_QUADRILATERAL;
      } else {
#if CAMAL_VERSION > 500
         CMLTriAdvance tri_mesher(&geom_eval);
#else
         CMLTriMesher tri_mesher(&geom_eval);
#endif

         success = tri_mesher.set_boundary_mesh(bdy_coords.size() / 3,
               &bdy_coords[0], (int) loop_sizes.size(), &loop_sizes[0],
               &loops[0]);
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

         std::cout << "Meshed surface " << cmel->get_gentity_id(gentity)
               << " with " << new_points << " new vertices and " << num_tris
               << " triangles." << std::endl;

         // get the generated mesh
         //bdy_coords.resize(3*(bdy_verts.size() + new_points));
         bdy_coords.resize(3 * new_points);
         connect.resize(3 * num_tris);
         success = tri_mesher.get_mesh(new_points, &bdy_coords[0], num_tris,
               &connect[0]);
         if (!success) {
            std::cerr << "Failed get generated mesh" << std::endl;
            return success;
         }

         etop = iMesh_TRIANGLE;
      }
   } else if (geom_eval.get_dimension() == 3) {

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
         if (iBase_SUCCESS != result)
            return false;
      }
#if 0
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
#endif
   }

   // put new mesh back into interface
   success = cmel->create_vertices_elements(gentity, bdy_verts, &bdy_coords[3
         * bdy_verts.size()], new_points - bdy_verts.size(), connect, etop,
         new_entities);
   if (!success)
      std::cerr << "Failed create generated mesh" << std::endl;

   return success;
}

bool CAMAL_mesh_trimmed_surface(CMEL * cmel, iBase_EntityHandle surface,
      double mesh_size, std::vector<iBase_EntityHandle> &new_entities,
      std::vector<double> trimmingBoundary, const bool quadMesh) {
   std::vector<iBase_EntityHandle> bdy_verts;// these should be empty now
   std::vector<double> bdy_coords;
   std::vector<int> connect;
   bool success = true; // have a positive attitude
   int new_points;
   iMesh_EntityTopology etop = iMesh_ALL_TOPOLOGIES;
   std::vector<double> coords;

   CAMALGeomEval geom_eval(cmel->geomIface, surface);
   geom_eval.set_mesh_size(mesh_size);
#if CAMAL_VERSION > 500
   CAMALSizeEval size_eval(mesh_size);
#endif
   // boundary is in a loop
   int numBoundPoints = trimmingBoundary.size() / 3;
   // the new boundary will be decided based on mesh size
   double totalLength = 0;
   std::vector<double> lengs;
   lengs.push_back(0);// start with 0
   for (int i = 0; i < numBoundPoints; i++) {
      int nextIndex = (i + 1) % numBoundPoints;
      double
            lenSeg =
                  DIST2( (&trimmingBoundary[3*i]), (&trimmingBoundary[3*nextIndex]) );
      totalLength += lenSeg;
      lengs.push_back(totalLength);
   }
   int numSegs = (int) (totalLength / mesh_size);

   double newMeshSize = totalLength / numSegs;
   // generate new points on the segments (in natural parametric space)
   // use a sort of bsearch to get the index in lengs
   // create some vertices /edges  in the mesh, which will be used for quad generation
   // we will not create them at the "boundary mesh" stage
   // although it would be easier to just have another gentity with a new edge as boundary
   // how to do that in igeom?

   // first node is the first one on trimming boundary (for sure on our surface)
   // std::vector <double> bdy_coords;
   int i = 0;
   for (i = 0; i < 3; i++)
      bdy_coords.push_back(trimmingBoundary[i]);
   // find a new point in the param space of the boundary curve
   double param = 0; // current parameter
   double sizeLocal = mesh_size;
   int ix = 0;
   while (param < totalLength - 3 * sizeLocal) {
      // latest point pushed is at index ix
      size_eval.size_at_point(bdy_coords[3 * ix], bdy_coords[3 * ix + 1],
            bdy_coords[3 * ix + 2], sizeLocal, 0);
      param += sizeLocal;
      // this is the par position of the next point, in natural coordinate
      double * pos = std::lower_bound(&lengs[0], &lengs[numBoundPoints], param);
      if (pos == &lengs[0])
         return false;
      pos = pos - 1; // get the previous position, it cannot be 0
      double extraLen = param - *pos;
      int index = pos - &lengs[0];
      int nextV = (index + 1) % numBoundPoints; // it can't be here, but hey...
      double direction[3] = { trimmingBoundary[3 * nextV] - trimmingBoundary[3
            * index], trimmingBoundary[3 * nextV + 1] - trimmingBoundary[3
            * index + 1], trimmingBoundary[3 * nextV + 2] - trimmingBoundary[3
            * index + 2] };
      NORMALIZE(direction);
      double pp[3];
      for (i = 0; i < 3; i++) {
         pp[i] = trimmingBoundary[3 * index + i] + direction[i] * extraLen;
      }
      // now see if this is closer than what we want it to be
      geom_eval.move_to_surface(pp[0], pp[1], pp[2]);
      for (i = 0; i < 3; i++)
         bdy_coords.push_back(pp[i]);
      ix++;

   }
   // the current index must be odd or even, less than 3, 4 segments left
   // to generate; it will be a little different
   int segRemain = (int) (totalLength - param) / sizeLocal;
   if ((bdy_coords.size() / 3 + segRemain) % 2 == 0 && quadMesh) {
      segRemain++;// make total number even
      std::cout << "nb segs so far: " << bdy_coords.size() / 3
            << " local size: " << sizeLocal << "\n";
      std::cout << "remains increased to " << segRemain << "\n";
      std::cout << "localSize reduced to " << sizeLocal << " \n";
      // 
   }
   sizeLocal = (totalLength - param) / segRemain;// this should be 3 or 4...
   for (int k = 1; k < segRemain; k++) {
      double parPosition = param + k * sizeLocal;
      double * pos = std::lower_bound(&lengs[0], &lengs[numBoundPoints],
            parPosition);
      if (pos == &lengs[0])
         return false;
      pos = pos - 1; // get the previous position, it cannot be 0
      double extraLen = parPosition - *pos;
      int index = pos - &lengs[0];
      int nextV = (index + 1) % numBoundPoints;
      double direction[3] = { trimmingBoundary[3 * nextV] - trimmingBoundary[3
            * index], trimmingBoundary[3 * nextV + 1] - trimmingBoundary[3
            * index + 1], trimmingBoundary[3 * nextV + 2] - trimmingBoundary[3
            * index + 2] };
      NORMALIZE(direction);
      double pp[3];
      for (i = 0; i < 3; i++) {
         pp[i] = trimmingBoundary[3 * index + i] + direction[i] * extraLen;
      }
      // now see if this is closer than what we want it to be
      geom_eval.move_to_surface(pp[0], pp[1], pp[2]);
      for (i = 0; i < 3; i++)
         bdy_coords.push_back(pp[i]);
   }
   std::cout << " boundary has " << bdy_coords.size() / 3 << " nodes\n";
   std::vector<int> loops;
   std::vector<int> loop_sizes;
   numSegs = bdy_coords.size() / 3;
   for (i = 0; i < numSegs; i++)
      loops.push_back(i);
   loop_sizes.push_back(numSegs);// one loop , external ....
   if (quadMesh) {
      // start copy
      //CAMALSizeEval size_eval(mesh_size);
      CMLPaver pave_mesher(&geom_eval, &size_eval);
      // the loops are established, get the mesh sets from each edge, and create the boundary loops arrays
      // mesh faces;

      // set only num_points_out -1 , because the last one is repeated
      success = pave_mesher.set_boundary_mesh(bdy_coords.size() / 3,
            &bdy_coords[0], (int) loop_sizes.size(), &loop_sizes[0], &loops[0]);
      if (!success) {
         std::cerr << "Failed setting boundary mesh" << std::endl;
         return false;
      }

      pave_mesher.set_sizing_function(CML::LINEAR_SIZING);

      // generate the mesh
      int num_quads;
      success = pave_mesher.generate_mesh(new_points, num_quads);
      if (!success) {
         std::cerr << "Failed generating mesh" << std::endl;
         //return success;
      }

      std::cout << "Meshed surface with " << new_points << " new vertices and "
            << num_quads << " quadrilaterals." << std::endl;

      // get the generated mesh
      //std::vector<double> coords;
      coords.resize(3 * new_points);
      //std::vector<double> new_coords;
      //new_coords.resize(3 * new_points);
      connect.resize(4 * num_quads);
      success = pave_mesher.get_mesh(new_points, &coords[0], num_quads,
            &connect[0]);
      if (!success) {
         std::cerr << "Failed to get generated mesh" << std::endl;
         //return success;
      }
      etop = iMesh_QUADRILATERAL;
   } else {
#if CAMAL_VERSION > 500
      CMLTriAdvance tri_mesher(&geom_eval);
#else
      CMLTriMesher tri_mesher(&geom_eval);
#endif

      success = tri_mesher.set_boundary_mesh(bdy_coords.size() / 3,
            &bdy_coords[0], (int) loop_sizes.size(), &loop_sizes[0], &loops[0]);
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

      std::cout << "Meshed surface "
            << " with " << new_points << " new vertices and " << num_tris
            << " triangles." << std::endl;

      // get the generated mesh
      //bdy_coords.resize(3*(bdy_verts.size() + new_points));
      coords.resize(3 * new_points);
      connect.resize(3 * num_tris);
      success = tri_mesher.get_mesh(new_points, &coords[0], num_tris,
            &connect[0]);
      if (!success) {
         std::cerr << "Failed get generated mesh" << std::endl;
         return success;
      }

      etop = iMesh_TRIANGLE;
   }
   // put new mesh back into interface
   success = cmel->create_vertices_elements(surface, bdy_verts, &coords[0],
         new_points, connect, etop, new_entities);
   return success;
}
