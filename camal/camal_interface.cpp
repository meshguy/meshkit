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

bool mesh_line_on_surface(CAMALGeomEval & geom_eval,
                          CAMALSizeEval & size_eval,
                          double mesh_size,
                          std::vector<double> & points,
                          bool periodic,
                          int force_evenify, // 0 no action, 1 force odd
                                             //  2 force even
                          std::vector<double> & mesh_points)
{

   int szLine = points.size()/3;
   double * point1 = &points[0];
   double * point2 = &points[3*(szLine-1)];

   if (periodic)
        point2 = point1;
   int numPoints = points.size() / 3;

   // first, mesh the grounding line
   double lineLength = 0;
   std::vector<double> lengsLine;
   lengsLine.push_back(0);// start with 0
   int i=0;
   for (i = 0; i < numPoints-1; i++) {

      double
            lenSeg =
                  DIST2( (&points[3*i]), (&points[3*(i+1)]) );
      lineLength += lenSeg;
      lengsLine.push_back(lineLength);
   }
   if (periodic)
   {
      double lenSeg = DIST2((&points[0]), (&points[3*numPoints-3]) );
      lineLength += lenSeg;
      lengsLine.push_back(lineLength);
   }
   int numSegsLine = (int) (lineLength / mesh_size);

   double newMeshSizeLine = lineLength / numSegsLine;
   // generate new points on the segments (in natural parametric space
   //                           on grounding line, between index 1 and 2)
   // use a sort of bsearch to get the index in lengsLine
   // create some vertices /edges  in the mesh, which will be used for quad generation
   // first node is the first one on line (for sure on our surface)
   // std::vector <double> bdy_coords;
   for (i = 0; i < 3; i++)
      mesh_points.push_back(point1[i]);
   // find a new point in the param space of the boundary curve
   double param = 0; // current parameter
   double sizeLocal = mesh_size;
   int ix = 0;
   while (param < lineLength - 3 * sizeLocal) {
      // latest point pushed is at index ix
      size_eval.size_at_point(mesh_points[3 * ix], mesh_points[3 * ix + 1],
            mesh_points[3 * ix + 2], sizeLocal, 0);
      param += sizeLocal;
      // this is the par position of the next point, in natural coordinate
      double * pos = std::lower_bound(&lengsLine[0], &lengsLine[szLine], param);
      if (pos == &lengsLine[0])
         return false;
      pos = pos - 1; // get the previous position, it cannot be 0
      double extraLen = param - *pos;
      int index = pos - &lengsLine[0];
      int nextV = (index + 1) ; // it can't be here, but hey...
      double direction[3] = { points[3 * nextV] - points[3 * index],
            points[3 * nextV + 1] - points[3 * index + 1],
            points[3 * nextV + 2] - points[3 * index + 2] };
      NORMALIZE(direction);
      double pp[3];
      for (i = 0; i < 3; i++) {
         pp[i] = points[3 * index + i] + direction[i] * extraLen;
      }
      // now see if this is closer than what we want it to be
      geom_eval.move_to_surface(pp[0], pp[1], pp[2]);
      for (i = 0; i < 3; i++)
         mesh_points.push_back(pp[i]);
      ix++;
   }
   // to generate; it will be a little different
   int segRemainLine = (int) (lineLength - param) / sizeLocal;
   int numSegSoFar = mesh_points.size()/3 -1;
   int resultSeg= (segRemainLine + numSegSoFar)%2;
   if (force_evenify)
   {
      if(1==force_evenify && 0==resultSeg) // force odd
         segRemainLine++;
      if (2==force_evenify && 1==resultSeg) // force even
         segRemainLine++;
   }
   sizeLocal = (lineLength - param) / segRemainLine;// this should be 3 or 4...
   for (int k = 1; k < segRemainLine; k++) {
      double parPosition = param + k * sizeLocal;
      double * pos = std::lower_bound(&lengsLine[0], &lengsLine[szLine],
            parPosition);
      if (pos == &lengsLine[0])
         return false;
      pos = pos - 1; // get the previous position, it cannot be 0
      double extraLen = parPosition - *pos;
      int index = pos - &lengsLine[0];
      int nextV = (index + 1) ;
      double direction[3] = { points[3 * nextV] - points[3 * index],
            points[3 * nextV + 1] - points[3 * index + 1],
            points[3 * nextV + 2] - points[3 * index + 2] };

      NORMALIZE(direction);
      double pp[3];
      for (i = 0; i < 3; i++) {
         pp[i] = points[3 * index + i] + direction[i] * extraLen;
      }
      // now see if this is closer than what we want it to be
      geom_eval.move_to_surface(pp[0], pp[1], pp[2]);
      for (i = 0; i < 3; i++)
         mesh_points.push_back(pp[i]);
   }
// also, add the last point of the line, point 2, if not periodic
   if (!periodic)
     for (i=0; i<3; i++)
        mesh_points.push_back(point2[i]);
   std::cout << " line has " << mesh_points.size() / 3 << " nodes\n";

   if (debug)
   {
      for (int i=0; i<mesh_points.size()/3; i++)
      {
         int i3 = i*3;
         std::cout<< i<< " " <<  mesh_points[i3]<< " " << mesh_points[i3+1]
               << " " << mesh_points[i3+2] << "\n";
      }
   }
   return true;// success

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

   success =  mesh_line_on_surface( geom_eval,
                              size_eval,
                              mesh_size,
                              trimmingBoundary,
                             /*bool periodic*/ true,
                            /* int force_evenify*/ 2 , // 0 no action, 1 force odd
                                                //  2 force even
                            bdy_coords);

   std::cout << " boundary has " << bdy_coords.size() / 3 << " nodes\n";
   std::vector<int> loops;
   std::vector<int> loop_sizes;
   int numSegs = bdy_coords.size() / 3;
   int i=0;
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
// this will create 2 loops, separated by the internalBoundary (grounding line)
// each will be meshed separately, but, of course, we start with the common line,
// the grounding line

bool CAMAL_mesh_trimmed_surface_with_grounding_line(CMEL * cmel, iBase_EntityHandle surface,
      double mesh_size, std::vector<iBase_EntityHandle> &new_entities,
      std::vector<double> trimmingBoundary, std::vector<double> internalBoundary,
      const bool quadMesh )
{
   int i = 0;
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
   // create 2 loops from one loop and one internal boundary
   // first, find the index in the loop of the start and end of internal boundary
   int szLine = internalBoundary.size()/3;
   double * point1 = &internalBoundary[0];
   double * point2 = &internalBoundary[3*(szLine-1)];

   int numBoundPoints = trimmingBoundary.size() / 3;

   int index1=-1, index2=-1;
   for (i=0; i<numBoundPoints; i++)
   {
      double * currentPoint = &(trimmingBoundary[3*i]);
      double dist1 = DIST2(point1, currentPoint);
      double dist2 = DIST2(point2, currentPoint);
      if (dist1<1.e-2)
         index1= i;
      if(dist2<1.e-2)
         index2=i;
   }
   if (index1<0 || index2<0)
   {
      // the grounding line is not what we want, we cannot separate
      std::cerr<<"can't use defined grounding line; it must connect outer boundary points\n";
      return false;
   }
   std::vector<double> gr_line_mesh;
   success =  mesh_line_on_surface( geom_eval,
                           size_eval,
                           mesh_size,
                           internalBoundary,
                          /*bool periodic*/ false,
                         /* int force_evenify*/ 0 , // 0 no action, 1 force odd
                                             //  2 force even
                          gr_line_mesh);
   if (!success)
      return success;
   int num_seg_ground_line = gr_line_mesh.size()/3 -1;
   int oddGroundLine = num_seg_ground_line%2; // 1 means odd, 0 means even

   /* if index1<index2, loop 1 will be grounding line, index2 to end, 0 - index1
                     loop 2 will be  reverse grounding line, index1, index2


    *   if index1>index2, then loop1 will be grounding line, index2, index1
    *                          loop2 will be reverse grounding line, index1, end, 0, index2
    *

    so, break the external loop in 2, from small index to large index

    */

   int largerIndex = (index1>index2)? index1 : index2;
   int smallerIndex = (index1<index2) ? index1 : index2;
   // first loop is smallerIndex to largerIndex
   // second loop is largerIndex, end, 0 to smallerIndex
   std::vector<double> boundary1;

   for (i=smallerIndex; i<=largerIndex; i++)
   {
      int i3=i*3;
      for (int j=0; j<3; j++)
      {
         boundary1.push_back(trimmingBoundary[i3+j]);
      }
   }
   // the other loop goes over end
   std::vector<double> boundary2;
   for (i=largerIndex; i<numBoundPoints; i++)
   {
      int i3=i*3;
      for (int j=0; j<3; j++)
      {
         boundary2.push_back(trimmingBoundary[i3+j]);
      }
   }
   for (i=0; i<=smallerIndex; i++)
   {
      int i3=i*3;
      for (int j=0; j<3; j++)
      {
         boundary2.push_back(trimmingBoundary[i3+j]);
      }
   }
   // now mesh them
   std::vector<double> mesh_points1;
   // if odd ground, force odd, if even, force even: 2-oddGroundLine
   success =  mesh_line_on_surface( geom_eval,
                              size_eval,
                              mesh_size,
                              boundary1,
                             /*bool periodic*/ false,
                            /* int force_evenify*/ 2-oddGroundLine, // 0 no action, 1 force odd
                                                //  2 force even
                            mesh_points1);
   if (!success)
       return success;

   std::vector<double> mesh_points2;
   success =  mesh_line_on_surface( geom_eval,
                                 size_eval,
                                 mesh_size,
                                 boundary2,
                                /*bool periodic*/ false,
                               /* int force_evenify*/ 2-oddGroundLine, // 0 no action, 1 force odd
                                                   //  2 force even
                               mesh_points2);
   if (!success)
       return success;
   // now, we have 2 faces, separated by the grounding line
   // each face has 2 edges, grounding line and a boundary edge
   // together, they form the boundary for camal paver.
   std::vector<double> bdy_coords1;
   std::vector<double> bdy_coords2; // reverse ground line
   for (i=0; i<num_seg_ground_line+1; i++)
   {
      int i3=i*3;
      int irev3= 3*(num_seg_ground_line-i);
      for (int j=0; j<3; j++)
      {
         bdy_coords1.push_back(gr_line_mesh[i3+j]);
         bdy_coords2.push_back(gr_line_mesh[irev3+j]);
      }
   }
   int sz1 = (int) mesh_points1.size()/3;
   int sz2 = (int) mesh_points2.size()/3;
   if (index1<index2)
   {
      // add the boundary 2 to loop1, and boundary 2 to loop 1 (without start and end points)
      for (i=1; i<sz2-1; i++)
      {
         int i3 = i*3;
         for (int j=0; j<3; j++)
         {
            bdy_coords1.push_back(mesh_points2[i3+j]);
         }
      }
      for (i=1; i<sz1-1; i++)
      {
         int i3 = i*3;
         for (int j=0; j<3; j++)
         {
            bdy_coords2.push_back(mesh_points1[i3+j]);
         }
      }

   }
   else
   {
      // add the boundary 2 to loop2, and boundary 1 to loop 1 (without start and end points)
      for (i=1; i<sz1-1; i++)
      {
         int i3 = i*3;
         for (int j=0; j<3; j++)
         {
            bdy_coords1.push_back(mesh_points1[i3+j]);
         }
      }
      for (i=1; i<sz2-1; i++)
      {
         int i3 = i*3;
         for (int j=0; j<3; j++)
         {
            bdy_coords2.push_back(mesh_points2[i3+j]);
         }
      }
   }


   std::cout << " boundary 1 has " << bdy_coords1.size() / 3 << " nodes\n";
   std::vector<int> loops;
   std::vector<int> loop_sizes;
   int numSegs = bdy_coords1.size() / 3;
   for (i = 0; i < numSegs; i++)
      loops.push_back(i);
   loop_sizes.push_back(numSegs);// one loop , external ....
   // only quad mesh

   CMLPaver * pave_mesher = new CMLPaver(&geom_eval, &size_eval);
      // the loops are established, get the mesh sets from each edge, and create the boundary loops arrays
      // mesh faces;

      // set only num_points_out -1 , because the last one is repeated
   success = pave_mesher->set_boundary_mesh(bdy_coords1.size() / 3,
            &bdy_coords1[0], (int) loop_sizes.size(), &loop_sizes[0], &loops[0]);
   if (!success) {
      std::cerr << "Failed setting boundary mesh" << std::endl;
      return false;
   }

   pave_mesher->set_sizing_function(CML::LINEAR_SIZING);

   // generate the mesh
   int num_quads;
   success = pave_mesher->generate_mesh(new_points, num_quads);
   if (!success) {
      std::cerr << "Failed generating mesh" << std::endl;
      //return success;
   }

   std::cout << "Meshed surface 1 with " << new_points << " new vertices and "
         << num_quads << " quadrilaterals." << std::endl;

   // get the generated mesh
   //std::vector<double> coords;
   coords.resize(3 * new_points);
   //std::vector<double> new_coords;
   //new_coords.resize(3 * new_points);
   connect.resize(4 * num_quads);
   success = pave_mesher->get_mesh(new_points, &coords[0], num_quads,
         &connect[0]);
   if (!success) {
      std::cerr << "Failed to get generated mesh" << std::endl;
      //return success;
   }
   std::cout << " boundary 2 has " << bdy_coords2.size() / 3 << " nodes\n";
   loops.clear();
   loop_sizes.clear();
   numSegs = bdy_coords2.size() / 3;
   for (i = 0; i < numSegs; i++)
      loops.push_back(i);
   loop_sizes.push_back(numSegs);// one loop , external ....
   // only quad mesh

   delete pave_mesher;
   pave_mesher= new CMLPaver(&geom_eval, &size_eval);
   success = pave_mesher->set_boundary_mesh(bdy_coords2.size() / 3,
            &bdy_coords2[0], (int) loop_sizes.size(), &loop_sizes[0], &loops[0]);
   if (!success) {
      std::cerr << "Failed setting boundary mesh" << std::endl;
      return false;
   }

   pave_mesher->set_sizing_function(CML::LINEAR_SIZING);

   // generate the mesh
   int num_quads2;
   int new_points2=0;
   success = pave_mesher->generate_mesh(new_points2, num_quads2);
   if (!success) {
      std::cerr << "Failed generating mesh" << std::endl;
      //return success;
   }

   // the first new points were the points
   std::cout << "Meshed surface 2 with " << new_points2 << " vertices and "
         << num_quads2 << " quadrilaterals." << std::endl;

   // get the generated mesh
   // now, everything will be on one surface, but it should be in 2 sets, eventually
   std::vector<double> coords2;
   coords2.resize(3 * new_points2);
   //std::vector<double> new_coords;
   //new_coords.resize(3 * new_points);
   std::vector<int> connect2;
   connect2.resize(4 * num_quads2);
   success = pave_mesher->get_mesh(new_points2, &coords2[0], num_quads2,
         &connect2[0]);
   if (!success) {
      std::cerr << "Failed to get generated mesh" << std::endl;
      //return success;
   }

   delete pave_mesher;

   etop = iMesh_QUADRILATERAL;
   // put all arrays together
   // the first nodes are repeated in the second list, but reversed
   // there are num_seg_ground_line+1 nodes on the grounding line
   //    0, 1, 2, ..., num_seg_ground_line, ..., num_points-1
   // second node array
   //  num_seg_ground_line, ..., 1, 0, [----------------------], ..., num_points2-1,
   int num_points = new_points+new_points2-num_seg_ground_line-1;
   coords.resize(3*(new_points+new_points2-num_seg_ground_line-1));
   for (i=3*(num_seg_ground_line+1); i<new_points2*3; i++)
   {
      coords[3*new_points+i-3*(num_seg_ground_line+1)]=coords2[i];
   }
   connect.resize( (num_quads+num_quads2)*4);
   // the connectivity of the second mesher is affected to the first nodes
   for (i=0; i<num_quads2*4; i++)
   {
      int indexInArr2 = connect2[i];
      int indexInBiggerArray = 0;
      if (indexInArr2<num_seg_ground_line+1)
         indexInBiggerArray = num_seg_ground_line-indexInArr2;
      else
         // else we just translate
         indexInBiggerArray = new_points-num_seg_ground_line-1+indexInArr2;
      connect[4*num_quads+i]=indexInBiggerArray;
   }


   // put new mesh back into interface
   success = cmel->create_vertices_elements(surface, bdy_verts, &coords[0],
         num_points, connect, etop, new_entities);
   return success;

}
