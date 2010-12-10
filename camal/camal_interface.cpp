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
#include "CMLSurfMapper.hpp"
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

bool point_at_param(std::vector<double> & points,
                    std::vector<double> & lengsLine,
                    double param,
                    double * pp, int & index)
{
     
      // this is the par position of the next point, in natural coordinate
      int szLine = lengsLine.size();
      double * pos = std::lower_bound(&lengsLine[0], &lengsLine[szLine], param);
      if (pos == &lengsLine[0])
         return false;
      pos = pos - 1; // get the previous position, it cannot be 0
      double extraLen = param - *pos;
      index = pos - &lengsLine[0];
      int nextV = (index + 1) ; // it can't be here, but hey...
      double direction[3] = { points[3 * nextV] - points[3 * index],
            points[3 * nextV + 1] - points[3 * index + 1],
            points[3 * nextV + 2] - points[3 * index + 2] };
      NORMALIZE(direction);
      //double pp[3];
      for (int i = 0; i < 3; i++) {
         pp[i] = points[3 * index + i] + direction[i] * extraLen;
      }
      return true;
}
// this method will mesh a polyline with the given mesh count
bool mesh_line_on_surface_with_meshcount(CAMALGeomEval & geom_eval,
                          CAMALSizeEval & size_eval,
                          std::vector<double> & points,
                          int meshCount,
                          std::vector<double> & mesh_points)
{
   // the mesh size will be just "informative", but we will try to keep it proportional to the
   // meshsize evaluator. It will be tough
   // this method will be used as a precursor to camal mapper
   // mesh count is defined as number of mesh line elements in the line
   // so it is number of mesh points - 1
   int szLine = points.size()/3;
   double * point1 = &points[0];
   double * point2 = &points[3*(szLine-1)];

   if (meshCount < 2)
   {
      mesh_points.push_back(point1[0]);
      mesh_points.push_back(point1[1]);
      mesh_points.push_back(point1[2]);
      mesh_points.push_back(point2[0]);
      mesh_points.push_back(point2[1]);
      mesh_points.push_back(point2[2]);
      return true; // nothing to do; we will just keep the start and end point, as mesh points
   }
   int numPoints = points.size() / 3;

   // first, build the "parametric" space;
   // the parameter is the length along the polyline
   // there are enough points in the polyline to have a curvature;
   // the mesh size will vary linearly from start to end, and the ratio will be the same
   //
   //
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

   double startMeshSize=-1;
   size_eval.size_at_point(point1[0], point1[1], point1[2], startMeshSize, 0);
   double endMeshSize=-1;
   size_eval.size_at_point(point2[0], point2[1], point2[2], endMeshSize, 0);
   double ratio = endMeshSize/startMeshSize;
   // so the lengths will vary from 1*start, (1+..)*start, ..., ratio*start
   // (1+b*index)*start
   // so total length is (lenStart+lenEnd)*meshCount/2 = lineLength
   // so,                 lenStart*(1+ratio)*meshCount/2 = lineLength
   double lenStart = 2*lineLength/(1+ratio)/meshCount;
   // this could be very different from mesh size

   //
   std::cout<<"Mesh size at start of line: " << lenStart << ", compared to field one:" <<
         startMeshSize << "\n";
   // each mesh length is increased linearly, with the formula
   //     meshCurrentLength = lenStart*(1 + index/(meshCount-1)*(ratio-1))
   // so the if index is meshCount -1, length is 1+ratio-1 = ratio (*lenStart)
   // find the mesh point;


   // generate new points on the segments (in natural parametric space
   //                           on line)
   // use a sort of bsearch to get the index in lengsLine
   // create some vertices /edges  in the mesh, which will be used for quad generation
   // first node is the first one on line (for sure on our surface)
   // std::vector <double> bdy_coords;
   for (i = 0; i < 3; i++)
      mesh_points.push_back(point1[i]);
   // find a new point in the param space of the polyline (boundary curve)
   double param = 0; // current parameter
   double sizeLocal = lenStart;
   for (int ix=0; ix<meshCount-1; ix++)
   {
      // latest point pushed is at index ix
      sizeLocal = lenStart*(1 + ix/(meshCount-1)*(ratio-1));
      param+=sizeLocal;
      // this is the par position of the mesh point, in natural coordinate
      double pp[3];
      int index=-1;
      bool success = point_at_param( points,
                     lengsLine,
                     param,
                     pp, index);
      if (!success)
         return success;
      // now see if this is closer than what we want it to be
      geom_eval.move_to_surface(pp[0], pp[1], pp[2]);
      for (i = 0; i < 3; i++)
         mesh_points.push_back(pp[i]);
      // now compute the
   }
   // now add the last point
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

// this method is meshing a polyline using the approx size evaluator
// the mesh count is resulting

bool mesh_line_on_surface(CAMALGeomEval & geom_eval,
                          CAMALSizeEval & size_eval,
                          double mesh_size,
                          std::vector<double> & points,
                          bool periodic,
                          int force_evenify, // 0 no action, 1 force odd
                                             //  2 force even
                          std::vector<double> & mesh_points)
{

   bool success = true;
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
      double pp[3];
      // this is the par position of the next point, in natural coordinate
      int index=-1;
      success = point_at_param( points,
                     lengsLine,
                     param,
                     pp, index);
      if (!success)
         return success;
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
      double pp[3];
      int index=-1;
      success = point_at_param( points,
                     lengsLine,
                     parPosition,
                     pp, index);
      
      if (!success)
         return success;
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


void FindScalingMatrix(double * A, double * b1, double * b2, double * c1, double * c2)
{
   /*
     the equations are A * b1 = c1
                       A * b2 = c2
   these are 4 unknowns, 4 equations
   A = [a11 a12; a21 a22]
     a11*b1[0] + a12*b1[1] = c1[0]
     a11*b2[0] + a12*b2[1] = c2[0]
     delta = (b1[0]*b2[1] - b2[0]*b1[1])

second set of equations (for second row of A)
     a21*b1[0] + a22*b1[1] = c1[1]
     a21*b2[0] + a22*b2[1] = c2[1]
     delta = (b1[0]*b2[1] - b2[0]*b1[1])
    
   */
   double delta = (b1[0]*b2[1] - b2[0]*b1[1]);
   // if delta == 0, we have a problem, maybe we should treat it
   // it means origin, and b1 and b2 are collinear; unlikely but possible

   A[0]=A[1] = A[2] = A[3] = 1;
   if (delta == 0)
   {
     std::cout<<"problem with scaling \n";
     return;
   }
   A[0] =  (c1[0]*b2[1] - c2[0]*b1[1])/delta;
   A[1] = -(c1[0]*b2[0] - c2[0]*b1[0])/delta;

   A[2] =  (c1[1]*b2[1] - c2[1]*b1[1])/delta;
   A[3] = -(c1[1]*b2[0] - c2[1]*b1[0])/delta;
}
/* so here, we cut a surface described with 2 boundaries, one short, one longer, in
 * 2 other faces
 *   input:
 *      |--------------\
 *     /                \
 *     b                bt
 *     |                 |
 *     \-----------------|
 *   output:
 *      |-e1--|--------\
 *     /     /          \
 *     | 1  ie     2     \
 *     |     |           |
 *     \-e2--\-----------|
 *
 *     the resulting face 1 will have width, new internal boundary will be approx
 *       width distance from the short boundary
 *
*/
bool trim_boundary_width(CAMALGeomEval & geom_eval,
      std::vector<double> & baseLine, // (IN) the grounding line
      std::vector<double> & boundaryToTrim, // (IN) the rest of the boundary, ccw
      double width,  // (IN) the desired width
      std::vector<double> & edge1, // (OUT)
      std::vector<double> & edge2,
      std::vector<double> & boundaryRemain,
      std::vector<double> & internalEdge)
{
   // first, march along the boundary, width distance, to get e1 and e2
   // after we get the points, construct e1, e2, and boundary remain
   // most work is ie, internal edge, which is complicated
   bool success = true;
   double lineLength = 0;
   int numPoints = boundaryToTrim.size()/3;// this is the boundary to trim
   std::vector<double> lengsLine;
   lengsLine.push_back(0);// start with 0
   int i=0;
   for (i = 0; i < numPoints-1; i++) {

      double
            lenSeg =
                  DIST2( (&boundaryToTrim[3*i]), (&boundaryToTrim[3*(i+1)]) );
      lineLength += lenSeg;
      lengsLine.push_back(lineLength);
   }
   // edge 1 is at the end of boundaryRemain, and edge2 is at start
   if (2*width >= lineLength)
   {
      std::cout<<" not enough length, too wide \n";
      return false;
   }
   double pp[3];// a new point 
   double param = width; // for edge 2 end
   int index=-1;
   success = point_at_param( boundaryToTrim,
                 lengsLine,
                 param,
                 pp, index);
   if (index<0 || !success)
      return false;
   geom_eval.move_to_surface(pp[0], pp[1], pp[2]);
   edge2.clear();// just to be sure
   for (i=0; i<=index; i++)
   {
      int i3 = 3*i;
      for (int k=0; k<3; k++)
         edge2.push_back(boundaryToTrim[i3+k]);
   }

   boundaryRemain.clear(); // just to be sure
   // add also pp point
   for (i=0; i<3; i++)
   {
      edge2.push_back(pp[i]);
      boundaryRemain.push_back(pp[i]);
   }
   // we should march along the baseline, to get some more points for the internal edge
   int index2=-1;
   param = lineLength-width;
   double p1[3];
   success = point_at_param( boundaryToTrim,
                 lengsLine,
                 param,
                 p1, index2);
   geom_eval.move_to_surface(pp[0], pp[1], pp[2]);
   if (!success || index2<0)
      return false;
   for (i=index+1; i<=index2; i++)
   {
      int i3 = 3*i;
      for (int k=0; k<3; k++)
         boundaryRemain.push_back(boundaryToTrim[i3+k]);
   }
   // add pp point to end of boundary remain and to start of edge 1
   edge1.clear(); //just to be sure


   // reverse edge 1, to be oriented the same as edge 2 (from the boundary)
   for (i=numPoints-1; i>=index2+1; i--)
   {
      int i3 = 3*i;
      for (int k=0; k<3; k++)
         edge1.push_back(boundaryToTrim[i3+k]);
   }
   internalEdge.clear();
   // add the last point to edge1 and boundaryRemain
   for (i=0; i<3; i++)
   {
      edge1.push_back(p1[i]);
      boundaryRemain.push_back(p1[i]);
      internalEdge.push_back(p1[i]);
   }
   // so now the trick is to create more points on the internal edge

   // march along the baseline, and create some more points at a "width" distance
   // the most important thing is that we do not want self intersections on the surface
   // we should just project to 2d xy plane, and intersect back the surface of interest
   // we will solve a simple problem in 2d, scaling

   // we will say that in 2d, a linear transformation takes the baseline to another line
   // we assume that everything is "nice"

   // the transformation is a matrix A(2,2), applied over x, y coordinates of baseline
   // find matrix A, such that A*b = c, for b the ends of baseline, and c the pp and p1
   int numPtsInternal = baseLine.size()/3;
   double A[4];
   FindScalingMatrix (A, &(baseLine[0]), &(baseLine[numPtsInternal*3-3]), p1, pp);
   // now, we know for each point on baseline, the formula 
   // for corresponding point is A*b = c

   // compute some new z, from intersection?
   double z = (baseLine[2]+baseLine[numPtsInternal+3-1])/2;
   for (i=1; i<numPtsInternal-1; i++)
   {
      double x = baseLine[3*i];
      double y = baseLine[3*i+1];
      double x1 = A[0]*x+A[1]*y;
      double y1 = A[2]*x+A[3]*y;
      // an approx point would be on surface
      double z1 = z-2000;// just to be sure it is below surface, then the first intersection will be fine
      // maybe here we want an actual intersection with a ray piercing
      geom_eval.pierce_surface_with_ray(x1, y1, z1, 0, 0, 1);
      internalEdge.push_back(x1);
      internalEdge.push_back(y1);
      internalEdge.push_back(z1);
   }

   // the last point is pp, on internal edge
   for (i=0; i<3; i++)
   {
      internalEdge.push_back(pp[i]);
   }


   return true;// success
}
// this will create 2 loops, separated by the internalBoundary (grounding line)
// each will be meshed separately, but, of course, we start with the common line,
// the grounding line

// add edges to create a boundary in ccw fashion 
// used for mapper and for paver, too
void   add_oriented( std::vector<double> &  bm1, 
        std::vector<double> edge, int reversed)
{
   int numPointsToAdd = edge.size()/3-1; // either 
   for (int i=0; i<numPointsToAdd; i++)
   {
      int index = i*3;
      if (reversed)
         index = numPointsToAdd*3 - index;
      for (int k=0; k<3; k++)
         bm1.push_back( edge[index+k]);
   }
   return ; 
}

bool CAMAL_mesh_trimmed_surface_with_grounding_line(CMEL * cmel, iBase_EntityHandle surface,
      double mesh_size, std::vector<iBase_EntityHandle> &new_entities,
      std::vector<double> trimmingBoundary, std::vector<double> internalBoundary,
      double widthLeft, double widthRight, const bool quadMesh )
{
   // first, we will create a new area around the grounding line
   //  which will be meshed using mapped mesher from Camal.
   //
   //  create new "lines", parallel to the grounding line, on the left and right of it
   // they should be at a distance approx width from the grounding line

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

   std::vector<double> boundary2; // from index 1 towards index 2
   // this is no the left of gr line

   for (i=index1; i>-1; i++) // i>-1 always true , will break at some point
   {
      int i3=(i%numBoundPoints)*3;
      for (int j=0; j<3; j++)
      {
         boundary2.push_back(trimmingBoundary[i3+j]);
      }
      if (i%numBoundPoints == index2)
           break;
   }
   // the other loop goes from index2 to index 1
   // on the right of gr line
   std::vector<double> boundary1;
   for (i=index2; i>-1; i++) // i>-1 always true , will break at some point
   {
      int i3=(i%numBoundPoints)*3;
      for (int j=0; j<3; j++)
      {
         boundary1.push_back(trimmingBoundary[i3+j]);
      }
      if (i%numBoundPoints == index1)
           break;
   }
   // trim boundary 1 or 2, with the width, and create another polyline on the surface
   //  we assume that the grounding line is almost perpendicular to the boundaries

// boundary 1 is to the right of grounding line, boundary 2 is to the left

   // now  break the surface represented by boundary 1 and grounding line
   std::vector<double>  edge11; // (OUT)
   std::vector<double>  edge12;
   std::vector<double>  boundaryRemain1;
   std::vector<double>  internalEdge1;
   std::vector<double> & baseLine = internalBoundary; // just rename it
   success  = trim_boundary_width(geom_eval,
                              baseLine, // (IN) the grounding line
                              boundary1, // (IN) the rest of the boundary, ccw
                              widthLeft,  // (IN) the desired width
  /*  std::vector<double> & */  edge11, // (OUT)
  /*  std::vector<double> & */  edge12,
  /*  std::vector<double> & */  boundaryRemain1,
  /*  std::vector<double> & */  internalEdge1);
   if (!success)
      return success;
   // now  break the surface represented by boundary 2 and grounding line reversed
   std::vector<double> baseLineReversed;
   int baseSize= baseLine.size()/3;
   for (i=0; i<baseSize; i++)
   {
      int index = (baseSize-1 -i)*3;
      for (int k=0; k<3; k++)
         baseLineReversed.push_back(baseLine[index+k]);
   }
   std::vector<double>  edge21; // (OUT)
   std::vector<double>  edge22;
   std::vector<double>  boundaryRemain2;
   std::vector<double>  internalEdge2;
   success  = trim_boundary_width(geom_eval,
                              baseLineReversed, // (IN) the grounding line
                              boundary2, // (IN) the rest of the boundary, ccw
                              widthRight,  // (IN) the desired width
  /*  std::vector<double> & */  edge21, // (OUT)
  /*  std::vector<double> & */  edge22,
  /*  std::vector<double> & */  boundaryRemain2,
  /*  std::vector<double> & */  internalEdge2);
   if (!success)
      return success;

   // now use the mapper to mesh surfaces close to the grounding line
   CMLSurfMapper * cmlMapper = new CMLSurfMapper (&geom_eval);
  // first mesh the boundaries, edeg11, edge12, internalEdge1
   std::vector<double> mesh_points11;
   success =  mesh_line_on_surface( geom_eval,
                              size_eval,
                              mesh_size,
                              edge11,
                             /*bool periodic*/ false,
                            /* int force_evenify*/ 0, // 0 no action
                            mesh_points11);
   if (!success)
       return success;

// this method will mesh a polyline with the given mesh count
   std::vector<double>  mp_edge12;
   success =  mesh_line_on_surface_with_meshcount( geom_eval,
                           size_eval,
                           edge12, 
                          /*int meshCount*/ mesh_points11.size()/3-1,
                           mp_edge12 );
   if (!success)
       return success;

// this method will mesh a polyline with the given mesh count
   std::vector<double>  mp_Internal1;
   success =  mesh_line_on_surface_with_meshcount( geom_eval,
                           size_eval,
                           internalEdge1, 
                          /*int meshCount*/ gr_line_mesh.size()/3-1,
                           mp_Internal1 );
   if (!success)
       return success;

   // we need to set the boundary in ccw fashion, for face1:
    //!@}

    //! @name Mesh input

    //!@{
    //! \brief Supply the surface mesh describing the surface boundary.
    //!
    //! \param num_points_i The number of points that make up the first side of
    //! the four sides in the boundary.
    //! The 3rd side has the same number of points.
    //! \param num_points_j The number of points that make up the second side
    //! of the four sides in the boundary.
    //! The 4th side has the same number of points.
    //! \param points An array of points (array size = 3 * num_points_in.)<br>
    //! The first three array values are the x, y and z coordinates of
    //! the first point. The next three are for the second point, then the
    //! third, etc.  The number of points is
    //! 2*[(num_points_i-1)+(num_points_j-1)]. The points are ordered
    //! counter-clockwise around the boundary when viewed from a positive
    //! distance along the surface normal..
    //! \param point_ids An \a optional array whose values are the user's
    //! identifiers for the points (array size = num_points_in.)<br>
    //! Messages generated by the mapper use this array to translate from
    //! internal point identifiers to user identifiers.
    //!
    //! \return \a true if successful, \a false otherwise.
   std::vector<double> bm1;
   // add edges
   add_oriented(bm1, gr_line_mesh, 0); // not reversed
   add_oriented(bm1, mp_edge12, 0); // not reversed
   add_oriented(bm1, mp_Internal1, 1); //  reversed
   add_oriented(bm1, mesh_points11, 1); //  reversed

   success= cmlMapper->set_boundary_mesh(/*int num_points_i*/ gr_line_mesh.size()/3,
         /*int num_points_j*/mp_edge12.size()/3,
         /*double *const points*/(double*)(&bm1[0])
         /*, int *const point_ids = NULL*/ );
   // 
   if (!success)
      return success;

   int numPM1, numQM1;
   success =  cmlMapper->generate_mesh(numPM1, numQM1);

   std::vector<double> pm1;
   pm1.resize(3*numPM1);
   std::vector<int> connM1;
   connM1.resize(4*numQM1);
   success =  cmlMapper->get_mesh(numPM1, (double *)&pm1[0], numQM1, (int *)&connM1[0]);

   if (!success)
      return success;

   delete cmlMapper;
   // second mapped mesh
   // now use the mapper to mesh surfaces close to the grounding line
   cmlMapper = new CMLSurfMapper (&geom_eval);
  // first mesh the boundaries, edge22, edge21, internalEdge2
   /*
    * success  = trim_boundary_width(geom_eval,
                              baseLineReversed, // (IN) the grounding line
                              boundary2, // (IN) the rest of the boundary, ccw
                              widthRight,  // (IN) the desired width
     std::vector<double> &   edge21, // (OUT)
     std::vector<double> &   edge22,
     std::vector<double> &   boundaryRemain2,
     std::vector<double> &    internalEdge2);
    *
    */
   std::vector<double> mesh_points21;
   success =  mesh_line_on_surface( geom_eval,
                              size_eval,
                              mesh_size,
                              edge21,
                             /*bool periodic*/ false,
                            /* int force_evenify*/ 0, // 0 no action
                            mesh_points21);
   if (!success)
       return success;

// this method will mesh a polyline with the given mesh count
   std::vector<double>  mp_edge22;
   success =  mesh_line_on_surface_with_meshcount( geom_eval,
                           size_eval,
                           edge22,
                          /*int meshCount*/ mesh_points21.size()/3-1,
                           mp_edge22 );
   if (!success)
       return success;

// this method will mesh a polyline with the given mesh count
   std::vector<double>  mp_Internal2;
   success =  mesh_line_on_surface_with_meshcount( geom_eval,
                           size_eval,
                           internalEdge2,
                          /*int meshCount*/ gr_line_mesh.size()/3-1, // reversed has the same
                           mp_Internal2 );
   if (!success)
       return success;

   std::vector<double> bm2;
   // add edges
   add_oriented(bm2, gr_line_mesh, 1); //   reversed
   add_oriented(bm2, mp_edge22, 0); // not reversed
   add_oriented(bm2, mp_Internal2, 1); //  reversed
   add_oriented(bm2, mesh_points21, 1); //  not reversed

   success= cmlMapper->set_boundary_mesh(/*int num_points_i*/ gr_line_mesh.size()/3,
         /*int num_points_j*/mp_edge22.size()/3,
         /*double *const points*/(double*)(&bm2[0])
         /*, int *const point_ids = NULL*/ );
   //
   if (!success)
       return success;

   int numPM2, numQM2;
   success =  cmlMapper->generate_mesh(numPM2, numQM2);

   std::vector<double> pm2;
   pm2.resize(3*numPM2);
   std::vector<int> connM2;
   connM2.resize(4*numQM2);
   success =  cmlMapper->get_mesh(numPM2, (double *)&pm2[0], numQM2, (int *)&connM2[0]);

   if (!success)
      return success;

   delete cmlMapper;
   std::vector<double> mesh_points1;
   // if odd ground, force odd, if even, force even: 2-oddGroundLine
   success =  mesh_line_on_surface( geom_eval,
                              size_eval,
                              mesh_size,
                              boundaryRemain1,
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
                                 boundaryRemain2,
                                /*bool periodic*/ false,
                               /* int force_evenify*/ 2-oddGroundLine, // 0 no action, 1 force odd
                                                   //  2 force even
                               mesh_points2);
   if (!success)
       return success;
   // now, we have 2 faces left, thet will be meshed with paver, separated by the grounding line
   // mapped regions
   // each face has 2 edges, grounding line and a boundary edge
   // together, they form the boundary for camal paver.
   std::vector<double> bdy_coords1;
   std::vector<double> bdy_coords2; //

   add_oriented(bdy_coords1, mp_Internal1, 0); //   not reversed
   add_oriented(bdy_coords2, mp_Internal2, 0); //   not reversed

   add_oriented(bdy_coords1, mesh_points1, 0);
   add_oriented(bdy_coords2, mesh_points2, 0);


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
   // now, everything will be on one surface, but it should be in 4 sets, eventually
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
   /*
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

   */
   // append everything to coords and connect
   // put new mesh back into interface  numPM2, numQM2
   //coords.resize(3*(new_points+new_points2+numPM1+numPM2) );
   //connect.resize (4*(num_quads+num_quads2+numQM1+numQM2));

   int correction = new_points;
   coords.insert(coords.end(), coords2.begin(), coords2.end());
   // pm2;
   for (i=0; i<connect2.size(); i++)
      connect.push_back(correction+connect2[i]);
   //pm2.resize(3*numPM2);
   //std::vector<int> connM2;
   coords.insert(coords.end(), pm1.begin(), pm1.end());
   correction+=new_points2;
   for (i=0; i<connM1.size(); i++)
      connect.push_back(correction+connM1[i]);

   coords.insert(coords.end(), pm2.begin(), pm2.end());
   correction+=numPM1;
   for (i=0; i<connM2.size(); i++)
      connect.push_back(correction+connM2[i]);
   int num_points = new_points+new_points2+numPM1+numPM2;
   success = cmel->create_vertices_elements(surface, bdy_verts, &coords[0],
         num_points, connect, etop, new_entities);
   return success;

}
