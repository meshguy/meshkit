/**
 * \file ComputeThickness.cpp
 *
 * \brief from 2 surfaces, top and bottom, that occupy a different xy domain,
 * compute a thickness between those 2 sheets
 *
 *  typical scenario: take a simplified top and simplified bed, load in iGeom
 *   augment with the thickness tag; this could be eventually used in a mesh size for camal
 *   as an indication (how many layers do we really need?)
 *
 */

#include "iMesh.h"
#include "iGeom.h"
#include "iRel.h"

#include "stdlib.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "moab/CartVect.hpp"
#include <vector>

#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

//bool debug_surf_eval = false;

bool ptIsOutside (double * pos, std::vector<double> & polCoords)
{
   int sizePolygon = polCoords.size()/3;
   if (sizePolygon < 1)
      return false;
   double ang = 0;
   for (int i = 0; i<sizePolygon; i++)
   {
      moab::CartVect v1(polCoords[3*i]-pos[0], polCoords[3*i+1]-pos[1], 0);
      int nexti = (i+1)%sizePolygon;
      moab::CartVect v2(polCoords[3*nexti]-pos[0], polCoords[3*nexti+1]-pos[1], 0);
      if (v1.length_squared()>0. && v2.length_squared()>0)
      {
         double ang1 = angle(v1,v2);
         moab::CartVect cr=v1*v2;
         if (cr[2] > 0)
           ang +=ang1;
         if (cr[2] < 0)
           ang -=ang1;
         if ( cr[2] ==0)
            return false; // actually on boundary
         
      }
   }
   if (ang >3) // PI
      return false;

   return true;
}

int main(int argc, char *argv[]) {
   // Check command line arg
   std::string bottom_filename;
   std::string top_filename;
   std::string out_mesh;
   std::string polygon_filename;
   bool smooth = true;

   if (argc < 4) {
      std::cout << "Usage: " << argv[0]
            << " <bottom_filename> <top_filename> <out_mesh> [-t <trimming polygon>]  [-n] "
            << std::endl
            << "  -t <trim_polygon_xy_file> = file with 2d trimming polygon, oriented ccw "
            << std::endl
            << "  -n regular facetting for top surface (default smooth faceting) \n"

      << std::endl;

      return 0;
   } else {
      bottom_filename = argv[1];
      top_filename = argv[2];
      out_mesh = argv[3];
      int argno = 4;
      while (argno < argc) {
         if (!strcmp(argv[argno], "-t")) {
            argno++;
            polygon_filename = argv[argno];
            argno++;
         } else if (!strcmp(argv[argno], "-n")) {
            argno++;
            smooth = false;
         } else {
            std::cerr << "Unrecognized option: " << argv[argno] << std::endl;
            return 1;
         }
      }
   }
   clock_t start_time = clock();
   int err = 0;
   // read initial mesh (triangular surface bottom)
   iMesh_Instance mesh1;
   iMesh_newMesh(0, &mesh1, &err, 0);
   assert(iBase_SUCCESS == err);

   iGeom_Instance geom2;
   iGeom_newGeom(0, &geom2, &err, 0);
   assert(iBase_SUCCESS == err);

   // read bottom mesh
   // we use Smooth MOAB:
   char * options = NULL; // "SMOOTH;";
   iBase_EntitySetHandle root_set;
   iMesh_getRootSet(mesh1, &root_set, &err);
   ERRORR("Couldn't get root set.", 1);

   // read  mesh

   iMesh_load(mesh1, root_set, bottom_filename.c_str(), options, &err,
         bottom_filename.size(), 0);
   if (iBase_SUCCESS != err) {
      std::cerr << "ERROR : can not load a mesh from " << bottom_filename
            << std::endl;
      return 1;
   }
   clock_t load_time1 = clock();

   char * opts2 = "SMOOTH;";
   if (smooth)
     iGeom_load(geom2, top_filename.c_str(), opts2, &err, top_filename.size(), 8);
   else
     iGeom_load(geom2, top_filename.c_str(), 0, &err, top_filename.size(), 0);
      
   if (iBase_SUCCESS != err) {
      std::cerr << "ERROR : can not load a geometry from " << top_filename
            << std::endl;
      return 1;
   }

   clock_t load_time2 = clock();// load the smooth file

   std::vector<double> polCoords;
   double direction[3];// normalized
   bool success;
   // first trim the first surface using z direction
   // this will be used by camel to get a new boundary / loop for the first surface
   // it will create a new boundary for surface 1
   if (!polygon_filename.empty()) {
      // check number of surfaces in the model
      // if more than 1, this is not really supported
      iBase_EntityHandle *these_gents = NULL;
      int these_gents_size = 0, these_gents_alloc = 0;
      // get all entities of this dimension
      iGeom_getEntities(geom2, 0, 2, &these_gents, &these_gents_alloc,
            &these_gents_size, &err);
      if (iBase_SUCCESS != err) {
         std::cerr << "Trouble getting gentities of dimension 2 \n";
      }
      if (these_gents_size != 1)
         std::cerr << " Only one surface should be trimmed/meshed at a time \n";

      // read the file with the polygon user data
      std::ifstream datafile(polygon_filename.c_str(), std::ifstream::in);
      if (!datafile) {
         std::cout << "can't read file\n";
         return 1;
      }
      //
      char temp[100];

      double gridSize;
      datafile.getline(temp, 100);// first line

      // get direction and mesh size along polygon segments, from file
      sscanf(temp, " %lf %lf %lf %lf ", direction, direction + 1,
            direction + 2, &gridSize);
      //NORMALIZE(direction);// just to be sure

      //std::vector<double> xs, ys, zs;
      while (!datafile.eof()) {
         datafile.getline(temp, 100);
         //int id = 0;
         double x, y, z;
         int nr = sscanf(temp, "%lf %lf %lf", &x, &y, &z);
         if (nr == 3) {
            polCoords.push_back(x);
            polCoords.push_back(y);
            polCoords.push_back(z);
         }
      }
      int sizePolygon = polCoords.size() / 3;
      if (sizePolygon < 3) {
         std::cerr << " Not enough points in the polygon" << std::endl;
         return 1;
      }
      /*success = cmel.trimSurface(polygon_filename.c_str(),
       polygon_filename.size());
       if (!success)
       std::cerr << "Problems trimming the surface." << std::endl;*/
   }
   /* success = cmel.mesh_geometry(mesh_size, mesh_intervals, force_intervals,
    true);
    if (!success)
    std::cerr << "Problems meshing." << std::endl;*/

   // add another tag, with the thickness, on the
   // write the mesh
   // get all the vertices, and shoot rays from their position -1000, then compute thickness
   // first get all nodes from bottom mesh
   iBase_EntityHandle *verts = NULL;
   int verts_alloc = 0;
   int numNodes = 0;
   iMesh_getEntities(mesh1, root_set, iBase_VERTEX, iMesh_POINT, &verts,
         &verts_alloc, &numNodes, &err);
   ERRORR("failed to get vertices.", 1);

   /* get the coordinates in one array */

   int vert_coords_alloc = 0;
   double * xyz = 0; // not allocated
   int vertex_coord_size = 0;

   iMesh_getVtxArrCoords(mesh1, verts, numNodes, iBase_INTERLEAVED, &xyz,
         &vert_coords_alloc, &vertex_coord_size, &err);
   ERRORR("failed to get vertex coordinates of entities in getMeshData.", 1);

   // then, go to shoot rays
   iBase_TagHandle elev_tag_handle;
   const char * tagName1 = "Thickness";
   iMesh_createTag(mesh1, tagName1,
   /*  size ? */1, iBase_DOUBLE, &elev_tag_handle, &err, strlen(tagName1));
   ERRORR("failed to create tag.", 1);
   double * dArr = new double[numNodes];
   // then, go to shoot rays
   iBase_TagHandle situation_tag_handle;
   const char * tagName2 = "Situation"; // < 0 not intersected 
   // 1 less than 200
   // 2 more than 200, but floating 
   // floating if ztop/thickness < 1 - roIce/roWater
   // 3 more than 200, ground
   iMesh_createTag(mesh1, tagName2,
   /*  size ? */1, iBase_DOUBLE, &situation_tag_handle, &err, strlen(tagName2));
   ERRORR("failed to create second tag.", 1);
   double * dArr2 = new double[numNodes];
   int j = 0;
   int numRaysIntersected =0 ;
   // double factorFloating = (1.-937./1026.);
   for (j = 0; j < numNodes; j++) {
      //dArr[j] = xyz[j * 3 + 2];
      dArr[j] = 0; // no intersection situation, marked by 0
      dArr2[j] = 0; // not computed
      // for a point, see if it is inside the polygon, with winding number

      double pos[3] = {xyz[j*3], xyz[3*j+1], xyz[3*j+2]-1000}; // subtract 1000 from shoot ray pt
      if (ptIsOutside (pos, polCoords))
         continue;
      iBase_EntityHandle * intersect_entity_handles = NULL;
      int intersect_entity_handles_allocated = 0, intersect_entity_handles_size = 0;
      double * intersect_coords = NULL;
      int intersect_coords_allocated =0 ,  intersect_coords_size = 0;
       double * param_coords = NULL;
      int param_coords_allocated = 0, param_coords_size =0;
      iGeom_getPntRayIntsct( geom2,
            pos[0], pos[1], pos[2],
            0, 0, 1,
            &intersect_entity_handles, &intersect_entity_handles_allocated,
            &intersect_entity_handles_size, iBase_INTERLEAVED,
            &intersect_coords, &intersect_coords_allocated, &intersect_coords_size,
            &param_coords, &param_coords_allocated, &param_coords_size,
            &err );
      // get the first coordinate
      if (err != 0 || intersect_entity_handles_size ==0 || param_coords == NULL)
         continue;
      numRaysIntersected++;
      // consider only the first intersection point
      dArr[j] = param_coords[0] - 1000; // the first intersection only
      // this is the thickness
      // intersect point has z = intersect_coords[0]
      if (dArr[j] < 200 ) 
         dArr2[j] = 1; // computed but less than 200
      else
      {
         dArr2[j] = 3; // ground
         // decide if it is floating  
         double z = xyz[3*j+2];
         if (z<0) // below sea level, could float
         {
            double floating = dArr[j]*910+z*1026;
            if (floating < 0)
               dArr2[j] = 2;
         }
      }
         
      free(intersect_entity_handles);
      free(intersect_coords);
      free(param_coords);

   }

   iMesh_setDblArrData(mesh1,
          /*in const iBase_EntityHandle* */ verts,
          /*in const int */  numNodes,
          /*in const iBase_TagHandle*/ elev_tag_handle,
          /*in const double* */  dArr,
          /*in const int */ numNodes,
          /*out int * */ &err);

   ERRORR("failed to set tag values.", 1);

   iMesh_setDblArrData(mesh1,
          /*in const iBase_EntityHandle* */ verts,
          /*in const int */  numNodes,
          /*in const iBase_TagHandle*/ situation_tag_handle,
          /*in const double* */  dArr2,
          /*in const int */ numNodes,
          /*out int * */ &err);

   ERRORR("failed to set tag2 values.", 1);
   clock_t compute_time = clock();

   assert(iBase_SUCCESS == err);
   iMesh_save(mesh1, root_set, out_mesh.c_str(), 0, &err, out_mesh.length(), 0);
   if (iBase_SUCCESS != err) {
      std::cerr << "ERROR saving mesh to " << out_mesh << std::endl;
      return 1;
   }
   clock_t out_time = clock();
   std::cout << "Total time is " << (double) (out_time - start_time)
         / CLOCKS_PER_SEC << " s\n  load bottom : " << (double) (load_time1
         - start_time) / CLOCKS_PER_SEC << " s\n  load top : "
         << (double) (load_time2 - load_time1) / CLOCKS_PER_SEC
         << " s\n  compute time : " << (double) (compute_time - load_time2)
         / CLOCKS_PER_SEC << " s\n  write time : " << (double) (out_time
         - compute_time) / CLOCKS_PER_SEC << std::endl;
   std::cout << "num rays intersected : " << numRaysIntersected << "\n";  

   free(xyz);
   delete [] dArr;
   delete [] dArr2;
   free (verts);
   return !success;
}

