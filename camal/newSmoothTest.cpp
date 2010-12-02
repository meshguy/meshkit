/**
 * \file smooth_test.cpp
 *
 * \brief main for camel driver for Paver, which is using smooth evaluator 
 *
 */
#include "iMesh.h"
#include "iGeom.h"
#include "iRel.h"

#include "camel.hpp"

#include <iostream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

extern bool debug;

//bool debug_surf_eval = false;

int main(int argc, char *argv[]) {
   // Check command line arg
   std::string geom_filename;
   std::string out_mesh_filename;
   std::string polygon_filename;
   std::string grounding_line_filename;
   double mesh_size = -1.0;
   double angle = 135.; // feature decider
   int mesh_intervals = -1;
   bool force_intervals = false;
   std::cout<<"Command line:\n";
   for (int argi=0; argi<argc; argi++)
      std::cout<< argv[argi] << " ";
   std::cout<<"\n";

   double width = -1; // grounding line

   if (argc < 4) {
      std::cout << "Usage: " << argv[0]
            << " <input_filename> <out_mesh_filename> [-s <uniform_size>] [-i <uniform_int>] [-f] "
            << std::endl << "  -s <uniform_size> = mesh with this size"
            << std::endl
            << "  -i <uniform_int> = mesh curves with this # intervals"
            << std::endl
            << "  -f = force these size/interval settings even if geometry has interval settings"
            << std::endl << "  -d = print debugging info" << std::endl
            << "  -a <angle> = feature angle decider" << std::endl
            << "  -t <trim_polygon_xy_file> = file with 2d trimming polygon, oriented ccw " << std::endl 
            << "  -g <grounding_line_file> = file with 2d line \n"
            << " -w <length> = width of the corridor along the grounding line \n"

      << std::endl;

      return 0;
   } else {
      geom_filename = argv[1];
      out_mesh_filename = argv[2];
      int argno = 3;
      while (argno < argc) {
         if (!strcmp(argv[argno], "-s")) {
            argno++;
            sscanf(argv[argno], "%lf", &mesh_size);
            argno++;
         } else if (!strcmp(argv[argno], "-a")) {
            argno++;
            sscanf(argv[argno], "%lf", &angle);
            argno++;
         } else if (!strcmp(argv[argno], "-i")) {
            argno++;
            sscanf(argv[argno], "%d", &mesh_intervals);
            argno++;
         } else if (!strcmp(argv[argno], "-f")) {
            argno++;
            force_intervals = true;
         } else if (!strcmp(argv[argno], "-d")) {
            argno++;
            debug = true;
         } else if (!strcmp(argv[argno], "-t")) {
            argno++;
            polygon_filename = argv[argno];
            argno++;
         } else if (!strcmp(argv[argno], "-g")) {
            argno++;
            grounding_line_filename = argv[argno];
            argno++;

         } else if (!strcmp(argv[argno], "-w")) {
            argno++;
            sscanf(argv[argno], "%lf", &width);
            argno++;

         }else {
            std::cerr << "Unrecognized option: " << argv[argno] << std::endl;
            return 1;
         }
      }
   }
   clock_t start_time = clock();
   int err = 0;
   // read initial mesh (triangular surface of one ice sheet)
   iGeom_Instance geom;
   iGeom_newGeom(0, &geom, &err, 0);
   assert(iBase_SUCCESS == err);

   iMesh_Instance mesh;
   iMesh_newMesh(0, &mesh, &err, 0);
   assert(iBase_SUCCESS == err);

   iRel_Instance relate;
   iRel_newRel(0, &relate, &err, 0);
   assert(iBase_SUCCESS == err);

   // create an association pair
   iRel_PairHandle classification;
   iRel_createPair(relate, geom, 0, iRel_IGEOM_IFACE, mesh, 1,
         iRel_IMESH_IFACE, &classification, &err);
   if (iBase_SUCCESS != err) {
      std::cerr << "ERROR : can not create an assoc pair." << std::endl;
      return 1;
   }

   // read geometry
   // we use Smooth MOAB:
   char * options = "SMOOTH;";
   iGeom_load(geom, geom_filename.c_str(), options, &err, geom_filename.size(),
         8);
   if (iBase_SUCCESS != err) {
      std::cerr << "ERROR : can not load a geometry from " << geom_filename
            << std::endl;
      return 1;
   }
   clock_t load_time = clock();

   // mesh the geometry
   CMEL cmel(geom, mesh, relate, classification);
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
      iGeom_getEntities(geom, 0, 2, &these_gents, &these_gents_alloc,
            &these_gents_size, &err);
      if (iBase_SUCCESS != err) {
         std::cerr << "Trouble getting gentities of dimension 2 \n";
      }
      if (these_gents_size != 1)
         std::cerr << " Only one surface should be trimmed/meshed at a time \n";
      success = cmel.trimSurface(polygon_filename.c_str(),
            polygon_filename.size());
      if (!success)
         std::cerr << "Problems trimming the surface." << std::endl;
      if (!grounding_line_filename.empty())
      {
         // we will try to respect the grounding line, as an interior 
         // line to the trimmed surface; this will add another boundary,
         // interior to the trimmed surface
         success = cmel.grounding_line(grounding_line_filename.c_str(),
            grounding_line_filename.size(), width);
         if (!success)
            std::cerr<< " problems adding new interior boundaries for grounding line \n";
      }
   }
   success = cmel.mesh_geometry(mesh_size, mesh_intervals, force_intervals,
         true);
   if (!success)
      std::cerr << "Problems meshing." << std::endl;

   // write the mesh
   clock_t mesh_time = clock();
   iBase_EntitySetHandle root;
   iMesh_getRootSet(mesh, &root, &err);
   assert(iBase_SUCCESS == err);
   iMesh_save(mesh, root, out_mesh_filename.c_str(), 0, &err,
         out_mesh_filename.length(), 0);
   if (iBase_SUCCESS != err) {
      std::cerr << "ERROR saving mesh to " << out_mesh_filename << std::endl;
      return 1;
   }
   clock_t out_time = clock();
   std::cout << "Total time is " << (double) (out_time - start_time)
         / CLOCKS_PER_SEC << " s\n  load time : " << (double) (load_time
         - start_time) / CLOCKS_PER_SEC << " s\n  mesh time : "
         << (double) (mesh_time - load_time) / CLOCKS_PER_SEC
         << " s\n  write time : " << (double) (out_time - mesh_time)
         / CLOCKS_PER_SEC << std::endl;

   return !success;
}

#if 0
// remesh the surface with paver , using smooth evaluator down deep
CamalPaveDriver kmlPave(MBI, rootSet, MBO, angle);
bool success = kmlPave.remesh(mesh_size, mesh_intervals, force_intervals);
if (!success)
std::cerr << "Problems meshing." << std::endl;
#endif

