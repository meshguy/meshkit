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

extern bool debug ;

//bool debug_surf_eval = false;

int main(int argc, char *argv[]) {
	// Check command line arg
	std::string geom_filename;
	std::string out_mesh_filename;
	double mesh_size = -1.0;
	double angle = 135.; // feature decider
	int mesh_intervals = -1;
	bool force_intervals = false;

	if (argc < 4) {
		std::cout << "Usage: " << argv[0]
				<< " <input_filename> <out_mesh_filename> [-s <uniform_size>] [-i <uniform_int>] [-f] "
				<< std::endl << "  -s <uniform_size> = mesh with this size"
				<< std::endl
				<< "  -i <uniform_int> = mesh curves with this # intervals"
				<< std::endl
				<< "  -f = force these size/interval settings even if geometry has interval settings" << std::endl
				<< "  -d = print debugging info" << std::endl
				<< "  -a <angle> = feature angle decider"  << std::endl
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
#ifdef USE_CGM
			} else if (!strcmp(argv[argno], "-cgm")) {
				argno++;
				use_cgm = true;
#endif
			}else {
				std::cerr << "Unrecognized option: " << argv[argno]
						<< std::endl;
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
  iRel_createPair( relate, geom, 0, iRel_IGEOM_IFACE,
                                  mesh, 1, iRel_IMESH_IFACE,
                                  &classification, &err );
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR : can not create an assoc pair." << std::endl;
    return 1;
  }


    // read geometry
  // we use Smooth MOAB:
  char * options = "SMOOTH;";
  iGeom_load( geom, geom_filename.c_str(), options, &err, geom_filename.size(), 8 );
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR : can not load a geometry from " << geom_filename << std
::endl;
    return 1;
  }
  clock_t load_time = clock();

  // mesh the geometry
  CMEL cmel(geom, mesh, relate, classification);
  bool success = cmel.mesh_geometry(mesh_size, mesh_intervals, force_intervals, true);
  if (!success)
  std::cerr << "Problems meshing." << std::endl;

    // write the mesh
  clock_t mesh_time = clock();
  iBase_EntitySetHandle root;
  iMesh_getRootSet( mesh, &root, &err );
  assert(iBase_SUCCESS == err);
  iMesh_save(mesh, root, out_mesh_filename.c_str(), 0, &err, out_mesh_filename.length(),
 0);
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR saving mesh to " << out_mesh_filename << std::endl;
    return 1;
  }
  clock_t out_time = clock();
  std::cout << "Total time is "
       << (double) (out_time - start_time)/CLOCKS_PER_SEC
       << " s\n  load time : "
       << (double) (load_time - start_time)/CLOCKS_PER_SEC
       << " s\n  mesh time : "
       << (double) (mesh_time - load_time)/CLOCKS_PER_SEC
       << " s\n  write time : "
       << (double) (out_time - mesh_time)/CLOCKS_PER_SEC
       << std::endl;

  return !success;
}

#if 0
	// remesh the surface with paver , using smooth evaluator down deep
	CamalPaveDriver kmlPave(MBI, rootSet, MBO, angle);
	bool success = kmlPave.remesh(mesh_size, mesh_intervals, force_intervals);
	if (!success)
		std::cerr << "Problems meshing." << std::endl;
#endif

