/**
 * \file smooth_test.cpp
 *
 * \brief main for camel driver for Paver, which is using smooth evaluator 
 *
 */
// go directly to MOAB
//#include "iMesh.h"
//#include "iGeom.h"

#include "MBCore.hpp"
#include "CamalPaveDriver.hpp"

#include <iostream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <string.h>

bool debug = false;

bool debug_surf_eval = false;
bool use_cgm = false;

int main(int argc, char *argv[]) {
	// Check command line arg
	std::string mesh_filename;
	std::string out_mesh_filename;
	double mesh_size = -1.0;
	double angle = 135.; // feature decider
	int mesh_intervals = -1;
	bool force_intervals = false;

	if (argc < 3) {
		std::cout << "Usage: " << argv[0]
				<< " <mesh_filename> <out_mesh_filename> [-s <uniform_size>] [-i <uniform_int>] [-f] "
				<< std::endl << "  -s <uniform_size> = mesh with this size"
				<< std::endl
				<< "  -i <uniform_int> = mesh curves with this # intervals"
				<< std::endl
				<< "  -f = force these size/interval settings even if geometry has interval settings" << std::endl
				<< "  -d = print debugging info" << std::endl
#ifdef USE_CGM
				<< " -cgm = use Cholla from CGM (default use evaluator from MOAB) " << std::endl
#endif
				<< "  -a <angle> = feature angle decider"  << std::endl
		<< std::endl;

		return 0;
	} else {
		mesh_filename = argv[1];
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
				debug_surf_eval = true;
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

	// initialize mesh interface instances
	MBErrorCode error;
	MBCore moab;
	MBInterface * MBI = &moab;

	MBCore mbOut;
	MBInterface * MBO = &mbOut;

	// read initial mesh (triangular surface of one ice sheet)

	error = MBI->load_mesh(mesh_filename.c_str(), NULL, 0);
	if (error != MB_SUCCESS)
	    return 1;

	// get the root set of the initial mesh, to pass it along; in general, it "could" be just a set

	MBEntityHandle rootSet = MBI->get_root_set(); // 0;// in MOAB this is the root set
	//assert(iBase_SUCCESS == err);
	// remesh the surface with paver , using smooth evaluator down deep
	CamalPaveDriver kmlPave(MBI, rootSet, MBO, angle);
	bool success = kmlPave.remesh(mesh_size, mesh_intervals, force_intervals);
	if (!success)
		std::cerr << "Problems meshing." << std::endl;


	MBO->write_mesh(out_mesh_filename.c_str());

	return !success;
}
