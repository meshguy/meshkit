/**
 * \file driverAF.cpp \test
 *
 * \brief Launch AF2DfltTriangleMeshOp on a given input file from the user
 *
 * use with
 *    driverAF <geo_file> <mesh_size> <outfile>
 */
// C++
#include <cstddef>
#include <iostream>
#include <string>

// MeshKit
#include "meshkit/MKCore.hpp"
#include "meshkit/AF2DfltTriangleMeshOp.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

// MeshKit testing
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model
#ifdef HAVE_ACIS
std::string geomExt = ".sat";
#elif HAVE_OCC
std::string geomExt = ".stp";
#endif

using namespace MeshKit;


// This variable is at global scope because calling deleteAll on
// the MKCore geometry instance appears to cause memory inconsistencies
// with later use of the geometry instance
MKCore* mk = NULL;

// These variables are at global scope because they affect all of
// the tests.  There is no way provided to save only some of the
// meshes
bool saveMesh = false;
std::string geofile = "";
std::string outfile = "";

int main(int argc, char **argv)
{
  // This variable is defined and used in main because a new MKCore
  // instance cannot be easily constructed after another MKCore
  // instance is deleted; there are problems with a tag left behind in
  // iGeom.
  mk = new MeshKit::MKCore();

  double msize = 1.0;

  if (argc <= 1)
  {
    std::cout << " use with driverAF <geo_file> <mesh_size> <outfile> \n" ;
    return 0;
  }
  if (argc >= 2)
  {
    geofile = argv[1];
  }
  if (argc >= 3)
    msize = atof(argv[2]);
  else
    msize = 1.0;

  if (argc>= 4)
  {
    outfile = argv[3];
    saveMesh = true;
  }

  std::cout <<" execute: " << argv[0] <<" " << geofile << " " << msize << "\n";
  mk->load_geometry(geofile.c_str());
  MEntVector newSurfs;
  mk->get_entities_by_dimension(2, newSurfs);

  // Construct the AF2DfltTriangleMeshOp on the surfaces
  mk->construct_meshop("AF2DfltTriangleMeshOp", newSurfs);
  SizingFunction* sfPtr = new SizingFunction(mk, -1, msize); // support only fixed mesh size

  for (unsigned int si = 0u; si < newSurfs.size(); ++si)
  {
    newSurfs[si]->sizing_function_index(sfPtr->core_index());
  }
  mk->setup_and_execute();

  // report the number of triangles
  moab::Range tris;
  moab::ErrorCode rval =
     mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() << " tris generated." << std::endl;
  // remove the MeshOp and the rest of the graph
  if (saveMesh)
  {
    std::cout << "save mesh file " << outfile << "\n";
    rval = mk->moab_instance()->write_file(outfile.c_str());
  }
  mk->clear_graph();

  delete mk;

}

