/*
 * test_AF2meshop.cpp
 *
 *  create a test/tool for meshop, to run different input files with various sizes
 *
 *  example of use
 *  test_AF2meshop <input_geo_file> <out_mesh> <size> <dbg_level>
 */

// C++
#include <cstddef>
#include <iostream>
#include <string>

// MeshKit
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include "meshkit/ModelEnt.hpp"

// MeshKit testing
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model


using namespace MeshKit;

int main(int argc, char **argv)
{
  // This variable is defined and used in main because a new MKCore
  // instance cannot be easily constructed after another MKCore
  // instance is deleted; there are problems with a tag left behind in
  // iGeom.
  MKCore * mk = new MeshKit::MKCore();

  std::string file_name, meshop, out_mesh;
  double size = 1.0;
  int debug=0;

  if (argc < 5)
  {
    std::cout <<" usage: " << argv[0] << " <geo_file> <mesh_op> <size> <out_file> [debug_level] \n ";
    return 1;
  }
  if (argc >= 5)
  {
    file_name = argv[1];
    meshop = argv[2];
    size = atof(argv[3]);
    out_mesh = argv[4];
  }

  if (argc>=6)
  {
    debug = atoi(argv[5]);
  }

  mk->load_geometry(file_name.c_str());
  // get the surfaces
  int dim_primary_dim=2; // could change for tet mesher, for example
  // this tool will work only for surf meshers ...
  MEntVector primary_ents;
  mk->get_entities_by_dimension(dim_primary_dim, primary_ents);

  // Construct the AF2DfltTriangleMeshOp on the surfaces
  MeshOp * meshOp = mk->construct_meshop(meshop.c_str(), primary_ents);

  SizingFunction* sfPtr = new SizingFunction(mk, -1, size);
  for (unsigned int si = 0u; si < primary_ents.size(); ++si)
  {
    primary_ents[si]->sizing_function_index(sfPtr->core_index());
  }

  meshOp->set_debug_verbosity(debug);

  mk->setup_and_execute();

  delete sfPtr;

  // save output
  mk->moab_instance()->write_file(out_mesh.c_str());

  delete mk;
  return 0;
}




