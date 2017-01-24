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
#include "meshkit/AF2DfltTriangleMeshOp.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include "meshkit/ModelEnt.hpp"

// MeshKit testing
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model

#if HAVE_OCC
std::string geomExt = ".stp";
#else
std::string geomExt = ".facet";
#define HAVE_FACET
#endif

using namespace MeshKit;

int main(int argc, char **argv)
{
  // This variable is defined and used in main because a new MKCore
  // instance cannot be easily constructed after another MKCore
  // instance is deleted; there are problems with a tag left behind in
  // iGeom.
  MKCore * mk = new MeshKit::MKCore();

  std::string file_name = TestDir + "/" + "squaresurf" + geomExt;


  if (argc > 1)
  {
    file_name = argv[1];
  }
  mk->load_geometry(file_name.c_str());
  // get the surfaces
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);

  // Construct the AF2DfltTriangleMeshOp on the surfaces
  MeshOp * AFM = mk->construct_meshop("AF2DfltTriangleMeshOp", surfs);
  std::string out_file = "";
  if (argc > 2)
  {
    out_file = argv[2];
  }
  double size = 5.;
  if (argc > 3)
  {
    size = atof(argv[3]);
  }
  SizingFunction* sfPtr = new SizingFunction(mk, -1, size);
  for (unsigned int si = 0u; si < surfs.size(); ++si)
  {
    surfs[si]->sizing_function_index(sfPtr->core_index());
  }
  if (argc>4 ) // enable debugging
  {
    int debugLevel=0;
    debugLevel = atoi(argv[4]);
    AFM->set_debug_verbosity(debugLevel);
  }
  mk->setup_and_execute();

  delete sfPtr;

  if (out_file.length()>1)
  {
    // save output
     mk->moab_instance()->write_file(out_file.c_str());
  }
  delete mk;
  return 0;
}




