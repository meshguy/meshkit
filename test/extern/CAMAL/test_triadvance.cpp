/** \file test_triadvance.cpp
 *
 * Test the CAMALTriAdvance for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/CAMALTriAdvance.hpp"
// for now, need to #include a vertexmesher so that its statics get initialized (namely, the mesher
// gets registered)
#include "meshkit/VertexMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

typedef bool (*fptr)(ModelEnt*);

int main(int argc, char **argv) 
{
  
    // start up MK and load the geometry
  MKCore mk;
  std::string file_name = TestDir + "/holysurf.sat";
  if (2 <= argc) file_name = argv[1];
  mk.load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs;
  mk.get_entities_by_dimension(2, surfs);
  CAMALTriAdvance *tm = (CAMALTriAdvance*) mk.construct_meshop("CAMALTriAdvance", surfs);

    // make a sizing function and set it on the surface
  SizingFunction esize(&mk, -1, 0.5);
  surfs[0]->sizing_function_index(esize.core_index());
  
    // mesh the surface, by calling execute
  mk.setup_and_execute();

    // report the number of tris
  moab::Range tris;
  moab::ErrorCode rval = mk.moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() << " tris generated." << std::endl;

  if (3 <= argc) {
    try {
        // output mesh
      mk.save_mesh(argv[2]);
    }
    catch (Error err) {
      std::cout << "Error occurred: " << err.what() << std::endl;
    }
  }
}

  
