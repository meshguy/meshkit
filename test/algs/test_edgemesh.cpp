/** \file test_edgemesh.cpp
 *
 * Test the EdgeMesher for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

int main(int argc, char **argv) 
{
    // start up MK and load the geometry
  MKCore mk;
  std::string file_name = TestDir + "/holysurf.sat";
  mk.load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs, curves, loops;
  mk.get_entities_by_dimension(2, surfs);
  CHECK_EQUAL(1, (int)surfs.size());

    // test getting loops
  std::vector<int> senses, loop_sizes;
  surfs[0]->boundary(0, loops, &senses, &loop_sizes);
  CHECK_EQUAL(4, (int)loop_sizes.size());
  
    // add up the loop sizes, should add to 8
  unsigned int l = 0;
  for (unsigned int i = 0; i < loop_sizes.size(); i++) l += loop_sizes[i];
  CHECK_EQUAL(8, (int)l);
  
    // make an edge mesher
  mk.get_entities_by_dimension(1, curves);
  EdgeMesher *em = (EdgeMesher*) mk.construct_meshop("EdgeMesher", curves);

    // make a sizing function
  SizingFunction esize(&mk, -1, 0.25);
  
    // mesh the edges, by calling execute
  mk.setup_and_execute();
}

  
