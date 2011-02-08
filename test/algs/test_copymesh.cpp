/** \file test_copymesh.cpp
 *
 * Test CopyMesh
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

#define DEFAULT_TEST_FILE "cube.cub"

int load_and_copymove(const char *input_filename,
		      const char *output_filename,
		      double *x);

int main(int argc, char **argv) 
{
  // check command line arg
  std::string input_filename;
  const char *output_filename = NULL;
  double x[3];

  if (argc ==6) {
    input_filename = argv[1];
    output_filename = argv[5];
    x[0] = atof(argv[2]);
    x[1] = atof(argv[3]);
    x[2] = atof(argv[4]);
  }

  else {
    std::cout << "Usage: " << argv[0] << "<input_mesh_filename> {x} {y} {z} <output_mesh_filename>" << std::endl;

    std::cout << std::endl;
    if (argc != 1) return 1;
    std::cout << "No file specified.  Defaulting to: " << DEFAULT_TEST_FILE << std::endl;
    std::string file_name = TestDir + "/" + DEFAULT_TEST_FILE;
    input_filename += TestDir;
    input_filename += "/";
    input_filename += DEFAULT_TEST_FILE;
    x[0] = 10.0;
    x[1] = 0.0;
    x[2] = 0.0;
  }
  
  if (load_and_copymove(input_filename.c_str(), output_filename, x))
    return 1;
  
  return 0;
}

int load_and_copymove(const char *input_filename,
		      const char *output_filename, double *x)
{
  // registration of this  mesh scheme
  // I can't register it with static register function in CopyMesh
  // should be changed.
  moab::EntityType CopyMesh_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBHEX};
  iBase_EntityType CopyMesh_mtp = iBase_REGION;
  int success = MKCore::register_meshop("CopyMesh", &CopyMesh_mtp, 1,
					CopyMesh_tps, 3, CopyMesh::factory,
					MeshOp::canmesh_region);

  // start up MK and load the geometry
  MKCore mk;
  mk.load_mesh(input_filename);

  // populate mesh to relate geometry entities and mesh sets
  mk.populate_mesh();

  // get the hexes
  MEntVector vols;
  mk.get_entities_by_dimension(3, vols);

  CopyMesh *cm = (CopyMesh*) mk.construct_meshop("CopyMesh", vols);
  cm->set_name("copy_move_mesh");

  // some entity tag types are always copy or expand
  cm->expand_sets().add_tag("MATERIAL_SET");
  cm->expand_sets().add_tag("DIRICHLET_SET");
  cm->expand_sets().add_tag("NEUMANN_SET");

  // put them in the graph
  mk.get_graph().addArc(mk.root_node()->get_node(), cm->get_node());
  mk.get_graph().addArc(cm->get_node(), mk.leaf_node()->get_node());
 
  // mesh embedded boundary mesh, by calling execute
  mk.setup_and_execute();

  // save if o/p file specified
  if(output_filename != NULL){
    mk.save_mesh(output_filename, NULL);
  }

  std::cout << "CopyMesh_test is successfully finished." << std::endl;
  return 0;
}

  
