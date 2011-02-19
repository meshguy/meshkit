/** \file test_copymesh.cpp
 *
 * Test CopyMesh
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

#define DEFAULT_TEST_FILE "cube.cub"

int load_and_extrude(const char *input_filename,
                     const char *output_filename,
                     const Vector<3> &dx);

int main(int argc, char **argv)
{
  // check command line arg
  std::string input_filename;
  const char *output_filename = NULL;
  Vector<3> x;

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

  if (load_and_extrude(input_filename.c_str(), output_filename, x))
    return 1;

  return 0;
}

int load_and_extrude(const char *input_filename,
                     const char *output_filename,
                     const Vector<3> &dx)
{
  // start up MK and load the geometry
  MKCore mk;
  mk.load_mesh(input_filename);

  // populate mesh to relate geometry entities and mesh sets
  mk.populate_mesh();

  // get the faces
  MEntVector faces;
  mk.get_entities_by_dimension(2, faces);

  ExtrudeMesh *em = (ExtrudeMesh*) mk.construct_meshop("ExtrudeMesh", faces);
  em->set_name("extrude_mesh");

  em->set_transform(Extrude::Translate(dx, 5));

  // put them in the graph
  mk.get_graph().addArc(mk.root_node()->get_node(), em->get_node());
  mk.get_graph().addArc(em->get_node(), mk.leaf_node()->get_node());

  // mesh embedded boundary mesh, by calling execute
  mk.setup_and_execute();

  // save if o/p file specified
  if(output_filename != NULL){
    mk.save_mesh(output_filename, NULL);
  }

  std::cout << "ExtrudeMesh_test is successfully finished." << std::endl;
  return 0;
}
