/** \file test_solidmesher.cpp \test
 *
 * Test the SolidMesher for a basic example.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/SolidSurfaceMesher.hpp"
#include "meshkit/SolidCurveMesher.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

#ifdef HAVE_ACIS
std::string extension = ".sat";
#elif HAVE_OCC
std::string extension = ".stp";
#endif


void read_cube_tris_test();

int main(int argc, char **argv) 
{
  
  // start up MK and load the geometry
  mk = new MKCore();

  std::string filename = "cube" ;

  filename = TestDir + "/" + filename + extension;

  mk->load_geometry(&filename[0], NULL, 0, 0, 0, false, false);

  MEntVector surfs;
  mk->get_entities_by_dimension(2,surfs);
  SolidSurfaceMesher *ssm;

  ssm = (SolidSurfaceMesher*) mk->construct_meshop("SolidSurfaceMesher", surfs);

  double facet_tol = 1e-04, geom_resabs = 1e-06;
  ssm->set_mesh_params(facet_tol, geom_resabs);

  mk->setup();
  mk->execute();



  int num_fail = 0;
  num_fail += RUN_TEST(read_cube_tris_test);

#if HAVE_OCC
  return 0;
#else
  return num_fail;
#endif
}


//Tests
// NOTE: all tests should be performed using the iMesh interface as that is where our faceted
//       information lives. It can, however, be compared to the iGeom information if desired.


void read_cube_tris_test()
{

  std::vector<iMesh::EntityHandle> tris;
  mk->imesh_instance()->getEntities(0, iBase_ALL_TYPES, iMesh_TRIANGLE, tris);

  //For a cube, there should be exactly 2 triangles per face
  int num_of_tris = tris.size();
  CHECK_EQUAL(12, num_of_tris);

}

