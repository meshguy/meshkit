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
void read_cube_curves_test();
void read_cube_surfs_test();

int count_topo_for_dim( int dim, iMesh::EntityTopology topo);

int main(int argc, char **argv) 
{
  
  // start up MK and load the geometry
  mk = new MKCore();

  std::string filename = "cube" ;

  filename = TestDir + "/" + filename + extension;

  mk->load_geometry(&filename[0]);

  MEntVector surfs;
  mk->get_entities_by_dimension(2,surfs);
  SolidSurfaceMesher *ssm;

  ssm = (SolidSurfaceMesher*) mk->construct_meshop("SolidSurfaceMesher", surfs);

  double facet_tol = 1e-04, geom_resabs = 1e-06;
  ssm->set_mesh_params(facet_tol, geom_resabs);

  mk->setup();
  mk->execute();


  //RUN TESTS
  int num_fail = 0;
  num_fail += RUN_TEST(read_cube_tris_test);
  num_fail += RUN_TEST(read_cube_curves_test);
  num_fail += RUN_TEST(read_cube_surfs_test);

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
  
  // get the number of tris in the surface model ents
  int num_of_tris = count_topo_for_dim( 2, iMesh_TRIANGLE);

  //For a cube, there should be exactly 2 triangles per face, totaling 12 for the cube.
  CHECK_EQUAL(12, num_of_tris);

}

void read_cube_curves_test()
{
  
  MEntVector ents;
  
  // there should be 12 curve (dim==1) ModelEnts for a cube
  mk->get_entities_by_dimension(1,ents);

  int num_of_curves = ents.size();

  CHECK_EQUAL(12, num_of_curves);
}

void read_cube_surfs_test()
{
  
  MEntVector ents;
  
  // there should be 6 surf (dim==1) ModelEnts for a cube
  mk->get_entities_by_dimension(2,ents);

  int num_of_surfs = ents.size();

  CHECK_EQUAL(6, num_of_surfs);
}

int count_topo_for_dim( int dim, iMesh::EntityTopology topo)
{

  MEntVector ents;
  mk->get_entities_by_dimension(dim, ents);

  int counter = 0;

  for(unsigned int i = 0; i < ents.size(); i++)
    {
      int temp = 0;
      iMesh::EntitySetHandle set_handle = IBSH(ents[i]->mesh_handle());
      mk->imesh_instance()->getNumOfTopo(set_handle, topo, temp);

      counter += temp;
      temp = 0;
    }

  return counter;
}
