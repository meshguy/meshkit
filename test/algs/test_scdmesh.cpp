//---------------------------------------------------------------------------//
// \file test/algs/test_scdmesh.cpp
// \author Stuart R. Slattery
// \date Friday February 4 18:4:54 2011
// \brief Unit test for SCDMesh MeshOp
//---------------------------------------------------------------------------//

#include "meshkit/MKCore.hpp"
#include "meshkit/SCDMesh.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/VertexMesher.hpp"
#include "meshkit/EdgeMesher.hpp"
#include <vector>

using namespace MeshKit;

#include "TestUtil.hpp"

#ifdef HAVE_ACIS
#define DEFAULT_TEST_FILE_1 "brick.sat"
#define DEFAULT_TEST_FILE_2 "three_bricks.sat"
#elif defined(HAVE_OCC)
#define DEFAULT_TEST_FILE_1 "sphere.stp"
#define DEFAULT_TEST_FILE_2 "sphere.stp"
#endif

//---------------------------------------------------------------------------//
// brief function definitions and core instance declaration
MKCore *mk;

void scd_test_1();
void scd_test_2();
void scd_test_3();

//---------------------------------------------------------------------------//
// main function
int main(int argc, char **argv)
{
  mk = new MKCore();

  int num_fail = 0;

  num_fail += RUN_TEST(scd_test_1);
  num_fail += RUN_TEST(scd_test_2);
  num_fail += RUN_TEST(scd_test_3);

  delete mk;
  return num_fail;
}

//---------------------------------------------------------------------------//
// Test 1:
// *: Full mesh representation
// *: Cartesian bounding box
// *: Coarse/fine grid sizing
// *: All volumes meshed with a single grid
// *: Should produce the same mesh as test 2
void scd_test_1()
{
  // load the test geometry
  std::string scd_geom = TestDir + "/" + DEFAULT_TEST_FILE_1;
  mk->load_geometry(scd_geom.c_str());

  // get the volumes
  MEntVector vols;
  mk->get_entities_by_dimension(3, vols);

  // make an SCD mesh instance with all volumes as separate model entities
  SCDMesh *scdmesh = (SCDMesh*) mk->construct_meshop("SCDMesh", vols);

  // provide the SCD mesh parameters for a cartesian grid
  scdmesh->set_interface_scheme(SCDMesh::full);
  scdmesh->set_grid_scheme(SCDMesh::cfMesh);
  scdmesh->set_axis_scheme(SCDMesh::cartesian);
  scdmesh->set_geometry_scheme(SCDMesh::all);

  // i direction parameters
  int ci_size = 5;
  std::vector<int> fine_i(ci_size);
  fine_i[0] = 2;
  fine_i[1] = 2;
  fine_i[2] = 10;
  fine_i[3] = 2;
  fine_i[4] = 2;
  scdmesh->set_coarse_i_grid(ci_size);
  scdmesh->set_fine_i_grid(fine_i);

  // j direction parameters
  int cj_size = 5;
  std::vector<int> fine_j(cj_size);
  fine_j[0] = 2;
  fine_j[1] = 2;
  fine_j[2] = 10;
  fine_j[3] = 2;
  fine_j[4] = 2;
  scdmesh->set_coarse_j_grid(cj_size);
  scdmesh->set_fine_j_grid(fine_j);

  // k direction parameters
  int ck_size = 5;
  std::vector<int> fine_k(ck_size);
  fine_k[0] = 2;
  fine_k[1] = 2;
  fine_k[2] = 10;
  fine_k[3] = 2;
  fine_k[4] = 2;
  scdmesh->set_coarse_k_grid(ck_size);
  scdmesh->set_fine_k_grid(fine_k);

  // execute and create the structured grid
  mk->setup_and_execute();

  // write the mesh to a file
  mk->save_mesh("SCDmesh1.vtk");

  // free memory
  delete scdmesh;
}

//---------------------------------------------------------------------------//
// Test 2
// *: Light-weight ScdInterface mesh representation
// *: Cartesian bounding box
// *: Coarse/fine grid sizing
// *: All volumes meshed with a single grid
// *: Should produce the same mesh as test 1
void scd_test_2()
{
  // load the test geometry
  std::string scd_geom = TestDir + "/" + DEFAULT_TEST_FILE_1;
  mk->load_geometry(scd_geom.c_str());

  // get the volumes
  MEntVector vols;
  mk->get_entities_by_dimension(3, vols);

  // make an SCD mesh instance with all volumes as separate model entities
  SCDMesh *scdmesh = (SCDMesh*) mk->construct_meshop("SCDMesh", vols);

  // provide the SCD mesh parameters for a cartesian grid
  scdmesh->set_interface_scheme(SCDMesh::scd);
  scdmesh->set_grid_scheme(SCDMesh::cfMesh);
  scdmesh->set_axis_scheme(SCDMesh::cartesian);
  scdmesh->set_geometry_scheme(SCDMesh::all);

  // i direction parameters
  int ci_size = 5;
  std::vector<int> fine_i(ci_size);
  fine_i[0] = 2;
  fine_i[1] = 2;
  fine_i[2] = 10;
  fine_i[3] = 2;
  fine_i[4] = 2;
  scdmesh->set_coarse_i_grid(ci_size);
  scdmesh->set_fine_i_grid(fine_i);

  // j direction parameters
  int cj_size = 5;
  std::vector<int> fine_j(cj_size);
  fine_j[0] = 2;
  fine_j[1] = 2;
  fine_j[2] = 10;
  fine_j[3] = 2;
  fine_j[4] = 2;
  scdmesh->set_coarse_j_grid(cj_size);
  scdmesh->set_fine_j_grid(fine_j);

  // k direction parameters
  int ck_size = 5;
  std::vector<int> fine_k(ck_size);
  fine_k[0] = 2;
  fine_k[1] = 2;
  fine_k[2] = 10;
  fine_k[3] = 2;
  fine_k[4] = 2;
  scdmesh->set_coarse_k_grid(ck_size);
  scdmesh->set_fine_k_grid(fine_k);

  // execute and create the structured grid
  mk->setup_and_execute();

  // write the mesh to a file
  mk->save_mesh("SCDmesh2.vtk");

  // free memory
  delete scdmesh;
}

//---------------------------------------------------------------------------//
// Test 3
// *: Full mesh representation
// *: Cartesian bounding box
// *: Coarse/fine grid sizing
// *: Individual volumes meshed with their own grid
// *: Should produce the same mesh as test 1
void scd_test_3()
{
  // load the test geometry
  std::string scd_geom = TestDir + "/" + DEFAULT_TEST_FILE_2;
  mk->load_geometry(scd_geom.c_str());

  // get the volumes
  MEntVector vols;
  mk->get_entities_by_dimension(3, vols);

  // make an SCD mesh instance with all volumes as separate model entities
  SCDMesh *scdmesh = (SCDMesh*) mk->construct_meshop("SCDMesh", vols);

  // provide the SCD mesh parameters for a cartesian grid
  scdmesh->set_interface_scheme(SCDMesh::full);
  scdmesh->set_grid_scheme(SCDMesh::cfMesh);
  scdmesh->set_axis_scheme(SCDMesh::cartesian);
  scdmesh->set_geometry_scheme(SCDMesh::individual);

  // i direction parameters
  int ci_size = 5;
  std::vector<int> fine_i(ci_size);
  fine_i[0] = 2;
  fine_i[1] = 2;
  fine_i[2] = 10;
  fine_i[3] = 2;
  fine_i[4] = 2;
  scdmesh->set_coarse_i_grid(ci_size);
  scdmesh->set_fine_i_grid(fine_i);

  // j direction parameters
  int cj_size = 5;
  std::vector<int> fine_j(cj_size);
  fine_j[0] = 2;
  fine_j[1] = 2;
  fine_j[2] = 10;
  fine_j[3] = 2;
  fine_j[4] = 2;
  scdmesh->set_coarse_j_grid(cj_size);
  scdmesh->set_fine_j_grid(fine_j);

  // k direction parameters
  int ck_size = 5;
  std::vector<int> fine_k(ck_size);
  fine_k[0] = 2;
  fine_k[1] = 2;
  fine_k[2] = 10;
  fine_k[3] = 2;
  fine_k[4] = 2;
  scdmesh->set_coarse_k_grid(ck_size);
  scdmesh->set_fine_k_grid(fine_k);

  // execute and create the structured grid
  mk->setup_and_execute();

  // write the mesh to a file
  mk->save_mesh("SCDmesh3.vtk");

  // free memory
  delete scdmesh;
}


//---------------------------------------------------------------------------//
// end test_scdmesh.cpp
//---------------------------------------------------------------------------//
