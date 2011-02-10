//---------------------------------------------------------------------------//
// \file test/algs/test_scdmesh.cpp
// \author Stuart R. Slattery
// Friday February 4 18:4:54 2011
// \brief Unit test for SCDMesh MeshOp
//---------------------------------------------------------------------------//

#include "meshkit/MKCore.hpp"
#include "meshkit/SCDMesh.hpp"
#include "meshkit/ModelEnt.hpp"

// for now, need to #include a vertexmesher so that its statics get initialized (namely, the mesher
// gets registered)
#include "meshkit/VertexMesher.hpp"
#include "meshkit/EdgeMesher.hpp"
#include <vector>

using namespace MeshKit;

#include "TestUtil.hpp"

#ifdef HAVE_ACIS
#define DEFAULT_TEST_FILE "sphere.sat"
#elif defined(HAVE_OCC)
#define DEFAULT_TEST_FILE "sphere.stp"
#endif
//---------------------------------------------------------------------------//

int main(int argc, char **argv)
{
  // Create an instance of the MeshKit core object
  MKCore mk_scdt;
  std::string scd_geom = TestDir + "/" + DEFAULT_TEST_FILE;
  mk_scdt.load_geometry(scd_geom.c_str());

  // get the volumes
  MEntVector vols;
  mk_scdt.get_entities_by_dimension(3, vols);
  CHECK_EQUAL(1, (int)vols.size());

  // make an SCD mesh instance
  SCDMesh *scdmesh = (SCDMesh*) mk_scdt.construct_meshop("SCDMesh", vols);

  // provide the SCD mesh parameters for a cartesian grid
  // i need to update some of this to use the sizing functions
  scdmesh->set_interface_scheme(SCDMesh::full);
  scdmesh->set_grid_scheme(SCDMesh::cfMesh);
  scdmesh->set_axis_scheme(SCDMesh::cartesian);

  int ci_size = 5;
  std::vector<int> fine_i (ci_size);
  fine_i[0] = 2;
  fine_i[1] = 2;
  fine_i[2] = 10;
  fine_i[3] = 2;
  fine_i[4] = 2;
  scdmesh->set_coarse_i_grid(ci_size);
  scdmesh->set_fine_i_grid(fine_i);

  int cj_size = 5;
  std::vector<int> fine_j (cj_size);
  fine_j[0] = 2;
  fine_j[1] = 2;
  fine_j[2] = 10;
  fine_j[3] = 2;
  fine_j[4] = 2;
  scdmesh->set_coarse_j_grid(cj_size);
  scdmesh->set_fine_j_grid(fine_j);

  int ck_size = 5;
  std::vector<int> fine_k (ck_size);
  fine_k[0] = 2;
  fine_k[1] = 2;
  fine_k[2] = 10;
  fine_k[3] = 2;
  fine_k[4] = 2;
  scdmesh->set_coarse_k_grid(ck_size);
  scdmesh->set_fine_k_grid(fine_k);

  // execute and create the SCD mesh
  mk_scdt.setup_and_execute();

  // write the mesh to a file
  scdmesh->export_mesh("SCDmesh.vtk");


} // end main()

//---------------------------------------------------------------------------//
// end test_scdmesh.cpp
//---------------------------------------------------------------------------//
