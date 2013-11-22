/** \file test_tfimapping.cpp \test
 *
 * Test the TFIMapping for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/EdgeMesher.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

void test_TFImapping();
void test_TFImappingcubit();

int main(int argc, char **argv)
{

  // start up MK and load the geometry
  int num_fail = 0;

  num_fail += RUN_TEST(test_TFImapping);

  //num_fail += RUN_TEST(test_TFImappingcubit);

  return num_fail;

}

void test_TFImappingcubit()
{
  mk = new MKCore();

  std::string file_name = TestDir + "/SquareWithEdgesMeshed.cub";
  mk->load_geometry_mesh(file_name.c_str(), file_name.c_str());

  //check the number of geometrical edges
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  ModelEnt *this_surf = (*surfs.rbegin());

  this_surf->get_adjacencies(1, curves);

  CHECK_EQUAL(4, (int)curves.size());

  //check the number of mesh line segments
  moab::Range edges;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 1,
      edges);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(40, (int)edges.size());

  //now, do the TFIMapping
  TFIMapping *tm = (TFIMapping*) mk->construct_meshop("TFIMapping", surfs);
  mk->setup_and_execute();

  //check the number of quads
  moab::Range faces;
  rval = mk->moab_instance()->get_entities_by_dimension(0, 2, faces);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(100, (int)faces.size());

  mk->save_mesh("TFIMappingFromCubit.vtk");

  delete tm;

  delete mk;

}

void test_TFImapping()
{
  mk = new MKCore();


#if HAVE_OCC
  std::string file_name_geo = TestDir + "/square.stp";
  std::string file_name_edgeMeshOnly=TestDir + "/squareEdge.h5m";
  mk->load_geometry_mesh(file_name_geo.c_str(), file_name_edgeMeshOnly.c_str());
#else
  std::string file_name = TestDir + "/SquareWithOneEdgeMeshed.cub";
  mk->load_geometry_mesh(file_name.c_str(), file_name.c_str());
#endif
  // get the surface
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  CHECK_EQUAL(1, (int)surfs.size());

  // make an edge mesher
  mk->get_entities_by_dimension(1, curves);
  //test there are 4 edges bounding the surface
  CHECK_EQUAL(4, (int)curves.size());

  SizingFunction * esize = new SizingFunction(mk, 6, -1);
  surfs[0]->sizing_function_index(esize->core_index());

  /*TFIMapping *tm = (TFIMapping*) */mk->construct_meshop("TFIMapping", surfs);
  mk->setup_and_execute();

  //check whether we got the right number of quads after TFIMapping
  moab::Range faces;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 2,
      faces);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
#if HAVE_OCC
  CHECK_EQUAL(100, (int)faces.size());
#else
  CHECK_EQUAL(60, (int)faces.size());
#endif
  //output the mesh to vtk file
  mk->save_mesh("TFIMapping.vtk");
  //delete mk->vertex_mesher();
  delete mk;
}
