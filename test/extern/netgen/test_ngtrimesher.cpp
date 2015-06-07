#include <cstddef>
#include <cstdio>

namespace nglib {
#include "nglib.h"
}

#include "meshkit/MKCore.hpp"
#include "meshkit/NGTriMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

#include "TestUtil.hpp"

void testSTL();
void testOCC();
void mesh_baseballfield();
void mesh_brick();
void mesh_holycyl();
void mesh_holysurf();
void mesh_test(std::string testDirFile);

MeshKit::MKCore* mk;

int main ()
{
  
  int num_fail = 0;
  
  // Test 1
  num_fail += RUN_TEST(testSTL);

  // Test 2
  num_fail += RUN_TEST(testOCC);

  mk = new MeshKit::MKCore;

  // Test 3
  num_fail += RUN_TEST(mesh_brick);
  
  // Test 4
  num_fail += RUN_TEST(mesh_holysurf);
  
  // Test 5
  num_fail += RUN_TEST(mesh_holycyl);
  
  // Test 6
  num_fail += RUN_TEST(mesh_baseballfield);
  
  return num_fail;
}

void testSTL()
{
  nglib::Ng_Init();

  std::string file_name = TestDir + "/" + "brick.stl";
  nglib::Ng_STL_Geometry* stlGeom =
    nglib::Ng_STL_LoadGeometry(file_name.c_str());
  nglib::Ng_STL_InitSTLGeometry(stlGeom);
  nglib::Ng_Meshing_Parameters* meshParams = new nglib::Ng_Meshing_Parameters();
  meshParams->maxh = 0.1;

  nglib::Ng_Mesh* netGenMesh = nglib::Ng_NewMesh();

  nglib::Ng_Result result =
    nglib::Ng_STL_MakeEdges(stlGeom, netGenMesh, meshParams);
  printf("Make Edges Result: %d\n", result);

  result =
    nglib::Ng_STL_GenerateSurfaceMesh(stlGeom, netGenMesh, meshParams);

  printf("Result: %d\n", result);
  printf("Number of Points: %d\n", nglib::Ng_GetNP(netGenMesh));
  printf("Number of Surface Elements: %d\n", nglib::Ng_GetNSE(netGenMesh));

  int lmnt[3];
  nglib::Ng_GetSurfaceElement(netGenMesh, 1, lmnt);
  printf("First Element: %d %d %d\n", lmnt[0], lmnt[1], lmnt[2]);

  nglib::Ng_DeleteMesh(netGenMesh);
  delete stlGeom;
  delete meshParams;

  nglib::Ng_Exit();
}

void testOCC() {

  nglib::Ng_Init();

  std::string file_name = TestDir + "/" + "brick.stp";
  nglib::Ng_OCC_Geometry* occGeom = nglib::Ng_OCC_Load_STEP(file_name.c_str());
  printf("Finished reading geometry\n");
  nglib::Ng_Meshing_Parameters* meshParams = new nglib::Ng_Meshing_Parameters();
  meshParams->maxh = 0.1;

  nglib::Ng_Result result;

  nglib::Ng_Mesh* occNetGenMesh = nglib::Ng_NewMesh();
  printf("Allocated mesh\n");

  result =
    nglib::Ng_OCC_GenerateEdgeMesh(occGeom, occNetGenMesh, meshParams);
  printf("OCC Generate Edge Mesh Result: %d\n", result);
  int numPnts = nglib::Ng_GetNP(occNetGenMesh);
  printf("OCC Number of Points: %d\n", numPnts);
  double pnt[3];
  for (int i = 0; i < numPnts; ++i)
  {
    nglib::Ng_GetPoint(occNetGenMesh, i + 1, pnt);
    printf("Point %d: %f %f %f\n", i + 1, pnt[0], pnt[1], pnt[2]);
  }

  result =
    nglib::Ng_OCC_GenerateSurfaceMesh(occGeom, occNetGenMesh, meshParams);

  printf("OCC Result: %d\n", result);
  printf("OCC Number of Points: %d\n", nglib::Ng_GetNP(occNetGenMesh));
  printf("OCC Number of Surface Elements: %d\n",
    nglib::Ng_GetNSE(occNetGenMesh));

  int lmnt[3];
  nglib::Ng_GetSurfaceElement(occNetGenMesh, 1, lmnt);
  printf("OCC First Element: %d %d %d\n", lmnt[0], lmnt[1], lmnt[2]);

  nglib::Ng_DeleteMesh(occNetGenMesh);
  nglib::Ng_OCC_DeleteGeometry(occGeom);
  delete meshParams;

  nglib::Ng_Exit();
}

void mesh_baseballfield()
{
  mesh_test("baseballfield.stp");
}

void mesh_brick()
{
  mesh_test("brick.stp");
}

void mesh_holysurf()
{
  mesh_test("holysurf.stp");
}

void mesh_holycyl()
{
  mesh_test("holycyl.stp");
}

void mesh_test(std::string testDirFile)
{
  std::string file_name = TestDir + "/" + testDirFile;
  mk->load_geometry(file_name.c_str());

  // get the surfaces
  MeshKit::MEntVector dum, surfaces;
  mk->get_entities_by_dimension(2, dum);
  surfaces.push_back(*dum.rbegin());
  mk->construct_meshop("NGTriMesher", surfaces);
  std::cout << surfaces.size() << " surfaces to mesh." << std::endl;
 
  // make a sizing function and set it on the surface
  MeshKit::SizingFunction esize(mk, -1, 0.25);
  surfaces[0]->sizing_function_index(esize.core_index());

  // mesh the surface, by calling execute
  mk->setup_and_execute();

  // report the number of triangles
  moab::Range triangles;
  moab::ErrorCode rval =
      mk->moab_instance()->get_entities_by_dimension(0, 2, triangles);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << triangles.size() << " triangles generated." << std::endl;

  mk->clear_graph();
  mk->delete_all();
}
