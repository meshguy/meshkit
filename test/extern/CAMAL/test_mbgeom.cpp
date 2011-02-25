/** \file test_mbgeom.cpp
 *
 * Test the mesh based geometry in meshkit
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/FBiGeom.hpp"
#include "meshkit/CAMALTriAdvance.hpp"
#include "meshkit/CAMALPaver.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

MKCore * mk;
bool save_mesh = false;
#include "TestUtil.hpp"

void meshTri(std::string , double);
void meshQuad(std::string , double);

void meshTriShell()
{
  meshTri("shell.h5m", 1.2);
}
void meshTriPB()
{
  meshTri("PB.h5m", 150.);
}
void meshQuadShell()
{
  meshQuad("shell.h5m", 1.5);
}
void meshQuadPB()
{
  meshQuad("PB.h5m", 100.);
}

int main(int argc, char **argv) 
{
    // start up a new iGeom engine, based on moab
  /*moab::Interface * mb = new moab::Core();
  iGeom * fbiGeom = new FBiGeom(mb); // true means smooth
  MKCore mk(fbiGeom, mb); // iMesh, iRel, will be constructed*/

  if (argc ==2 ) save_mesh = true;
  mk = new MKCore;

  int num_fail = 0;

  num_fail += RUN_TEST(meshTriShell);
  //num_fail += RUN_TEST(meshTriPB);
  //num_fail += RUN_TEST(meshQuadShell);
  //num_fail += RUN_TEST(meshQuadPB);

  return num_fail;
}
void meshTri(std::string model, double size)
{
  std::string file_name = TestDir + "/" + model;
  FBiGeom * fbiGeom = new FBiGeom(mk, true); // true for smooth, false for linear

  // this will do the reading of the moab db in memory

  //this also will do heavy stuff, like smoothing
  // we do not want to do it in the constructor
  // this should also populate ModelEnts in MKCore
  fbiGeom->load(file_name.c_str());

  moab::Range tris;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 2, tris);

  int nbInitial = tris.size();
  tris.clear();
  // initial number of triangles
  std::cout << nbInitial << " initial triangles in the model" << std::endl;
  //fbiGeom->Init();

  // get the surface model ents
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);

    // create surface mesher for them
  CAMALTriAdvance *cp = (CAMALTriAdvance*) mk->construct_meshop("CAMALTriAdvance", surfs);

    // size the mesh
  SizingFunction *sf = new SizingFunction(mk, -1 , size);

    // now mesh them
  mk->setup_and_execute();

  // report the number of triangles generated
  // put it in a new set in moab

  rval = mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() - nbInitial << " triangles generated." << std::endl;

  if (save_mesh) {
        // output mesh only for the first surface
    std::string outfile = model + std::string(".vtk");
    moab::EntityHandle out_set = surfs[0]->mesh_handle();
    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL, &out_set, 1);
    MBERRCHK(rval, mk->moab_instance());
  }

  mk->clear_graph();
}

void meshQuad(std::string model, double size)
{
  std::string file_name = TestDir + "/" + model;
  FBiGeom * fbiGeom = new FBiGeom(mk, true); // true for smooth, false for linear

  // this will do the reading of the moab db in memory

  //this also will do heavy stuff, like smoothing
  // we do not want to do it in the constructor
  // this should also populate ModelEnts in MKCore
  fbiGeom->load(file_name.c_str());

  moab::Range tris;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 2, tris);

  int nbInitial = tris.size();
  tris.clear();
  // initial number of triangles
  std::cout << nbInitial << " initial triangles in the model" << std::endl;
  //fbiGeom->Init();

  // get the surface model ents
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);

    // create surface mesher for them
  CAMALPaver *tm = (CAMALPaver*) mk->construct_meshop("CAMALPaver", surfs);
    // size the mesh
  SizingFunction *sf = new SizingFunction(mk, -1 , size);

    // now mesh them
  mk->setup_and_execute();

  // report the number of triangles generated
  // put it in a new set in moab

  rval = mk->moab_instance()->get_entities_by_type(0, moab::MBQUAD, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() << " quads generated." << std::endl;

  if (save_mesh) {
        // output mesh only for the first surface
    std::string outfile = model + std::string(".vtk");
    moab::EntityHandle out_set = surfs[0]->mesh_handle();
    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL, &out_set, 1);
    MBERRCHK(rval, mk->moab_instance());
  }

  mk->clear_graph();
}

  
