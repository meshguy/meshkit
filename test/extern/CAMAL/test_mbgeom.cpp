/** \file test_mbgeom.cpp
 *
 * Test the mesh based geometry in meshkit
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/FBiGeom.hpp"
#include "meshkit/CAMALTriAdvance.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

MKCore * mk;
#include "TestUtil.hpp"

void partBed();
void shell2();

int main(int argc, char **argv) 
{
    // start up a new iGeom engine, based on moab
  /*moab::Interface * mb = new moab::Core();
  iGeom * fbiGeom = new FBiGeom(mb); // true means smooth
  MKCore mk(fbiGeom, mb); // iMesh, iRel, will be constructed*/

  mk = new MKCore;

  int num_fail = 0;

  num_fail += RUN_TEST(shell2);
  num_fail += RUN_TEST(partBed);

  return num_fail;
}
void partBed()
{
  std::string file_name = TestDir + "/" + "PB.h5m";
  FBiGeom * fbiGeom = new FBiGeom(mk, true); // true for smooth, false for linear

  // this will do the reading of the moab db in memory

  //this also will do heavy stuff, like smoothing
  // we do not want to do it in the constructor
  // this should also populate ModelEnts in MKCore
  fbiGeom->load(file_name.c_str());

  //fbiGeom->Init();

  // get the surface model ents
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);

    // create surface mesher for them
  CAMALTriAdvance *cp = (CAMALTriAdvance*) mk->construct_meshop("CAMALTriAdvance", surfs);

    // size the mesh
  SizingFunction *sf = new SizingFunction(mk, -1, 200.0);

    // now mesh them
  mk->setup_and_execute();

}

void shell2()
{
  std::string file_name = TestDir + "/" + "shell.h5m";
  FBiGeom * fbiGeom = new FBiGeom(mk, true); // true for smooth, false for linear

  // this will do the reading of the moab db in memory

  //this also will do heavy stuff, like smoothing
  // we do not want to do it in the constructor
  // this should also populate ModelEnts in MKCore
  fbiGeom->load(file_name.c_str());

  //fbiGeom->Init();

  // get the surface model ents
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);

    // create surface mesher for them
  CAMALTriAdvance *cp = (CAMALTriAdvance*) mk->construct_meshop("CAMALTriAdvance", surfs);

    // size the mesh
  SizingFunction *sf = new SizingFunction(mk, -1, 1.2);

    // now mesh them
  mk->setup_and_execute();

}

  
