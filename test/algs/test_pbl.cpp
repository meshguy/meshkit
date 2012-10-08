/** \file test_pbl.cpp
 *
 * Test Pbl
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Pbl.hpp"
#include "TestUtil.hpp"
#include "meshkit/CAMALTetMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/CopyGeom.hpp"

using namespace MeshKit;

MKCore *mk;

void test_Pbl_default(int argc, char **argv);

int main(int argc, char *argv[])
{
  mk = new MKCore();
  test_Pbl_default(argc, argv);
  delete mk;
  return 0;
}

void test_Pbl_default(int argc, char **argv)
{
  // create a model entity vector for construting Pbl meshop, note that NO model entities are required for Pbl meshop.
  MEntVector volso;

  // construct the meshop and set name
  Pbl *ag = (Pbl*) mk->construct_meshop("Pbl", volso);
   ag->set_name("Pbl");

  // setup input/output Pbl files for creating the 'Reactor Assembly' geometry
  ag->PrepareIO(argc, argv, TestDir);
  ag->setup_this();
  ag->execute_this();

  delete ag;
}



