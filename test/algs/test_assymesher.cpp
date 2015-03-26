/** \file test_assymesher.cpp \test
 *
 * Test AssyMesher
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "TestUtil.hpp"
#include "meshkit/CAMALTetMesher.hpp"
#include "meshkit/AssyMesher.hpp"
using namespace MeshKit;

MKCore *mk;

void test_assymesher_default(int argc, char **argv);

int main(int argc, char *argv[])
{
  mk = new MKCore();
  test_assymesher_default(argc, argv);
  delete mk;
  return 0;
}

void test_assymesher_default(int argc, char **argv)
{
  // Create a model entity vector for construting assymesher meshop.
  // No model entities are required for assymesher meshop, so the vector
  // remains empty.
  MEntVector volso;

  // construct the meshop and set name
  AssyMesher *am = (AssyMesher*) mk->construct_meshop("AssyMesher", volso);
  am->set_name("assymesher");

  // setup input/output assymesher files for meshing the
  // 'Reactor Assembly' geometry
  am->PrepareIO(argc, argv, TestDir);
  mk->setup_and_execute();

  //  mk->save_geometry("t.sat");
  // TODO: mesh using camal and parallel mesher

  delete am;
}



