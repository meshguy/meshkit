/** \file test_assygen.cpp \test
 *
 * Test AssyGen
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
  // create a model entity vector for construting assygen meshop, note that NO model entities are required for assygen meshop.
  MEntVector volso;

  // construct the meshop and set name
  AssyMesher *am = (AssyMesher*) mk->construct_meshop("AssyMesher", volso);
  am->set_name("assygen");

  // setup input/output assygen files for creating the 'Reactor Assembly' geometry
  am->PrepareIO(argc, argv, TestDir);
  mk->setup_and_execute();

  //  mk->save_geometry("t.sat");
  // TODO: mesh using camal and parallel mesher

  delete am;
}



