/** \file test_assygen.cpp \test
 *
 * Test AssyGen
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/AssyGen.hpp"
#include "TestUtil.hpp"

using namespace MeshKit;

MKCore *mk;

void test_assygen_default(int argc, char **argv);

int main(int argc, char *argv[])
{
  mk = new MKCore();
  test_assygen_default(argc, argv);
  delete mk;
  return 0;
}

void test_assygen_default(int argc, char **argv)
{
  // create a model entity vector for construting assygen meshop, note that NO model entities are required for assygen meshop.
  MEntVector volso;

  // construct the meshop and set name
  AssyGen *ag = (AssyGen*) mk->construct_meshop("AssyGen", volso);
  ag->set_name("assygen");

  // setup input/output assygen files for creating the 'Reactor Assembly' geometry
  ag->PrepareIO(argc, argv, TestDir);
  ag->setup_this();
  ag->execute_this();
}



