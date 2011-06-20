/** \file test_assygen.cpp
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
  MEntVector vols;

  // construct the meshop and set name
  AssyGen *ag = (AssyGen*) mk->construct_meshop("AssyGen", vols);
  ag->set_name("assygen");

  // setup input/output assygen files for creating the 'Reactor Assembly' geometry
  ag->PrepareIO(argc, argv, TestDir);

  // put them in the graph
  mk->get_graph().addArc(mk->root_node()->get_node(), ag->get_node());
  mk->get_graph().addArc(ag->get_node(), mk->leaf_node()->get_node());

  // execute the assygen graph
  mk->setup_and_execute();

  delete ag;
}



