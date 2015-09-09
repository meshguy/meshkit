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
#include "meshkit/CAMALTetMesher.hpp"
#include "meshkit/CopyGeom.hpp"

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

//  // now populate model ents to get the geometry created
//  mk->populate_model_ents(0, -1, -1);

//  MEntVector vols;
//  mk->get_entities_by_dimension(3, vols);

//  CopyGeom *cg = (CopyGeom*) mk->construct_meshop("CopyGeom", vols);
//  cg->set_name("copy_move_geom");

//  // set the location
//  Vector<3> dx; dx[0] = 23.5; dx[1] = 0; dx[2] = 0;
//  cg->set_location(dx);
//  cg->setup_this();
//  cg->execute_this();

//  // merge and imprint the o/p geometry
//  std::vector<iBase_EntityHandle> entities_out;
//  mk->igeom_instance()->getEntities(0, iBase_REGION, entities_out);
//  double dTol = 1.0e-2;
//  mk->igeom_instance()->mergeEnts(&entities_out[0],entities_out.size(), dTol);

//  //  mk->save_geometry("t.sat");
//  // TODO: mesh using camal and parallel mesher

//  delete cg;
  delete ag;
}



