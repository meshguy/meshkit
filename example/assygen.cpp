/*!
\example assygen.cpp

\section assygen_cpp_title Assembly Gen

\subsection assygen_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning Currently requires multiple setup/execute cycles.

\subsection assygen_cpp_goal Goal

\subsection assygen_cpp_cw Code Walkthrough

\subsection assygen_cpp_in Input
\image html assygen.in.jpg
There is no input.

\subsection assygen_cpp_out Output
\image html assygen.out.jpg

\subsection assygen_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/AssyGen.hpp"

#include "meshkit/CAMALTetMesher.hpp"
#include "meshkit/CopyGeom.hpp"

using namespace MeshKit;

MKCore *mk;
bool save_mesh = true;

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

  if(save_mesh)
	  ag->save_mesh("assygen.in.exo");
  ag->setup_this();
  ag->execute_this();

  // now populate model ents to get the geometry created
  mk->populate_model_ents(0, -1, -1);

  MEntVector vols;
  mk->get_entities_by_dimension(3, vols);

  CopyGeom *cg = (CopyGeom*) mk->construct_meshop("CopyGeom", vols);
  cg->set_name("copy_move_geom");

  // set the location
  Vector<3> dx; dx[0] = 23.5; dx[1] = 0; dx[2] = 0;
  cg->set_location(dx);

  cg->setup_this();
  cg->execute_this();

  // merge and imprint the o/p geometry
  std::vector<iBase_EntityHandle> entities_out;
  mk->igeom_instance()->getEntities(0, iBase_REGION, entities_out);
  double dTol = 1.0e-2;
  mk->igeom_instance()->mergeEnts(&entities_out[0],entities_out.size(), dTol);

  //  mk->save_geometry("t.sat");
  // TODO: mesh using camal and parallel mesher

  if(save_mesh)
    mk->save_mesh("assygen.out.exo");

  delete cg;
  delete ag;
}



