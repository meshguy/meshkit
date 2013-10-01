/*!
\example example_copygeom.cpp

\section example_CopyGeom_cpp_title CopyGeom Example

\subsection example_CopyGeom_cpp_in Input
A geometry file (brick) from data folder, location to move dx.
\subsection example_CopyGeom_cpp_out Output
An example to load and move geometry file.
Geometric entities are merged and imprinted after copy/move operation using the iGeom instance.
\subsection example_CopyGeom_cpp_inf Misc. Information
\date 9-30-2013
\bug <placeholder>
\warning <placeholder>

\subsection example_CopyGeom_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/CopyGeom.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#ifdef HAVE_ACIS
#define DEFAULT_TEST_FILE "brick.sat"
#elif defined(HAVE_OCC)
#define DEFAULT_TEST_FILE "brick.stp"
#endif

MKCore *mk;

int main(int argc, char **argv)
{
  mk = new MKCore();
  std::string filename = std::string (MESH_DIR) + "/" + DEFAULT_TEST_FILE;
  mk->load_geometry(filename.c_str());

  MEntVector vols;
  mk->get_entities_by_dimension(3, vols);

  CopyGeom *cg = (CopyGeom*) mk->construct_meshop("CopyGeom", vols);
  cg->set_name("copy_move_geom");

  // set the location
  Vector<3> dx; dx[0] = 1; dx[1] = 0; dx[2] = 0;
  cg->set_location(dx);

  // put them in the graph
  mk->get_graph().addArc(mk->root_node()->get_node(), cg->get_node());
  mk->get_graph().addArc(cg->get_node(), mk->leaf_node()->get_node());

  // execute the CopyGeom graph
  mk->setup_and_execute();

  // merge and imprint the o/p geometry
  std::vector<iBase_EntityHandle> entities_out;
  mk->igeom_instance()->getEntities(0, iBase_REGION, entities_out);
  double dTol = 1.0e-2;
  mk->igeom_instance()->mergeEnts(&entities_out[0],entities_out.size(), dTol);
  mk->igeom_instance()->imprintEnts(&entities_out[0], entities_out.size());

  // delete the cg instance
  delete cg;
  return 0;
}


