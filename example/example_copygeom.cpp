/*!
\example example_copygeom.cpp

\section example_CopyGeom_cpp_title <pretty-name-of-this-file>

\subsection example_CopyGeom_cpp_in Input
\image html CopyGeom.in.jpg "(description of image)"

\subsection example_CopyGeom_cpp_out Output
\image html CopyGeom.out.jpg "(description of image)"

\subsection example_CopyGeom_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
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

void test_load_and_copymove();

int main(int argc, char **argv)
{
  mk = new MKCore();
  int num_fail = 0;

test_load_and_copymove();

  delete mk;
  return num_fail;
}

void test_load_and_copymove()
{
  std::string filename = TestDir + "/" + DEFAULT_TEST_FILE;
  mk->load_geometry(filename.c_str());

  MEntVector vols;
  mk->get_entities_by_dimension(3, vols);

  CopyGeom *cg = (CopyGeom*) mk->construct_meshop("CopyGeom", vols);
  cg->set_name("copy_move_geom");

// some entity tag types are always copy or expand
//  cg->expand_sets().add_tag("MATERIAL_SET");
//  cg->expand_sets().add_tag("DIRICHLET_SET");
//  cg->expand_sets().add_tag("NEUMANN_SET");

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
}



