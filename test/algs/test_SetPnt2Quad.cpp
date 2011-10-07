#include <iostream>

#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <cstring>

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/TriangleMesher.hpp"

#include "meshkit/QslimMesher.hpp"
#include "meshkit/QslimOptions.hpp"


using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk;


int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename = 0;
  //const char *outfile = 0;

  std::string fstr=TestDir+"/pts.h5m";
  filename = fstr.c_str();
  char * opts  = (char *)("pc");
  int direction = 3; // default
  double fretting = 1.e+38;// huge, do not eliminate anything in general

  mk = new MKCore();
  mk->load_mesh(filename);
  MEntVector selection, dum;
  mk->get_entities_by_dimension(2, dum);
  selection.push_back(*dum.rbegin());// push just the last one retrieved from core

  TriangleMesher *tm = (TriangleMesher*) mk->construct_meshop("TriangleMesher", selection);

  tm->set_options(opts, direction, fretting);

  // done with triangulation, now qslim it:
  QslimOptions options;
  options.face_target = 4500;
  options.will_constrain_boundaries = true;
  options.boundary_constraint_weight = 1000;
  options.height_fields = 1;

  // use the same model ents as from triangle
  QslimMesher *qm = (QslimMesher*) mk->construct_meshop("QslimMesher", selection);

  qm->set_options(options);

  // mk.add_arc(tm, qm )
  mk->get_graph().addArc(tm->get_node(), qm->get_node());
  mk->setup_and_execute();

  mk->clear_graph();
  mk->delete_model_entities();
// clear model ents
  selection.clear();
  //
  mk->populate_model_ents(-1, 0, -1);

  mk->get_entities_by_dimension(2, selection);
  mk->construct_meshop("MBGeomOp", selection);

  mk->setup_and_execute();


  return 0;
}
