#include <iostream>
#include <fstream>

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

#include "meshkit/MBGeomOp.hpp"
#include "meshkit/MBSplitOp.hpp"


using namespace MeshKit;

#include "TestUtil.hpp"
#include "ReadPolyLine.hpp"

MKCore *mk;


int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename = 0;
  //const char *outfile = 0;

  std::string fstr=TestDir+"/pts.h5m";
  filename = fstr.c_str();
  std::string polyFile = TestDir +"/polyPB.txt";
  const char * file_poly = polyFile.c_str();

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

  // selection.clear();
  MBGeomOp * geo = (MBGeomOp *) mk->construct_meshop("MBGeomOp", selection);

  // maybe we should look at the moab set from the previous node in the
  // graph, to decide what set to geometrize
  mk->get_graph().addArc(qm->get_node(), geo->get_node());

  // now add a split operation
  // we should add it to the graph
  MEntVector noEntsYet;
  MBSplitOp * splitOp = (MBSplitOp*) mk->construct_meshop("MBSplitOp", noEntsYet);

  // make sure that we split the face we want...
  mk->get_graph().addArc(geo->get_node(), splitOp->get_node());
  // how???

  std::vector<double> xyz;
  double dirSplit[3];

  int rc = ReadPolyLineFromFile(file_poly, dirSplit, xyz);
  if (rc !=0)
  {
    std::cout<<" can't read from polyline file\n";
    return rc;
  }
  int sizePolygon = (int)xyz.size()/3;
  if (sizePolygon < 3) {
    std::cerr << " Not enough points in the polygon" << std::endl;
    return 1;
  }

  // we know we will split a face with id 1; there should be only one so far
  //  !!!
  splitOp->set_options( /* int globalId*/ 1, dirSplit[0], dirSplit[1],
      dirSplit[2], /* int closed*/ 1);

  for (int k=0 ; k<sizePolygon; k++)
      splitOp->add_points(xyz[3*k], xyz[3*k+1], xyz[3*k+2]);

  mk->setup_and_execute();

/*  mk->clear_graph();
  mk->delete_model_entities();
// clear model ents
  selection.clear();
  //
  mk->populate_model_ents(-1, 0, -1);

  mk->get_entities_by_dimension(2, selection);
  mk->construct_meshop("MBGeomOp", selection);*/

  //mk->setup_and_execute();

  mk->save_mesh("oo.h5m");

  return 0;
}
