/*!
\example example_SetPnt2Quad.cpp

\section example_SetPnt2Quad_cpp_title <pretty-name-of-this-file>

\subsection example_SetPnt2Quad_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection example_SetPnt2Quad_cpp_goal Goal

\subsection example_SetPnt2Quad_cpp_cw Code Walkthrough

\subsection example_SetPnt2Quad_cpp_in Input
\image html example_SetPnt2Quad.in.jpg
There is no input.

\subsection example_SetPnt2Quad_cpp_out Output
\image html example_SetPnt2Quad.out.jpg

\subsection example_SetPnt2Quad_cpp_src Source Code
*/

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

#include "meshkit/CAMALPaver.hpp"
#include "meshkit/SizingFunction.hpp"


using namespace MeshKit;


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
  options.boundary_constraint_weight = 100;
  options.height_fields = 1;
  options.create_range = 1;

  // use the same model ents as from triangle
  QslimMesher *qm = (QslimMesher*) mk->construct_meshop("QslimMesher", selection);

  qm->set_options(options);

  // mk.add_arc(tm, qm )
  //mk->get_graph().addArc(tm->get_node(), qm->get_node());
  // insert_node(GraphNode *inserted, GraphNode *before, GraphNode *after= NULL);
  mk->insert_node( qm, mk->leaf_node(), tm);

  // selection.clear();
  MBGeomOp * geo = (MBGeomOp *) mk->construct_meshop("MBGeomOp", selection);

  // maybe we should look at the moab set from the previous node in the
  // graph, to decide what set to geometrize
  // mk->get_graph().addArc(qm->get_node(), geo->get_node());
  // insert_node(GraphNode *inserted, GraphNode *before, GraphNode *after= NULL);
  mk->insert_node(geo, mk->leaf_node(), qm);

  // now add a split operation
  // we should add it to the graph
  MEntVector noEntsYet;
  MBSplitOp * splitOp = (MBSplitOp*) mk->construct_meshop("MBSplitOp", noEntsYet);

  // make sure that we split the face we want...
  //mk->get_graph().addArc(geo->get_node(), splitOp->get_node());
  //insert_node(GraphNode *inserted, GraphNode *before, GraphNode *after= NULL);
  mk->insert_node(splitOp, mk->leaf_node(), geo);
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
      dirSplit[2], /* int closed*/ 1, /*min_dot*/ 0.8);

  for (int k=0 ; k<sizePolygon; k++)
      splitOp->add_points(xyz[3*k], xyz[3*k+1], xyz[3*k+2]);

  // add an operation for CAMAL Paver mesher...
  MEntVector surfs;// empty so far
  CAMALPaver * camalPaver = (CAMALPaver*) mk->construct_meshop("CAMALPaver", surfs);
  // this will make sure that camal will be executed after split
  // if the model ents sel is empty, trigger a special case...
  // size the mesh
  double mesh_size = 50;
  /*SizingFunction *sf =*/ new SizingFunction(mk, -1, mesh_size);
  // convention: will add the latest SizingFunction available to the model ents
  // before the "setup" in execute
  /*{
    for (unsigned int i = 0; i < surfs.size(); i++)
      surfs[i]->sizing_function_index(sf->core_index());
  }*/

  //mk->get_graph().addArc(splitOp->get_node(), camalPaver->get_node());
  // insert_node(GraphNode *inserted, GraphNode *before, GraphNode *after= NULL);
  mk->insert_node(camalPaver, mk->leaf_node(), splitOp);

  mk->setup_and_execute();

  MEntSelection & mentSel = camalPaver->me_selection();
  MEntVector ments;
  for (MEntSelection::iterator sit = mentSel.begin(); sit != mentSel.end(); sit++) {
        // make a me, for convenience
      ModelEnt *me = (*sit).first;
      ments.push_back(me);
  }

  mk->save_mesh_from_model_ents("test.h5m", ments);

  return 0;
}
