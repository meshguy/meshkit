#include <string>
#include <iostream>
#include <cassert>
#include <cstring>

#include "Tri2Quad.h"
#include "QuadCleanUp.h"
#include "circumcenter.h"
#include "MeshRefine2D.h"

using namespace Jaal;

int main(int argc, char **argv)
{
  Tri2Quads t2quad;

  int err;
  if( argc != 2) {
      cout << "Usage: Executable trimeshfile " << endl;
      return 1;
  }
  string fname = argv[1];

  Mesh *qm = new Mesh;
  qm->readFromFile( argv[1] );

  // Preprocessing Steps ...
  qm->saveAs( "inmesh");
  QuadCleanUp qClean(qm);
  cout << "Input Mesh: " << endl;
  cout << "      #Nodes     : " << qm->getSize(0) << endl;
  cout << "      #Triangles : " << qm->getSize(2) << endl;
  qClean.getVertexFaceDegrees();

  // Core Stuff ...
  t2quad.getQuadMesh( qm, 1);

  // Post Processing ...
  cout << "Info: Storing QuadMesh in qmesh.dat" << endl;
  cout << "Output Mesh: " << endl;
  cout << "      #Nodes   : " << qm->getSize(0) << endl;
  cout << "      #Quads   : " << qm->getSize(2) << endl;
  qClean.getVertexFaceDegrees();

  Jaal::laplacian_smoothing( qm, 100, 1);
  qClean.search_and_remove_doublets();

  qm->makeConsistentlyOriented();
// qm->reverse();
  qm->saveAs( "qmesh0");

  iMesh_Instance qmesh = 0;
  qm->toMOAB( qmesh );

  // call mesquite on ...
  qm->shapeOptimize();

  string name1 = "quad0.vtk";
  int namelen  = strlen( name1.c_str() );
  iMesh_save(qmesh, 0, name1.c_str(), NULL, &err, namelen, 0);







/*
  qClean.advancing_front_cleanup();
  qm->saveAs( "clean");

  exit(0);


  iMesh_Instance qmesh = 0;
  qm->saveAs("quad0");
  exit(0);

  qClean.cleanup_boundary();
// qClean.search_flat_quads();
   Jaal::laplacian_smoothing(qm, 500);
   qClean.remove_bridges();
   qClean.remove_diamonds(1, 1, 1);
   qClean.remove_doublets();
   qClean.remove_doublets();
   qClean.remove_bridges();
   qClean.remove_diamonds(1, 1, 1);
   qClean.remove_doublets();
// qm->search_boundary();
// qClean.search_bridges();
// qm->set_strip_markers();
// qm->getAspectRatio();
// cout << " Surface Area " << qm->getSurfaceArea() << endl;
// qm->search_boundary();
// qClean.search_yrings(); 
// qm->check_convexity();
//  qClean.search_diamonds(); 
   qm->setWavefront(2);
   qm->saveAs("quad1");
   exit(0);

   qClean.remove_doublets();
   qClean.remove_diamonds();

//  qClean.search_doublets();
//  qm->getVertexFaceDegrees();

  iMesh_dtor( trimesh, &err);
  iMesh_dtor( quadmesh, &err);

*/

}
