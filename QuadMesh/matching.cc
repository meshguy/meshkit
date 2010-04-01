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
  int err;
  if( argc != 2) {
      cout << "Usage: Executable trimeshfile " << endl;
      return 1;
  }
  string fname = argv[1];
  Mesh *mesh = new Mesh;
  mesh->readData( fname );

/*
  CentroidRefine2D mrefine(mesh);
  mrefine.setNumOfIterations(3);
  mrefine.execute();
  exit(0);
*/

  cout << "Edge Flipping " << endl;


  EdgeFlip edgeflip(mesh);
  edgeflip.execute();
  mesh->saveAs("ref");
  exit(0);

  cout << " Smoothing " << endl;

  laplacian_smoothing( mesh, 10 );

  mesh->saveAs("ref");


  exit(0);
  string vtkfile;

  iMesh_Instance trimesh;
  iMesh_newMesh(NULL, &trimesh, &err, 0);
  Jaal::readMeshData( trimesh, fname );

/*
  cout << " Refinement " << endl;
//Refine2D14  meshrefine;
//CentroidRefine2D  meshrefine;
//LongestEdgeRefine2D  meshrefine;
  ObtuseRefine2D  meshrefine;
  meshrefine.setBoundarySplitFlag(0);
  meshrefine.setMesh( trimesh );

  for( int i = 0; i < 5; i++) {
       ostringstream oss;
       oss << "trimesh" << i << ".vtk";
       vtkfile = oss.str();
       cout << vtkfile << endl;
       int namelen0  = strlen( vtkfile.c_str() );
       iMesh_save(trimesh, 0, vtkfile.c_str(), NULL, &err, namelen0, 0);
       meshrefine.execute();
  }
  Jaal::getVertexFaceDegrees(trimesh);
  Jaal::mesh_shape_optimization( trimesh );
  */

  Tri2Quads t2quad;
  t2quad.getQuadMesh(trimesh, 0, 1);

/*
  string name1 = "quad0.vtk";
  int namelen  = strlen( name1.c_str() );
  iMesh_save(trimesh, 0, name1.c_str(), NULL, &err, namelen, 0);
*/

//cout << "Before Cleanup .... " << endl;
//Jaal::getVertexFaceDegrees(trimesh);

  Mesh *qm = new Mesh;
  qm->fromMOAB( trimesh );

  iMesh_Instance qmesh = 0;

  QuadCleanUp qClean(qm);
// qClean.cleanup_boundary();
// qClean.search_flat_quads();
   Jaal::laplacian_smoothing(qm, 500);
   qClean.remove_doublets();
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

/*
   qClean.remove_doublets();
   qClean.remove_diamonds();

//  qClean.search_doublets();
//  qm->getVertexFaceDegrees();

  iMesh_dtor( trimesh, &err);
  iMesh_dtor( quadmesh, &err);
*/

}
