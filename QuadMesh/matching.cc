#include <string>
#include <iostream>
#include <cassert>
#include <cstring>

#include "Tri2Quad.h"
#include "QuadCleanUp.h"

using namespace Jaal;

int main(int argc, char **argv)
{
  int err;
  if( argc != 2) {
      cout << "Usage: Executable trimeshfile " << endl;
      return 1;
  }
  string fname = argv[1];

  iMesh_Instance trimesh;
  iMesh_newMesh(NULL, &trimesh, &err, 0);

  Jaal::readMeshData( trimesh, fname );

  Jaal::getVertexFaceDegrees(trimesh);

  Jaal::mesh_shape_optimization( trimesh );

  string name0 = "trimesh.vtk";
  int namelen0  = strlen( name0.c_str() );
  iMesh_save(trimesh, 0, name0.c_str(), NULL, &err, namelen0, 0);

  Tri2Quads t2quad;
  t2quad.getQuadMesh(trimesh, 0, 1);

  string name1 = "quad0.vtk";
  int namelen  = strlen( name1.c_str() );
  iMesh_save(trimesh, 0, name1.c_str(), NULL, &err, namelen, 0);

  string name2 = "quad1.vtk";
  namelen  = strlen( name2.c_str() );
  iMesh_save(trimesh, 0, name2.c_str(), NULL, &err, namelen, 0);

  cout << "Before Cleanup .... " << endl;
  Jaal::getVertexFaceDegrees(trimesh);

  Mesh *qm = new Mesh;
  qm->fromMOAB( trimesh );
  qm->saveAs("quad0");

  QuadCleanUp qClean(qm);
// qClean.cleanup_boundary();
// qClean.search_flat_quads();
// Jaal::laplacian_smoothing(qm, 500);

   qClean.remove_doublets();
   qClean.remove_diamonds(1, 1, 0);
   qClean.remove_doublets();

   /*
   qClean.remove_diamonds();
   qClean.remove_doublets();
   qClean.remove_doublets();
   */

// qm->set_strip_markers();
// qClean.search_bridges(); 
//   qm->setWavefront(2);
//   qm->getAspectRatio();
//   cout << " Surface Area " << qm->getSurfaceArea() << endl;

  qm->check_convexity();
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
