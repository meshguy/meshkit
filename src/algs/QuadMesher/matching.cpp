#include <iostream>
#include <string>

#include "Tri2Quad.hpp"
#include "JaalMoabConverter.hpp"
#include "StopWatch.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

Mesh* tri_quad_conversion (Mesh *trimesh)
{
  Tri2Quads t2quad;

  cout << "Input: Triangle Mesh " << endl;
  cout << "# Nodes " << trimesh->getSize(0) << endl;
  cout << "# Faces " << trimesh->getSize(2) << endl;

  Mesh *quadmesh = t2quad.getQuadMesh(trimesh, 1);

  cout << "Input: Quad Mesh " << endl;
  cout << "# Nodes " << quadmesh->getSize(0) << endl;
  cout << "# Faces " << quadmesh->getSize(2) << endl;

  return quadmesh;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  if( argc != 3) {
      cout << "Usage: Executable infile outfile" << endl;
      return 1;
  }

  Mesh *trimesh = new Mesh;
  trimesh->readFromFile( argv[1] );

  Mesh *quadmesh = NULL;

  StopWatch watch;
  watch.start();
  quadmesh = tri_quad_conversion( trimesh );
  watch.stop();
  cout << "Info: Tri-Quad conversion time (sec) : " << watch.getSeconds() << endl;

  quadmesh->saveAs( argv[2] );

  delete quadmesh;
  delete trimesh;
}

///////////////////////////////////////////////////////////////////////////////
