#include <iostream>
/*
#include <cassert>
#include <cstring>
*/
#include <string>

#include <meshkit/Tri2Quad.h>

using namespace Jaal;

int main(int argc, char **argv)
{
  Tri2Quads t2quad;

  if( argc != 3) {
      cout << "Usage: Executable infile outfile" << endl;
      return 1;
  }

  Mesh *trimesh = new Mesh;
  trimesh->readFromFile( argv[1] );

  cout << "Input: Triangle Mesh " << endl;
  cout << "# Nodes " << trimesh->getSize(0) << endl;
  cout << "# Faces " << trimesh->getSize(2) << endl;

  Mesh *quadmesh = t2quad.getQuadMesh(trimesh, 1);

  cout << "Input: Quad Mesh " << endl;
  cout << "# Nodes " << quadmesh->getSize(0) << endl;
  cout << "# Faces " << quadmesh->getSize(2) << endl;

  quadmesh->saveAs( argv[2] );

  delete quadmesh;
  delete trimesh;
}
