#include <iostream>
#include <string>

#include <meshkit/Tri2Quad.hpp>
#include <meshkit/JaalMoabConverter.hpp>

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

Mesh* tri_quad_conversion (iMesh_Instance imesh)
{
  JaalMoabConverter meshconverter;
  Mesh *trimesh  = meshconverter.fromMOAB(imesh);
  Mesh* quadmesh = tri_quad_conversion ( trimesh );

  meshconverter.toMOAB(quadmesh, imesh);
  delete trimesh;
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

#ifdef USE_MOAB
  iMesh_Instance imesh  = 0;
  JaalMoabConverter meshconverter;
  meshconverter.toMOAB(trimesh, imesh);
  quadmesh = tri_quad_conversion( imesh );
#else
   quadmesh = tri_quad_conversion( trimesh );
#endif

  delete quadmesh;
  delete trimesh;
}

///////////////////////////////////////////////////////////////////////////////
