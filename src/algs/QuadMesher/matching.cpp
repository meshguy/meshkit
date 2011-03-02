#include <iostream>
#include <string>

#include "Tri2Quad.hpp"
#include "JaalMoabConverter.hpp"

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

#ifdef HAVE_MOAB
Mesh* tri_quad_conversion (iMesh_Instance imesh)
{
  JaalMoabConverter meshconverter;
  Mesh *trimesh  = meshconverter.fromMOAB(imesh);
  Mesh* quadmesh = tri_quad_conversion ( trimesh );

  meshconverter.toMOAB(quadmesh, imesh);
  delete trimesh;
  return quadmesh;
}
#endif

///////////////////////////////////////////////////////////////////////////////

#ifdef MAIN_QUADMESH
int main(int argc, char **argv)
{
  if( argc != 3) {
      cout << "Usage: Executable infile outfile" << endl;
      return 1;
  }

  MKCore mk;
  mk.load_mesh( argv[1] );

  exit(0);
  

  Mesh *trimesh = new Mesh;
  trimesh->readFromFile( argv[1] );

  Mesh *quadmesh = NULL;

#ifdef HAVE_MOAB
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
#endif

///////////////////////////////////////////////////////////////////////////////
