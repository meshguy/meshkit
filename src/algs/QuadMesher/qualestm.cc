#include <string>
#include <iostream>
#include <cassert>
#include <cstring>

#include "QuadCleanUp.hpp"
#include "circumcenter.hpp"

using namespace Jaal;

int main(int argc, char **argv)
{
  int err;
  if( argc != 2) {
      cout << "Usage: Executable Quad Mesh file " << endl;
      return 1;
  }
  string fname = argv[1];

  Mesh *qm = new Mesh;
  qm->readFromFile( argv[1] );

  QuadCleanUp qClean(qm);

  //qClean.getVertexFaceDegrees();

  cout << " Input Mesh    : " << endl;
  cout << "      #Nodes   : " << qm->getSize(0) << endl;
  cout << "      #Quads   : " << qm->getSize(2) << endl;
 // qClean.getVertexFaceDegrees();

  qClean.search_boundary_singlets();
  qClean.search_interior_doublets();
  //qClean.search_bridges();
  qClean.search_diamonds();

}
