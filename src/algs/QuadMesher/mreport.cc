#include <string>
#include <iostream>
#include <cassert>
#include <cstring>

#include "Mesh.hpp"
#include "QuadCleanUp.hpp"

using namespace Jaal;

int main(int argc, char **argv)
{
  string infilename;
  assert( argc == 2 ) ;
 infilename = argv[1];
  
  if( infilename.empty() ) {
      cout <<"Warning: No input file specified " << endl;
      return 1 ;
  }

  Jaal::Mesh *mesh = new Jaal::Mesh;
  mesh->readFromFile( infilename );
  
  QuadCleanUp qClean( mesh );

  cout <<"# of concave faces : " << mesh->countConcaveCells() << endl;

  vector<Doublet>  doublets = qClean.search_interior_doublets();
  vector<Singlet>  singlets = qClean.search_boundary_singlets();

  qClean.getVertexFaceDegrees();

  



}
