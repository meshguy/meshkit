#include <iostream>
#include <fstream>

using namespace std;

namespace nglib {
#include <nglib.h>
}

#include <occgeom.hpp>

using namespace nglib;
using namespace netgen;

namespace netgen 
{
extern int OCCGenerateMesh( OCCGeometry &geom, Mesh *&mesh, int statr, int end, char *optstr);
extern void WriteNeutralFormat( const Mesh &mesh, const CSGeometry &geom, const string &f);
}

///////////////////////////////////////////////////////////////////////////////

void WriteMesh( const Mesh &mesh, const string &filename)
{
/*
   CSGeometry *geom = NULL;
   WriteNeutralFormat( mesh, *geom, filename);
*/
   mesh.Save(filename);
}

///////////////////////////////////////////////////////////////////////////////

int main (int argc, char ** argv)
{
   if (argc < 3)
   {
      cerr << "Usage : executable inputfile(*.brep)  outfile(*.unv) " << endl;
      return 1;
   }

   Ng_Init();

   Mesh *mesh = new Mesh();

   int np, ne;

   OCCGeometry *geom = LoadOCC_BREP( argv[1]);
   exit(0);
   if(!geom)
   {
      cout << "Error reading in BRep File: " << argv[1] << endl;
      return 1;
   }
   cout << "Successfully loaded BRep  File: " << argv[1] << endl;

   OCCGenerateMesh( *geom, mesh, MESHCONST_ANALYSE, MESHCONST_MESHSURFACE, NULL);

   WriteMesh( *mesh, argv[2] );

   return 0;
}
