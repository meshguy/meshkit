#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>
#include <iomanip>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include "SimpleArray.h"
#include "ITAP_NetGen_Mesh.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////

int readGeometry(const string &filename, iGeom_Instance &geom)
{
    string engine_opt = ";engine=OCC";

    int err;
    iGeom_newGeom(engine_opt.c_str(), &geom, &err, engine_opt.length());
    iGeom_load(geom, &filename[0], 0, &err, filename.length(), 0);

    iBase_EntitySetHandle rootSet;
    iGeom_getRootSet(geom, &rootSet, &err);

    cout << "Model Contents " << endl;
    const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};

    int count;
    for (int i = 0; i <= 3; ++i)
    {
        iGeom_getNumOfType(geom, rootSet, i, &count, &err);
        std::cout << gtype[i] << count << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////

int main( int argc, char **argv)
{
  int err;
  if( argc != 2 ) {
      cout << " Usage: executable geomfile" << endl;
      return 1;
  }

  iGeom_Instance geom;
  readGeometry( argv[1], geom);

  Mesh mesh; 
  ITAP_NetGen_EdgeMesh emesher( geom, mesh);
  emesher.execute();

  exit(0);

  NetGen_SurfaceMeshing(geom, mesh);

  exit(0);

 //ITAP_NetGen_GenerateMesh( geom, bbox, mesh, 0, 1);

/*
  iMesh_Instance mesh;
  readMesh( argv[2], mesh);

  verify_surface_mesh( geom, mesh);

  iBase_EntitySetHandle rootSet;
  iMesh_getRootSet( mesh, &rootSet, &err);
 
  const char *outfile = argv[3];
  int namelen = strlen( outfile );
  char  *options = NULL;
  int    optlen  = 0;
  iMesh_save( mesh, rootSet, outfile , options, &err, namelen, optlen);
*/
}

