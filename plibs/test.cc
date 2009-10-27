#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

using namespace std;

//#include "SurfMesher.h"
//#include "VolMesher.h"

//int  example_volume_mesher( const string &);

int main()
{
    vector<int> a;

    for( int i = 0; i < 10; i++) 
         a.push_back( i );
    reverse( a.begin(), a.end() );

    for( int i = 0; i < 10; i++) 
         cout << a[i] << endl;

/*
    example_volume_mesher( "model.off");

    Ng_Init();
    Ng_Geometry_2D *splineGeometry = Ng_NewSplineGeometry2D();

    int numEdges = 10;
    int numNodes = 10;
  
    double dtheta = 2*M_PI/(double)numEdges;
    for( int i = 0; i < numNodes; i++) {
         double u = cos( i*dtheta);
         double v = sin( i*dtheta);
         Ng_AddPoint(splineGeometry, u, v);
    }

    for( int i = 0; i < numEdges; i++)  {
         int v1 = (i+0)%numNodes;
         int v2 = (i+1)%numNodes;
         Ng_AddLineSegment(splineGeometry, v1, v2);
    }

    Ng_Mesh *ngMesh;
    Ng_Meshing_Parameters  meshParams;

    Ng_GenerateMesh_2D( splineGeometry, &ngMesh, &meshParams);
    Ng_SaveMesh( ngMesh, "ngmesh.unv");
*/

}
