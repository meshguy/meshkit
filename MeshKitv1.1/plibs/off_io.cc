#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <assert.h>

using namespace std;

#include <iMesh.h>
#include "SimpleArray.hpp"

iMesh_Instance read_off_file( const string &filename )
{
   iMesh_Instance mesh;

   ifstream ifile( filename.c_str(), ios::in);
   if( ifile.fail() ) {
       cout << "Warning: Cann't open file " << filename << endl;
       return mesh;
   }

   char  *options = NULL;
   int    optlen  = 0;
   int    err, status, result;
   iMesh_newMesh( options, &mesh, &err, optlen);  

   string str;

   ifile >> str;
   if( str != "OFF") {
       cout << "Error: Invalid header: OFF must be first string " << endl;
       return mesh;
   }

   int numNodes, numEdges, numFaces;

   ifile >> numNodes >> numFaces >> numEdges;

   const char *tagname = "GLOBAL_ID";
   int namelen = strlen( tagname);

   iBase_TagHandle idtag;
   iMesh_createTag(mesh, tagname, 1, iBase_INTEGER, &idtag, &result, namelen);

   iBase_EntityHandle  newhandle;

   vector<iBase_EntityHandle> nodeHandles;
   nodeHandles.resize( numNodes );
   double x, y, z;
   int gid;
   for( int i = 0; i < numNodes; i++) 
   {
        ifile >> x >> y >> z;
        iMesh_createVtx(mesh, x, y, z, &newhandle, &err);
        gid = i;
        iMesh_setIntData( mesh, newhandle, idtag, gid, &err);
        nodeHandles[i] = newhandle;
   }
   int nnodes, n0, n1, n2;

   vector<iBase_EntityHandle> facehandles, eConnect;
   facehandles.resize( numFaces );
   for( int i = 0; i < numFaces; i++) 
   {
        ifile >> nnodes;
        eConnect.resize(nnodes);
        for(  int j = 0; j < nnodes; j++) { 
              ifile >> n0;
              eConnect[j] = nodeHandles[n0];
        }
        switch( nnodes) 
        {
            case 3:
                 iMesh_createEnt(mesh, iMesh_TRIANGLE, &eConnect[0], nnodes, &newhandle, &status, &err);
                 break;
            case 4:
                 iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &eConnect[0], nnodes, &newhandle, &status, &err);
                 break;
            default:
                 iMesh_createEnt(mesh, iMesh_POLYGON, &eConnect[0], nnodes, &newhandle, &status, &err);
                 break;
        }
   }
   return mesh;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

int write_off_file( iMesh_Instance &mesh, const string &filename )
{
   int err;
   ofstream ofile( filename.c_str(), ios::out);
   if( ofile.fail() ) {
       cout << "Warning: Cann't open file " << filename << endl;
       return 1;
   }
   iBase_EntitySetHandle rootSet;
   iMesh_getRootSet(mesh, &rootSet, &err);

   ofile << "OFF" << endl;

   int numNodes, numEdges, numFaces;

   iMesh_getNumOfType( mesh, rootSet, iBase_VERTEX, &numNodes, &err);
   iMesh_getNumOfType( mesh, rootSet, iBase_EDGE,   &numEdges, &err);
   iMesh_getNumOfType( mesh, rootSet, iBase_FACE,   &numFaces, &err);

   ofile << numNodes << " " <<  numFaces << " " << numEdges << endl;

   SimpleArray<iBase_EntityHandle> nodeHandles;
   iMesh_getEntities(mesh, rootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(nodeHandles), &err);

   double x, y, z;
   for( int i = 0; i < numNodes; i++) 
   {
      iMesh_getVtxCoord(mesh, nodeHandles[i], &x, &y, &z, &err);
      ofile << x << "  " << y << " " << z << endl;
   }
   ofile << endl;

   iBase_TagHandle  idtag;
   iMesh_getTagHandle( mesh, "GLOBAL_ID", &idtag, &err, strlen("GLOBAL_ID") );
   assert(!err);

   int gid;
   if( numFaces ) {
       SimpleArray<iBase_EntityHandle> faceHandles, facenodes;
       iMesh_getEntities(mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(faceHandles), &err);
       for( int i = 0; i < numFaces; i++) {
            iMesh_getEntAdj(mesh, faceHandles[i], iBase_VERTEX, ARRAY_INOUT(facenodes), &err);
            ofile << facenodes.size() << " ";
            for( int j = 0; j < facenodes.size() ; j++) {
                 iMesh_getIntData(mesh, facenodes[j], idtag, &gid, &err);
                 ofile << gid << " ";
            }
            ofile << endl;
       }
   }

   return 0;
           

}
