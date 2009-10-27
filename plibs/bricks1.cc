//////////////////////////////////////////////////////////////////////////////
// Objective:     Convert 8 nodes brick element into 20 and 27 nodes element.
// Developer:     Chaman Singh Verma
//                Argonne National Lab. Argonne. USA
// Date:          11th April, 2009.
// PI:            Dr. Tim Tautges

//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include "SimpleArray.h"

using namespace std;

iMesh_Instance mesh;
iBase_EntitySetHandle rootSet;
int nxCells, nyCells, nzCells;
double xLength, yLength, zLength;
double xOrigin, yOrigin, zOrigin;

////////////////////////////////////////////////////////////////////////////////

class Face
{
    typedef std::vector<iBase_EntityHandle> Container;
public:

    Face() { }

    Face(vector<iBase_EntityHandle> &v)
    {
        connect = v;
    }

    void setConnection(const vector<iBase_EntityHandle> &v)
    {
        connect = v;
    }

    iBase_EntityHandle getHash() const
    {
        return *min_element(connect.begin(), connect.end());
    }

    int getNumNodes() const
    {
        return connect.size();
    }

    iBase_EntityHandle getVertex(int i) const
    {
        return connect[i];
    }

    Container getConnection() const
    {
        return connect;
    }

    bool isSame(const Face &face)
    {
        size_t msize = face.getNumNodes();
        for (size_t i = 0; i < msize; i++)
        {
            if (std::find(connect.begin(), connect.end(), face.getVertex(i)) == connect.end())
                return 0;
        }

        Container mcon = face.getConnection();
        Container::iterator middle = std::find(mcon.begin(), mcon.end(), connect[0]);
        std::rotate(mcon.begin(), middle, mcon.end());

        // Possibly in the same order
        int stat = 1;
        for (size_t i = 0; i < msize; i++)
            if (connect[i] != mcon[i]) stat = 0;
        if (stat) return 1;

        std::reverse(mcon.begin(), mcon.end());
        middle = std::find(mcon.begin(), mcon.end(), connect[0]);
        std::rotate(mcon.begin(), middle, mcon.end());

        stat = 1;
        // Possibly in the reverse order ..
        for (size_t i = 0; i < msize; i++)
            if (connect[i] != mcon[i]) stat = 0;
        if (stat) return -1;

        return 0;
    }

private:
    vector<iBase_EntityHandle> connect;
};

///////////////////////////////////////////////////////////////////////////////

class EdgeMap
{
    typedef std::multimap<iBase_EntityHandle, iBase_EntityHandle> Container;

public:

    bool contains(iBase_EntityHandle &v1, iBase_EntityHandle &v2) const
    {
        iBase_EntityHandle vmin = min(v1, v2);
        iBase_EntityHandle vmax = max(v1, v2);

        Container::const_iterator it = container.find(vmin);

        if (it != container.end())
        {
            Container::const_iterator ibegin = container.lower_bound(vmin);
            Container::const_iterator iend = container.upper_bound(vmin);
            for (it = ibegin; it != iend; ++it)
                if (it->second == vmax) return 1;
        }
        return 0;
    }

    int add(iBase_EntityHandle &v1, iBase_EntityHandle &v2)
    {
        iBase_EntityHandle vmin = min(v1, v2);
        iBase_EntityHandle vmax = max(v1, v2);
        container.insert(make_pair(vmin, vmax));
    }

    int clear()
    {
        container.clear();
    }
    size_t getSize() const { return container.size(); }

private:
    Container container;
};

///////////////////////////////////////////////////////////////////////////////

class FaceMap
{
    typedef std::multimap<iBase_EntityHandle, Face> Container;

public:

    bool contains(const Face &checkface) const
    {
        iBase_EntityHandle vhash = checkface.getHash();

        Container::const_iterator it = container.find(vhash);

        if (it != container.end())
        {
            Container::const_iterator ibegin = container.lower_bound(vhash);
            Container::const_iterator iend = container.upper_bound(vhash);
            for (it = ibegin; it != iend; ++it)
            {
                Face existface = it->second;
                int status = existface.isSame(checkface);
                if (status != 0) return status;
            }
        }
        return 0;
    }

    int add(const Face &face)
    {
        iBase_EntityHandle vhash = face.getHash();
        container.insert(make_pair(vhash, face));
    }

    int clear()
    {
        container.clear();
    }
    size_t getSize() const { return container.size(); }

private:
    Container container;
};

///////////////////////////////////////////////////////////////////////////////

int addEdges(const vector<iBase_EntityHandle> &hexConnect, EdgeMap &edgemap)
{
    int err, status, isPresent;
    iBase_EntityHandle eConnect[2], edgehandle;

    // Edge 1
    eConnect[0] = hexConnect[0];
    eConnect[1] = hexConnect[1];

    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 2
    eConnect[0] = hexConnect[1];
    eConnect[1] = hexConnect[2];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 3
    eConnect[0] = hexConnect[2];
    eConnect[1] = hexConnect[3];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 4
    eConnect[0] = hexConnect[3];
    eConnect[1] = hexConnect[0];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 5
    eConnect[0] = hexConnect[4];
    eConnect[1] = hexConnect[5];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 6
    eConnect[0] = hexConnect[5];
    eConnect[1] = hexConnect[6];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 7
    eConnect[0] = hexConnect[6];
    eConnect[1] = hexConnect[7];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 8
    eConnect[0] = hexConnect[7];
    eConnect[1] = hexConnect[4];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 9
    eConnect[0] = hexConnect[0];
    eConnect[1] = hexConnect[4];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 10
    eConnect[0] = hexConnect[1];
    eConnect[1] = hexConnect[5];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 11
    eConnect[0] = hexConnect[3];
    eConnect[1] = hexConnect[7];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }

    // Edge 12
    eConnect[0] = hexConnect[2];
    eConnect[1] = hexConnect[6];
    if (!edgemap.contains(eConnect[0], eConnect[1]))
    {
        edgemap.add(eConnect[0], eConnect[1]);
        iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &edgehandle, &status, &err);
    }
}

///////////////////////////////////////////////////////////////////////////////

int addFaces(const vector<iBase_EntityHandle> &hexConnect, FaceMap &facemap)
{
    int err, status, isPresent;
    Face face;
    iBase_EntityHandle facehandle;
    vector<iBase_EntityHandle> fConnect(4);

    // Face 1:  Left Side
    fConnect[0] = hexConnect[0];
    fConnect[1] = hexConnect[4];
    fConnect[2] = hexConnect[7];
    fConnect[3] = hexConnect[3];
    face.setConnection(fConnect);
    if (!facemap.contains(face))
    {
        facemap.add(face);
        iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &fConnect[0], 4, &facehandle, &status, &err);
    }

    // Face 2:  Right side
    fConnect[0] = hexConnect[1];
    fConnect[1] = hexConnect[2];
    fConnect[2] = hexConnect[6];
    fConnect[3] = hexConnect[5];
    face.setConnection(fConnect);
    if (!facemap.contains(face))
    {
        facemap.add(face);
        iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &fConnect[0], 4, &facehandle, &status, &err);
    }

    // Face 3:  Bottom side
    fConnect[0] = hexConnect[0];
    fConnect[1] = hexConnect[1];
    fConnect[2] = hexConnect[5];
    fConnect[3] = hexConnect[4];
    face.setConnection(fConnect);
    if (!facemap.contains(face))
    {
        facemap.add(face);
        iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &fConnect[0], 4, &facehandle, &status, &err);
    }

    // Face 4:  Top  side
    fConnect[0] = hexConnect[3];
    fConnect[1] = hexConnect[7];
    fConnect[2] = hexConnect[6];
    fConnect[3] = hexConnect[2];
    face.setConnection(fConnect);
    if (!facemap.contains(face))
    {
        facemap.add(face);
        iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &fConnect[0], 4, &facehandle, &status, &err);
    }

    // Face 5:  Back Side
    fConnect[0] = hexConnect[0];
    fConnect[1] = hexConnect[3];
    fConnect[2] = hexConnect[2];
    fConnect[3] = hexConnect[1];
    face.setConnection(fConnect);
    if (!facemap.contains(face))
    {
        facemap.add(face);
        iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &fConnect[0], 4, &facehandle, &status, &err);
    }

    // Face 6:  Front Side
    fConnect[0] = hexConnect[4];
    fConnect[1] = hexConnect[5];
    fConnect[2] = hexConnect[6];
    fConnect[3] = hexConnect[7];
    face.setConnection(fConnect);
    if (!facemap.contains(face))
    {
        facemap.add(face);
        iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &fConnect[0], 4, &facehandle, &status, &err);
    }

}

///////////////////////////////////////////////////////////////////////////////
int  read_brick8_nodes( char *filename)
{
    ifstream infile(filename, ios::in);

    if( infile.fail() ) {
        cout << "Error: cann't read file " << filename << endl;
        return 1;
    }

    EdgeMap edgemap;
    FaceMap facemap;

    string s;
    infile >> s;  assert(s == "#Nodes");

    int numNodes;
    infile >> numNodes;
   
    int err;

    vector<iBase_EntityHandle> nodeHandles;
    nodeHandles.resize(numNodes);

    double x, y, z;

    iBase_EntityHandle newhandle;
    int id, index = 0;
    for (int i = 0; i < numNodes; i++) {
         infile >> id >> x >> y >> z;
         iMesh_createVtx(mesh, x, y, z, &newhandle, &err);
         nodeHandles[index++] = newhandle;
    }

    infile >> s;  assert( s == "#Hex");
    int numCells;
    infile >> numCells;
    
    vector<iBase_EntityHandle> cellHandles;
    cellHandles.resize(numCells);

    vector<iBase_EntityHandle> connect(8);
    int status, n0, n1, n2, n3, n4, n5, n6, n7;
    index = 0;
    for (int i = 0; i < numCells; i++)
    {
       infile >> n0 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7;
       connect[0] = nodeHandles[n0];
       connect[1] = nodeHandles[n1];
       connect[2] = nodeHandles[n3];
       connect[3] = nodeHandles[n2];
       connect[4] = nodeHandles[n4];
       connect[5] = nodeHandles[n5];
       connect[6] = nodeHandles[n7];
       connect[7] = nodeHandles[n6];
       iMesh_createEnt(mesh, iMesh_HEXAHEDRON, &connect[0], 8, &newhandle, &status, &err);
       cellHandles[index++] = newhandle;
       addEdges(connect, edgemap);
       addFaces(connect, facemap);
    }
}

int simple_brick8_nodes()
{
    EdgeMap edgemap;
    FaceMap facemap;

    int err;
    int nxNodes = nxCells + 1;
    int nyNodes = nyCells + 1;
    int nzNodes = nzCells + 1;
    int numNodes = nxNodes * nyNodes*nzNodes;

    vector<iBase_EntityHandle> nodeHandles;
    nodeHandles.resize(numNodes);

    double x, y, z;
    double dx = xLength / (double) nxCells;
    double dy = yLength / (double) nyCells;
    double dz = zLength / (double) nzCells;

    iBase_EntityHandle newhandle;
    int index = 0;
    for (int k = 0; k < nzNodes; k++)
    {
        z = zOrigin + k*dz;
        for (int j = 0; j < nyNodes; j++)
        {
            y = yOrigin + j*dy;
            for (int i = 0; i < nxNodes; i++)
            {
                x = xOrigin + i*dx;
                iMesh_createVtx(mesh, x, y, z, &newhandle, &err);
                nodeHandles[index++] = newhandle;
            }
        }
    }

    int numCells = nxCells * nyCells*nzCells;
    vector<iBase_EntityHandle> cellHandles;
    cellHandles.resize(numCells);

    vector<iBase_EntityHandle> connect(8);
    int offset, status;
    index = 0;
    for (int k = 0; k < nzCells; k++)
    {
        for (int j = 0; j < nyCells; j++)
        {
            for (int i = 0; i < nxCells; i++)
            {
                offset = k * nxNodes * nyNodes + j * nxNodes + i;
                connect[0] = nodeHandles[offset];
                connect[1] = nodeHandles[offset + 1];
                connect[2] = nodeHandles[offset + nxNodes + 1];
                connect[3] = nodeHandles[offset + nxNodes];

                offset = (k + 1) * nxNodes * nyNodes + j * nxNodes + i;
                connect[4] = nodeHandles[offset];
                connect[5] = nodeHandles[offset + 1];
                connect[6] = nodeHandles[offset + nxNodes + 1];
                connect[7] = nodeHandles[offset + nxNodes];
                iMesh_createEnt(mesh, iMesh_HEXAHEDRON, &connect[0], 8, &newhandle, &status, &err);
                cellHandles[index++] = newhandle;
                addEdges(connect, edgemap);
                addFaces(connect, facemap);
            }
        }
    }

    char *options = NULL;
    int optlen = 0;
    const char *outfile = "model0.vtk";
    int namelen = strlen(outfile);
    iMesh_save(mesh, rootSet, outfile, options, &err, namelen, optlen);
}
/////////////////////////////////////////////////////////////////////////////////
iBase_EntityHandle searchEntity( SimpleArray<iBase_EntityHandle> &elementHandles, 
                                 SimpleArray<iBase_EntityHandle> &vertexHandles)
{
   int err, nodematch;
   SimpleArray<iBase_EntityHandle> elemNodes;

   std::sort( vertexHandles.begin(), vertexHandles.end() );

   for( int i = 0; i < elementHandles.size(); i++) 
   {
        elemNodes.clear();
        iMesh_getEntAdj(mesh, elementHandles[i], iBase_VERTEX, ARRAY_INOUT(elemNodes), &err);
         if( elemNodes.size() == vertexHandles.size() ) 
         {
             std::sort( elemNodes.begin(), elemNodes.end() );
             nodematch = 1;
             for( int j = 0; j < elemNodes.size(); j++) 
                  if( elemNodes[j] != vertexHandles[j] ) nodematch = 0;
             if( nodematch ) return elementHandles[i];
         }
   }

   cout << "Fatal Error: No Entity match in array " << endl;
   exit(0);
}

/////////////////////////////////////////////////////////////////////////////////

void get_brick20_nodes(iBase_EntityHandle cellHandle, iBase_TagHandle hordertag, 
                       vector<iBase_EntityHandle> &brick20nodes)
{
   int err;

   SimpleArray<iBase_EntityHandle> hexVertexHandles;
   iMesh_getEntAdj(mesh, cellHandle,  iBase_VERTEX, ARRAY_INOUT(hexVertexHandles), &err);

   brick20nodes.resize(20);

   brick20nodes[0]  =  hexVertexHandles[0];
   brick20nodes[2]  =  hexVertexHandles[1];
   brick20nodes[5]  =  hexVertexHandles[3];
   brick20nodes[7]  =  hexVertexHandles[2];
   brick20nodes[12] =  hexVertexHandles[4];
   brick20nodes[14] =  hexVertexHandles[5];
   brick20nodes[17] =  hexVertexHandles[7];
   brick20nodes[19] =  hexVertexHandles[6];

   // Now search vertices on the edges;
   iBase_EntityHandle edgehandle, vhandle;
   SimpleArray<iBase_EntityHandle> edgeHandles;
   iMesh_getEntAdj(mesh, cellHandle, iBase_EDGE, ARRAY_INOUT(edgeHandles), &err);

   SimpleArray<iBase_EntityHandle>  eConnect;
   eConnect.resize(2);

   eConnect[0] = hexVertexHandles[0];
   eConnect[1] = hexVertexHandles[1];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[1] =  vhandle;

   eConnect[0] = hexVertexHandles[3];
   eConnect[1] = hexVertexHandles[0];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[3] =  vhandle;

   eConnect[0] = hexVertexHandles[1];
   eConnect[1] = hexVertexHandles[2];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[4] =  vhandle;

   eConnect[0] = hexVertexHandles[2];
   eConnect[1] = hexVertexHandles[3];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[6] =  vhandle;

   eConnect[0] = hexVertexHandles[0];
   eConnect[1] = hexVertexHandles[4];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[8] =  vhandle;

   eConnect[0] = hexVertexHandles[1];
   eConnect[1] = hexVertexHandles[5];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[9] =  vhandle;

   eConnect[0] = hexVertexHandles[3];
   eConnect[1] = hexVertexHandles[7];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[10] =  vhandle;

   eConnect[0] = hexVertexHandles[2];
   eConnect[1] = hexVertexHandles[6];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[11] =  vhandle;

   eConnect[0] = hexVertexHandles[4];
   eConnect[1] = hexVertexHandles[5];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[13] =  vhandle;

   eConnect[0] = hexVertexHandles[7];
   eConnect[1] = hexVertexHandles[4];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[15] =  vhandle;
   
   eConnect[0] = hexVertexHandles[5];
   eConnect[1] = hexVertexHandles[6];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[16] =  vhandle;

   eConnect[0] = hexVertexHandles[6];
   eConnect[1] = hexVertexHandles[7];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick20nodes[18] =  vhandle;
}

////////////////////////////////////////////////////////////////////////////////

int gen_brick20nodes()
{
    int err, result;

    // Generate high order node on each edge of the mesh.
    SimpleArray<iBase_EntityHandle> edgeHandles, edgenodes;
    iMesh_getEntities(mesh, rootSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(edgeHandles), &err);

    const char *tagname = "HO_POINTS";
    int namelen = strlen( tagname);

    iBase_TagHandle tag;
    iMesh_createTag(mesh, tagname, 1, iBase_ENTITY_HANDLE, &tag, &result, namelen); 

    double x0, y0, z0;
    double x1, y1, z1;
    double xm, ym, zm;

    iBase_EntitySetHandle vertexSet;
    iMesh_createEntSet( mesh, 0, &vertexSet, &err);
 
    iBase_EntityHandle newhandle;
    int numEdges = edgeHandles.size(); 
    for( int i = 0; i < numEdges; i++) {
         iMesh_getEntAdj( mesh, edgeHandles[i], iBase_VERTEX, ARRAY_INOUT(edgenodes), &err);
         iMesh_getVtxCoord( mesh, edgenodes[0], &x0, &y0, &z0, &err);
         iMesh_getVtxCoord( mesh, edgenodes[1], &x1, &y1, &z1, &err);
         xm = 0.5*(x0+x1);
         ym = 0.5*(y0+y1);
         zm = 0.5*(z0+z1);
         iMesh_createVtx(mesh, xm, ym, zm, &newhandle, &err);
         iMesh_setEHData( mesh, edgeHandles[i], tag, newhandle, &err);
    }

    // All nodes generated. Now collect and order the nodes for storage.
   
    SimpleArray<iBase_EntityHandle> meshnodes;
    iMesh_getEntities( mesh, rootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(meshnodes), &err);

    iBase_TagHandle  idtag;
    iMesh_getTagHandle( mesh, "GLOBAL_ID", &idtag, &err, strlen("GLOBAL_ID") );
    assert(!err);

    int gid;
    int numNodes = meshnodes.size();

    const char *outfile = "hgelems20.dat";
    ofstream ofile( outfile, ios::out);

    ofile << "#Nodes " << numNodes << endl;
    for( int i = 0; i < numNodes; i++) 
    {
        iMesh_getVtxCoord( mesh, meshnodes[i], &x0, &y0, &z0, &err);
        iMesh_setIntData( mesh,  meshnodes[i], idtag, i,  &err);
        ofile << x0 << " " << y0 << " " << z0 << endl;
    }

    SimpleArray<iBase_EntityHandle> cellHandles;
    iMesh_getEntities(mesh, rootSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(cellHandles), &err);

    ofile << "#Hex " << cellHandles.size() << endl;

    vector<iBase_EntityHandle> brick20nodes;
    for( int i = 0; i < cellHandles.size(); i++) {
         get_brick20_nodes(cellHandles[i], tag, brick20nodes);
         
         for( int j = 0; j < 20; j++) {
              iMesh_getIntData( mesh, brick20nodes[j], idtag, &gid, &err);
              ofile << gid << " ";
         }
         ofile << endl;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void get_brick27_nodes( iBase_EntityHandle cellHandle, iBase_TagHandle hordertag, 
                        vector<iBase_EntityHandle> &brick27nodes)
{
   int err;

   SimpleArray<iBase_EntityHandle> hexVertexHandles;
   iMesh_getEntAdj(mesh, cellHandle,  iBase_VERTEX, ARRAY_INOUT(hexVertexHandles), &err);

   brick27nodes.resize(27);

   brick27nodes[0]  =  hexVertexHandles[0];
   brick27nodes[2]  =  hexVertexHandles[1];
   brick27nodes[6]  =  hexVertexHandles[3];
   brick27nodes[8]  =  hexVertexHandles[2];
   brick27nodes[18] =  hexVertexHandles[4];
   brick27nodes[20] =  hexVertexHandles[5];
   brick27nodes[24] =  hexVertexHandles[7];
   brick27nodes[26] =  hexVertexHandles[6];

   // Now search vertices on the edges;
   iBase_EntityHandle edgehandle, vhandle;
   SimpleArray<iBase_EntityHandle> edgeHandles;
   iMesh_getEntAdj(mesh, cellHandle, iBase_EDGE, ARRAY_INOUT(edgeHandles), &err);

   SimpleArray<iBase_EntityHandle>  eConnect;
   eConnect.resize(2);

   eConnect[0] = hexVertexHandles[0];
   eConnect[1] = hexVertexHandles[1];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[1] =  vhandle;

   eConnect[0] = hexVertexHandles[3];
   eConnect[1] = hexVertexHandles[0];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[3] =  vhandle;

   eConnect[0] = hexVertexHandles[1];
   eConnect[1] = hexVertexHandles[2];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[5] =  vhandle;

   eConnect[0] = hexVertexHandles[2];
   eConnect[1] = hexVertexHandles[3];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[7] =  vhandle;

   eConnect[0] = hexVertexHandles[0];
   eConnect[1] = hexVertexHandles[4];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[9] =  vhandle;

   eConnect[0] = hexVertexHandles[1];
   eConnect[1] = hexVertexHandles[5];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[11] =  vhandle;

   eConnect[0] = hexVertexHandles[3];
   eConnect[1] = hexVertexHandles[7];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[15] =  vhandle;

   eConnect[0] = hexVertexHandles[2];
   eConnect[1] = hexVertexHandles[6];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[17] =  vhandle;

   eConnect[0] = hexVertexHandles[4];
   eConnect[1] = hexVertexHandles[5];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[19] =  vhandle;

   eConnect[0] = hexVertexHandles[7];
   eConnect[1] = hexVertexHandles[4];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[21] =  vhandle;
   
   eConnect[0] = hexVertexHandles[5];
   eConnect[1] = hexVertexHandles[6];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[23] =  vhandle;

   eConnect[0] = hexVertexHandles[6];
   eConnect[1] = hexVertexHandles[7];
   edgehandle  = searchEntity(edgeHandles, eConnect);
   iMesh_getEHData( mesh, edgehandle, hordertag, &vhandle, &err);
   brick27nodes[25] =  vhandle;

   // Now search vertices on the face;
   iBase_EntityHandle facehandle;
   SimpleArray<iBase_EntityHandle> faceHandles;
   iMesh_getEntAdj(mesh, cellHandle, iBase_FACE, ARRAY_INOUT(faceHandles), &err);

   SimpleArray<iBase_EntityHandle> fConnect;
   fConnect.resize(4);

   fConnect[0] = hexVertexHandles[0];
   fConnect[1] = hexVertexHandles[1];
   fConnect[2] = hexVertexHandles[2];
   fConnect[3] = hexVertexHandles[3];
   facehandle  = searchEntity(faceHandles, fConnect);
   iMesh_getEHData( mesh, facehandle, hordertag, &vhandle, &err);
   brick27nodes[4] =  vhandle;

   fConnect[0] = hexVertexHandles[4];
   fConnect[1] = hexVertexHandles[5];
   fConnect[2] = hexVertexHandles[6];
   fConnect[3] = hexVertexHandles[7];
   facehandle  = searchEntity(faceHandles, fConnect);
   iMesh_getEHData( mesh, facehandle, hordertag, &vhandle, &err);
   brick27nodes[22] =  vhandle;

   fConnect[0] = hexVertexHandles[0];
   fConnect[1] = hexVertexHandles[4];
   fConnect[2] = hexVertexHandles[7];
   fConnect[3] = hexVertexHandles[3];
   facehandle  = searchEntity(faceHandles, fConnect);
   iMesh_getEHData( mesh, facehandle, hordertag, &vhandle, &err);
   brick27nodes[12] =  vhandle;

   fConnect[0] = hexVertexHandles[1];
   fConnect[1] = hexVertexHandles[2];
   fConnect[2] = hexVertexHandles[6];
   fConnect[3] = hexVertexHandles[5];
   facehandle  = searchEntity(faceHandles, fConnect);
   iMesh_getEHData( mesh, facehandle, hordertag, &vhandle, &err);
   brick27nodes[14] =  vhandle;

   fConnect[0] = hexVertexHandles[0];
   fConnect[1] = hexVertexHandles[1];
   fConnect[2] = hexVertexHandles[5];
   fConnect[3] = hexVertexHandles[4];
   facehandle  = searchEntity(faceHandles, fConnect);
   iMesh_getEHData( mesh, facehandle, hordertag, &vhandle, &err);
   brick27nodes[10] =  vhandle;

   fConnect[0] = hexVertexHandles[3];
   fConnect[1] = hexVertexHandles[7];
   fConnect[2] = hexVertexHandles[6];
   fConnect[3] = hexVertexHandles[2];
   facehandle  = searchEntity(faceHandles, fConnect);
   iMesh_getEHData( mesh, facehandle, hordertag, &vhandle, &err);
   brick27nodes[16] =  vhandle;

   iMesh_getEHData( mesh, cellHandle, hordertag, &vhandle, &err);
   brick27nodes[13] =  vhandle;
}

////////////////////////////////////////////////////////////////////////////////

int gen_brick27nodes()
{
    int err, result;
    double x0, y0, z0;
    double xsum, ysum, zsum;
    double xm, ym, zm;

    const char *tagname = "HO_POINTS";
    int namelen = strlen( tagname);

    iBase_TagHandle tag;
    iMesh_createTag(mesh, tagname, 1, iBase_ENTITY_HANDLE, &tag, &result, namelen); 

    SimpleArray<iBase_EntityHandle> cellHandles;
    SimpleArray<iBase_EntityHandle> faceHandles;
    SimpleArray<iBase_EntityHandle> edgeHandles;

    iBase_EntityHandle vhandle, ehandle, fhandle;

    iMesh_getEntities(mesh, rootSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(faceHandles), &err);

    int numFaces = faceHandles.size(); 
    for( int i = 0; i < numFaces; i++) {
         edgeHandles.clear();
         iMesh_getEntAdj( mesh, faceHandles[i], iBase_EDGE, ARRAY_INOUT(edgeHandles), &err);
         xsum = 0.0;
         ysum = 0.0;
         zsum = 0.0;
         int  numEdges = edgeHandles.size();
         for( int j = 0; j < numEdges; j++) {
              iMesh_getEHData( mesh, edgeHandles[j], tag, &vhandle, &err);
              iMesh_getVtxCoord( mesh, vhandle, &x0, &y0, &z0, &err);
              xsum += x0;
              ysum += y0;
              zsum += z0;
         }
         xm = xsum/(double)numEdges;
         ym = ysum/(double)numEdges;
         zm = zsum/(double)numEdges;
         iMesh_createVtx(mesh, xm, ym, zm, &vhandle, &err);
         iMesh_setEHData( mesh, faceHandles[i], tag, vhandle, &err);
    }

    iMesh_getEntities(mesh, rootSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(cellHandles), &err);

    int numCells = cellHandles.size(); 
    for( int i = 0; i < numCells; i++) {
         faceHandles.clear();
         iMesh_getEntAdj( mesh, cellHandles[i], iBase_FACE, ARRAY_INOUT(faceHandles), &err);
         xsum = 0.0;
         ysum = 0.0;
         zsum = 0.0;
         int numFaces = faceHandles.size();
         for( int j = 0; j < numFaces; j++) {
              iMesh_getEHData( mesh, faceHandles[j], tag, &vhandle, &err);
              iMesh_getVtxCoord( mesh, vhandle, &x0, &y0, &z0, &err);
              xsum += x0;
              ysum += y0;
              zsum += z0;
         }
         xm = xsum/(double)numFaces;
         ym = ysum/(double)numFaces;
         zm = zsum/(double)numFaces;
         iMesh_createVtx(mesh, xm, ym, zm, &vhandle, &err);
         iMesh_setEHData( mesh, cellHandles[i], tag, vhandle, &err);
    }

    SimpleArray<iBase_EntityHandle> meshnodes;
    iMesh_getEntities( mesh, rootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(meshnodes), &err);

    //
    // Presently MOAB Vtk Writer doesn't write only nodes elements properly, therefore
    // I am writing coordinates in "OFF" format first and then converting into "VTK" format
    // using external package. This step is only for visualization purpose.
    //  
    iBase_TagHandle  idtag;
    iMesh_getTagHandle( mesh, "GLOBAL_ID", &idtag, &err, strlen("GLOBAL_ID") );
    assert(!err);
    
    int numNodes = meshnodes.size();
    const char *outfile = "homodel27.dat";
    ofstream ofile( outfile, ios::out);

    ofile << "#Nodes " << numNodes << endl;

    for( int i = 0; i < numNodes; i++) {
        iMesh_getVtxCoord( mesh, meshnodes[i], &x0, &y0, &z0, &err);
        iMesh_setIntData( mesh,  meshnodes[i], idtag, i,  &err);
        ofile << x0 << " " << y0 << " " << z0 << endl;
    }

    ofile << "#Hex " << cellHandles.size() << endl;

    int gid;
    vector<iBase_EntityHandle> brick27nodes;

    for( int i = 0; i < cellHandles.size(); i++) {
         get_brick27_nodes(cellHandles[i], tag, brick27nodes);
         for( int j = 0; j < 27; j++) {
              iMesh_getIntData( mesh, brick27nodes[j], idtag, &gid, &err);
              ofile << gid << " ";
         }
         ofile << endl;
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    char *options = NULL;
    int optlen = 0;
    int err;
    iMesh_newMesh(options, &mesh, &err, optlen);
    iMesh_getRootSet(mesh, &rootSet, &err);

/*
//  Generate Simple Hex Mesh ....
    nxCells = 10;   nyCells = 1;   nzCells = 1;
    xOrigin = 0.0; yOrigin = 0.0; zOrigin = 0.0;
    xLength = 1.0; yLength = 1.0; zLength = 1.0;
    simple_brick8_nodes();
*/
    
    if( argc != 2 ) {
        cout << "Usage: 8 nodes hex elements file " << endl;
        return 1;
    }
    read_brick8_nodes( argv[1] );

//
//  Given 8 node bricks, generate 20 nodes or 27 nodes bricks. For 27 nodes, first
//  generate 20 nodes brick element.
//
    gen_brick20nodes();
    gen_brick27nodes();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////


