#ifndef MY_MESH_H
#define MY_MESH_H

#include <vector>
using namespace std;

#define REAL double
#define ANSI_DECLARATORS

extern "C"
{
   double exactinit();
}

#include <string.h>

struct Vertex
{
   Vertex() { onBoundary = 0; }
   int    id;
   bool   onBoundary;
   bool   onCorner;
   double uCoord, uvCoords[2], xyzCoords[3];
};

struct Edge
{
  Vertex* connect[2];
  int  geomEdgeID;
};

struct Face
{
    int  getNumNodes() const { return connect.size(); }
    Vertex* getVertex(int i) const { return connect[i]; }

    vector<Vertex*> connect;
    int geomFaceID;
};

struct Cell
{
    int  getNumNodes() const { return connect.size(); }
    Vertex*  getVertex(int i) const { return connect[i]; }

    vector<Vertex*> connect;
    int geomCellID;
};

#endif

