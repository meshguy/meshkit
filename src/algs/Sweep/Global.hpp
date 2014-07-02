//-----------------------------------C++-------------------------------------//
// File: src/algs/Global.hpp
// Wednesday February 11 10:50 2011
// Brief: HarmonicMap class definition: do the harmonic mapping for the surface
//        mesh 
//---------------------------------------------------------------------------//


#ifndef GLOBAL1_HPP
#define GLOBAL1_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>
#include "meshkit/SimpleArray.hpp"

#include <iGeom.h>
#include "meshkit/Matrix.hpp"
#include <vector>
#include <set>
#include <list>


using namespace std;
typedef MeshKit::Vector<3> Vector3D;
typedef MeshKit::Vector<2> Vector2D;
typedef MeshKit::Matrix<3, 3> Matrix3D;

//===========================================================================//
  /*!
   * \class Global
   * \brief define the data structure for Sweeping
   * 
   */
//===========================================================================//
namespace MeshKit
{
static const double dist_tolerance = 1.0e-1;
static const double eps = 1.0e-5;
struct Vertex {
  Vertex()
  {
    onBoundary = 0;
  }
  int id;
  int index;
  bool onBoundary;
  Vector3D xyz;
  Vector2D uv;
  bool onCorner;
  iBase_EntityHandle gVertexHandle;
};
struct Edge {
  int getNumNodes() const
  {
    return connect.size();
  }
  Vertex* getVertex(int i) const
  {
    return connect[i];
  }
  vector<Vertex*> connect;//could be 1 or 2
  iBase_EntityHandle gEdgeHandle;
  int id;
  int index;
  int edge_type;//-1  corner, 0  side, 1  end, -2  reversal
  double e;
};

struct Face {
  int getNumNodes() const
  {
    return connect.size();
  }
  Vertex* getVertex(int i) const
  {
    return connect[i];
  }
  int index;
  vector<Vertex*> connect;
  vector<Edge*> connEdges;
  vector<vector<int> > vertexloops;
  vector<vector<int> > edgeloops;
  iBase_EntityHandle gFaceHandle;
  int src_tgt_link;//0--source, 1--target, 2--linking
};

}

#endif
