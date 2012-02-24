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

#include <iGeom.h>
#include "meshkit/Matrix.hpp"
#include <set>
#include <vector>
#include <set>
#include <list>

using namespace std;


//===========================================================================//
  /*!
   * \class Global
   * \brief define the data structure for Sweeping
   * 
   */
//===========================================================================//
namespace MeshKit
{

struct Vertex {
  Vertex()
  {
    onBoundary = 0;
  }
  int id;
  int index;
  bool onBoundary;
  double uCoord;
  Vector<3> xyz;
  iBase_EntityHandle gVertexHandle;
};
struct Edge {
  Vertex* connect[2];
  int EdgeID;
  int index;
  bool onBoundary;
  iBase_EntityHandle gEdgeHandle;
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
  vector<Vertex*> connect;
  vector<Edge*> connEdges;
  int FaceID;
  int index;
  iBase_EntityHandle gFaceHandle;
};

}

#endif
