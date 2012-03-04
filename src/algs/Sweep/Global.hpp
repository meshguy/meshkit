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
#include <vector>
#include <set>
#include <list>

using namespace std;
typedef MeshKit::Vector<3> Vector3D;
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

struct Vertex {
  Vertex()
  {
    onBoundary = 0;
  }
  int id;
  int index;
  bool onBoundary;
  Vector3D xyz;
  iBase_EntityHandle gVertexHandle;
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
  iBase_EntityHandle gFaceHandle;
};

}

#endif
