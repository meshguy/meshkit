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
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iMesh.h>
#include <iGeom.h>
#include <set>
#include <iRel.h>
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

struct Point2D
{
	double pu;
	double pv;	
};
struct Point3D
{
	double px;
	double py;
	double pz;
};
struct Vertex
{
   	Vertex() { onBoundary = 0; }
   	int    id;
	int index;
   	bool   onBoundary;
   	bool   onCorner;
   	double uCoord, uvCoords[2], xyzCoords[3], pre[2];
	iBase_EntityHandle gVertexHandle;
};
struct Edge
{
  	Vertex* connect[2];
	int  EdgeID;
	int index;
	bool onBoundary;
	iBase_EntityHandle gEdgeHandle;
};

struct Face
{
    	int  getNumNodes() const { return connect.size(); }
    	Vertex* getVertex(int i) const { return connect[i]; }
    	vector<Vertex*> connect;
	vector<Edge*> connEdges;
    	int FaceID;
	int index;
	iBase_EntityHandle gFaceHandle;
};

struct Matrix2D
{
	int index;
	int index_x;
	int index_y;
	double value[4];
};

#endif
