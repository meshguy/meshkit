//-----------------------------------C++-------------------------------------//
// File: src/algs/OneToOneSwept.hpp
// Wednesday February 11 10:50 2011
// Brief: OneToOneSwept class definition: generate the all-hexahedral mesh by
//        general sweeping 
//---------------------------------------------------------------------------//


#ifndef __ONETOONESWEPT_H
#define __ONETOONESWEPT_H

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
#include <set>
#include <iRel.h>
#include <vector>

#include "meshkit/MeshScheme.hpp"

#include "SimpleArray.hpp"

using namespace std;

namespace MeshKit
{
//===========================================================================//
  /*!
   * \class OneToOneSwept
   * \brief Generate the all-hexahedral mesh by sweeping
   * 
   * OneToOneSwept generates the#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp" all-hexahedral mesh by sweeping the source mesh
   * to the target surface
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
   	double uCoord, uvCoords[2], xyzCoords[3];
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
    	int FaceID;
	int index;
	iBase_EntityHandle gFaceHandle;
};

class OneToOneSwept :  public MeshScheme
{	
public:

	OneToOneSwept(MKCore *mk_core, const MEntVector &me_vec);
	~OneToOneSwept();
	virtual void setup_this();
	virtual void execute_this();

	//make an instance of the EdgeMesher class
	static MeshOp *factory(MKCore *mkcore, const MEntVector &me_vec);
	void mesh_types(std::vector<moab::EntityType> &tps);
	


private:
	//Build Association function: try to build the association between the geometry and mesh
	void buildAssociation(iGeom_Instance &geom, iMesh_Instance &mesh, iRel_Instance &assoc, iRel_PairHandle &rel);

	//specify the source surface for OneToOneSwept class
	void SetSourceSurface();
	//specify the target surface for OneToOneSwept class
	void SetTargetSurface();
	
	//function for all-quad meshing on the target surface
	int TargetSurfProjection();

	//function for obtaining the parametric coordinates from x,y,z coordinates
	int getUVCoords(iBase_EntityHandle gFaceHandle, Point3D pts3, Point2D &pts2);

	//function for obtaining the x,y,z coordinates from parametric coordinates
	int getXYZCoords(iBase_EntityHandle gFaceHandle, Point3D &pts3, double uv[2]);

private://private member variable
	iBase_EntityHandle sourceSurface;
	iBase_EntityHandle targetSurface;
	std::vector<Edge> gsEdgeList;  //geometrical edges on the source surfaces
	std::vector<Edge> gtEdgeList;  //geometrical edges on the target surfaces
	std::vector<Vertex> NodeList;  //mesh nodes on the source surface
	map<int, int> edgePairs;  //store the relationship between the  edge id on the source surface and target surface
	map<int, int> cornerPairs; ////store the relationship between the corner vertex id on the source surface and target surface	
	std::vector<Vertex> gVertexList; //geometrical vertex
	std::vector<Face> gLinkFaceList;  //geometrical face list for linking surface
	std::vector<Vertex> TVertexList;
	std::vector<Edge> TEdgeList;
	std::vector<Face> TFaceList;
	std::vector<Edge> EdgeList;
	std::vector<Face> FaceList;	

	
	
};

}

#endif


