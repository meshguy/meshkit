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
#include <iMesh.h>
#include <set>
#include <iRel.h>
#include <vector>
#include <map>


#include "SimpleArray.h"

using namespace std;

class OneToOneSwept
{	
public:
	//public function
	OneToOneSwept(iGeom_Instance &geometry, iMesh_Instance &Mesh, iRel_Instance &association, iRel_RelationHandle &irel);
	~OneToOneSwept();
	void SurfaceSpecifying();
	void Execute();
	int SaveMesh(const char *FileName);

public:
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

private:
	//member variables
	iGeom_Instance geom;
	iMesh_Instance mesh;
	iRel_Instance assoc;
	iRel_RelationHandle rel;
	iBase_TagHandle  geom_id_tag, mesh_id_tag;
	iBase_EntityHandle sourceSurface;
	iBase_EntityHandle targetSurface;
	iBase_EntityHandle sourceRegion;
	iBase_EntityHandle targetRegion;
	iBase_EntitySetHandle geom_root_set, mesh_root_set;
	iBase_EntitySetHandle volumeSet;
	iBase_EntityHandle volEntity;
	vector<Vertex> NodeList;
	vector<Edge> EdgeList;
	vector<Face> FaceList;
	vector<Vertex> TVertexList;
	vector<Edge> TEdgeList;
	vector<Face> TFaceList;
	vector<Edge> gLinkSides;  //geometrical edges for linking sides between source and target
	vector<Edge> gsEdgeList;  //geometrical edges on the source surfaces
	vector<Edge> gtEdgeList;  //geometrical edges on the target surfaces
	vector<Vertex> gVertexList; //geometrical vertex
	vector<Face> gLinkFaceList;  //geometrical face list for linking surface
	int numLayers;
	map<int, int> edgePairs;  //store the relationship between the  edge id on the source surface and target surface
	map<int, int> cornerPairs; ////store the relationship between the corner vertex id on the source surface and target surface
	
	

private:
	//private functions
	void buildAssociation(iGeom_Instance &geom, iMesh_Instance &mesh, iRel_Instance &assoc, iRel_RelationHandle &rel);
	int TargetSurfProjection();
	int InnerLayerMeshing();
	int LinkSurfMeshing(vector<vector <Vertex> > &linkVertexList);
	int InnerNodesProjection(vector<vector <Vertex> > &linkVertexList);
	int CreateElements(vector<vector <Vertex> > &linkVertexList);
	double parametricTFI2D(double r, double s, double pt_0s, double pt_1s, double pt_r0, double pt_r1);
	double linear_interpolation(double r, double x0, double x1);	
	int getUVCoords(iBase_EntityHandle gFaceHandle, Point3D pts3, Point2D &pts2);
	int getXYZCoords(iBase_EntityHandle gFaceHandle, Point3D &pts3, double uv[2]);
	

	int FindCorners();
	int getList();
	int isEdgeBoundary(iBase_EntityHandle gEdgeHandle);
	int isVertexCorner(iBase_EntityHandle gNodeHandle, iBase_EntityHandle gFaceHandle);
	int is_SamePoint(double ptx1, double pty1, double ptz1, double ptx2, double pty2, double ptz2);
	int getLinkingMesh();
	
	int MeshSetting();
};


#endif

