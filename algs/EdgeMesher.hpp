#ifndef MESHKIT_EDGE_MESHER_HPP
#define MESHKIT_EDGE_MESHER_HPP

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

//#include <boost/array.hpp>

#include "SimpleArray.hpp"

struct Point3D
{
	double px;
	double py;
	double pz;	
};

using namespace std;

class EdgeMesher
{
public:
	
	enum SchemeType {equalMesh=0, biasMesh, dualMesh, curvatureMesh};

	
public:
	EdgeMesher(iGeom_Instance &geom, iBase_EntityHandle *EdgeHandle, iMesh_Instance &Mesh, iRel_Instance &association, iRel_PairHandle *irel);
	~EdgeMesher();
	double getLength() const;
	double getLength(double ustart, double uend) const;
	void SetScheme(int SetOption);
	void SetStepSize(double StepSize);
	const char *getSchemeName();
	double getStepSize();
	void ShowCoorData();
	
	int Execute();
	int SaveFile(const char *FileName);
	
private:
	//member variables
	iGeom_Instance geometry;
	iMesh_Instance mesh;
	iRel_Instance assoc;
	iRel_PairHandle rel;
	
	//int SchemeOption;
	int NumEdges;
	double stepSize;
	SchemeType SchemeOption;
	const char *SelectedShemeName;

	iBase_EntityHandle gEdgeHandle;

	double umin, umax;
	bool periodic;
	int NumVertex;
	
	//private functions
	void EdgeMesh();
	vector<double> EqualMeshing();
	vector<double> BiasMeshing();
	vector<double> DualBiasMeshing();
	vector<double> CurvatureMeshing();
	//void CreateAssociation();
	Point3D getXYZCoords(double u) const;
	double getUCoord(double ustart, double dist, double uguess) const;
	void DivideIntoMore(Point3D p0, Point3D pMid, Point3D p1, double u0, double u1, double uMid, int &index, vector<double> &nodes, vector<double> &URecord);
	bool ErrorCalculate(Point3D p0, Point3D p1, Point3D pMid);
	void RapidSorting(vector<double> &nodes, vector<double> &URecord, int left, int right);
	void QuickSorting(vector<double> &nodes, vector<double> &URecord, int count);
	void GetRelatedEntitySet(iBase_EntityHandle edgeHandle, iBase_EntitySetHandle &entitySet);
	void get_related_entityset(iBase_EntitySetHandle &mesh_entityset);

};


#endif
