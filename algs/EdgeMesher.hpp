#ifndef __EDGEMESHER_H
#define __EDGEMESHER_H

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

#include "meshkit/MeshScheme.hpp"

//#include <boost/array.hpp>

#include "SimpleArray.hpp"

namespace MeshKit
{


struct Point3D
{
	double px;
	double py;
	double pz;	
};

using namespace std;

class EdgeMesher : public MeshScheme
{
public:
	
	enum EdgeSchemeType={equalMesh=0, biasMesh, dualMesh, curvatureMesh};

	
public:
	//EdgeMesher(iGeom_Instance &geom, iBase_EntityHandle *EdgeHandle, iMesh_Instance &Mesh, iRel_Instance &association, iRel_PairHandle *irel);
	EdgeMesher(MKCore *mk_core, const MEVector &me_vec);
	void mesh_types(std::vector<EntityType> &tps);
	virtual void setup_this();
	virtual void execute_this();
	static MeshOp *factory(MKCore *mkcore, const MEVector &me_vec);


	~EdgeMesher();
	
private:
	EdgeSchemeType schemeType;	
	
	//member variables
	iGeom_Instance geometry;
	iMesh_Instance mesh;
	iRel_Instance assoc;
	iRel_PairHandle rel;

	iBase_EntityHandle gEdgeHandle;

	Point3D getXYZCoords(double u) const;
	double getUCoord(double ustart, double dist, double uguess) const;
	void DivideIntoMore(Point3D p0, Point3D pMid, Point3D p1, double u0, double u1, double uMid, int &index, vector<double> &nodes, vector<double> &URecord);
	bool ErrorCalculate(Point3D p0, Point3D p1, Point3D pMid);
	void RapidSorting(vector<double> &nodes, vector<double> &URecord, int left, int right);
	void QuickSorting(vector<double> &nodes, vector<double> &URecord, int count);
	void GetRelatedEntitySet(iBase_EntityHandle edgeHandle, iBase_EntitySetHandle &entitySet);
	void get_related_entityset(iBase_EntitySetHandle &mesh_entityset);


	void EqualMeshing(ModelEnt *ent, int num_edges, std::vector<double> &coords);
	void BiasMeshing(ModelEnt *ent, int num_edges, std::vector<double> &coords);
	void DualBiasMeshing(ModelEnt *ent, int &num_edges, std::vector<double> &coords);
	void CurvatureMeshing(ModelEnt *ent, int &num_edges, std::vector<double> &coords);

	double measure(iGeom::EntityHandle ent, double ustart, double uend) const;

};

inline void EdgeMesher::set_edge_scheme(EdgeMesher::EdgeSchemeType scheme)
{
	schemeType = scheme;
}

inline EdgeMesher::EdgeSchemeType EdgeMesher::get_edge_scheme() const
{
	return schemeType;
}

}

#endif
