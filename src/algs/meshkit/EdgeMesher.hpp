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

#include "iGeom.hh"
#include <set>
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
	
	enum EdgeSchemeType {equalMesh=0, biasMesh, dualMesh, curvatureMesh};

	
public:
	//EdgeMesher(iGeom_Instance &geom, iBase_EntityHandle *EdgeHandle, iMesh_Instance &Mesh, iRel_Instance &association, iRel_PairHandle *irel);
	EdgeMesher(MKCore *mk_core, const MEntVector &me_vec);
	void mesh_types(std::vector<moab::EntityType> &tps);
	virtual void setup_this();
	virtual void execute_this();
	static MeshOp *factory(MKCore *mkcore, const MEntVector &me_vec);
	EdgeSchemeType get_edge_scheme() const;
	 void set_edge_scheme(EdgeSchemeType scheme);


	~EdgeMesher();
	
private:
	EdgeSchemeType schemeType;	
	
	//member variables
	iGeom_Instance geometry;
	iMesh_Instance mesh;
	iRel_Instance assoc;
	iRel_PairHandle rel;

	iBase_EntityHandle gEdgeHandle;

	Point3D getXYZCoords(ModelEnt *ent, double u) const;
	double getUCoord(ModelEnt *ent, double ustart, double dist, double uguess, double umin, double umax) const;
	void DivideIntoMore(ModelEnt *ent, Point3D p0, Point3D pMid, Point3D p1, double u0, double u1, double uMid, int &index, vector<double> &nodes, vector<double> &URecord);
	bool ErrorCalculate(ModelEnt *ent, Point3D p0, Point3D p1, Point3D pMid);
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
