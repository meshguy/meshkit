#ifndef __SURFPROHARMONICMAP_HPP
#define __SURFPROHARMONICMAP_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/iRel.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SimpleArray.hpp"
#include "Global.hpp"
#include "HarmonicMapper.hpp"
#include <set>
#include <vector>
#include <map>


using namespace std;

namespace MeshKit {

class SurfProHarmonicMap
{	
public:
	//public function
	SurfProHarmonicMap(MKCore* core, iBase_EntityHandle s, iBase_EntityHandle t, iBase_EntityHandle v);
	~SurfProHarmonicMap();
	void match();
	void projection();
	void setMeshData(vector<Vertex> &s, vector<Vertex> &t, vector<Face> &f);
	void getMeshData(vector<Vertex> &v);
private:
	void preprocessing();
	void addFaceToList(iBase_EntityHandle entity, Face& f, int index, bool is_set_int);
	void GetGeomLoops(Face surf, vector<vector<int> > &loops_vertex, vector<vector<int> > &loops_edge);
	void postProcessGeomLoops(Face& surf);
	void getFacets();
	void MapFacetGeom(Face f_surf, vector<Vertex> facet_node, std::map<int, int> &map_data, vector<list<int> > &geom_facet, int size_facet_v);
	//boundary distribution of the outmost boundary on the unit disk
	void boundaryDistribution();
	//distribute facet vertices of outmost boundary on an unit disk
	void boundaryUnitDisk(Face f, vector<double> dist, vector<Vertex> &facet_v, std::vector<std::list<int> > geom_facet);
    //establish the correspondence of boundary vertices on the outmost boundary
    void LocateBoundaryNodesTarget(); 
	//Distribute the interior nodes on the unit disk
	void test();
	//compute the weight for each edge
	void ComputeWeight();
	//add the facet edge to the list
	void addEdgeToList(int a, int b, int &count, vector<set<int> > &edge_connect, vector<Edge> &f_edges, vector<Vertex> &facet_v);
	//compute the edge weight
	void computeEdgeWeight(vector<Edge> &f_edges, vector<Face> &f, vector<Vertex> f_v);
	//add extra stuff to the facet mesh
	void addExtra(Face f, vector<Vertex> &facet_v, vector<Edge> &facet_e, vector<Face> &facet_tri, vector<list<int> > geom_facet, int size_facet_v);
	//find a specific triangle where a point is located
	int findFacetTri(vector<Face> &facet_tri, vector<Vector3D> nrml, Vector3D xyz, Vector3D &uvw);
	int findFacetTri(vector<Face> &facet_tri, Vector2D uv, Vector3D &uvw);
	void prjPtsToTri(Face tri, Vector3D pts, Vector3D nrml, Vector3D &xyz);
	bool ComputeBarycentric(Vector3D a, Vector3D b, Vector3D c, Vector3D xyz, Vector3D &uvw);
	bool ComputeBarycentric(Vector2D a, Vector2D b, Vector2D c, Vector2D xy, Vector3D &uvw);
	void computeNormalTri(Face &tri, Vector3D& nrml, Face surf);
	void cleanup();
	void adjustVtxEdges(Face &f);
private:
	//member variables
	MKCore* mk_core;

	Face source;
	Face target;
    iBase_EntityHandle volume;
	vector<Face> link;
	vector<Edge> edges;
	vector<Vertex> vertices;
	iGeom *igeom_instance;
	iMesh *imesh_instance;
	iRel *irel_instance;
	iRel::PairHandle *irel_pair;
	iBase_TagHandle global_geom_tag;
	iBase_TagHandle global_mesh_tag;
	iBase_TagHandle harmonic_surf_pro, facet_mesh_tag;
	iGeom::Error g_err;
	iMesh::Error m_err;
	iRel::Error r_err;
	vector<Vertex> src_facet_v;
	vector<Vertex> tgt_facet_v;
	vector<Face> src_facet_tri;
	vector<Face> tgt_facet_tri;
	vector<Edge> src_facet_e;
	vector<Edge> tgt_facet_e;
	vector<set<int> > adj_src, adj_tgt;
	int size_src_v, size_src_e, size_src_f;
	int size_tgt_v, size_tgt_e, size_tgt_f;
	//map between geometry and facet
	std::map<int, int> src_facet_geom, tgt_facet_geom;
	std::vector<std::list<int> > src_geom_facet, tgt_geom_facet;

	vector<Vertex> quad_mesh_src, quad_mesh_tgt;
	
	vector<Face> facelist;
	
};

}
#endif

