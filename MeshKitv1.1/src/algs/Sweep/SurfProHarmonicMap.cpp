#include "SurfProHarmonicMap.hpp"
#include <iostream>
#include <math.h>
#include <map>

namespace MeshKit {

SurfProHarmonicMap::SurfProHarmonicMap(MKCore* core, iBase_EntityHandle s, iBase_EntityHandle t, iBase_EntityHandle v)
{
  mk_core = core;
  source.gFaceHandle = s;
  target.gFaceHandle = t;
  volume = v;
  
  irel_pair = mk_core->irel_pair();
  irel_instance = mk_core->irel_instance();
  igeom_instance = mk_core->igeom_instance();
  imesh_instance = mk_core->imesh_instance();

  g_err = igeom_instance->getTagHandle("GLOBAL_ID", global_geom_tag);
  IBERRCHK(g_err, "Trouble get the geom_id_tag for 'GLOBAL_ID'.");
  m_err = imesh_instance->getTagHandle("GLOBAL_ID", global_mesh_tag);
  IBERRCHK(m_err, "Trouble get the mesh_id_tag for 'GLOBAL_ID'.");
  g_err = igeom_instance->createTag("HARMONIC_SURF_PROJECTION", 1, iBase_INTEGER, harmonic_surf_pro);
  IBERRCHK(g_err, "Trouble create the tag handle.");
  m_err = imesh_instance->createTag("HARMONIC_FACET_MESH", 1, iBase_INTEGER, facet_mesh_tag);
  IBERRCHK(m_err, "Trouble create the tag handle.");
  preprocessing();
  getFacets();
}

void SurfProHarmonicMap::projection()
{
	//compute the normal vector
	vector<Vector3D> normal_src, normal_tgt;
	normal_src.resize(src_facet_tri.size());
	normal_tgt.resize(tgt_facet_tri.size());
	for (unsigned int i = 0; i < src_facet_tri.size(); i++){
		computeNormalTri(src_facet_tri[i], normal_src[i], source);
		//std::cout << "source index = " << i << "\tnrml = {" << normal_src[i][0] << "," << normal_src[i][1] << "," << normal_src[i][2] << "}\n";
	}
	for (unsigned int i = 0; i < tgt_facet_tri.size(); i++){
		computeNormalTri(tgt_facet_tri[i], normal_tgt[i], target);
		//std::cout << "target index = " << i << "\tnrml = {" << normal_tgt[i][0] << "," << normal_tgt[i][1] << "," << normal_tgt[i][2] << "}\n";
	}

	//loop over the interior nodes of quad mesh on the source surface and compute barycentric coordinates in the physical source surface
	for (vector<Vertex>::iterator it = quad_mesh_src.begin(); it != quad_mesh_src.end(); it++){
		if (it->onBoundary)
			continue;
        //std::cout << "index=" << it->index  << "\tsrc quad node xyz = {" << it->xyz[0] << "," << it->xyz[1] << "," << it->xyz[2] << "}\t";
        Vector3D bary_coord_src, bary_coord_tgt;
		int tri_index = findFacetTri(src_facet_tri, normal_src, it->xyz, bary_coord_src);
		//compute the positions of all the quad mesh nodes on the unit disk
		it->uv[0] = 0.0;
		it->uv[1] = 0.0;
		for (int i = 0; i < 3; i++){
			it->uv[0] += bary_coord_src[i]*src_facet_tri[tri_index].connect[i]->uv[0];
			it->uv[1] += bary_coord_src[i]*src_facet_tri[tri_index].connect[i]->uv[1];
		}
        //std::cout << "\t\tsrc uv={" << it->uv[0] << "," << it->uv[1] << "}  facet index=" << tri_index << "\tv1={" << src_facet_tri[tri_index].connect[0]->xyz[0] << "," << src_facet_tri[tri_index].connect[0]->xyz[1] << "," << src_facet_tri[tri_index].connect[0]->xyz[2] << "} v2={" << src_facet_tri[tri_index].connect[1]->xyz[0] << "," << src_facet_tri[tri_index].connect[1]->xyz[1] << "," << src_facet_tri[tri_index].connect[1]->xyz[2] << "} v3={" << src_facet_tri[tri_index].connect[2]->xyz[0] << "," << src_facet_tri[tri_index].connect[2]->xyz[1] << "," << src_facet_tri[tri_index].connect[2]->xyz[2] << "}\n"; 
		//compute the barycentric coordinates of those quad mesh nodes on the unit disk of target surface
		tri_index = findFacetTri(tgt_facet_tri, it->uv, bary_coord_tgt);
		Vector3D xyz(0.0);	
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				xyz[j] += bary_coord_tgt[i]*tgt_facet_tri[tri_index].connect[i]->xyz[j];
		

		g_err = igeom_instance->getEntClosestPtTrimmed(target.gFaceHandle, xyz[0], xyz[1], xyz[2], quad_mesh_tgt[it->index].xyz[0], quad_mesh_tgt[it->index].xyz[1], quad_mesh_tgt[it->index].xyz[2]);
		IBERRCHK(g_err, "Trouble get the trimmed closed point positions on the target surface!");

         //std::cout << "\t\ttgt uvw={" << bary_coord_tgt[0] << "," << bary_coord_tgt[1] << ","  << bary_coord_tgt[2] << "}  facet index=" << tri_index << "\tv1={" << tgt_facet_tri[tri_index].connect[0]->xyz[0] << "," << tgt_facet_tri[tri_index].connect[0]->xyz[1] << "," << tgt_facet_tri[tri_index].connect[0]->xyz[2] << "} v2={" << tgt_facet_tri[tri_index].connect[1]->xyz[0] << "," << tgt_facet_tri[tri_index].connect[1]->xyz[1] << "," << tgt_facet_tri[tri_index].connect[1]->xyz[2] << "} v3={" << tgt_facet_tri[tri_index].connect[2]->xyz[0] << "," << tgt_facet_tri[tri_index].connect[2]->xyz[1] << "," << tgt_facet_tri[tri_index].connect[2]->xyz[2] << "}\n";
         //std::cout << "\t\tprojected pts xyz = {" << quad_mesh_tgt[it->index].xyz[0] << "," << quad_mesh_tgt[it->index].xyz[1] << "," << quad_mesh_tgt[it->index].xyz[2] << "}\ttmp = {" << xyz[0] << "," << xyz[1] << "," << xyz[2] << "}\n";
	}
	/*
	vector<iBase_EntityHandle> nodes(quad_mesh_tgt.size());
	for (unsigned int i = 0; i < quad_mesh_src.size(); i++){
		if (quad_mesh_src[i].onBoundary)
			continue;
		//create vtx
		m_err = imesh_instance->createVtx(quad_mesh_tgt[i].xyz[0], quad_mesh_tgt[i].xyz[1], quad_mesh_tgt[i].xyz[2], nodes[i]);
		IBERRCHK(m_err, "Trouble create vtx ent!");
	}

	
	for (unsigned int i = 0; i < facelist.size(); i++){
		bool is_boundary = true;
		for (int j = 0; j < 4; j++)
			is_boundary = is_boundary && !facelist[i].connect[j]->onBoundary;
		if (is_boundary){
			iBase_EntityHandle quads;
			vector<iBase_EntityHandle> tmp;
			tmp.push_back(nodes[facelist[i].connect[0]->index]);
			tmp.push_back(nodes[facelist[i].connect[1]->index]);
			tmp.push_back(nodes[facelist[i].connect[2]->index]);
			tmp.push_back(nodes[facelist[i].connect[3]->index]);
			m_err = imesh_instance->createEnt(iMesh_QUADRILATERAL, &tmp[0], 4, quads);
			IBERRCHK(m_err, "Trouble create an face entity!");						
		}
	}
	*/
	cleanup();
}

bool SurfProHarmonicMap::ComputeBarycentric(Vector3D a, Vector3D b, Vector3D c, Vector3D xyz, Vector3D &uvw)
{
	//Compute vectors;
	Vector3D v0 = b - a;
	Vector3D v1 = c - a;
	Vector3D v2 = xyz - a;
	//compute dot products
	double dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2];
	double dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
	double dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2];
	double dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
	double dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	//compute barycentric coordinates
	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	uvw[1] = (dot11 * dot02 - dot01 * dot12) * invDenom;
	uvw[2] = (dot00 * dot12 - dot01 * dot02) * invDenom;
	uvw[0] = 1.0 - uvw[1] - uvw[2];

	//check if the point is in the triangle
	return (uvw[0]>=0)&&(uvw[1]>=0)&&(uvw[0]+uvw[1]<=1);
}

void SurfProHarmonicMap::test()
{
 	vector<iBase_EntityHandle> test_nodes(src_facet_v.size()), test_edges(src_facet_e.size()), test_tri(src_facet_tri.size());
	for (vector<Vertex>::iterator it = src_facet_v.begin(); it != src_facet_v.end(); it++){
		m_err = imesh_instance->createVtx(it->uv[0], it->uv[1], 200.0, test_nodes[it->index]);
		IBERRCHK(m_err, "Trouble create a test vertex!");
	}
	for (vector<Face>::iterator it = src_facet_tri.begin(); it != src_facet_tri.end(); it++){
		vector<iBase_EntityHandle> tmp;
		tmp.push_back(test_nodes[it->connect[0]->index]);
		tmp.push_back(test_nodes[it->connect[1]->index]);
		tmp.push_back(test_nodes[it->connect[2]->index]);
		m_err = imesh_instance->createEnt(iMesh_TRIANGLE, &tmp[0], tmp.size(), test_tri[it->index]);
		IBERRCHK(m_err, "Trouble create a triangle mesh entity!");
	}
	test_nodes.clear();
	test_edges.clear();
	test_tri.clear();
	test_nodes.resize(tgt_facet_v.size());
	test_tri.resize(tgt_facet_tri.size());

	for (vector<Vertex>::iterator it = tgt_facet_v.begin(); it != tgt_facet_v.end(); it++){
		m_err = imesh_instance->createVtx(it->uv[0], it->uv[1], -200.0, test_nodes[it->index]);
		IBERRCHK(m_err, "Trouble create a test vertex!");
		//std::cout << "unit disk vertex index = " << it->index << "\tuv = {" << it->uv[0] << ", " << it->uv[1] << "}\tis_boundary = " << it->onBoundary << std::endl;
	}
	for (vector<Face>::iterator it = tgt_facet_tri.begin(); it != tgt_facet_tri.end(); it++){
		vector<iBase_EntityHandle> tmp;
		tmp.push_back(test_nodes[it->connect[0]->index]);
		tmp.push_back(test_nodes[it->connect[1]->index]);
		tmp.push_back(test_nodes[it->connect[2]->index]);
		m_err = imesh_instance->createEnt(iMesh_TRIANGLE, &tmp[0], tmp.size(), test_tri[it->index]);
		IBERRCHK(m_err, "Trouble create a triangle mesh entity!");
	}	
}

void SurfProHarmonicMap::cleanup(){
	vector<iBase_EntityHandle> ents_to_del;
	for (vector<Edge>::iterator it = src_facet_e.begin(); it != src_facet_e.end(); it++)
		ents_to_del.push_back(it->gEdgeHandle);
	for (vector<Edge>::iterator it = tgt_facet_e.begin(); it != tgt_facet_e.end(); it++)
		ents_to_del.push_back(it->gEdgeHandle);
	for (vector<Face>::iterator it = src_facet_tri.begin(); it != src_facet_tri.end(); it++)
		ents_to_del.push_back(it->gFaceHandle);
	for (vector<Face>::iterator it = tgt_facet_tri.begin(); it != tgt_facet_tri.end(); it++)
		ents_to_del.push_back(it->gFaceHandle);
	for (vector<Vertex>::iterator it = src_facet_v.begin(); it != src_facet_v.end(); it++)
		ents_to_del.push_back(it->gVertexHandle);
	for (vector<Vertex>::iterator it = tgt_facet_v.begin(); it != tgt_facet_v.end(); it++)
		ents_to_del.push_back(it->gVertexHandle);
	m_err = imesh_instance->rmvArrTag(&ents_to_del[0], ents_to_del.size(), facet_mesh_tag);
	IBERRCHK(m_err, "Trouble remove the tag data!");
	m_err = imesh_instance->deleteEntArr(&ents_to_del[0], ents_to_del.size());
	IBERRCHK(m_err, "Trouble delete entities!");
	g_err = igeom_instance->destroyTag(harmonic_surf_pro, true);
	IBERRCHK(g_err, "Trouble delete a tag!");
	m_err = imesh_instance->destroyTag(facet_mesh_tag, true);
	IBERRCHK(g_err, "Trouble delete a tag!");
		
}

int SurfProHarmonicMap::findFacetTri(vector<Face> &facet_tri, Vector2D uv, Vector3D &uvw){
	for (unsigned int i = 0; i < facet_tri.size(); i++){
		if (ComputeBarycentric(facet_tri[i].connect[0]->uv, facet_tri[i].connect[1]->uv, facet_tri[i].connect[2]->uv, uv, uvw))
			return i;
	}

	return -1;
}

bool SurfProHarmonicMap::ComputeBarycentric(Vector2D a, Vector2D b, Vector2D c, Vector2D xy, Vector3D &uvw)
{
	//Compute vectors;
	Vector2D v0 = b - a;
	Vector2D v1 = c - a;
	Vector2D v2 = xy - a;
	//compute dot products
	double dot00 = v0[0]*v0[0] + v0[1]*v0[1];
	double dot01 = v0[0]*v1[0] + v0[1]*v1[1];
	double dot02 = v0[0]*v2[0] + v0[1]*v2[1];
	double dot11 = v1[0]*v1[0] + v1[1]*v1[1];
	double dot12 = v1[0]*v2[0] + v1[1]*v2[1];
	//compute barycentric coordinates
	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	uvw[1] = (dot11 * dot02 - dot01 * dot12) * invDenom;
	uvw[2] = (dot00 * dot12 - dot01 * dot02) * invDenom;
	uvw[0] = 1.0 - uvw[1] - uvw[2];

	//check if the point is in the triangle
	return (uvw[0]>=0)&&(uvw[1]>=0)&&(uvw[0]+uvw[1]<=1);
}
//compute the normal of a triangle
void SurfProHarmonicMap::computeNormalTri(Face &tri, Vector3D& nrml, Face surf)
{
	Vector3D a = tri.connect[0]->xyz, b = tri.connect[1]->xyz, c = tri.connect[2]->xyz;
	Vector3D va = b - a, vb = c - b;
	nrml[0] = va[1]*vb[2] - va[2]*vb[1];
	nrml[1] = va[2]*vb[0] - va[0]*vb[2];
	nrml[2] = va[0]*vb[1] - va[1]*vb[0];
	double len = sqrt(nrml[0]*nrml[0] + nrml[1]*nrml[1] + nrml[2]*nrml[2]);
	for (int i = 0; i < 3; i++)
		nrml[i] /= len;

	Vector3D surf_nrml(0.0);
	g_err = igeom_instance->getEntNrmlXYZ(surf.gFaceHandle, tri.connect[0]->xyz[0], tri.connect[0]->xyz[1], tri.connect[0]->xyz[2], surf_nrml[0], surf_nrml[1], surf_nrml[2]);
	IBERRCHK(g_err, "Trouble get the surface normal at a given position!");
	double dotproduct = surf_nrml[0]*nrml[0] + surf_nrml[1]*nrml[1] + surf_nrml[2]*nrml[2];
	if (dotproduct < 0){
		for (int j = 0; j < 3; j++)
			nrml[j] = -1.0*nrml[j];
		Vertex *tmp_vtx = tri.connect[1];
		tri.connect[1] = tri.connect[2];
		tri.connect[2] = tmp_vtx;
	}
}
//find a specific triangle where the point is located
int SurfProHarmonicMap::findFacetTri(vector<Face> &facet_tri, vector<Vector3D> nrml, Vector3D xyz, Vector3D &uvw)
{
	for (unsigned int i = 0; i < facet_tri.size(); i++){
		Vector3D prj_pts;
		prjPtsToTri(facet_tri[i], xyz, nrml[i], prj_pts);
		
		//bool test_bool = ComputeBarycentric(facet_tri[i].connect[0]->xyz, facet_tri[i].connect[1]->xyz, facet_tri[i].connect[2]->xyz, prj_pts, uvw);
		//std::cout << "\t\t\tindex = " << i << "\tbool = "<< test_bool << "\tbary uvw = {" << uvw[0] << "," << uvw[1] << "," << uvw[2] << "}\n";

		if (ComputeBarycentric(facet_tri[i].connect[0]->xyz, facet_tri[i].connect[1]->xyz, facet_tri[i].connect[2]->xyz, prj_pts, uvw)){
			//std::cout << "projected pts = {" << prj_pts[0] << "," << prj_pts[1] << "," << prj_pts[2] << "}\n";
			//double x = uvw[0]*facet_tri[i].connect[0]->xyz[0]+uvw[1]*facet_tri[i].connect[1]->xyz[0]+uvw[2]*facet_tri[i].connect[2]->xyz[0];
			//double y = uvw[0]*facet_tri[i].connect[0]->xyz[1]+uvw[1]*facet_tri[i].connect[1]->xyz[1]+uvw[2]*facet_tri[i].connect[2]->xyz[1];
			//double z = uvw[0]*facet_tri[i].connect[0]->xyz[2]+uvw[1]*facet_tri[i].connect[1]->xyz[2]+uvw[2]*facet_tri[i].connect[2]->xyz[2];
			//std::cout << "\t\t\tx = " << x << "  y = " << y << "  z = " << z << std::endl;
			return i;

		}
	}

	return -1;
}
//project the point onto a plane determined by triangle
void SurfProHarmonicMap::prjPtsToTri(Face tri, Vector3D pts, Vector3D nrml, Vector3D &xyz)
{
	Vector3D origin = tri.connect[0]->xyz;
	Vector3D v = pts - origin;
	double dist = v[0]*nrml[0] + v[1]*nrml[1] + v[2]*nrml[2];
	xyz = origin + v - dist*nrml;
}

void SurfProHarmonicMap::setMeshData(vector<Vertex> &s, vector<Vertex> &t, vector<Face> &f){
	quad_mesh_src.insert(quad_mesh_src.begin(), s.begin(), s.end());
	quad_mesh_tgt.insert(quad_mesh_tgt.begin(), t.begin(), t.end());
	facelist.insert(facelist.begin(), f.begin(), f.end());	
}

void SurfProHarmonicMap::getMeshData(vector<Vertex> &v){
	int index = -1;
	for (vector<Vertex>::iterator it = quad_mesh_tgt.begin(); it != quad_mesh_tgt.end(); it++){
		index++;
		if (it->onBoundary)
			continue;
		v[index].xyz[0] = it->xyz[0];
		v[index].xyz[1] = it->xyz[1];
		v[index].xyz[2] = it->xyz[2];
	}

}

void SurfProHarmonicMap::match()
{
	//outmost boundary will be distributed on the unit disk
	boundaryDistribution();
	ComputeWeight();
	//extra processing of interior loops

    HarmonicMapper hm_src(mk_core, src_facet_v, src_facet_tri, src_facet_e, adj_src);
	hm_src.execute();
	hm_src.getUV(src_facet_v);
	
	LocateBoundaryNodesTarget();

	HarmonicMapper hm_tgt(mk_core, tgt_facet_v, tgt_facet_tri, tgt_facet_e, adj_tgt);
	hm_tgt.execute();
	hm_tgt.getUV(tgt_facet_v);

	//test();
}

void SurfProHarmonicMap::LocateBoundaryNodesTarget()
{
	std::set<int> set_edges, set_vertices;
	for (unsigned int i = 0; i < target.connEdges.size(); i++)
		set_edges.insert(target.connEdges[i]->index);
	for (unsigned int i = 0; i < source.connect.size(); i++)
		set_vertices.insert(source.connect[i]->index);

    //match the vertices between the source and target surface.
    std::map<int,int> tgt_src_v, tgt_src_e;
    int count = -1;
    for (unsigned int mm = 0; mm < target.vertexloops.size(); mm++){
        for (unsigned int i = 0; i < target.vertexloops[mm].size(); i++){
		    count++;
            tgt_src_v[count] = -1;
		    bool is_found = false;
		    int edge_index = -1;
		    int pre_vertex = target.vertexloops[mm][i];
		    int next_vertex = -1;
		    while(!is_found){
			    vector<iBase_EntityHandle> adj;
			    //get the adjacent edges of a vertex
			    g_err = igeom_instance->getEntAdj(vertices[pre_vertex].gVertexHandle, iBase_EDGE, adj);
			    IBERRCHK(g_err, "Trouble get the adjacent edges around a vertex");
			    //find the edge on the linking surface			
			    if (edge_index != -1){
				    //get edges on the linking surfaces
				    set_edges.clear();
				    std::vector<iBase_EntityHandle> adj_faces;
				    g_err = igeom_instance->getEntAdj(edges[edge_index].gEdgeHandle, iBase_FACE, adj_faces);
				    IBERRCHK(g_err, "Trouble get the adjacent faces around an edge");
				    for (std::vector<iBase_EntityHandle>::iterator it = adj_faces.begin(); it != adj_faces.end(); it++){
					    vector<iBase_EntityHandle> adj_edges;
					    g_err = igeom_instance->getEntAdj(*it, iBase_EDGE, adj_edges);
					    IBERRCHK(g_err, "Trouble get the adjacent edges around a face!");
					    for (vector<iBase_EntityHandle>::iterator it_e = adj_edges.begin(); it_e != adj_edges.end(); it_e++){
						    int intdata = -1;
						    g_err = igeom_instance->getIntData(*it_e, harmonic_surf_pro, intdata);
						    IBERRCHK(g_err, "Trouble get the int data for an edge");
						    set_edges.insert(intdata);
					    }
				    }
			    }
			    //proceed to the next linking edge
			    is_found = true;
			    for (vector<iBase_EntityHandle>::iterator it = adj.begin(); it != adj.end(); it++){
				    g_err = igeom_instance->getIntData(*it, harmonic_surf_pro, edge_index);
				    IBERRCHK(g_err, "Trouble get the adjacent edges of a vertex!");
				    if (std::find(set_edges.begin(), set_edges.end(), edge_index)==set_edges.end()){//this is an edge on the linking surface.
					    is_found = is_found && false;
					    break;
				    }
				    else{
					    is_found = is_found && true;
					    continue;
				    }
			    }
			    //find the next vertex connected by the linking edge
			    if (edges[edge_index].connect[0]->index == pre_vertex)
				    next_vertex = edges[edge_index].connect[1]->index;
			    else
				    next_vertex = edges[edge_index].connect[0]->index;
			    //check whether next vertex is on the source surface or not
			    if (is_found)
				    break;
		    }
		    tgt_src_v[target.vertexloops[mm][i]] = next_vertex;
        }
    }
	//done with the mapping of vertices between the source and target surfaces;
    //test the vertex mapping between source and target
    count = -1;
    for (unsigned int mm = 0; mm < target.vertexloops.size(); mm++){
        for (unsigned int i = 0; i < target.vertexloops[mm].size(); i++){
			int tgt_edge_index = target.edgeloops[mm][i];
			int vtx_a = tgt_src_v[edges[tgt_edge_index].connect[0]->index], vtx_b = tgt_src_v[edges[tgt_edge_index].connect[1]->index];
			for (unsigned int j = 0; j < source.connEdges.size(); j++){
				int vtx_c = (source.connEdges[j]->connect[0])->index, vtx_d = (source.connEdges[j]->connect[1])->index;
				if (((vtx_c == vtx_a)&&(vtx_d == vtx_b))||((vtx_c == vtx_b)&&(vtx_d == vtx_a))){
					tgt_src_e[tgt_edge_index] = source.connEdges[j]->index;
					break;
				}
			} 
        }
    }
    
    count = -1;
    map<int, int> src_vtx_facet;
    vector<vector<int> > src_nodes_list;
    vector<vector<double> > src_nodes_u;
    src_nodes_list.resize(edges.size());
	src_nodes_u.resize(edges.size());
    for (unsigned int i = 0; i < source.edgeloops.size(); i++){

        for (unsigned int j = 0; j < source.edgeloops[i].size(); j++){
           count++;
		   double u;
           for (list<int>::iterator it = src_geom_facet[count+source.connect.size()].begin(); it != src_geom_facet[count+source.connect.size()].end(); it++){
                src_nodes_list[source.edgeloops[i][j]].push_back(*it);
                g_err = igeom_instance->getEntXYZtoU(edges[source.edgeloops[i][j]].gEdgeHandle, src_facet_v[*it].xyz[0], src_facet_v[*it].xyz[1], src_facet_v[*it].xyz[2], u);
                IBERRCHK(g_err, "Trouble get the u coordinates!");
                src_nodes_u[source.edgeloops[i][j]].push_back(u);
           }
        }
    }
	count = -1;
	for (unsigned int i = 0; i < source.vertexloops.size(); i++){
		for (unsigned int j = 0; j < source.vertexloops[i].size(); j++){
			count++;
			src_vtx_facet[source.vertexloops[i][j]] = *(src_geom_facet[count].begin());

		}
	}
    //match the boundary facet nodes between source and target
	count = -1;
    for (unsigned int i = 0; i < target.vertexloops.size(); i++){
        //loop i
		vector<iBase_EntityHandle> tmp_src, tmp_tgt;
		vector<double> dist_source, dist_target;
		for (unsigned int j = 0; j < target.edgeloops[i].size(); j++){
			tmp_src.push_back(edges[target.edgeloops[i][j]].gEdgeHandle);
			tmp_tgt.push_back(edges[tgt_src_e[target.edgeloops[i][j]]].gEdgeHandle);
		}
		dist_source.resize(tmp_src.size());
		g_err = igeom_instance->measure(&tmp_src[0], tmp_src.size(), &dist_source[0]);
		IBERRCHK(g_err, "Trouble get the measure of geometric edges!");
		
		dist_target.resize(tmp_tgt.size());
		g_err = igeom_instance->measure(&tmp_tgt[0], tmp_tgt.size(), &dist_target[0]);
		IBERRCHK(g_err, "Trouble get the measure of geometric edges!");     
		for (unsigned int j = 0; j < target.vertexloops[i].size(); j++){
			count++;            
			int current_tgt_vtx = target.vertexloops[i][j], current_src_vtx = tgt_src_v[current_tgt_vtx];
            int next_tgt_vtx = target.vertexloops[i][(j+1+target.vertexloops[i].size())%target.vertexloops[i].size()], next_src_vtx = tgt_src_v[next_tgt_vtx];
            int tgt_edge_index = target.edgeloops[i][j];
            int src_edge_index = tgt_src_e[tgt_edge_index];
            double u1, u2, u3, u4, u;
            g_err = igeom_instance->getEntXYZtoU(edges[tgt_edge_index].gEdgeHandle, vertices[current_tgt_vtx].xyz[0], vertices[current_tgt_vtx].xyz[1], vertices[current_tgt_vtx].xyz[2], u1);
            IBERRCHK(g_err, "Trouble get the parametric coordinate of a vertex on an edge!");
            g_err = igeom_instance->getEntXYZtoU(edges[tgt_edge_index].gEdgeHandle, vertices[next_tgt_vtx].xyz[0], vertices[next_tgt_vtx].xyz[1], vertices[next_tgt_vtx].xyz[2], u2);
            IBERRCHK(g_err, "Trouble get the parametric coordinate of a vertex on an edge!");
            g_err = igeom_instance->getEntXYZtoU(edges[src_edge_index].gEdgeHandle, vertices[current_src_vtx].xyz[0], vertices[current_src_vtx].xyz[1], vertices[current_src_vtx].xyz[2], u3);
            IBERRCHK(g_err, "Trouble get the parametric coordinate of a vertex on an edge!");
            g_err = igeom_instance->getEntXYZtoU(edges[src_edge_index].gEdgeHandle, vertices[next_src_vtx].xyz[0], vertices[next_src_vtx].xyz[1], vertices[next_src_vtx].xyz[2], u4);
            IBERRCHK(g_err, "Trouble get the parametric coordinate of a vertex on an edge!");
			int count_pts = 0;
			int tgt_vtx_index = *(tgt_geom_facet[count].begin());
			tgt_facet_v[tgt_vtx_index].onBoundary = true;
			tgt_facet_v[tgt_vtx_index].uv[0] = src_facet_v[src_vtx_facet[current_src_vtx]].uv[0];
			tgt_facet_v[tgt_vtx_index].uv[1] = src_facet_v[src_vtx_facet[current_src_vtx]].uv[1];

			std::cout << "--------\nu = {" << u3 << " ";
			for (int m = src_nodes_u[src_edge_index].size()-1; m >= 0; m--)
				std::cout << "[" << m << "]" << src_nodes_u[src_edge_index][m] << " ";
			std::cout << u4 << "}\n--------\n";
			for (list<int>::iterator it = tgt_geom_facet[count+target.connect.size()].begin(); it != tgt_geom_facet[count+target.connect.size()].end(); it++){
				count_pts++;
                g_err = igeom_instance->getEntXYZtoU(edges[tgt_edge_index].gEdgeHandle, tgt_facet_v[*it].xyz[0], tgt_facet_v[*it].xyz[1], tgt_facet_v[*it].xyz[2], u);
                IBERRCHK(g_err, "Trouble get the parametric coordinate of a vertex on an edge!");
                double u_src = (u4 - u3)*(u - u1)/(u2 - u1) + u3;
	
                tgt_facet_v[*it].onBoundary = true;
                //how to compute the corresponding uv coordinates of each nodes on the source H_1
                //loop over facet nodes
                //cases to be discussed: (1) no facet nodes (2) last facet nodes (3) first facet nodes
				int pre_facet_v = -1, next_facet_v = -1;
				double pre_facet_u = 0.0, next_facet_u = 0.0;
				if (src_nodes_u[src_edge_index].size() == 0){//no facet nodes
					//it will be a straight line between current_src_vtx and next_src_vtx
					pre_facet_v = src_vtx_facet[current_src_vtx];
					next_facet_v = src_vtx_facet[next_src_vtx];
					pre_facet_u = u3;
					next_facet_u = u4;
				}
				else{
					for (int m = src_nodes_u[src_edge_index].size()-1; m >= 0; m--){
						if (u4 > u3){
							if (m == 0)//it is located before the first item
								if (u_src >= src_nodes_u[src_edge_index][0]){
									next_facet_v = src_vtx_facet[next_src_vtx];
									next_facet_u = u4;
									pre_facet_v = src_nodes_list[src_edge_index][0];
									pre_facet_u = src_nodes_u[src_edge_index][0];
									break;		
								}
							if (m == (src_nodes_u[src_edge_index].size()-1)){//it is located after the last item
								int tmp_index = src_nodes_list[src_edge_index].size() -1;
								if (u_src <= src_nodes_u[src_edge_index][tmp_index]){
									next_facet_v = src_nodes_list[src_edge_index][tmp_index];
									next_facet_u = src_nodes_u[src_edge_index][tmp_index];
									pre_facet_v = src_vtx_facet[current_src_vtx];
									pre_facet_u = u3;
									break;
								}
							}				
							if (src_nodes_u[src_edge_index][m] >= u_src){//it is located between items
								if ((m+1)==src_nodes_list[src_edge_index].size()){
									pre_facet_v = src_vtx_facet[current_src_vtx];
									pre_facet_u = u3;
								}
								else{
									pre_facet_v = src_nodes_list[src_edge_index][m+1];
									pre_facet_u = src_nodes_u[src_edge_index][m+1];
								}
								
								next_facet_v = src_nodes_list[src_edge_index][m];
								next_facet_u = src_nodes_u[src_edge_index][m];
								break;
							}
							
						}
						else{
							if (m == 0)//it is located before the first item
								if (u_src <= src_nodes_u[src_edge_index][0]){
									next_facet_v = src_vtx_facet[next_src_vtx];
									next_facet_u = u4;
									pre_facet_v = src_nodes_list[src_edge_index][0];
									pre_facet_u = src_nodes_u[src_edge_index][0];
									break;		
								}
							if (m == (src_nodes_u[src_edge_index].size()-1)){//it is located after the last item
								int tmp_index = src_nodes_list[src_edge_index].size() -1;
								if (u_src >= src_nodes_u[src_edge_index][tmp_index]){
									next_facet_v = src_nodes_list[src_edge_index][tmp_index];
									next_facet_u = src_nodes_u[src_edge_index][tmp_index];
									pre_facet_v = src_vtx_facet[current_src_vtx];
									pre_facet_u = u3;
									break;
								}
							}
							if (src_nodes_u[src_edge_index][m] <= u_src){
								next_facet_v = src_nodes_list[src_edge_index][m];
								next_facet_u = src_nodes_u[src_edge_index][m];
								if ((m+1)==src_nodes_list[src_edge_index].size()){
									pre_facet_v = src_vtx_facet[current_src_vtx];
									pre_facet_u = u3;
								}
								else{
									pre_facet_v = src_nodes_list[src_edge_index][m+1];
									pre_facet_u = src_nodes_u[src_edge_index][m+1];
								}
								break;
							}
						}
					}
				}
				double alpha = (u_src - pre_facet_u)/(next_facet_u - pre_facet_u);				
				tgt_facet_v[*it].uv[0] = (1.0 - alpha)*src_facet_v[pre_facet_v].uv[0]+alpha*src_facet_v[next_facet_v].uv[0];
				tgt_facet_v[*it].uv[1] = (1.0 - alpha)*src_facet_v[pre_facet_v].uv[1]+alpha*src_facet_v[next_facet_v].uv[1];
			}
        }
    }
}

void SurfProHarmonicMap::ComputeWeight()
{
	//create the facet mesh on the source surface
	int count = -1;
	for (vector<Vertex>::iterator it = src_facet_v.begin(); it != src_facet_v.end(); it++){
		m_err = imesh_instance->createVtx(it->xyz[0], it->xyz[1], it->xyz[2], it->gVertexHandle);
		IBERRCHK(m_err, "Trouble create a mesh vertex!");
		count++;
		m_err = imesh_instance->setIntData(it->gVertexHandle, facet_mesh_tag, count);
		IBERRCHK(m_err, "Trouble set the int data for a facet vertex");
	}
	//create the facet edge entities on the source surface
	vector<set<int> > edge_connect;
	edge_connect.resize(src_facet_v.size());
	count = -1;
	for (vector<Face>::iterator it = src_facet_tri.begin(); it != src_facet_tri.end(); it++){
		int a = it->connect[0]->index, b = it->connect[1]->index, c = it->connect[2]->index;
		if (edge_connect[a].find(b)==edge_connect[a].end())
			addEdgeToList(a, b, count, edge_connect, src_facet_e, src_facet_v);
		if (edge_connect[a].find(c)==edge_connect[a].end())
			addEdgeToList(a, c, count, edge_connect, src_facet_e, src_facet_v);
		if (edge_connect[b].find(c)==edge_connect[b].end())
			addEdgeToList(b, c, count, edge_connect, src_facet_e, src_facet_v);
	}
	size_src_e = src_facet_e.size();
	//create the facet face entities on the source surface
	count = -1;
	for (vector<Face>::iterator it = src_facet_tri.begin(); it != src_facet_tri.end(); it++){
		vector<iBase_EntityHandle> nodes;
		for (unsigned int i = 0; i < it->connect.size(); i++)
			nodes.push_back((it->connect[i])->gVertexHandle);
		m_err = imesh_instance->createEnt(iMesh_TRIANGLE, &nodes[0], nodes.size(), it->gFaceHandle);
		IBERRCHK(m_err, "Trouble create a triangle mesh entity on the source surface!");
		count++;
		m_err = imesh_instance->setIntData(it->gFaceHandle, facet_mesh_tag, count);
		IBERRCHK(m_err, "Trouble set the int data for a facet vertex");
	}
	//add extra stuff in order to seal the interior loops
	addExtra(source, src_facet_v, src_facet_e, src_facet_tri, src_geom_facet, size_src_v);

	//get the adjacent edges of every vertex on the source surface
	adj_src.resize(src_facet_v.size());
	for (vector<Vertex>::iterator it = src_facet_v.begin(); it != src_facet_v.end(); it++){
		vector<iBase_EntityHandle> adj_edges;
		m_err = imesh_instance->getEntAdj(it->gVertexHandle, iBase_EDGE, adj_edges);
		IBERRCHK(m_err, "Trouble get the adjacent edges around a vertex!");
		int intdata = -1;
		for (vector<iBase_EntityHandle>::iterator it_e = adj_edges.begin(); it_e != adj_edges.end(); it_e++){
			m_err = imesh_instance->getIntData(*it_e, facet_mesh_tag, intdata);
			IBERRCHK(m_err, "Trouble get the int data for the adjacent edges!");
			adj_src[it->index].insert(intdata);
		}
	}
	//compute the weight for edges on the source surface
	computeEdgeWeight(src_facet_e, src_facet_tri, src_facet_v);

	//create the facet mesh on the target surface
	count = -1;
	for (vector<Vertex>::iterator it = tgt_facet_v.begin(); it != tgt_facet_v.end(); it++){
		m_err = imesh_instance->createVtx(it->xyz[0], it->xyz[1], it->xyz[2], it->gVertexHandle);
		IBERRCHK(m_err, "Trouble create a mesh vertex");
		count++;
		m_err = imesh_instance->setIntData(it->gVertexHandle, facet_mesh_tag, count);
		IBERRCHK(m_err, "Trouble set the int data for a facet vertex");
	}
	//get the facet edges on the target surface
	edge_connect.clear();
	edge_connect.resize(tgt_facet_v.size());
	count = -1;
	for (vector<Face>::iterator it = tgt_facet_tri.begin(); it != tgt_facet_tri.end(); it++){
		int a = it->connect[0]->index, b = it->connect[1]->index, c = it->connect[2]->index;
		if (edge_connect[a].find(b)==edge_connect[a].end())
			addEdgeToList(a, b, count, edge_connect, tgt_facet_e, tgt_facet_v);
		if (edge_connect[a].find(c)==edge_connect[a].end())
			addEdgeToList(a, c, count, edge_connect, tgt_facet_e, tgt_facet_v);
		if (edge_connect[b].find(c)==edge_connect[b].end())
			addEdgeToList(b, c, count, edge_connect, tgt_facet_e, tgt_facet_v);
	}
	size_tgt_e = tgt_facet_e.size();
	//get the facet faces on the target surface
	count = -1;
	for (vector<Face>::iterator it = tgt_facet_tri.begin(); it != tgt_facet_tri.end(); it++){
		vector<iBase_EntityHandle> nodes;
		for (unsigned int i = 0; i < it->connect.size(); i++)
			nodes.push_back((it->connect[i])->gVertexHandle);
		m_err = imesh_instance->createEnt(iMesh_TRIANGLE, &nodes[0], nodes.size(), it->gFaceHandle);
		IBERRCHK(m_err, "Trouble create a triangle mesh entity on the target surface!");
		count++;
		m_err = imesh_instance->setIntData(it->gFaceHandle, facet_mesh_tag, count);
		IBERRCHK(m_err, "Trouble set the int data for a facet vertex!");
	}
	//add extra stuff in order to seal the interior loops
	//addExtra(target, tgt_facet_v, tgt_facet_e, tgt_facet_tri, tgt_geom_facet, size_tgt_v);
	//get the adjacent edges of every vertex on the target surface
	adj_tgt.resize(tgt_facet_v.size());
	for (vector<Vertex>::iterator it = tgt_facet_v.begin(); it != tgt_facet_v.end(); it++){
		vector<iBase_EntityHandle> adj_edges;
		m_err = imesh_instance->getEntAdj(it->gVertexHandle, iBase_EDGE, adj_edges);
		IBERRCHK(m_err, "Trouble get the adjacent edges around a vertex!");
		int intdata = -1;
		for (vector<iBase_EntityHandle>::iterator it_e = adj_edges.begin(); it_e != adj_edges.end(); it_e++){
			m_err = imesh_instance->getIntData(*it_e, facet_mesh_tag, intdata);
			IBERRCHK(m_err, "Trouble get the int data for the adjacent edges!");
			adj_tgt[it->index].insert(intdata);
		}
	}
	//compute the edge weight on the target surface
	computeEdgeWeight(tgt_facet_e, tgt_facet_tri, tgt_facet_v);

}

void SurfProHarmonicMap::addExtra(Face f, vector<Vertex> &facet_v, vector<Edge> &facet_e, vector<Face> &facet_tri, vector<list<int> > geom_facet, int size_facet_v)
{
	//add extra nodes to the list
	if (f.vertexloops.size() > 1){
		int count = 0;
		for (unsigned int i = 1; i < f.vertexloops.size(); i++){
			double xyz[3] = {0.0, 0.0, 0.0};
			count += f.vertexloops[i-1].size();
			int vtx_center = size_facet_v+i-1;
			vector<int> pointlist;
			for (unsigned int j = 0; j < f.vertexloops[i].size(); j++){
				pointlist.push_back(*(geom_facet[count+j].begin()));
				for (int m = 0; m < 3; m++)
					xyz[m] += facet_v[*(geom_facet[count+j].begin())].xyz[m];
				for (list<int>::iterator it = geom_facet[count+j+f.connect.size()].begin(); it != geom_facet[count+j+f.connect.size()].end(); it++){
					pointlist.push_back(*it);
					for (int m = 0; m < 3; m++)
						xyz[m] += facet_v[*it].xyz[m];
				}
			}
			//done
			//add extra vertex
			for (int m = 0; m < 3; m++)
				 facet_v[vtx_center].xyz[m] = xyz[m]/double(pointlist.size());
			m_err = imesh_instance->createVtx(facet_v[vtx_center].xyz[0], facet_v[vtx_center].xyz[1], facet_v[vtx_center].xyz[2], facet_v[vtx_center].gVertexHandle);
			IBERRCHK(m_err, "Trouble create a facet vertex on the interior loops of source surface!");
			m_err = imesh_instance->setIntData(facet_v[vtx_center].gVertexHandle, facet_mesh_tag, facet_v[vtx_center].index);
			IBERRCHK(m_err, "Trouble set the int data for a facet vertex");
			facet_v[vtx_center].index = vtx_center;
			facet_v[vtx_center].onBoundary = false;
			//add extra facet edges and adjacent edges
			int ent_index = src_facet_e.size();
			vector<iBase_EntityHandle> tmp;			
			for (vector<int>::iterator it = pointlist.begin(); it != pointlist.end(); it++, ent_index++){
				//std::cout << "test ent_index = " << ent_index << std::endl;
				Edge tmp_edge;		
				tmp_edge.index = ent_index;
				tmp_edge.connect.resize(2);
				tmp_edge.connect[0] = &facet_v[*it];
				tmp_edge.connect[1] = &facet_v[vtx_center];
				tmp.clear();
				tmp.push_back(facet_v[*it].gVertexHandle);
				tmp.push_back(facet_v[vtx_center].gVertexHandle);
				m_err = imesh_instance->createEnt(iMesh_LINE_SEGMENT, &tmp[0], 2, tmp_edge.gEdgeHandle);
				IBERRCHK(m_err, "Trouble create a extra facet edge on the interior loops of source surface!");
				m_err = imesh_instance->setIntData(tmp_edge.gEdgeHandle, facet_mesh_tag, tmp_edge.index);
				IBERRCHK(m_err, "Trouble set the int data for a facet vertex");
				
				facet_e.push_back(tmp_edge);
			}
			//add extra facet faces
			ent_index = facet_tri.size();
			//src_facet_tri.resize(ent_index+pointlist.size());
			for (unsigned int j = 0; j < pointlist.size(); j++, ent_index++){
				Face tmp_f;
				int vtx_index[3];
				if (j == (pointlist.size()-1)){
					vtx_index[0] = pointlist[j];
					vtx_index[1] = pointlist[0];
				}
				else{
					vtx_index[0] = pointlist[j];
					vtx_index[1] = pointlist[j+1];
				}
				vtx_index[2] = vtx_center;
				tmp_f.connect.resize(3);
				tmp.clear();
				for (int m = 0; m < 3; m++){
					tmp_f.connect[m] = &facet_v[vtx_index[m]];
					tmp.push_back(facet_v[vtx_index[m]].gVertexHandle);
				}
				tmp_f.index = ent_index;
				m_err = imesh_instance->createEnt(iMesh_TRIANGLE, &tmp[0], tmp.size(), tmp_f.gFaceHandle);
				IBERRCHK(m_err, "Trouble create a extra triangle facet face entity!");
				m_err = imesh_instance->setIntData(tmp_f.gFaceHandle, facet_mesh_tag, tmp_f.index);
				IBERRCHK(m_err, "Trouble set the int data for a facet vertex");
				facet_tri.push_back(tmp_f);
			}			
		}
	}



}

void SurfProHarmonicMap::computeEdgeWeight(vector<Edge> &f_edges, vector<Face> &f, vector<Vertex> f_v)
{
	//compute the weight 
	for (vector<Edge>::iterator it = f_edges.begin(); it != f_edges.end(); it++){
		vector<iBase_EntityHandle> adj_faces;
		m_err = imesh_instance->getEntAdj(it->gEdgeHandle, iBase_FACE, adj_faces);
		IBERRCHK(m_err, "Trouble get the adjacent faces for an facet edge!");		
		assert(adj_faces.size()<=2);
		double weight = 0.0;
		int a = (it->connect[0])->index, b = (it->connect[1])->index;
		for (vector<iBase_EntityHandle>::iterator it_ent = adj_faces.begin(); it_ent != adj_faces.end(); it_ent++){
			int intdata = -1, c = -1;
			m_err = imesh_instance->getIntData(*it_ent, facet_mesh_tag, intdata);
			IBERRCHK(m_err, "Trouble get the int data for the adjacent face!");
			if (((f[intdata].connect[0]->index == a)&&(f[intdata].connect[1]->index==b))||((f[intdata].connect[1]->index == a)&&(f[intdata].connect[0]->index==b)))
				c = f[intdata].connect[2]->index;
			else if (((f[intdata].connect[0]->index == a)&&(f[intdata].connect[2]->index==b))||((f[intdata].connect[2]->index == a)&&(f[intdata].connect[0]->index==b)))
				c = f[intdata].connect[1]->index;
			else if (((f[intdata].connect[1]->index == a)&&(f[intdata].connect[2]->index==b))||((f[intdata].connect[2]->index == a)&&(f[intdata].connect[1]->index==b)))
				c = f[intdata].connect[0]->index;
			else
				continue;
			double vec1[3], vec2[3];
			for (int k = 0; k < 3; k++){
				vec1[k] = f_v[a].xyz[k] - f_v[c].xyz[k];
				vec2[k] = f_v[b].xyz[k] - f_v[c].xyz[k];
			}
			double cos_value = (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2])/sqrt((vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])*(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]));
			//std::cout << "ctan(angle) = " << cos_value/sqrt(1-pow(cos_value,2)) << std::endl;
			weight += cos_value/sqrt(1-pow(cos_value,2));
		}
		it->e = 0.5*weight;
	}

}

void SurfProHarmonicMap::addEdgeToList(int a, int b, int &count, vector<set<int> > &edge_connect, vector<Edge> &f_edges, vector<Vertex> &facet_v)
{
	count++;
	f_edges.resize(count+1);
	edge_connect[a].insert(b);
	edge_connect[b].insert(a);
	f_edges[count].connect.resize(2);
	f_edges[count].connect[0] = &facet_v[a];
	f_edges[count].connect[1] = &facet_v[b];
	f_edges[count].index = count;
	vector<iBase_EntityHandle> nodes;
	nodes.push_back(facet_v[a].gVertexHandle);
	nodes.push_back(facet_v[b].gVertexHandle);
	m_err = imesh_instance->createEnt(iMesh_LINE_SEGMENT, &nodes[0], nodes.size(), f_edges[count].gEdgeHandle);
	IBERRCHK(m_err, "Trouble create a facet line segment!");
	m_err = imesh_instance->setIntData(f_edges[count].gEdgeHandle, facet_mesh_tag, count);
	IBERRCHK(m_err, "Trouble set the int data for a facet edge!");
	
}

void SurfProHarmonicMap::boundaryDistribution()
{
	//loop over the geometry entities on the outmost boundary
	//1. get the total distance of the outmost boundary loop
	vector<iBase_EntityHandle> tmp;
	vector<double> dist_source, dist_target;
	double total_dist_source = 0.0, total_dist_target = 0.0;
	for (unsigned int i = 0; i < source.edgeloops[0].size(); i++)
		tmp.push_back(edges[source.edgeloops[0][i]].gEdgeHandle);
	dist_source.resize(tmp.size());
	g_err = igeom_instance->measure(&tmp[0], tmp.size(), &dist_source[0]);
	IBERRCHK(g_err, "Trouble get the measure of geometric edges!");
	tmp.clear();
	for (unsigned int i = 0; i < target.edgeloops[0].size(); i++)
		tmp.push_back(edges[target.edgeloops[0][i]].gEdgeHandle);
	dist_target.resize(tmp.size());
	g_err = igeom_instance->measure(&tmp[0], tmp.size(), &dist_target[0]);
	IBERRCHK(g_err, "Trouble get the measure of geometric edges!");
	
	//Initialization of u,v coordinates
	for (unsigned int i = 0; i < src_facet_v.size(); i++){
		src_facet_v[i].uv[0] = 0.0;
		src_facet_v[i].uv[1] = 0.0;
	}
	for (unsigned int i = 0; i < tgt_facet_v.size(); i++){
		tgt_facet_v[i].uv[0] = 0.0;
		tgt_facet_v[i].uv[1] = 0.0;
	}		

	//assign u,v coordinates for facet vertices on the outmost boundary onto the unit disk, both source and target surface
	for (vector<Vertex>::iterator it = src_facet_v.begin(); it != src_facet_v.end(); it++)
		it->onBoundary = false;
	for (vector<Vertex>::iterator it = tgt_facet_v.begin(); it != tgt_facet_v.end(); it++)
		it->onBoundary = false;
	boundaryUnitDisk(source, dist_source, src_facet_v, src_geom_facet);
}


void SurfProHarmonicMap::boundaryUnitDisk(Face f, vector<double> dist, vector<Vertex> &facet_v, std::vector<std::list<int> > geom_facet)
{
	double local_dist = 0.0, dist_sum = 0.0, len_curve = 0.0;
	int start_vertex = f.vertexloops[0][0];
	for (std::vector<double>::iterator it = dist.begin(); it != dist.end(); it++)
		dist_sum += *it;
	for (unsigned int i = 0; i < f.vertexloops[0].size(); i++){
		if (i > 0){
			len_curve += dist[i-1];		
			local_dist = len_curve;
		}
		//get the facet on a geometry vertex of a loop	
		facet_v[*(geom_facet[i].begin())].uv[0] = 100.0*cos(2*PI*local_dist/dist_sum);
		facet_v[*(geom_facet[i].begin())].uv[1] = 100.0*sin(2*PI*local_dist/dist_sum);
		facet_v[*(geom_facet[i].begin())].onBoundary = true;
		int pre_v = *(geom_facet[i].begin());
		//get the facet on a geometry edge of a loop
		for (std::list<int>::iterator it = geom_facet[i+f.connect.size()].begin(); it != geom_facet[i+f.connect.size()].end(); it++){
			local_dist += sqrt(pow(facet_v[*it].xyz[0]-facet_v[pre_v].xyz[0],2)+pow(facet_v[*it].xyz[1]-facet_v[pre_v].xyz[1],2)+pow(facet_v[*it].xyz[2]-facet_v[pre_v].xyz[2],2));
			pre_v = *it;
			facet_v[*it].uv[0] = 100.0*cos(2*PI*local_dist/dist_sum);
			facet_v[*it].uv[1] = 100.0*sin(2*PI*local_dist/dist_sum);
			facet_v[*it].onBoundary = true;
		}
	}
}

void SurfProHarmonicMap::getFacets()
{
	SimpleArray<double> src_coords, tgt_coords;
	SimpleArray<int> src_facets, tgt_facets;
	double dist_tolerance = 1.0e-1;
	int err;
	//facets on the source surface
	iGeom_getFacets(igeom_instance->instance(), source.gFaceHandle, dist_tolerance, ARRAY_INOUT(src_coords), ARRAY_INOUT(src_facets), &err);
	assert(!err);
	src_facet_v.resize(src_coords.size()/3+source.vertexloops.size()-1);
	src_facet_tri.resize(src_facets.size()/3);
	size_src_v = src_coords.size()/3;
	size_src_f = src_facets.size()/3;
	
	for (unsigned int i = 0; i < src_coords.size(); i += 3){
		src_facet_v[i/3].index = i/3;
		for (int k = 0; k < 3; k++)
			src_facet_v[i/3].xyz[k] = src_coords[i+k];
	}
	for (unsigned int i = 0; i < src_facets.size(); i += 3){
		src_facet_tri[i/3].index = i/3;
		src_facet_tri[i/3].connect.resize(3);
		for (int k = 0; k < 3; k++)
			src_facet_tri[i/3].connect[k] = &src_facet_v[src_facets[i+k]];
	}
	//facets on the target surface
	iGeom_getFacets(igeom_instance->instance(), target.gFaceHandle, dist_tolerance, ARRAY_INOUT(tgt_coords), ARRAY_INOUT(tgt_facets), &err);
	assert(!err);
	//tgt_facet_v.resize(tgt_coords.size()/3+target.vertexloops.size()-1);
	tgt_facet_v.resize(tgt_coords.size()/3);
	tgt_facet_tri.resize(tgt_facets.size()/3);
	size_tgt_v = tgt_coords.size()/3;
	size_tgt_f = tgt_facet_tri.size();
	for (unsigned int i = 0; i < tgt_coords.size(); i+=3){
		tgt_facet_v[i/3].index = i/3;
		for (int k = 0; k < 3; k++)
			tgt_facet_v[i/3].xyz[k] = tgt_coords[i+k];
		tgt_facet_v[i/3].uv[0] = 0.0;
		tgt_facet_v[i/3].uv[1] = 0.0;
	}
	for (unsigned int i = 0; i < tgt_facets.size(); i += 3){
		tgt_facet_tri[i/3].index = i/3;
		tgt_facet_tri[i/3].connect.resize(3);
		for (int k = 0; k < 3; k++)
			tgt_facet_tri[i/3].connect[k] = &tgt_facet_v[tgt_facets[i+k]];

	}

	//match the facet mesh with geometry(source surface and target surface)
	//map-->vertices, edges, faces on the source surface
	
	std::cout << "starting the source surface mapping\n";
	MapFacetGeom(source, src_facet_v, src_facet_geom, src_geom_facet, size_src_v);
	std::cout << "starting the target surface mapping\n";
	MapFacetGeom(target, tgt_facet_v, tgt_facet_geom, tgt_geom_facet, size_tgt_v);	
}

void SurfProHarmonicMap::MapFacetGeom(Face f_surf, vector<Vertex> facet_node, std::map<int, int> &map_data, vector<list<int> > &geom_facet, int size_facet_v)
{
	geom_facet.resize(f_surf.connect.size()+f_surf.connEdges.size()+1);
	//0===f_surf.connect.size()-1													store vertex info
	//f_surf.connect.size()===f_surf.connect.size()+f_surf.connEdges.size()-1		store edge info
	//f_surf.connect.size()+f_surf.connEdges.size()									store face info	
	for (unsigned int i = 0; i < size_facet_v; i++){
		facet_node[i].index = i;
		int is_on = false;

		map_data[i] = -1;		

		double proj_xyz[3];
		g_err = igeom_instance->getEntClosestPtTrimmed(f_surf.gFaceHandle, facet_node[i].xyz[0], facet_node[i].xyz[1], facet_node[i].xyz[2], proj_xyz[0], proj_xyz[1], proj_xyz[2]);
		IBERRCHK(g_err, "Trouble get the closest point on a geometry entity!");
		facet_node[i].xyz[0] = proj_xyz[0];
		facet_node[i].xyz[1] = proj_xyz[1];
		facet_node[i].xyz[2] = proj_xyz[2];



		//check whether the facet mesh is on vertices or not.
		for (unsigned int j = 0; j < f_surf.connect.size(); j++){			
			iGeom_isPositionOn(igeom_instance->instance(), f_surf.connect[j]->gVertexHandle, facet_node[i].xyz[0], facet_node[i].xyz[1], facet_node[i].xyz[2], &is_on);			
			if (is_on){
				map_data[i] = j;
				geom_facet[j].push_back(i);		
				break;
			}
		}
		if (is_on)
			continue;
		//check whether the facet mesh is on edges or not
		for (unsigned int j = 0; j < f_surf.connEdges.size(); j++){
			iGeom_isPositionOn(igeom_instance->instance(), f_surf.connEdges[j]->gEdgeHandle, facet_node[i].xyz[0], facet_node[i].xyz[1], facet_node[i].xyz[2], &is_on);			
			if (is_on){
				map_data[i] = f_surf.connect.size()+j;
				geom_facet[f_surf.connect.size()+j].push_back(i);				
				break;
			}
		}
		if (is_on)
			continue;
		//check whether the facet mesh is on the surface or not
		iGeom_isPositionOn(igeom_instance->instance(), f_surf.gFaceHandle, facet_node[i].xyz[0], facet_node[i].xyz[1], facet_node[i].xyz[2], &is_on);		
		if (is_on){
			map_data[i] = f_surf.connect.size()+f_surf.connEdges.size();
			geom_facet[f_surf.connect.size()+f_surf.connEdges.size()].push_back(i);
		}
		else{
			std::cout << "Something is wrong for the facet mesh!\txyz = {" << facet_node[i].xyz[0] << "," << facet_node[i].xyz[1] << "," << facet_node[i].xyz[2] << "}\n";
			
		}
	}
}

void SurfProHarmonicMap::preprocessing()
{
	vector<iBase_EntityHandle> v;
	vector<iBase_EntityHandle> e;
	vector<iBase_EntityHandle> f;
	g_err = igeom_instance->getEntAdj(volume, iBase_VERTEX, v);
	IBERRCHK(g_err, "Trouble get the adjacent vertices on a volume!");
	g_err = igeom_instance->getEntAdj(volume, iBase_EDGE, e);
	IBERRCHK(g_err, "Trouble get the adjacent edges on a volume!");
	g_err = igeom_instance->getEntAdj(volume, iBase_FACE, f);
	IBERRCHK(g_err, "Trouble get the adjacent faces on a volume!");
	vertices.resize(v.size());
	for (unsigned int i = 0; i < v.size(); i++){
		vertices[i].gVertexHandle = v[i];
		g_err = igeom_instance->getVtxCoord(v[i], vertices[i].xyz[0], vertices[i].xyz[1], vertices[i].xyz[2]);
		IBERRCHK(g_err, "Trouble get the xyz coordinates!");
		vertices[i].index = i;
		g_err = igeom_instance->getIntData(vertices[i].gVertexHandle, global_geom_tag, vertices[i].id);
		IBERRCHK(g_err, "Trouble get the int data for vertex entity");
		g_err = igeom_instance->setIntData(vertices[i].gVertexHandle, harmonic_surf_pro, i);
		IBERRCHK(g_err, "Trouble set the int data for vertex entity");	
	}
	edges.resize(e.size());
	for (unsigned int i = 0; i < e.size(); i++){
		edges[i].index = i;
		edges[i].gEdgeHandle = e[i];
		g_err = igeom_instance->getIntData(edges[i].gEdgeHandle, global_geom_tag, edges[i].id);
		IBERRCHK(g_err, "Trouble get the int data for geometrical edge entity!");
		vector<iBase_EntityHandle> adjs;
		g_err = igeom_instance->getEntAdj(e[i], iBase_VERTEX, adjs);
		IBERRCHK(g_err, "Trouble get the adjacent vertices of a geometrical edge!");
		for (unsigned int j = 0; j < adjs.size(); j++){
			int tmp;
			g_err = igeom_instance->getIntData(adjs[j], harmonic_surf_pro, tmp);
			IBERRCHK(g_err, "Trouble get the int data for a geometrical vertex!");
			edges[i].connect.push_back(&vertices[tmp]);
		}
		g_err = igeom_instance->setIntData(edges[i].gEdgeHandle, harmonic_surf_pro, i);
		IBERRCHK(g_err, "Trouble set the int data for a geometrical edge!");
	}
	link.resize(f.size());
	for (unsigned int i = 0; i < f.size(); i++){
		if (f[i] == source.gFaceHandle){
			addFaceToList(f[i], source, i, false);
			GetGeomLoops(source, source.vertexloops, source.edgeloops);
			postProcessGeomLoops(source);
			source.src_tgt_link = 0;
			link[i].index = -1;
			continue;
		}			
		if (f[i] == target.gFaceHandle){
			addFaceToList(f[i], target, i, false);
			GetGeomLoops(target, target.vertexloops, target.edgeloops);
			postProcessGeomLoops(target);
			target.src_tgt_link = 1;
			link[i].index = -2;
			continue;
		}
		link[i].src_tgt_link = 2;		
		addFaceToList(f[i], link[i], i, true);

	}
	
		
}

void SurfProHarmonicMap::addFaceToList(iBase_EntityHandle entity, Face& f, int index, bool is_set_int)
{
	f.gFaceHandle = entity;
	if (is_set_int){
		g_err = igeom_instance->setIntData(f.gFaceHandle, harmonic_surf_pro, index);
		IBERRCHK(g_err, "Trouble set the int data for geometrical face entity!");
		f.index = index;
	}
	
	//get the adjacent vertices on a face
	vector<iBase_EntityHandle> tmp;
	g_err = igeom_instance->getEntAdj(f.gFaceHandle, iBase_VERTEX, tmp);
	IBERRCHK(g_err, "Trouble get the adjacent vertices on a geometrical face!");
	for (unsigned int i = 0; i < tmp.size(); i++){
		int tmpint = -1;
		g_err = igeom_instance->getIntData(tmp[i], harmonic_surf_pro, tmpint);
		IBERRCHK(g_err, "Trouble get the int data for a vertex!");
		f.connect.push_back(&vertices[tmpint]);
	}
	//get the adjacent edges on a face
	tmp.clear();
	g_err = igeom_instance->getEntAdj(f.gFaceHandle, iBase_EDGE, tmp);
	IBERRCHK(g_err, "Trouble get the adjacent edges on a geometrical face!");
	for (unsigned int i = 0; i < tmp.size(); i++){
		int tmpint = -1;
		g_err = igeom_instance->getIntData(tmp[i], harmonic_surf_pro, tmpint);
		IBERRCHK(g_err, "Trouble get the int data for a vertex!");
		f.connEdges.push_back(&edges[tmpint]);
	}
	
}

//post process the geometric loops on the source/target surface
//problem not to be solved: matching the loops between source and target surfaces should be done!!!!! 
void SurfProHarmonicMap::postProcessGeomLoops(Face& surf)
{
	//make sure the outmost boundary loop is at [0];
	int outmost_index = -1;
	double volume = -1.0;
	if (surf.edgeloops.size()<=1)
		return;
	for (unsigned int j = 0; j < surf.edgeloops.size(); j++){
		std::vector<iBase_EntityHandle> entities;
		double mincorner[3] = {1.0e10,1.0e10,1.0e10}, maxcorner[3] = {-1.0e10,-1.0e10,-1.0e10};
		for (unsigned int i = 0; i < surf.edgeloops[j].size(); i++){
			double t_min[3], t_max[3];
			entities.clear();
			entities.push_back(edges[surf.edgeloops[j][i]].gEdgeHandle);
			g_err = igeom_instance->getArrBoundBox(&entities[0], entities.size(), iBase_StorageOrder_MIN, &t_min[0], &t_max[0]);
			IBERRCHK(g_err, "Trouble get the bound box for an array of entities");
			if (fabs(mincorner[0])==fabs(mincorner[1])&&fabs(mincorner[1])==fabs(mincorner[2])&&fabs(maxcorner[0])==fabs(maxcorner[1])&&fabs(maxcorner[1])==fabs(mincorner[2])){
				for (int m = 0; m < 3; m++){
					mincorner[m] = t_min[m];
					maxcorner[m] = t_max[m];
				}
			}
			else{
				if (t_min[0] < mincorner[0])
					mincorner[0] = t_min[0];
				if (t_min[1] < mincorner[1])				
					mincorner[1] = t_min[1];
				if (t_min[2] < mincorner[2])					
					mincorner[2] = t_min[2];
				if (t_max[0] > maxcorner[0])
					maxcorner[0] = t_max[0];
				if (t_max[1] > maxcorner[1])
					maxcorner[1] = t_max[1];
				if (t_max[2] > maxcorner[2])
					maxcorner[2] = t_max[2];				
			}			
		}
		//g_err = igeom_instance->getArrBoundBox(&entities[0], entities.size(), iBase_StorageOrder_MIN, &mincorner[0], &maxcorner[0]);
		//IBERRCHK(g_err, "Trouble get the bound box for an array of entities");
		double len_x = fabs(maxcorner[0]-mincorner[0]), len_y = fabs(maxcorner[1]-mincorner[1]), len_z = fabs(maxcorner[2]-mincorner[2]);
		if (fabs(maxcorner[0]-mincorner[0]) < 1.0e-5)
			len_x = 1.0;
		if (fabs(maxcorner[1]-mincorner[1]) < 1.0e-5)
			len_y = 1.0;
		if (fabs(maxcorner[2]-mincorner[2]) < 1.0e-5)
			len_z = 1.0;
		double tmp_volume = len_x*len_y*len_z;
		if ( tmp_volume > volume){
			volume = tmp_volume;
			outmost_index = j;
		}
	}
	//exchange loop[0] and loop[outmost_index]
	std::vector<int> tmp_loop1, tmp_loop2;
	for (unsigned int j = 0; j < surf.vertexloops[0].size(); j++){
		tmp_loop1.push_back(surf.vertexloops[0][j]);
		tmp_loop2.push_back(surf.edgeloops[0][j]);
	}
	surf.vertexloops[0].clear();
	surf.edgeloops[0].clear();
	for (unsigned int j = 0; j < surf.vertexloops[outmost_index].size(); j++){
		surf.vertexloops[0].push_back(surf.vertexloops[outmost_index][j]);
		surf.edgeloops[0].push_back(surf.edgeloops[outmost_index][j]);
	}
	surf.vertexloops[outmost_index].clear();
	surf.edgeloops[outmost_index].clear();
	for (unsigned int j = 0; j < tmp_loop1.size(); j++){
		surf.vertexloops[outmost_index].push_back(tmp_loop1[j]);
		surf.edgeloops[outmost_index].push_back(tmp_loop2[j]);
	}
	//adjust surf.connect and surf.connEdges based on the boundary loops
	int count = 0;
	std::vector<int> v, e;
	for (unsigned int i = 0; i < surf.vertexloops.size(); i++)
		for (unsigned int j = 0; j < surf.vertexloops[i].size(); j++)
			v.push_back(surf.vertexloops[i][j]);
	for (unsigned int i = 0; i < surf.edgeloops.size(); i++)
		for (unsigned int j = 0; j < surf.edgeloops[i].size(); j++)
			e.push_back(surf.edgeloops[i][j]);	
	
	for (unsigned int i = 0; i < v.size(); i++)
		surf.connect[i] = &vertices[v[i]];
	for (unsigned int i = 0; i < e.size(); i++)
		surf.connEdges[i] = &edges[e[i]];
	
}

//get the geometric loops on target surfaces
void SurfProHarmonicMap::GetGeomLoops(Face surf, vector<vector<int> > &loops_vertex, vector<vector<int> > &loops_edge)
{
	//collect edges' information
	set<int> edge_index;
	for (unsigned int i = 0; i < surf.connEdges.size(); i++)
		edge_index.insert(surf.connEdges[i]->index);

	int	CurrentNum = 0;
	int TotalNum = surf.connect.size();
	
	while(CurrentNum < TotalNum){
		int start_vertex = 0, next_vertex = -1, begin_vertex = 0;
		int start_edge = 0, next_edge = -1;
		
		//get the starting edge
		start_vertex = edges[*(edge_index.begin())].connect[0]->index;
		start_edge = *(edge_index.begin());

		//get the edge sense and make sure that the nodes are oriented correctly
		int edge_sense = 0;
		int first_vertex = start_vertex, second_vertex = edges[*(edge_index.begin())].connect[1]->index;
		g_err = igeom_instance->getEgVtxSense(edges[start_edge].gEdgeHandle, vertices[first_vertex].gVertexHandle, vertices[second_vertex].gVertexHandle, edge_sense);
		IBERRCHK(g_err, "Trouble get the edge sense.");
		int face_sense = 0;
		g_err = igeom_instance->getEgFcSense(edges[start_edge].gEdgeHandle,  surf.gFaceHandle, face_sense);
		IBERRCHK(g_err, "Trouble get the face sense.");
		if (face_sense*edge_sense < 0){
			start_vertex = edges[*(edge_index.begin())].connect[1]->index;
		}

		//initialization
		loops_vertex.resize(loops_vertex.size()+1);
		loops_edge.resize(loops_edge.size()+1);
		begin_vertex = start_vertex;
		//find the geometric loop		
		while(next_vertex != begin_vertex){
			loops_vertex[loops_vertex.size()-1].push_back(start_vertex);
			loops_edge[loops_edge.size()-1].push_back(start_edge);
			edge_index.erase(start_edge);
			if (start_vertex == edges[start_edge].connect[0]->index)
				next_vertex = edges[start_edge].connect[1]->index;
			else
				next_vertex = edges[start_edge].connect[0]->index;
			//find the adjacent geometric edge
			vector<iBase_EntityHandle> adj_edges;
			g_err = igeom_instance->getEntAdj(vertices[next_vertex].gVertexHandle, iBase_EDGE, adj_edges);
			IBERRCHK(g_err, "Trouble get the adjacent edges around a geometric vertex.");
			//remove unnecessary geometric edges
			assert(adj_edges.size()>=1);
			int intdata = -1;
			int next_edge_index = -1;
			for (int i = adj_edges.size()-1; i >= 0; i--){
				g_err = igeom_instance->getIntData(adj_edges[i], harmonic_surf_pro, intdata);
				if (edge_index.find(intdata)==edge_index.end())//remove that geometric edge
					continue;
				else{
					next_edge_index = intdata;
					break;
				}
			}
			if (next_edge_index == -1)
				break;
			next_edge = intdata;
			start_vertex = next_vertex;
			start_edge = next_edge;
		}
		CurrentNum += loops_vertex[loops_vertex.size()-1].size();
	}
}


SurfProHarmonicMap::~SurfProHarmonicMap()
{
}

}
