#include "meshkit/SubMapping.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/VertexMesher.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/SimpleArray.hpp"
#include "lp_lib.h"
#include "meshkit/MeshImprove.hpp"
#include "LPSolveClass.hpp"
#include "IsoLaplace.hpp"
#include "EquipotentialSmooth.hpp"
#ifdef HAVE_MESQUITE
#include "meshkit/MeshImprove.hpp"
#endif
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>



namespace MeshKit
{

//---------------------------------------------------------------------------//
//Entity Type initilization for SubMapping meshing
moab::EntityType SubMapping_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBQUAD, moab::MBMAXTYPE};
const moab::EntityType* SubMapping::output_types()
  { return SubMapping_tps; }

//---------------------------------------------------------------------------//
// construction function for SubMapping class
SubMapping::SubMapping(MKCore *mk_core, const MEntVector &me_vec) : MeshScheme(mk_core, me_vec)
{
	//buildAssociation();
}


//---------------------------------------------------------------------------//
// deconstruction function for SubMapping class
SubMapping::~SubMapping()
{

}

//---------------------------------------------------------------------------//
// setup function: define the size between the different layers
void SubMapping::setup_this()
{
	if (mentSelection.empty())
    	return;
    //compute the number of intervals for the associated ModelEnts, from the size set on them
    //the sizing function they point to, or a default sizing function
  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
      ModelEnt *me = mit -> first;
	  int dimension = me->dimension();
	  if (dimension != 2)
		 ECERRCHK(MK_FAILURE, "bad input for TFI Mapping, we can only mesh surfaces");
	  //first check whether the surface is meshed or not
      if (me->get_meshed_state() >= COMPLETE_MESH)
      	   continue;

	  
	  SizingFunction *sf = mk_core()->sizing_function(me->sizing_function_index());
	  if (!sf && me -> mesh_intervals() < 0 && me -> interval_firmness() == DEFAULT &&
        mk_core()->sizing_function(0))
        sf = mk_core()->sizing_function(0);

	  if (!sf && me -> mesh_intervals() < 0 && me -> interval_firmness() == DEFAULT){
         //no sizing set, just assume default #intervals as 20
         me->mesh_intervals(20);
         me->interval_firmness(DEFAULT);
      }
      else{
			//check # intervals first, then size, and just choose for now
      		if (sf->intervals() > 0)
      			throw Error(MK_INCOMPLETE_MESH_SPECIFICATION,  "Sizing function for edge should have edge size instead of intervals for submapping.");
			else if (sf->size()>0){
				size_low_bound = sf->size();
        		me->interval_firmness(SOFT);
			}
			else
        		throw Error(MK_INCOMPLETE_MESH_SPECIFICATION,  "Sizing function for edge had neither positive size nor positive intervals.");
      }	  
  }


}


//---------------------------------------------------------------------------//
// execute function: generate the all-hex mesh through sweeping from source 
// surface to target surface
void SubMapping::execute_this()
{
	std::vector<double> coords;
	std::vector<moab::EntityHandle> nodes;

	iBase_TagHandle global_tag, global_tag1;
	iGeom::Error g_err = mk_core()->igeom_instance()->getTagHandle("GLOBAL_ID", global_tag);
	IBERRCHK(g_err, "Trouble get the mesh entity set from geometric edges.");

	iMesh::Error m_err = mk_core()->igeom_instance()->getTagHandle("GLOBAL_ID", global_tag1);
	IBERRCHK(m_err, "Trouble get the mesh entity set from geometric edges.");
	
	for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  	{
    	ModelEnt *me = mit -> first;		
		//first check whether the surface is meshed or not
		if (me->get_meshed_state() >= COMPLETE_MESH)
			continue;

		//classify the vertex as END, CORNER, REVERSAL, SIDE
	  	VertexClassification(me);
      	//classify the edge as +I, -I, +J and -J
	  	EdgeClassification();
      //low bound for mesh size
		EdgeDiscretization(me);

		InteriorNodeInterpolation(me);
		MeshSmoothing(me);

		//ok, we are done!
		vertices_types.clear();
		nodes.clear();
		vertices.clear();
		edges.clear();
		interior_angle.clear();
		sorted_vertex_list.clear();
		sorted_node_list.clear();
		sorted_edge_list.clear();
		edges_types.clear();
		coordinate_i_j.clear();
		edge_size.clear();
		quads.clear();

		me->commit_mesh(mit->second, COMPLETE_MESH);
		mk_core()->save_mesh("test_cai.vtk");
	}

}

//setup the mesh size
void SubMapping::SetupMeshSize(double size)
{
	size_low_bound = size;
}


//classify the vertices as side, end, reversal and corner
//end 		--an interior angle of approximately 90 degrees
//side 		--an interior angle of approximately 180 degrees = 0
//corner	--an interior angle of approximately 270 degrees = -90
//reversal  --an interior angle of approximately 360 degrees = -180
void SubMapping::VertexClassification(ModelEnt *ent)
{
	//create a taghandle
	iGeom::Error g_err = mk_core()->igeom_instance()->getTagHandle("GEOM_SUBMAPPING", g_taghandle);
	if (g_err){
		g_err = mk_core()->igeom_instance()->createTag("GEOM_SUBMAPPING", 1, iBase_INTEGER, g_taghandle);	
		IBERRCHK(g_err, "Trouble create a taghandle.");
	}

	//get the geometric vertices from the surface
	vertices.clear();
	vector<iBase_EntityHandle> tmp;
	g_err = mk_core()->igeom_instance()->getEntAdj(ent->geom_handle(), iBase_VERTEX, tmp);
	IBERRCHK(g_err, "Trouble get the adjacent geometric vertices on a surface.");
	assert(tmp.size()>0);

	vertices.resize(tmp.size());
	for (unsigned int i = 0; i < tmp.size(); i++){
		vertices[i].index = i;
		vertices[i].gVertexHandle = tmp[i];

		g_err = mk_core()->igeom_instance()->getVtxCoord(vertices[i].gVertexHandle, vertices[i].xyz[0], vertices[i].xyz[1], vertices[i].xyz[2]);
		IBERRCHK(g_err, "Trouble get the coordinates from vertices.");

		g_err = mk_core()->igeom_instance()->setIntData(vertices[i].gVertexHandle, g_taghandle, i);
		IBERRCHK(g_err, "Trouble set the int data for the geometric vertex.");
	}

	vertices_types.resize(vertices.size());
	
	//get the geometric edges
	edges.clear();
	
	tmp.clear();
	g_err = mk_core()->igeom_instance()->getEntAdj(ent->geom_handle(), iBase_EDGE, tmp);
	IBERRCHK(g_err, "Trouble get the adjacent geometric vertices on a surface.");

	std::set<iBase_EntityHandle> edge_set;
	edges.resize(tmp.size());
	for (unsigned int i = 0; i < tmp.size(); i++){
		edges[i].connect.resize(2);
		edges[i].index = i;
		edges[i].gEdgeHandle = tmp[i];
		edge_set.insert(tmp[i]);

		g_err = mk_core()->igeom_instance()->setIntData(edges[i].gEdgeHandle, g_taghandle, i);
		IBERRCHK(g_err, "Trouble set the int data for the geometric edge.");

		vector<iBase_EntityHandle> adj_vertices;
		g_err = mk_core()->igeom_instance()->getEntAdj(edges[i].gEdgeHandle, iBase_VERTEX, adj_vertices);
		IBERRCHK(g_err, "Trouble get the adjacent vertices for the geometric edge.");

		assert(adj_vertices.size()<=2);
		for (unsigned int j = 0; j < adj_vertices.size(); j++){
			int index = -1;
			g_err = mk_core()->igeom_instance()->getIntData(adj_vertices[j], g_taghandle, index);
			IBERRCHK(g_err, "Trouble get the int data of geometric vertex.");

			edges[i].connect[j] = &vertices[index];
		}
	}

	//organize the vertices and edges on the boundaries
	VerEdgOrganize(edge_set, tmp, ent->geom_handle());

	//calculate the angle for vertices on the boundaries
	interior_angle.resize(vertices.size());
	GetAngle(ent->geom_handle(), interior_angle);

	//initial vertex classification,
	//need to use the linear programming method to produce the valid vertex classification 
	//if the initial vertex classification is not a valid vertex classification
	vertices_types.resize(vertices.size());
	for (unsigned int i = 0; i < interior_angle.size(); i++){
		if ((interior_angle[i] >= 45)&&(interior_angle[i] <= 130)){
			vertices_types[i] = END;
		}
		else if ((interior_angle[i] < 225)&&(interior_angle[i] > 130)){
			vertices_types[i] = SIDE;
		}
		else if ((interior_angle[i] > 225)&&(interior_angle[i] < 315)){
			vertices_types[i] = CORNER;
		}
		else{
			vertices_types[i] = REVERSAL;
		}
	}

	//extra step
	for (unsigned int m = 0; m < vertices_types.size(); m++){
			
		if (vertices_types[m] == REVERSAL)
			vertices_types[m] = SIDE;
	}




	int sum = 0;
	for (std::vector<VertexTypes>::iterator it = vertices_types.begin(); it != vertices_types.end(); it++)
		sum += *it;

	std::cout << "\n\n\nbefore, sum of vertex type is " << sum << std::endl;
	
	if (sum != 4){//check whether the initial vertex classification is valid or not
		//use the linear programming to produce a valid vertex classification
				
		VtxClassificationLP();
	}
	

	std::cout << "\n\n\nsum = " << sum << std::endl;
}

//use the lpsolve library to solve the linear programming
void SubMapping::VtxClassificationLP()
{
	LPSolveClass lp;
	vector<int> VtxType1(vertices_types.size()), VtxType2;
	int m = 0;
	for (std::vector<VertexTypes>::iterator it = vertices_types.begin(); it != vertices_types.end(); it++){
			
		if (*it == REVERSAL)
			vertices_types[m] = SIDE;
		
		VtxType1[m] = int(vertices_types[m]);		
		m++;
	}

	//setup the objective function: minimize the max value, first half variables are new variables, the latter half variables are old variables
	vector<double> b(2*vertices_types.size());
	for (unsigned int i = 0; i < vertices_types.size(); i++){
		b[i] = 1.0;
		b[vertices_types.size()+i] = 0.0;
	}
	lp.SetupObj(b, 0.0);
	
	//setup the equality constraint
	vector<vector<double> > left;
	vector<double> right;
	right.push_back(4.0);
	left.resize(1);
	left[0].resize(2*vertices_types.size());
	for (unsigned int i = 0; i < vertices_types.size(); i++){
		left[0][i] = 0.0;
		left[0][vertices_types.size()+i] = 1.0;
	}
	/*
	for (unsigned int i = 0; i < vertices_types.size(); i++){
		right.push_back(0.0);
		left[i+1].resize(2*vertices_types.size());
		for (unsigned int j = 0; j < vertices_types.size(); j++){
			left[i+1][j] = 0.0;
			left[i+1][vertices_types.size()+j] = 0.0;
			if ((i==j)&&(vertices_types[i]==END)){//fix the END vertex type
				left[i+1][vertices_types.size()+j] = 1.0;
				right[i+1] = 1.0;
			}
		}
	}
	*/
	lp.SetupEqu(left, right);

	//setup the inequality constraint
	left.clear();
	right.clear();
	left.resize(6*vertices_types.size());
	right.resize(6*vertices_types.size());
	for (unsigned int i = 0; i < vertices_types.size(); i++){
		right[i] = VtxType1[i];
		right[vertices_types.size()+i] = VtxType1[i];
		right[2*vertices_types.size()+i] = 1;
		right[3*vertices_types.size()+i] = 2;
		right[4*vertices_types.size()+i] = 1+VtxType1[i];
		right[5*vertices_types.size()+i] = 1-VtxType1[i];
	}
	for (unsigned int i = 0; i < vertices_types.size(); i++){
		left[i].resize(2*vertices_types.size());
		left[vertices_types.size()+i].resize(2*vertices_types.size());
		for (unsigned int j= 0; j < vertices_types.size(); j++){
			left[i][j] = -1.0;
			left[i][vertices_types.size()+j] = 1.0;
			left[vertices_types.size()+i][j] = -1.0;
			left[vertices_types.size()+i][vertices_types.size()+j] = -1.0;
		}
		//constrain the variable range
		left[2*vertices_types.size()+i].resize(2*vertices_types.size());
		left[3*vertices_types.size()+i].resize(2*vertices_types.size());
		left[4*vertices_types.size()+i].resize(2*vertices_types.size());
		left[5*vertices_types.size()+i].resize(2*vertices_types.size());
		for (unsigned int j = 0; j < 2*vertices_types.size(); j++){
			left[2*vertices_types.size()+i][j] = 0.0;
			left[3*vertices_types.size()+i][j] = 0.0;
			left[4*vertices_types.size()+i][j] = 0.0;
			left[5*vertices_types.size()+i][j] = 0.0;
		}
		left[2*vertices_types.size()+i][vertices_types.size()+i] = 1.0;
		left[3*vertices_types.size()+i][vertices_types.size()+i] = -1.0;

		//constrain the difference between new vertex type and initial vertex type
		left[4*vertices_types.size()+i][vertices_types.size()+i] = 1.0;
		left[5*vertices_types.size()+i][vertices_types.size()+i] = -1.0;
			
	}
	lp.SetupInEqu(left, right);

	lp.Execute();
	//get the solved variables
	vector<int> vtx_types;
	lp.GetVariables(vtx_types);
	for (unsigned int i = 0; i < vertices_types.size(); i++){
		switch(vtx_types[vertices_types.size()+i]){
			case 1:
				vertices_types[i] = END;
				break;
			case 0:
				vertices_types[i] = SIDE;
				break;
			case -1:
				vertices_types[i] = CORNER;
				break;
			case -2:
				vertices_types[i] = REVERSAL;
				break;
			default:
				break;
		}
	}

	//check whether sum = 4;
	int sum = 0;
	for (std::vector<VertexTypes>::iterator it = vertices_types.begin(); it != vertices_types.end(); it++)
		sum += *it;
	assert(sum == 4);
}

//calculate the angle for vertices on the boundaries
void SubMapping::GetAngle(iBase_EntityHandle surf, vector<double> &angle)
{
	angle.resize(sorted_vertex_list.size());
	for (unsigned int i = 0; i < sorted_vertex_list.size(); i++){
		//vertices = i,  previous edge = (i+sorted_edge_list.size()-1)%sorted_edge_list.size(), next edge = (i+sorted_edge_list.size()+1)%sorted_edge_list.size()

		//calculate the tangent vector at the previous edge
		double tangvector_p[3];
		int previous = (i+sorted_edge_list.size()-1)%sorted_edge_list.size();
		//calculate the tangent vector and edge sense with respect to two vertices
		iGeom::Error g_err = mk_core()->igeom_instance()->getEntTgntXYZ(edges[sorted_edge_list[previous]].gEdgeHandle, vertices[sorted_vertex_list[i]].xyz[0], vertices[sorted_vertex_list[i]].xyz[1], vertices[sorted_vertex_list[i]].xyz[2], tangvector_p[0], tangvector_p[1], tangvector_p[2]);
		IBERRCHK(g_err, "Trouble get the tangent vector at a vertex of an edge.");

		int sense = -10;
		g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[sorted_edge_list[previous]].gEdgeHandle, vertices[sorted_vertex_list[i]].gVertexHandle, vertices[sorted_vertex_list[(i+vertices.size()-1)%vertices.size()]].gVertexHandle, sense);
		IBERRCHK(g_err, "Trouble get edge sense with respect to two vertices.");

		int pre_edge_sense_face = -10;
		g_err = mk_core()->igeom_instance()->getEgFcSense(edges[sorted_edge_list[previous]].gEdgeHandle, surf, pre_edge_sense_face);
		IBERRCHK(g_err, "Trouble get edge sense with respect to the face.");
		std::cout << "-------------------\nvertex index = " << i << "\tprevious edge sense with respect to face = " << pre_edge_sense_face << std::endl;

		if (sense < 0){
			tangvector_p[0] = -1.0*tangvector_p[0];
			tangvector_p[1] = -1.0*tangvector_p[1];
			tangvector_p[2] = -1.0*tangvector_p[2];
		}
					
		//calculate the tangent vector at next edge
		double tangvector_n[3];
		int next = (i+sorted_edge_list.size())%sorted_edge_list.size();
		g_err = mk_core()->igeom_instance()->getEntTgntXYZ(edges[sorted_edge_list[next]].gEdgeHandle, vertices[sorted_vertex_list[i]].xyz[0], vertices[sorted_vertex_list[i]].xyz[1], vertices[sorted_vertex_list[i]].xyz[2], tangvector_n[0], tangvector_n[1], tangvector_n[2]);
		IBERRCHK(g_err, "Trouble get the tangent vector at a vertex of an edge.");
		g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[sorted_edge_list[next]].gEdgeHandle, vertices[sorted_vertex_list[i]].gVertexHandle, vertices[sorted_vertex_list[(i+1)%vertices.size()]].gVertexHandle, sense);
		IBERRCHK(g_err, "Trouble get edge sense with respect to two vertices.");
		
		int next_edge_sense_face = -10;
		g_err = mk_core()->igeom_instance()->getEgFcSense(edges[sorted_edge_list[next]].gEdgeHandle, surf, next_edge_sense_face);
		IBERRCHK(g_err, "Trouble get edge sense with respect to the face.");
		std::cout << "\t" << i << "\tnext edge sense with respect to face = " << next_edge_sense_face << std::endl;

		if (sense < 0){
			tangvector_n[0] = -1.0*tangvector_n[0];
			tangvector_n[1] = -1.0*tangvector_n[1];
			tangvector_n[2] = -1.0*tangvector_n[2];
		}

		
		//calculate the crossproduct
		double crossproduct[3];
		crossproduct[0] = tangvector_n[1]*tangvector_p[2]-tangvector_n[2]*tangvector_p[1];
		crossproduct[1] = tangvector_n[2]*tangvector_p[0]-tangvector_n[0]*tangvector_p[2];
		crossproduct[2] = tangvector_n[0]*tangvector_p[1]-tangvector_n[1]*tangvector_p[0];

		//calculate the dotproduct
		double dotproduct = tangvector_n[0]*tangvector_p[0] + tangvector_n[1]*tangvector_p[1] + tangvector_n[2]*tangvector_p[2];
		double theta_pi = acos(dotproduct/(sqrt(pow(tangvector_n[0],2)+pow(tangvector_n[1], 2)+pow(tangvector_n[2], 2))*sqrt(pow(tangvector_p[0],2)+pow(tangvector_p[1],2)+pow(tangvector_p[2],2))));
		double theta = theta_pi*180.0/3.1415926;

		//get the surface norm at a specific point
		double normal[3];
		g_err = mk_core()->igeom_instance()->getEntNrmlXYZ(surf, vertices[sorted_vertex_list[i]].xyz[0], vertices[sorted_vertex_list[i]].xyz[1], vertices[sorted_vertex_list[i]].xyz[2], normal[0], normal[1], normal[2]);
		IBERRCHK(g_err, "Trouble get the normal at a specific point on the surface.");
		dotproduct = crossproduct[0]*normal[0]+crossproduct[1]*normal[1]+crossproduct[2]*normal[2];
		
		if (dotproduct < 0)
			theta = 360.0 - theta;
		
		angle[sorted_vertex_list[i]] = theta;
	}

}

bool SubMapping::isCurved(int vtx_index, vector<double> u1, vector<double> u2, vector<double> u3, vector<double> u4, vector<vector<double> > tang_pre, vector<vector<double> > tang_next)
{
	int previous = (vtx_index+sorted_edge_list.size()-1)%sorted_edge_list.size(), next = (vtx_index+sorted_edge_list.size())%sorted_edge_list.size();
	double	dotproduct, vec[3], theta1, theta2, pts2[3], pts1[3], u, length;
	u = u1[vtx_index] + 0.5*(u2[vtx_index]-u1[vtx_index]);
	iGeom::Error g_err = mk_core()->igeom_instance()->getEntUtoXYZ(edges[sorted_edge_list[previous]].gEdgeHandle, u, pts1[0], pts1[1], pts1[2]);
	IBERRCHK(g_err, "Trouble get xyz coordinates from parametric coordinates on edge.");
	u = u1[vtx_index] + 0.501*(u2[vtx_index]-u1[vtx_index]);
	g_err = mk_core()->igeom_instance()->getEntUtoXYZ(edges[sorted_edge_list[previous]].gEdgeHandle, u, pts2[0], pts2[1], pts2[2]);
	IBERRCHK(g_err, "Trouble get xyz coordinates from parametric coordinates on edge.");
	for (int j = 0; j < 3; j++)
		vec[j] = pts2[j] - pts1[j];
	dotproduct = vec[0]*tang_pre[vtx_index][0] + vec[1]*tang_pre[vtx_index][1] + vec[2]*tang_pre[vtx_index][2];
	length = sqrt(pow(vec[0],2)+pow(vec[1], 2)+pow(vec[2], 2))*sqrt(pow(tang_pre[vtx_index][0],2)+pow(tang_pre[vtx_index][1],2)+pow(tang_pre[vtx_index][2],2));
	if (fabs(dotproduct) > length){
		if (dotproduct > 0)
			dotproduct = length;
		else
			dotproduct = -1.0*length;
	}
	theta1 = 180.0/3.1415926*acos(dotproduct/length);
	u = u3[vtx_index] + 0.5*(u4[vtx_index]-u3[vtx_index]);
	g_err = mk_core()->igeom_instance()->getEntUtoXYZ(edges[sorted_edge_list[next]].gEdgeHandle, u, pts1[0], pts1[1], pts1[2]);
	IBERRCHK(g_err, "Trouble get xyz coordinates from parametric coordinates on edge.");
	u = u3[vtx_index] + 0.501*(u4[vtx_index]-u3[vtx_index]);
	g_err = mk_core()->igeom_instance()->getEntUtoXYZ(edges[sorted_edge_list[next]].gEdgeHandle, u, pts2[0], pts2[1], pts2[2]);
	IBERRCHK(g_err, "Trouble get xyz coordinates from parametric coordinates on edge.");
	for (int j = 0; j < 3; j++)
		vec[j] = pts2[j] - pts1[j];
	dotproduct = vec[0]*tang_next[vtx_index][0] + vec[1]*tang_next[vtx_index][1] + vec[2]*tang_next[vtx_index][2];
	length = sqrt(pow(vec[0],2)+pow(vec[1], 2)+pow(vec[2], 2))*sqrt(pow(tang_next[vtx_index][0],2)+pow(tang_next[vtx_index][1],2)+pow(tang_next[vtx_index][2],2));
	if (fabs(dotproduct) > length){
		if (dotproduct > 0)
			dotproduct = length;
		else
			dotproduct = -1.0*length;
	}
	theta2 = 180.0/3.1415926*acos(dotproduct/length);
	
	if ((fabs(theta1) > 0.5)&&(fabs(theta2) > 0.5))
		return true;
	else
		return false;

}


//reorganize the vertices and edges on the boudnaries of geometry
void SubMapping::VerEdgOrganize(std::set<iBase_EntityHandle> edge_set, std::vector<iBase_EntityHandle> g_edge, iBase_EntityHandle surf)
{
	int first_edge_index = 0;
	int first_index = edges[0].connect[1]->index;
	int second_index = edges[0].connect[0]->index;
	int start_index = edges[0].connect[0]->index;

	int test_sense = -10;
	iGeom::Error g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[first_edge_index].gEdgeHandle, vertices[second_index].gVertexHandle, vertices[first_index].gVertexHandle, test_sense);
	IBERRCHK(g_err, "Trouble get the sense of edge with respect to two vertices.");

	int test_face_sense = -10;
	g_err = mk_core()->igeom_instance()->getEgFcSense(edges[first_edge_index].gEdgeHandle, surf, test_face_sense);
	IBERRCHK(g_err, "Trouble get the sense of edge with respect to two vertices.");
	std::cout << "edge sense = " << test_sense << "\tface sense = " << test_face_sense << "combined(multiply) = " << test_sense*test_face_sense << std::endl;	

	if ((test_sense*test_face_sense) < 0){
		int tmp = second_index;
		second_index = first_index;
		first_index = tmp;
		start_index = second_index;
	}
	

	sorted_vertex_list.push_back(start_index);
	
	sorted_edge_list.push_back(first_edge_index);

	while(start_index != first_index){
		sorted_vertex_list.push_back(first_index);

		vector<iBase_EntityHandle> adj_edges;
		g_err = mk_core()->igeom_instance()->getEntAdj(vertices[first_index].gVertexHandle, iBase_EDGE, adj_edges);
		IBERRCHK(g_err, "Trouble get the adjacent edges w.r.t a vertex.");

		int index = -1;
		for (unsigned int i = 0; i < adj_edges.size(); i++){
			if ((adj_edges[i] != edges[first_edge_index].gEdgeHandle)&&(edge_set.find(adj_edges[i]) != edge_set.end())){
				g_err = mk_core()->igeom_instance()->getIntData(adj_edges[i], g_taghandle, index);
				IBERRCHK(g_err, "Trouble get the int data for geometric edge.");
				break;
			}
			else
				continue;
		}
		
		//insert the new edge into sorted_edge_list
		if (index > -1){		
			
			first_edge_index = index;
			sorted_edge_list.push_back(first_edge_index);
			second_index = first_index;
			if (edges[first_edge_index].connect[0]->index == first_index)
				first_index = edges[first_edge_index].connect[1]->index;
			else
				first_index = edges[first_edge_index].connect[0]->index;
			
		}		
		else
			break;	
	}
}


//classify the boundary edges as -I, +I, -J, +J
void SubMapping::EdgeClassification()
{
	edges_types.resize(sorted_edge_list.size());
	
	//Find a starting END vertex
	start_index = 0;
	for (; start_index < sorted_vertex_list.size(); start_index++)
		if (vertices_types[sorted_vertex_list[start_index]]==END)
			break;

	//setting up the initial direction for i-j coordinate system(positive i and positive j direction)
	int pre_edge = (start_index + sorted_edge_list.size() -1)%sorted_edge_list.size(), next_edge = start_index;
		
	edges_types[sorted_edge_list[next_edge]] = POSI_J;
	edges_types[sorted_edge_list[pre_edge]] = NEG_I;
	EdgeTypes i_direction = NEG_I, j_direction = POSI_J;
	VertexTypes pre_vertex_type = END;

	int vertex_index = start_index;
	for (unsigned int i = 1; i < sorted_vertex_list.size()-1; i++){
		pre_edge = (start_index + i + sorted_edge_list.size() -1)%sorted_edge_list.size(); //index in the sorted_edge_list
		next_edge = (start_index + i)%sorted_edge_list.size();   //index in the sorted_edge_list
		vertex_index = (i + start_index)%sorted_vertex_list.size();

		switch(vertices_types[sorted_vertex_list[vertex_index]]){
			case END://switch from i to j or from j to i
				if ((edges_types[sorted_edge_list[pre_edge]] == POSI_I)||(edges_types[sorted_edge_list[pre_edge]] == NEG_I)){
						if (j_direction == POSI_J)
							if (pre_vertex_type == END)
								edges_types[sorted_edge_list[next_edge]] = NEG_J;
							else
								edges_types[sorted_edge_list[next_edge]] = POSI_J;
						else
							if (pre_vertex_type == END)
								edges_types[sorted_edge_list[next_edge]] = POSI_J;
							else
								edges_types[sorted_edge_list[next_edge]] = NEG_J;
						j_direction = edges_types[sorted_edge_list[next_edge]];					
				}
				else{//switch from j-direction to i-direction
					if (i_direction == POSI_I)
						if (pre_vertex_type == END)
							edges_types[sorted_edge_list[next_edge]] = NEG_I;
						else
							edges_types[sorted_edge_list[next_edge]] = POSI_I;
					else
						if (pre_vertex_type == END)
							edges_types[sorted_edge_list[next_edge]] = POSI_I;
						else
							edges_types[sorted_edge_list[next_edge]] = NEG_I;
					i_direction = edges_types[sorted_edge_list[next_edge]];
				}
				pre_vertex_type = END;				
				break;
			case SIDE://keep the same as previous
				edges_types[sorted_edge_list[next_edge]] = edges_types[sorted_edge_list[pre_edge]];
				break;
		
			case CORNER://switch from i to j or from j to i
				if ((edges_types[sorted_edge_list[pre_edge]] == POSI_I)||(edges_types[sorted_edge_list[pre_edge]] == NEG_I)){
					if (j_direction == POSI_J)
						if (pre_vertex_type == END)
							edges_types[sorted_edge_list[next_edge]] = POSI_J;
						else
							edges_types[sorted_edge_list[next_edge]] = NEG_J;
					else
						if (pre_vertex_type == END)
							edges_types[sorted_edge_list[next_edge]] = NEG_J;
						else
							edges_types[sorted_edge_list[next_edge]] = POSI_J;
					j_direction = edges_types[sorted_edge_list[next_edge]];						
				}
				else{
					if (i_direction == POSI_I)
						if (pre_vertex_type == END)
							edges_types[sorted_edge_list[next_edge]] = POSI_I;
						else
							edges_types[sorted_edge_list[next_edge]] = NEG_I;
					else
						if (pre_vertex_type == END)
							edges_types[sorted_edge_list[next_edge]] = NEG_I;
						else
							edges_types[sorted_edge_list[next_edge]] = POSI_I;
					i_direction = edges_types[sorted_edge_list[next_edge]];	
				}
				pre_vertex_type = CORNER;				
				break;

			case REVERSAL://need to consider this point later on
				if (edges_types[sorted_edge_list[pre_edge]] == POSI_I)
					edges_types[sorted_edge_list[next_edge]] = NEG_I;
				else if (edges_types[sorted_edge_list[pre_edge]] == NEG_I)
					edges_types[sorted_edge_list[next_edge]] = POSI_I;
				else if (edges_types[sorted_edge_list[pre_edge]] == POSI_J)
					edges_types[sorted_edge_list[next_edge]] = NEG_J;
				else
					edges_types[sorted_edge_list[next_edge]] = POSI_J;
				break;

			default://do nothing
				break;
		}
	}	
}

//assign the i,j coordinates for each node on the boundary
void SubMapping::EdgeDiscretization(ModelEnt *me)
{
	edge_size.resize(sorted_edge_list.size());
	for (unsigned int i = 0; i < sorted_edge_list.size(); i++){
		double measure;
		iGeom::Error g_err = mk_core()->igeom_instance()->measure(&(edges[sorted_edge_list[i]].gEdgeHandle), 1, &measure);
		IBERRCHK(g_err, "Trouble measure the boundary edges.");
		edge_size[sorted_edge_list[i]] = int(measure/size_low_bound);
	}

	//linear programming to get # of line segments for each boundary edge
	LPSolveClass lp;
	//setup the model for linear programming
	vector<vector<double> > coeffs;
	vector<double> b(edges.size(), 1.0);
	
	//setup the objective function for linear programming
	lp.SetupObj(b, 0.0);

	//setup the equality constraint for linear programming
	coeffs.resize(2);
	coeffs[0].resize(edges.size());
	coeffs[1].resize(edges.size());
	b.clear();
	b.resize(2, 0.0);
	for (unsigned int i = 0; i < edges.size(); i++){
		if ((edges_types[i] == POSI_I)||(edges_types[i] == NEG_I)){//positive i and negative i
			if (edges_types[i] == POSI_I)
				coeffs[0][i] = 1.0;
			else if (edges_types[i] == NEG_I)
				coeffs[0][i] = -1.0;
			coeffs[1][i] = 0.0;
		}
		else{//positive j and negative j
			if (edges_types[i] == POSI_J)
				coeffs[1][i] = 1.0;
			else if (edges_types[i] == NEG_J)
				coeffs[1][i] = -1.0;
			coeffs[0][i] = 0.0;
		}
	}
	lp.SetupEqu(coeffs, b);
	
	//setup the constant constraint there is already mesh on the geometric edge
	vector<int> num_line_segments(edges.size());
	for (unsigned int i = 0; i < edges.size(); i++){
		iBase_EntitySetHandle entityset;
		iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(edges[i].gEdgeHandle, 0, entityset);
		IBERRCHK(r_err, "Trouble get the entity set for geometric edge.");
		iMesh::Error m_err = mk_core()->imesh_instance()->getNumOfType(entityset, iBase_EDGE, num_line_segments[i]);
		IBERRCHK(m_err, "Trouble get # of line segments on the edge mesh.");

		if ((num_line_segments[i] < 1)&&(!m_err))
			num_line_segments[i] = -1;
	}
	lp.SetupConst(num_line_segments);

	//setup the inequality constraint for linear programming
	coeffs.clear();
	coeffs.resize(edges.size(), vector<double>(edges.size()));
	for (unsigned int i = 0; i < edges.size(); i++)//diagonal matrix
		coeffs[i][i] = -1.0;
	b.clear();
	b.resize(edges.size());
	int min_line_segment = 1.0e10;
	for (unsigned int i = 0; i < num_line_segments.size(); i++)
		if ((num_line_segments[i] < min_line_segment)&&(num_line_segments[i]>0))
			min_line_segment = num_line_segments[i];

	//preprocess the inequality constraints
	int sum_posi_I = 0, sum_neg_I = 0, sum_posi_J = 0, sum_neg_J = 0;
	vector<int> posi_i, neg_i, posi_j, neg_j;
	for (int i = 0; i < num_line_segments.size(); i++){
		if (num_line_segments[i] > 0){
			switch(edges_types[i]){
				case POSI_I:
					sum_posi_I += num_line_segments[i]; break;
				case NEG_I:
					sum_neg_I += num_line_segments[i]; break;
				case POSI_J:
					sum_posi_J += num_line_segments[i]; break;
				case NEG_J:
					sum_neg_J += num_line_segments[i]; break;
				default:
					break;
			}
		}
		else{
			switch(edges_types[i]){
				case POSI_I:
					posi_i.push_back(i); break;
				case NEG_I:
					neg_i.push_back(i); break;
				case POSI_J:
					posi_j.push_back(i); break;
				case NEG_J:
					neg_j.push_back(i); break;
				default:
					break;		
			}
		}
	}
	//preprocess the edge size
	int tmp_sum_POS_I = sum_posi_I, tmp_sum_POS_J = sum_posi_J, tmp_sum_NEG_I = sum_neg_I, tmp_sum_NEG_J = sum_neg_J;
	for (unsigned int i = 0; i < posi_i.size(); i++)
		tmp_sum_POS_I += edge_size[posi_i[i]];
	for (unsigned int i = 0; i < posi_j.size(); i++)
		tmp_sum_POS_J += edge_size[posi_j[i]];
	for (unsigned int i = 0; i < neg_i.size(); i++)
		tmp_sum_NEG_I += edge_size[neg_i[i]];
	for (unsigned int i = 0; i < neg_j.size(); i++)
		tmp_sum_NEG_J += edge_size[neg_j[i]];
	
	//process i-direction ----new
	if (tmp_sum_POS_I == tmp_sum_NEG_I)
	{}//do nothing, it is always true
	else if (tmp_sum_POS_I > tmp_sum_NEG_I){//adjust the NEG_I
		if (neg_i.size() > 0){
			for (unsigned int i = 0; i < neg_i.size(); i++)
				edge_size[neg_i[i]] = int(double(tmp_sum_POS_I-sum_neg_I)*double(edge_size[neg_i[i]])/double(tmp_sum_NEG_I-sum_neg_I));
		}
		else{//neg_i.size == 0
			if (posi_i.size() == 0){
				std::cout << "Constraint check fails in i-direction_a\n";exit(1);
			}
			else{
				if (sum_posi_I >= tmp_sum_NEG_I){//here, tmp_sum_NEG_I = sum_neg_I
					std::cout << "Constraint check fails in i-direction_b\n";exit(1);
				}
				else{//sum_posi_I < tmp_sum_NEG_I, adjust the edge posi_i
					for (unsigned int i = 0; i < posi_i.size(); i++){
						edge_size[posi_i[i]] = int(double(tmp_sum_NEG_I-sum_posi_I)*double(edge_size[posi_i[i]])/double(tmp_sum_POS_I-sum_posi_I));
					}
				}
			}
		}
	}
	else{//tmp_sum_POS_I < tmp_sum_NEG_I, //adjust the POSS_I
		if (posi_i.size() > 0){
			for (unsigned int i = 0; i < posi_i.size(); i++)
				edge_size[posi_i[i]] = int(double(tmp_sum_NEG_I-sum_posi_I)*double(edge_size[posi_i[i]])/double(tmp_sum_POS_I-sum_posi_I));
		}
		else{//posi_i.size() == 0
			//further constraints, 
			if (neg_i.size() == 0){
				std::cout << "Constraint check fails in i-direction_c\n";exit(1);
			}
			else{
				if (sum_neg_I >= tmp_sum_POS_I){//here tmp_sum_POS_I = sum_posi_I
					std::cout << "Constraint check fails in i-direction_d\n";exit(1);
				}
				else{//sum_neg_I < tmp_sum_POS_I, adjust the edge neg_i
					for (unsigned int i = 0; i < neg_i.size(); i++){
						edge_size[neg_i[i]] = int(double(tmp_sum_POS_I-sum_neg_I)*double(edge_size[neg_i[i]])/double(tmp_sum_NEG_I-sum_neg_I));
					}
				}
			}
		}
	}


	
	//process j-direction ----new
	if (tmp_sum_POS_J == tmp_sum_NEG_J)
	{}//do nothing, it is always true
	else if (tmp_sum_POS_J > tmp_sum_NEG_J){//adjust the NEG_J
		if (neg_j.size() > 0){
			for (unsigned int i = 0; i < neg_j.size(); i++)
				edge_size[neg_j[i]] = int(double(tmp_sum_POS_J-sum_neg_J)*double(edge_size[neg_j[i]])/double(tmp_sum_NEG_J-sum_neg_J));
		}
		else{//neg_j.size == 0
			if (posi_j.size() == 0){
				std::cout << "Constraint check fails in j-direction_a\n";exit(1);
			}
			else{
				if (sum_posi_J >= tmp_sum_NEG_J){//here, tmp_sum_NEG_J = sum_neg_J
					std::cout << "Constraint check fails in j-direction_b\n";exit(1);
				}
				else{//sum_posi_J < tmp_sum_NEG_J, adjust the edge posi_j
					for (unsigned int i = 0; i < posi_j.size(); i++){
						edge_size[posi_j[i]] = int(double(tmp_sum_NEG_J-sum_posi_J)*double(edge_size[posi_j[i]])/double(tmp_sum_POS_J-sum_posi_J));
					}
				}
			}
		}
	}
	else{//tmp_sum_POS_J < tmp_sum_NEG_J, //adjust the POS_I
		if (posi_j.size() > 0){
			for (unsigned int i = 0; i < posi_j.size(); i++)
				edge_size[posi_j[i]] = int(double(tmp_sum_NEG_J-sum_posi_J)*double(edge_size[posi_j[i]])/double(tmp_sum_POS_J-sum_posi_J));
		}
		else{//posi_j.size() == 0
			//further constraints, 
			if (neg_j.size() == 0){
				std::cout << "Constraint check fails in j-direction_c\n";exit(1);
			}
			else{
				if (sum_neg_J >= tmp_sum_POS_J){//here tmp_sum_POS_J = sum_posi_J
					std::cout << "Constraint check fails in j-direction_d\n";exit(1);
				}
				else{//sum_neg_J < tmp_sum_POS_J, adjust the edge neg_j
					for (unsigned int i = 0; i < neg_j.size(); i++){
						edge_size[neg_j[i]] = int(double(tmp_sum_POS_J-sum_neg_J)*double(edge_size[neg_j[i]])/double(tmp_sum_NEG_J-sum_neg_J));
					}
				}
			}
		}
	}



	

	//process i-direction
	if ((sum_posi_I == sum_neg_I)&&(sum_posi_I != 0)){
		if ((posi_i.size() != 0)&&(neg_i.size() != 0))
		{}//these constraints should work
		else if ((posi_i.size() == 0)&&(neg_i.size() == 0))
		{}//these constraints should work
		else if ((posi_i.size() != 0)&&(neg_i.size() == 0))//this constraints won't work
		{std::cout << "Constraint check fails in i-direction3\n";exit(1);}
		else//these constraints won't work
		{std::cout << "Constraint check fails in i-direction4\n";exit(1);}
	}
	else{
		if (sum_posi_I > sum_neg_I){
			if (sum_neg_I == 0){//low bound for NEG_I edges should be less than {sum_posi_I+edge_size[POSI_I]}
				//this work	
			}
			else{//sum_neg_I > 0
				if ((int(posi_i.size()) == 0)&&(int(neg_i.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in i-direction5\n";exit(1);}
				else if ((int(posi_i.size()) == 0)&&(int(neg_i.size()) > 0)){//sum_neg_I + edge_size[NEG_I] <= sum_posi_I
					//this works
				}
				else if ((int(posi_i.size()) > 0)&&(int(neg_i.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in i-direction6\n";exit(1);}
				else{ //((int(posi_i.size()) > 0)&&(int(neg_i.size()) > 0))
					//these constraints should work.
				}
			}
		}
		else{//sum_posi_I < sum_neg_I
			if (sum_posi_I == 0){//low bound for POSI_I edges should be less than {sum_neg_I+edge_size[NEG_I]}
				//this works	
			}
			else{//sum_posi_I > 0
				if ((int(posi_i.size()) == 0)&&(int(neg_i.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in i-direction7\n";exit(1);}
				else if ((int(neg_i.size()) == 0)&&(int(posi_i.size()) > 0)){//sum_posi_I + edge_size[POSI_I] <= sum_neg_I
					//this works
				}
				else if ((int(neg_i.size()) > 0)&&(int(posi_i.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in i-direction8\n";exit(1);}
				else{ //((int(posi_i.size()) > 0)&&(int(neg_i.size()) > 0))
					//these constraints should work.
				}
			}
		}
	}
	//process j direction
	if ((sum_posi_J == sum_neg_J)&&(sum_posi_J != 0)){
		if ((posi_j.size() != 0)&&(neg_j.size() != 0))//these constraints should work
		{}
		else if ((posi_j.size() == 0)&&(neg_j.size() == 0))
		{}//these constraints should work
		else if ((posi_j.size() != 0)&&(neg_j.size() == 0))//this constraints won't work
		{std::cout << "Constraint check fails in j-direction3\n";exit(1);}
		else//these constraints won't work
		{std::cout << "Constraint check fails in j-direction4\n";exit(1);}
	}
	else{//sum_posi_J != sum_neg_J
		if (sum_posi_J > sum_neg_J){
			if (sum_neg_J == 0){//low bound for NEG_J edges should be less than {sum_posi_J+edge_size[POSI_J]}
				//this works	
			}
			else{//sum_neg_J > 0
				if ((int(posi_j.size()) == 0)&&(int(neg_j.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in j-direction5\n";exit(1);}
				else if ((int(posi_j.size()) == 0)&&(int(neg_j.size()) > 0)){//sum_neg_J + edge_size[NEG_J] <= sum_posi_J
					//this works
				}
				else if ((int(posi_j.size()) > 0)&&(int(neg_j.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in j-direction6\n";exit(1);}
				else{ //((int(posi_j.size()) > 0)&&(int(neg_j.size()) > 0))
					//these constraints should work.
				}
			}
		}
		else{//sum_posi_J < sum_neg_J
			if (sum_posi_J == 0){//low bound for POSI_J edges should be less than {sum_neg_J+edge_size[NEG_J]}
				//this works		
			}
			else{//sum_posi_J > 0
				if ((int(posi_j.size()) == 0)&&(int(neg_j.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in j-direction7\n";exit(1);}
				else if ((int(neg_j.size()) == 0)&&(int(posi_j.size()) > 0)){//sum_posi_J + edge_size[POSI_J] <= sum_neg_J
					//this works
				}
				else if ((int(neg_j.size()) > 0)&&(int(posi_j.size()) == 0))//these constraints won't work
				{std::cout << "Constraint check fails in j-direction8\n";exit(1);}
				else{ //((int(posi_j.size()) > 0)&&(int(neg_j.size()) > 0))
					//these constraints should work.
				}
			}
		}
	}
	
	for (unsigned int i = 0; i < edges.size(); i++)
		if (min_line_segment < 1.0e9){
			if ((num_line_segments[i] < edge_size[i])&&(num_line_segments[i] > 0)){//avoid the conflicts between const constraints and inequality constraints;
				b[i] = -1.0*num_line_segments[i];
			}
			else{
				//if (min_line_segment < edge_size[i])
				//	b[i] = -1.0*min_line_segment;
				//else
				b[i] = -1.0*edge_size[i];
			}
		}
		else
			b[i] = -1.0*edge_size[i];
	
	lp.SetupInEqu(coeffs, b);

	lp.Execute();
	lp.GetVariables(edge_size);

	//do the vertex mesher
	//get the mesh node on the geometric vertex
	for (unsigned int i = 0; i < sorted_vertex_list.size(); i++){
		iBase_EntitySetHandle entityset;
		iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(vertices[sorted_vertex_list[i]].gVertexHandle, 0, entityset);
		IBERRCHK(r_err, "Trouble get the entity set for edge from geometric vertex.");

		vector<iBase_EntityHandle> points;
		points.clear();
		iMesh::Error m_err = mk_core()->imesh_instance()->getEntities(entityset, iBase_VERTEX, iMesh_POINT, points);
		IBERRCHK(m_err, "Trouble get the node entities from mesh entity sets.");

		if (points.size()==0){
			points.resize(1);
			//there is no mesh nodes on the geometric vertices, it need to create a mesh node on the geometric vertex
			m_err = mk_core()->imesh_instance()->createVtx(vertices[sorted_vertex_list[i]].xyz[0], vertices[sorted_vertex_list[i]].xyz[1], vertices[sorted_vertex_list[i]].xyz[2], points[0]);
			IBERRCHK(m_err, "Trouble create the mesh nodes on the geometric vertex.");
			m_err = mk_core()->imesh_instance()->addEntToSet(points[0], entityset);
			IBERRCHK(m_err, "Trouble add the mesh node entity to entity set.");
		}
	}

	//do the edge mesher
	MEntVector curves;
	me->get_adjacencies(1, curves);
	assert(curves.size()==edges.size());

	//call the edge mesher to discretize the boundary edges
	for (unsigned int i = 0; i < edges.size(); i++){
		//initial size functon for edges, get the number of edges and assign it to the edge

		//do the edge mesher 
		SizingFunction esize(mk_core(), edge_size[i], -1);
  		me->sizing_function_index(esize.core_index());

		//detect the edge on the surface
		MEntVector edge_curve;
		edge_curve.resize(1);
		
		for (unsigned int j = 0; j < curves.size(); j++){
			int index_id;
			iGeom::Error g_err = mk_core()->igeom_instance()->getIntData(curves[j]->geom_handle(), g_taghandle, index_id);
			IBERRCHK(g_err, "Trouble get the int data for geometric edges.");

			if (index_id == (int)i){
				edge_curve[0] = curves[j];
				break;
			}
		}
		
		//do the edge mesher
		EdgeMesher *em = (EdgeMesher*) mk_core()->construct_meshop("EdgeMesher", edge_curve);
		//mk_core()->setup_and_execute();
		em->setup_this();
		em->execute_this();	
		////done with the meshing for edge i on the surface
	}

	//extract the mesh info
	iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("MeshSubMapping", m_taghandle);
	if (m_err){
		m_err = mk_core()->imesh_instance()->createTag("MeshSubMapping", 1, iBase_INTEGER, m_taghandle);	
		IBERRCHK(m_err, "Trouble create a taghandle.");
	}

	int index = 0;	
	for (unsigned int i = 0; i < sorted_vertex_list.size(); i++){
		int next_edge = sorted_edge_list[i%sorted_edge_list.size()];
		int pre_edge = sorted_edge_list[(i+sorted_edge_list.size()-1)%sorted_edge_list.size()];
		int next_vertex = sorted_vertex_list[(i+1)%sorted_vertex_list.size()];
		iBase_EntitySetHandle entityset;

		//get the mesh node on the geometric vertex
		iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(vertices[sorted_vertex_list[i]].gVertexHandle, 0, entityset);
		IBERRCHK(r_err, "Trouble get the entity set for edge from geometric vertex.");

		vector<iBase_EntityHandle> points;
		points.clear();
		m_err = mk_core()->imesh_instance()->getEntities(entityset, iBase_VERTEX, iMesh_POINT, points);
		IBERRCHK(m_err, "Trouble get the node entities from mesh entity sets.");

		assert(points.size()==1);
		index++;
		nodes.resize(index);
		nodes[index-1].index = index-1;
		nodes[index-1].onBoundary = false;
		nodes[index-1].onCorner = true;
		geom_mesh_node[sorted_vertex_list[i]] = index -1;
		mesh_geom_vertex[index-1] = sorted_vertex_list[i];
		nodes[index-1].gVertexHandle = points[0];
		m_err = mk_core()->imesh_instance()->setIntData(nodes[index-1].gVertexHandle, m_taghandle, index-1);
		IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
		for (int j = 0; j < 3; j++)
			nodes[index-1].xyz[j] = vertices[sorted_vertex_list[i]].xyz[j];

		//assign the i-j coordinates for boundary nodes
		if (i == 0){
			nodes[index-1].uv[0] = 0.0;
			nodes[index-1].uv[1] = 0.0;
		}
		else{
			AssignIJCoords(nodes[index-1].uv[0], nodes[index-1].uv[1], edges_types[pre_edge], index -1);
		}

		//get the mesh nodes on the geometric edges
		r_err = mk_core()->irel_pair()->getEntSetRelation(edges[next_edge].gEdgeHandle, 0, entityset);
		IBERRCHK(r_err, "Trouble get the entity set for edge from geometric vertex.");

		points.clear();
		m_err = mk_core()->imesh_instance()->getEntities(entityset, iBase_VERTEX, iMesh_POINT, points);
		IBERRCHK(m_err, "Trouble get the node entities from mesh entity sets.");

		int sense = 0;
		iGeom::Error g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[next_edge].gEdgeHandle, vertices[sorted_vertex_list[i]].gVertexHandle, vertices[next_vertex].gVertexHandle, sense);
		IBERRCHK(g_err, "Trouble get the edge sense with respect to two vertices.");

		unsigned int j;
		if (sense < 0) 
			j = points.size()-1;
		else
			j = 0; 
		for (; ((j < points.size())&&(j>=0));){						
			index++;
			nodes.resize(index);
			nodes[index-1].index = index-1;
			nodes[index-1].gVertexHandle = points[j];
			nodes[index-1].onBoundary = true;
			nodes[index-1].onCorner = false;
			m_err = mk_core()->imesh_instance()->setIntData(nodes[index-1].gVertexHandle, m_taghandle, index-1);
			IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
			m_err = mk_core()->imesh_instance()->getVtxCoord(nodes[index-1].gVertexHandle, nodes[index-1].xyz[0], nodes[index-1].xyz[1], nodes[index-1].xyz[2]);
			IBERRCHK(m_err, "Trouble get the coordinates from mesh nodes.");
			AssignIJCoords(nodes[index-1].uv[0], nodes[index-1].uv[1], edges_types[next_edge], index -1);
			if (sense < 0) j--;
			else j++;
		}
	}
}

//Assign the i j coordinates for boundary nodes
void SubMapping::AssignIJCoords(double &u, double &v, EdgeTypes type, int index){
	switch(type){
		case POSI_I:
			nodes[index].uv[0] = nodes[index-1].uv[0]+1;
			nodes[index].uv[1] = nodes[index-1].uv[1];
			break;
		case NEG_I:
			nodes[index].uv[0] = nodes[index-1].uv[0]-1;
			nodes[index].uv[1] = nodes[index-1].uv[1];
			break;
		case POSI_J:
			nodes[index].uv[0] = nodes[index-1].uv[0];
			nodes[index].uv[1] = nodes[index-1].uv[1]+1;
			break;
		case NEG_J:
			nodes[index].uv[0] = nodes[index-1].uv[0];
			nodes[index].uv[1] = nodes[index-1].uv[1]-1;
			break;
		default:
			break;
	}	
}

//subdivide the surface 
void SubMapping::InteriorNodeInterpolation(ModelEnt *me)
{
	iBase_EntitySetHandle entityset;
	iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(me->geom_handle(), 0, entityset);
	IBERRCHK(r_err, "Trouble get the entity set for edge from geometric vertex.");

	int i_max = -1.0e5, i_min = 1.0e5, j_max = -1.0e5, j_min = 1.0e5;
	//calculate the bounding i-j: i_max, i_min, j_max, j_min, these, these extreme points should be on the geometric vertices
	for (unsigned int i = 0; i < sorted_vertex_list.size(); i++){
		int index = geom_mesh_node[sorted_vertex_list[i]];
		if (nodes[index].uv[0] < i_min)
			i_min = int(nodes[index].uv[0]);
		if (nodes[index].uv[0] > i_max)
			i_max = int(nodes[index].uv[0]);
		if (nodes[index].uv[1] < j_min)
			j_min = int(nodes[index].uv[1]);
		if (nodes[index].uv[1] > j_max)
			j_max = int(nodes[index].uv[1]);
	}
	vector<Vertex> boundary_ij;
	for (unsigned int i = 0; i < nodes.size(); i++)
		if (nodes[i].onCorner){
			boundary_ij.resize(boundary_ij.size()+1);
			boundary_ij[boundary_ij.size()-1] = nodes[i];
	}

	//i_min, i_max, j_min and j_max are the boundaries
	int index = nodes.size();
	int num_node_boundary = index;
	for (int i = i_min+1; i < i_max; i++){
		for (int j = j_min+1; j < j_max; j++){
			double j_coord[2] = {j_min, j_max}, i_coord[2] = {i_min, i_max};
			int j_index[2] = {-1, -1}, i_index[2] = {-1, -1};
			for (int k = 0; k < num_node_boundary; k++){
				if ((int(nodes[k].uv[0]) == i)&&(int(nodes[k].uv[1]) == j)){
					//this point is a boundary node, i_index and j_index are always -1
					for (int m = 0; m < 2; m++){
						i_index[m] = -1;
						j_index[m] = -1;
					}
					break;
				}
				else if ((int(nodes[k].uv[0]) == i)&&(int(nodes[k].uv[1]) != j)){
					//this point has vertical nodes
					if ((nodes[k].uv[1] > j)&&(nodes[k].uv[1] <= j_coord[1])){//top point
						j_index[1] = k;
						j_coord[1] = nodes[k].uv[1];
					}
					else if ((nodes[k].uv[1] < j)&&(nodes[k].uv[1] >= j_coord[0])){//bottom point
						j_index[0] = k;
						j_coord[0] = nodes[k].uv[1];
					}
					else
						continue;
				}
				else if ((int(nodes[k].uv[0]) != i)&&(int(nodes[k].uv[1]) == j)){
					//this point has horizontal nodes
					if ((nodes[k].uv[0] > i)&&(nodes[k].uv[0] <= i_coord[1])){
						i_index[1] = k;
						i_coord[1] = nodes[k].uv[0];
					}
					else if ((nodes[k].uv[0] < i)&&(nodes[k].uv[0] >= i_coord[0])){
						i_index[0] = k;
						i_coord[0] = nodes[k].uv[0];
					}
					else
						continue;
				}
				else//this point is not interesting for me
					continue;
			}
			if ((i_index[0]!= -1)&&(i_index[1] != -1)&&(j_index[0] != -1)&&(j_index[1] != -1)){
				//interpolate the interior point
				double xyz[3];
				double weight[2] = {fabs(i-nodes[i_index[1]].uv[0]), fabs(j-nodes[j_index[1]].uv[1])};
				if (fabs(j-nodes[j_index[0]].uv[1]) < fabs(j-nodes[j_index[1]].uv[1]))
					weight[1] = double(fabs(j-nodes[j_index[0]].uv[1]));
				if (fabs(i-nodes[i_index[0]].uv[0]) < fabs(i-nodes[i_index[1]].uv[0]))
					weight[0] = fabs(i-nodes[i_index[0]].uv[0]);
				double total = weight[0]+weight[1];				
				double i_weight = weight[1]/total, j_weight = weight[0]/total;
				
				for (int k = 0; k < 3; k++){
					double firstpart = nodes[j_index[1]].xyz[k]*fabs(j-nodes[j_index[0]].uv[1])/fabs(nodes[j_index[1]].uv[1]-nodes[j_index[0]].uv[1]);
					double secondpart = nodes[j_index[0]].xyz[k]*fabs(nodes[j_index[1]].uv[1]-j)/fabs(nodes[j_index[1]].uv[1]-nodes[j_index[0]].uv[1]);
					double j_value = firstpart + secondpart;
					double thirdpart = nodes[i_index[1]].xyz[k]*fabs(i-nodes[i_index[0]].uv[0])/fabs(nodes[i_index[1]].uv[0]-nodes[i_index[0]].uv[0]);
					double fourthpart = nodes[i_index[0]].xyz[k]*fabs(nodes[i_index[1]].uv[0]-i)/fabs(nodes[i_index[1]].uv[0]-nodes[i_index[0]].uv[0]);
					double i_value = thirdpart + fourthpart;
					xyz[k] = j_weight*firstpart+j_weight*secondpart+i_weight*thirdpart+i_weight*fourthpart;
				}
				iGeom::Error g_err = mk_core()->igeom_instance()->getEntClosestPt(me->geom_handle(), xyz[0], xyz[1], xyz[2], xyz[0], xyz[1], xyz[2]);
				IBERRCHK(g_err, "Trouble get the closest point on the surface.");
				if (isOutSideSurf(boundary_ij, i, j)){
				//if ((num_i[0]%2==1)||(num_i[1]%2==1)||(num_j[0]%2==1)||(num_j[1]%2==1))			
					//create the vertex entity
					index++;
					nodes.resize(index);
					nodes[index-1].index = index - 1;
					for (int k = 0; k < 3; k++)
						nodes[index-1].xyz[k] = xyz[k];
					nodes[index-1].onBoundary = false;
					nodes[index-1].onCorner = false;
					nodes[index-1].uv[0] = i;
					nodes[index-1].uv[1] = j;	
				}
			}
		}
	}
	

	std::cout << "number of boundary nodes is " << num_node_boundary << std::endl;
	//create the interior node entities
	vector<iBase_EntityHandle> m_nodes;
	m_nodes.resize(nodes.size()-num_node_boundary);
	for (unsigned int i = num_node_boundary; i < nodes.size(); i++){
		iMesh::Error m_err = mk_core()->imesh_instance()->createVtx(nodes[i].xyz[0], nodes[i].xyz[1], nodes[i].xyz[2], m_nodes[i-num_node_boundary]);
		IBERRCHK(m_err, "Trouble create the node entity.");
		nodes[i].gVertexHandle = m_nodes[i-num_node_boundary];
		m_err = mk_core()->imesh_instance()->addEntToSet(m_nodes[i-num_node_boundary], entityset);
		IBERRCHK(m_err, "Trouble add the mesh node into the entityset.");
		m_err = mk_core()->imesh_instance()->setIntData(m_nodes[i-num_node_boundary], m_taghandle, i);
		IBERRCHK(m_err, "Trouble set the int data.");
	}
	//ok, we are done with the node entities
	
	//trick here, create an array of int data to quickly locate the node index
	vector<int> node_index((j_max-j_min+1)*(i_max-i_min+1),-1);//here, invalid position is -1
	for (unsigned int i = 0; i < nodes.size(); i++){
		node_index[(nodes[i].uv[1]-j_min)*(i_max-i_min+1)+nodes[i].uv[0]-i_min] = i;
	}

	int face_sense;
	iGeom::Error g_err = mk_core()->igeom_instance()->getEgFcSense(edges[sorted_edge_list[0]].gEdgeHandle, me->geom_handle(), face_sense);
	IBERRCHK(g_err, "Trouble get the edge sense with respect to a geometrical face.");

	int edge_sense;
	g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[sorted_edge_list[0]].gEdgeHandle, vertices[sorted_vertex_list[0]].gVertexHandle, vertices[sorted_vertex_list[1]].gVertexHandle, edge_sense);
	IBERRCHK(g_err, "Trouble get the sense of edge with respect to two vertices.");

	int test_sense = face_sense * edge_sense;

	index = 0;
	//create the quadrilateral elements on the surface
	for (int j = j_min; j < j_max; j++){
		for (int i = i_min; i < i_max; i++){
			
			if (node_index[(j-j_min)*(i_max-i_min+1)+i-i_min] != -1){//there exists such a mesh node
				//check the other three vertices whether they exist or not
				if ((j+1)>j_max)
					continue;
				if ((i+1)>i_max)
					continue;
				if ((node_index[(j+1-j_min)*(i_max-i_min+1)+i-i_min] != -1)&&(node_index[(j-j_min)*(i_max-i_min+1)+i+1-i_min] != -1)&&(node_index[(j+1-j_min)*(i_max-i_min+1)+i+1-i_min] != -1)){
					index++;
					quads.resize(index);
					quads[index-1].connect.resize(4);
					if (test_sense < 0){
						quads[index-1].connect[0] = &nodes[node_index[(j-j_min)*(i_max-i_min+1)+i-i_min]];
						quads[index-1].connect[3] = &nodes[node_index[(j+1-j_min)*(i_max-i_min+1)+i-i_min]];
						quads[index-1].connect[2] = &nodes[node_index[(j+1-j_min)*(i_max-i_min+1)+i+1-i_min]];
						quads[index-1].connect[1] = &nodes[node_index[(j-j_min)*(i_max-i_min+1)+i+1-i_min]];
					}
					else{
						quads[index-1].connect[0] = &nodes[node_index[(j-j_min)*(i_max-i_min+1)+i-i_min]];
						quads[index-1].connect[1] = &nodes[node_index[(j+1-j_min)*(i_max-i_min+1)+i-i_min]];
						quads[index-1].connect[2] = &nodes[node_index[(j+1-j_min)*(i_max-i_min+1)+i+1-i_min]];
						quads[index-1].connect[3] = &nodes[node_index[(j-j_min)*(i_max-i_min+1)+i+1-i_min]];
					}					
				}
			}
		}
	}
	std::cout << "test quad size = " << quads.size() << std::endl;

	vector<iBase_EntityHandle> m_faces;
	m_faces.resize(quads.size());
	for (unsigned int i = 0; i < quads.size(); i++){
		vector<iBase_EntityHandle> n_nodes;
		n_nodes.resize(4);
		for (int k = 0; k < 4; k++)
			n_nodes[k] = nodes[quads[i].connect[k]->index].gVertexHandle;

		iMesh::Error m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &n_nodes[0], 4, m_faces[i]);
		IBERRCHK(m_err, "Trouble create the quadrilateral element.");
		quads[i].gFaceHandle = m_faces[i];
		m_err = mk_core()->imesh_instance()->addEntToSet(m_faces[i], entityset);
		IBERRCHK(m_err, "Trouble add the mesh node into the entityset.");
		m_err = mk_core()->imesh_instance()->setIntData(m_faces[i], m_taghandle, i);
		IBERRCHK(m_err, "Trouble set the int data.");

		n_nodes.clear();
	}

	std::cout << "Quad size = " << quads.size() << std::endl;	
	std::cout << "imax = " << i_max << "\timin = " << i_min << "\tjmax = " << j_max << "\tjmin = " << j_min << std::endl;
}

//smooth the structured mesh
void SubMapping::MeshSmoothing(ModelEnt *ent)
{
	std::vector<std::set<int> > AdjElements;
	std::vector<std::vector<int> > Quads(quads.size(), vector<int>(4));
	std::vector<std::vector<double> > coords(nodes.size(), vector<double>(3));
	std::vector<bool> isBnd(nodes.size(), false);
	std::vector<std::vector<int> > connect(nodes.size(), std::vector<int>(9));
	
	AdjElements.resize(nodes.size());
	for (unsigned int i = 0; i < quads.size(); i++){
	    std::vector<iBase_EntityHandle> adjEnts;
	    adjEnts.clear();
	    iMesh::Error m_err = mk_core()->imesh_instance()->getEntAdj(quads[i].gFaceHandle, iBase_VERTEX, adjEnts);
	    IBERRCHK(m_err, "Trouble get the adjacent nodes wrt a quad.");

	    assert(adjEnts.size()==4);

	    for (unsigned int j = 0; j < adjEnts.size(); j++){
			int index_id = -1;
			m_err = mk_core()->imesh_instance()->getIntData(adjEnts[j], m_taghandle, index_id);
			IBERRCHK(m_err, "Trouble get int data for nodes.");
			Quads[i][j] = index_id;
	    }
	}

	for (unsigned int i = 0; i < nodes.size(); i++){
		if (nodes[i].onBoundary || nodes[i].onCorner)
			isBnd[i] = true;
		for (int j = 0; j < 3; j++)
			coords[i][j] = nodes[i].xyz[j];
		if (!isBnd[i]){
			vector<iBase_EntityHandle> adjEnts;
			iMesh::Error m_err = mk_core()->imesh_instance()->getEntAdj(nodes[i].gVertexHandle, iBase_FACE, adjEnts);
			IBERRCHK(m_err, "Trouble get the adjacent quads wrt a node.");
			for (unsigned int j = 0; j < adjEnts.size(); j++){
				int index_id = -1;
				m_err = mk_core()->imesh_instance()->getIntData(adjEnts[j], m_taghandle, index_id);
				IBERRCHK(m_err, "Trouble get int data for quads.");
				AdjElements[i].insert(index_id);
			}

			//process the connect info
			//there are 4 adjacent quadrilateral elements around node i	
			assert(AdjElements[i].size() == 4);
			std::set<int>::iterator it = AdjElements[i].begin();
			int st_index[4];
			//process 4 quad elements
			int j = -1;
			for (; it != AdjElements[i].end(); it++){
				j++;				
				if (int(i) == Quads[*it][0])
					st_index[j] = 0;
				else if (int(i) == Quads[*it][1])
					st_index[j] = 1;
				else if (int(i) == Quads[*it][2])
					st_index[j] = 2;
				else
					st_index[j] = 3;
			}
			it = AdjElements[i].begin();
			connect[i][2] = Quads[*it][(st_index[0]+3)%4];
			connect[i][8] = Quads[*it][(st_index[0]+1)%4];
			connect[i][1] = Quads[*it][(st_index[0]+2)%4];
			//finish processing the quad 1
			std::set<int>::iterator it1 = AdjElements[i].begin();
			it1++;
			for (j = 1; j < 4; j++, it1++){
				if (connect[i][8] == Quads[*it1][(st_index[j]+1)%4]){
					connect[i][7] = Quads[*it1][(st_index[j]+2)%4];
					connect[i][6] = Quads[*it1][(st_index[j]+3)%4];
					break;
				}
				else if (connect[i][8] == Quads[*it1][(st_index[j]+3)%4]){
					connect[i][7] = Quads[*it1][(st_index[j]+2)%4];
					connect[i][6] = Quads[*it1][(st_index[j]+1)%4];					
					break;
				}
				else
					continue;
			}
			//finish processing the quad 2
			std::set<int>::iterator it2 = AdjElements[i].begin();
			it2++;
			for (j=1; it2 != AdjElements[i].end(); it2++,j++){
				if (connect[i][2] == Quads[*it2][(st_index[j]+1)%4]){
					connect[i][3] = Quads[*it2][(st_index[j]+2)%4];
					connect[i][4] = Quads[*it2][(st_index[j]+3)%4];
					break;
				}
				else if (connect[i][2] == Quads[*it2][(st_index[j]+3)%4]){
					connect[i][3] = Quads[*it2][(st_index[j]+2)%4];
					connect[i][4] = Quads[*it2][(st_index[j]+1)%4];					
					break;
				}
				else
					continue;
			}
			//finish processing the quad 4;
			std::set<int>::iterator it3 = AdjElements[i].begin();
			it3++;
			for (j=1; it3 != AdjElements[i].end(); it3++,j++){
				if ((it3 != it1)&&(it3 != it2)){
					connect[i][5] = Quads[*it2][(st_index[j]+2)%4];	
				}
				else
					continue;
			}	
		}
	}

	mk_core()->save_mesh("BeforeEquipotential.vtk");
	
	EquipotentialSmooth smoother;	
		
	smoother.SetupData(AdjElements, Quads, coords, isBnd, connect);
	smoother.Execute();
	
	std::vector<std::vector<double> > coors(nodes.size(), vector<double>(3));
	smoother.GetCoords(coors);

	//update the new position for nodes
	for (unsigned int i = 0; i < nodes.size(); i++){
		if (!isBnd[i]){
			double tmp_coord[3] = {coords[i][0], coords[i][1], coords[i][2]};
			
			iGeom::Error g_err = mk_core()->igeom_instance()->getEntClosestPt(ent->geom_handle(), coords[i][0], coords[i][1], coords[i][2], tmp_coord[0], tmp_coord[1], tmp_coord[2]);
			IBERRCHK(g_err, "Trouble get a closest point on the linking surface.");
			iMesh::Error m_err = mk_core()->imesh_instance()->setVtxCoord(nodes[i].gVertexHandle, tmp_coord[0], tmp_coord[1], tmp_coord[2]);
			IBERRCHK(m_err, "Trouble set the new coordinates for nodes.");
		}
	}
	

#ifdef HAVE_MESQUITE
	iBase_EntitySetHandle entityset;
	iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(ent->geom_handle(), 0, entityset);
	IBERRCHK(r_err, "Trouble get the entity set for edge from geometric vertex.");
	MeshImprove shapesmooth(mk_core(), true, true, false, false, false);
    shapesmooth.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);

#endif
  vector<iBase_EntityHandle> entities;
  for (unsigned int i = 0; i < nodes.size(); i++)
      entities.push_back(nodes[i].gVertexHandle);
  iMesh::Error m_err = mk_core()->imesh_instance()->rmvArrTag(&entities[0], entities.size(), m_taghandle);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  entities.clear();
  for (unsigned int i = 0; i < quads.size(); i++)
      entities.push_back(quads[i].gFaceHandle);
  m_err = mk_core()->imesh_instance()->rmvArrTag(&entities[0], entities.size(), m_taghandle);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  entities.clear();
  for (unsigned int i = 0; i < vertices.size(); i++)
	  entities.push_back(vertices[i].gVertexHandle);
  iGeom::Error g_err = mk_core()->igeom_instance()->rmvArrTag(&entities[0], entities.size(), g_taghandle);
  IBERRCHK(g_err, "Trouble remove the tag values from an array of entities.");
  entities.clear();
  for (unsigned int i = 0; i < edges.size(); i++)
	  entities.push_back(edges[i].gEdgeHandle);
  g_err = mk_core()->igeom_instance()->rmvArrTag(&entities[0], entities.size(), g_taghandle);
  IBERRCHK(g_err, "Trouble remove the tag values from an array of entities.");

	
}

//check whether a point is outside the surface or not
bool SubMapping::isOutSideSurf(vector<Vertex> corner, int i, int j)
{//test whether the point (i,j) is inside the polygon or outside
   //use the wind number
   unsigned int m;
   double angle=0;
   double p1[2],p2[2];

   for (m = 0; m < corner.size(); m++) {
      p1[0] = corner[m].uv[0] - i;
      p1[1] = corner[m].uv[1] - j;
      p2[0] = corner[(m+1)%corner.size()].uv[0] - i;
      p2[1] = corner[(m+1)%corner.size()].uv[1] - j;
      angle += Angle2D(p1[0],p1[1],p2[0],p2[1]);
   }

   if (fabs(angle) < PI)
      return false;
   else
      return true;
}

/*
   Return the angle between two vectors on a plane
   The angle is from vector 1 to vector 2, positive anticlockwise
   The result is between -pi -> pi
*/
double SubMapping::Angle2D(double x1, double y1, double x2, double y2)
{
   double dtheta,theta1,theta2;

   theta1 = atan2(y1,x1);
   theta2 = atan2(y2,x2);
   dtheta = theta2 - theta1;
   while (dtheta > PI)
      dtheta -= 2.0*PI;
   while (dtheta < -PI)
      dtheta += 2.0*PI;

   return dtheta;
}

}

