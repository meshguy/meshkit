#include "meshkit/TFIMapping.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "moab/ReadUtilIface.hpp"
#include <vector>
#include <iostream>
#include <math.h>
#include <map>


namespace MeshKit
{

//---------------------------------------------------------------------------//
//Entity Type initilization for TFIMapping meshing
moab::EntityType TFIMapping_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBQUAD, moab::MBMAXTYPE};
const moab::EntityType* TFIMapping::output_types()
  { return TFIMapping_tps; }

//---------------------------------------------------------------------------//
// construction function for TFIMapping class
TFIMapping::TFIMapping(MKCore *mk_core, const MEntVector &me_vec) : MeshScheme(mk_core, me_vec)
{
	//buildAssociation();
}

//---------------------------------------------------------------------------//
// deconstruction function for TFIMapping class
TFIMapping::~TFIMapping()
{

}

//---------------------------------------------------------------------------//
// setup function: 
void TFIMapping::setup_this()
{
	
}

//---------------------------------------------------------------------------//
// execute function: generate the all-quad mesh through the TFI mapping
void TFIMapping::execute_this()
{
	std::vector<double> coords;
	MEntVector edges;
	
	//loop over the surfaces
	for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  	{
    		ModelEnt *me = mit -> first;
		
		//first check whether the surface is meshed or not
		if (me->get_meshed_state() >= COMPLETE_MESH)
			continue;
		

		//test number of edges bounding the surface
		std::vector<int> senses, edge_sizes;
		me->boundary(0, edges, &senses, &edge_sizes);
		assert(4==(int)edges.size());

		

		SurfMapping(me);
		
      		//ok, we are done, commit to ME
    		me->commit_mesh(mit->second, COMPLETE_MESH);	
  	}

}



//****************************************************************************//
// function   : SurfMapping
// Date       : Mar 3, 2011
// Description: Generate the mesh on the linking surface by using TFI
//***************************************************************************//
int TFIMapping::SurfMapping(ModelEnt *ent)
{
	std::vector<iBase_EntityHandle> edges, gNode;

	//get the four geometrical edges
	iGeom::Error g_err = ent->igeom_instance()->getEntAdj(ent->geom_handle(), iBase_EDGE, edges);
	IBERRCHK(g_err, "Trouble get the adjacent geometric edges on a surface.");

	//get the four corners on a surface
	g_err = ent->igeom_instance()->getEntAdj(ent->geom_handle(), iBase_VERTEX, gNode);
	IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a surface.");


	iBase_TagHandle taghandle;
	g_err = ent->igeom_instance()->createTag("TFIMapping", 1, iBase_INTEGER, taghandle);
	IBERRCHK(g_err, "Trouble create the taghandle for the surface.");
	
	//it has already been checked that there are 4 edges bounding a surface in the execute function
	std::vector<iBase_EntitySetHandle> EdgeMeshSets(4);
	std::vector<int>  num(edges.size());
	//get the edge mesh entity sets for 4 edges
	for (unsigned int i = 0; i < edges.size(); i++)
	{
		iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(edges[i], 0, EdgeMeshSets[i]);
		IBERRCHK(r_err, "Trouble get@ the mesh entity set of geometrical edge.");

		std::vector<iBase_EntityHandle> mEdges;
		mEdges.clear();
		iMesh::Error m_err = mk_core()->imesh_instance()->getEntities(EdgeMeshSets[i], iBase_EDGE, iMesh_LINE_SEGMENT, mEdges);
		IBERRCHK(m_err, "Trouble get the line segments in the edge mesh.");
		
		num[i] = mEdges.size();
		//iMesh::Error m_err = mk_core()->imesh_instance()->getNumOfTopo(EdgeMeshSets[i], iMesh_LINE_SEGMENT, num[i]);
		//IBERRCHK(m_err, "Trouble get the number of line segments in the edge mesh.");
		
		//set the new int value for edges and nodes
		g_err = ent->igeom_instance()->setIntData(edges[i], taghandle, i);
		IBERRCHK(g_err, "Trouble create the taghandle for the surface.");

		g_err = ent->igeom_instance()->setIntData(gNode[i], taghandle, i);
		IBERRCHK(g_err, "Trouble create the taghandle for the surface.");
	}

	//find the corresponding edge and node relationship: based on the internal constraints, there are two edges with equal number of line segments
	int num_i=0, num_j=0;
	std::map<int, int> edgePair, nodePair, adjEdgesOfNode, oppNodeOfNode;
	num_i = num[0];

	std::vector<iBase_EntityHandle> nNodes, nEdges;
	int index_a, index_b;
	//nodes' index for edge 0, 
	//node_index_a denotes one node in edge 0, node_index_b denotes the corresponding node on the opposite edge
	int node_index_a, node_index_b;	
	//pick up the edge 0, find its two end nodes	
	g_err = ent->igeom_instance()->getEntAdj(edges[0], iBase_VERTEX, nNodes);
	IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
	
	//find the int for edge 0
	g_err = ent->igeom_instance()->getIntData(edges[0], taghandle, index_a);
	IBERRCHK(g_err, "Trouble get the int data on edge 0.");
	
	//pick up nodes[0], find its adjacent edges
	g_err = ent->igeom_instance()->getEntAdj(nNodes[0], iBase_EDGE, nEdges);
	IBERRCHK(g_err, "Trouble get the adjacent geometric edges on a node.");

	g_err = ent->igeom_instance()->getIntData(nNodes[0], taghandle, node_index_a);
	IBERRCHK(g_err, "Trouble get the int data fro node 0.");	

	//find the index for one of the adjacent edges 
	g_err = ent->igeom_instance()->getIntData(nEdges[0], taghandle, index_b);
	IBERRCHK(g_err, "Trouble get the int data on edge 0.");	
	if (index_a != index_b)
	{
		//choose nEdge[0] to find the opposite node to nNodes[0]
		adjEdgesOfNode[0] = index_b;
		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;		
		g_err = ent->igeom_instance()->getEntAdj(nEdges[0], iBase_VERTEX, tnodes);
		IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
		
		g_err = ent->igeom_instance()->getIntData(tnodes[0], taghandle, node_index_b);
		IBERRCHK(g_err, "Trouble get the int data of node.");
		if (node_index_b != node_index_a)
			oppNodeOfNode[node_index_a] = node_index_b;			
			//oppNodeOfNode[node_index_a] = node_index_b;
		else
		{				
			g_err = ent->igeom_instance()->getIntData(tnodes[1], taghandle, node_index_b);
			IBERRCHK(g_err, "Trouble get the int data of node.");
			//oppNodeOfNode.insert( pair<int, int>(node_index_a, node_index_b));			
			oppNodeOfNode[node_index_a] = node_index_b;
		}	
	}
	else
	{//choose nEdges[1]
		g_err = ent->igeom_instance()->getIntData(nEdges[1], taghandle, index_b);
		IBERRCHK(g_err, "Trouble get the int data on edge 0.");
		//adjEdgesOfNode.insert ( pair<int,int>(0,index_b) );		
		adjEdgesOfNode[0] = index_b;

		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;		
		g_err = ent->igeom_instance()->getEntAdj(nEdges[1], iBase_VERTEX, tnodes);
		IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
		
		g_err = ent->igeom_instance()->getIntData(tnodes[0], taghandle, node_index_b);
		IBERRCHK(g_err, "Trouble get the int data of node.");
		if (node_index_b != node_index_a)
			//oppNodeOfNode.insert( pair<int, int>(node_index_a, node_index_b));
			oppNodeOfNode[node_index_a] = node_index_b;
		else
		{				
			g_err = ent->igeom_instance()->getIntData(tnodes[1], taghandle, node_index_b);
			IBERRCHK(g_err, "Trouble get the int data of node.");
			//oppNodeOfNode.insert( pair<int, int>(node_index_a, node_index_b));
			oppNodeOfNode[node_index_a] = node_index_b;
		}	
	}

	int node_index_c, node_index_d;

	//pick up nodes[1] find its adjacent edges
	nEdges.clear();	
	g_err = ent->igeom_instance()->getEntAdj(nNodes[1], iBase_EDGE, nEdges);
	IBERRCHK(g_err, "Trouble get the adjacent geometric edges on a node.");

	g_err = ent->igeom_instance()->getIntData(nNodes[1], taghandle, node_index_c);
	IBERRCHK(g_err, "Trouble get the int data fro node 0.");	

	//find the index for one of the adjacent edges 
	g_err = ent->igeom_instance()->getIntData(nEdges[0], taghandle, index_b);
	IBERRCHK(g_err, "Trouble get the int data on edge 0.");	
	if (index_a != index_b)
	{
		//choose nEdge[0] to find the opposite node to nNodes[0]
		adjEdgesOfNode[1] = index_b;
		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;		
		g_err = ent->igeom_instance()->getEntAdj(nEdges[0], iBase_VERTEX, tnodes);
		IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
		
		g_err = ent->igeom_instance()->getIntData(tnodes[0], taghandle, node_index_d);
		IBERRCHK(g_err, "Trouble get the int data of node.");
		if (node_index_d != node_index_c)
			oppNodeOfNode[node_index_c] = node_index_d;
		else
		{				
			g_err = ent->igeom_instance()->getIntData(tnodes[1], taghandle, node_index_d);
			IBERRCHK(g_err, "Trouble get the int data of node.");
			oppNodeOfNode[node_index_c] = node_index_d;
		}	
	}
	else
	{//choose nEdges[1]
		g_err = ent->igeom_instance()->getIntData(nEdges[1], taghandle, index_b);
		IBERRCHK(g_err, "Trouble get the int data on edge 0.");
		adjEdgesOfNode[1] = index_b;

		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;		
		g_err = ent->igeom_instance()->getEntAdj(nEdges[1], iBase_VERTEX, tnodes);
		IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
		
		g_err = ent->igeom_instance()->getIntData(tnodes[0], taghandle, node_index_d);
		IBERRCHK(g_err, "Trouble get the int data of node.");
		if (node_index_d != node_index_c)
			oppNodeOfNode[node_index_c] = node_index_d;
		else
		{				
			g_err = ent->igeom_instance()->getIntData(tnodes[1], taghandle, node_index_d);
			IBERRCHK(g_err, "Trouble get the int data of node.");
			oppNodeOfNode[node_index_c] = node_index_d;
		}	
	}

	num_j = num[adjEdgesOfNode[0]];

	//find the opposite geometric edge
	for (int i = 0; i < 4; i++)
		if ((i != 0) && (i != adjEdgesOfNode[0]) && (i != adjEdgesOfNode[1]))
			edgePair[0] = i;

	//finish all the geometry initialization
	
	//initialize the mesh
	std::vector<iBase_EntityHandle> NodeList_i, NodeList_j, NodeList_ii, NodeList_jj;
	std::vector<iBase_EntityHandle> corner(4);	
	//get mesh node lists from four edges
	//note: if we use getEntEntArrRelation, nodes and edges will be returned. It will return all the mesh entities related to mesh entity
	iBase_EntitySetHandle TmpSet;

	iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(edges[0], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical edge 0.");	
	iMesh::Error m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, NodeList_i);	
	IBERRCHK(m_err, "Trouble get the mesh node list from the geometrical edge 0.");

	r_err = mk_core()->irel_pair()->getEntSetRelation(edges[edgePair[0]], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical edge opposite to edge 0.");
	m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, NodeList_ii);	
	IBERRCHK(m_err, "Trouble get the mesh node list from the geometrical edge opposite to edge 0.");

	r_err = mk_core()->irel_pair()->getEntSetRelation(edges[adjEdgesOfNode[0]], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical edge adjacent to node 0.");
	m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, NodeList_j);	
	IBERRCHK(m_err, "Trouble get the mesh node list from the geometrical edge adjacent to node 0.");

	r_err = mk_core()->irel_pair()->getEntSetRelation(edges[adjEdgesOfNode[1]], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical edge adjacent to node 1.");
	m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, NodeList_jj);	
	IBERRCHK(m_err, "Trouble get the mesh node list from the geometrical edge adjacent to node 1.");
	

	//get mesh nodes from four corners
	std::vector<iBase_EntityHandle> tmpNodeHandle;
	r_err = mk_core()->irel_pair()->getEntSetRelation(gNode[node_index_a], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical corner 0.");
	m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, tmpNodeHandle);	
	IBERRCHK(m_err, "Trouble get the mesh node from the geometrical corner 0.");
	assert(tmpNodeHandle.size()==1);	
	corner[0] = tmpNodeHandle[0];

	r_err = mk_core()->irel_pair()->getEntSetRelation(gNode[node_index_c], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical corner 1.");
	tmpNodeHandle.clear();
	m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, tmpNodeHandle);	
	IBERRCHK(m_err, "Trouble get the mesh node from the geometrical corner 1.");
	assert(tmpNodeHandle.size()==1);	
	corner[1] = tmpNodeHandle[0];

	r_err = mk_core()->irel_pair()->getEntSetRelation(gNode[oppNodeOfNode[node_index_a]], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical corner 2.");
	tmpNodeHandle.clear();
	m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, tmpNodeHandle);	
	IBERRCHK(m_err, "Trouble get the mesh node from the geometrical corner 2.");
	assert(tmpNodeHandle.size()==1);	
	corner[2] = tmpNodeHandle[0];

	r_err = mk_core()->irel_pair()->getEntSetRelation(gNode[oppNodeOfNode[node_index_c]], 0, TmpSet);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical corner 3.");
	tmpNodeHandle.clear();
	m_err = mk_core()->imesh_instance()->getEntities(TmpSet, iBase_VERTEX, iMesh_POINT, tmpNodeHandle);	
	IBERRCHK(m_err, "Trouble get the mesh node from the geometrical corner 3.");
	assert(tmpNodeHandle.size()==1);	
	corner[3] = tmpNodeHandle[0];
	
	//get the sense of vertex pair with respect to edge
	std::vector<int> senses(4);
	g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[0], gNode[node_index_a], gNode[node_index_c], senses[0]);
	IBERRCHK(g_err, "Trouble get the sense of vertex pair with respect to edge 0.");

	g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[adjEdgesOfNode[0]], gNode[node_index_a], gNode[oppNodeOfNode[node_index_a]], senses[1]);
	IBERRCHK(g_err, "Trouble get the sense of vertex pair with respect to adjacent edge of node 0.");

	g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[adjEdgesOfNode[1]], gNode[node_index_c], gNode[oppNodeOfNode[node_index_c]], senses[2]);
	IBERRCHK(g_err, "Trouble get the sense of vertex pair with respect to adjacent edge of node 0.");

	g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[edgePair[0]], gNode[oppNodeOfNode[node_index_a]], gNode[oppNodeOfNode[node_index_c]], senses[3]);
	IBERRCHK(g_err, "Trouble get the sense of vertex pair with respect to adjacent edge of node 0.");
	
	//reorgainize the node list based on its sense
	std::vector<iBase_EntityHandle> List_i(NodeList_i.size()+2), List_j(NodeList_j.size()+2), List_ii(NodeList_ii.size()+2), List_jj(NodeList_jj.size()+2);
	List_i[0] = corner[0];
	List_ii[0] = corner[2];
	List_j[0] = corner[0];
	List_jj[0] = corner[1];

	List_i[List_i.size()-1] = corner[1];//node List for edge 0
	List_ii[List_ii.size()-1] = corner[3];
	List_j[List_j.size()-1] = corner[2];
	List_jj[List_jj.size()-1] = corner[3];	
	for (unsigned int i = 0; i < NodeList_i.size(); i++)
	{//list i and list ii should have the same number of nodes
		if (senses[0] == 1)
			List_i[i+1] = NodeList_i[i];
		else
			List_i[i+1] = NodeList_i[NodeList_i.size()-1-i];
		if (senses[3] == 1)
			List_ii[i+1] = NodeList_ii[i];
		else
			List_ii[i+1] = NodeList_ii[NodeList_ii.size()-1-i];
	}
	for (unsigned int i = 0; i < NodeList_j.size(); i++)
	{//list j and list jj should have the same number of nodes
		if (senses[1] == 1)
			List_j[i+1] = NodeList_j[i];
		else
			List_j[i+1] = NodeList_j[NodeList_j.size()-1-i];
		if (senses[2] == 1)
			List_jj[i+1] = NodeList_jj[i];
		else
			List_jj[i+1] = NodeList_jj[NodeList_jj.size()-1-i];
	}

	//ok, done with all the initializations

	double umax, vmax, umin, vmin;
	g_err =  ent->igeom_instance()->getEntUVRange(ent->geom_handle(), umin, vmin, umax, vmax);
  	IBERRCHK(g_err, "Trouble get parameter uv range for face.");
	
	//calculate the interior nodes based on TFI
	std::vector<iBase_EntityHandle> InteriorNodes((NodeList_j.size())*(NodeList_i.size()));
	for (unsigned int j = 1; j < (List_j.size()-1); j++)
	{
			//                             Node 2 (List_ii)
			//		  		 |		 
			//	Node 0 (List_j)---------Node------------Node 1 (List_jj)	
			//		  	 	 |		 
			//		               Node 3 (List_i)
		double coords[3];
		double uv_0[2], uv_1[2], uv_2[2], uv_3[2];
		
		g_err = mk_core()->imesh_instance()->getVtxCoord(List_j[j], coords[0], coords[1], coords[2]);
		IBERRCHK(g_err, "Trouble get the xyz coordinates for node 0.");
		g_err = mk_core()->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), coords[0], coords[1], coords[2], uv_0[0], uv_0[1]);
		IBERRCHK(g_err, "Trouble get the uv parametric coordinates from xyz coordinates for node 0.");

		g_err = mk_core()->imesh_instance()->getVtxCoord(List_jj[j], coords[0], coords[1], coords[2]);
		IBERRCHK(g_err, "Trouble get the xyz coordinates for node 1.");
		g_err = mk_core()->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), coords[0], coords[1], coords[2], uv_1[0], uv_1[1]);
		IBERRCHK(g_err, "Trouble get the uv parametric coordinates from xyz coordinates for node 1.");
		for (unsigned int i = 1; i < (List_i.size()-1); i++)
		{

			g_err = mk_core()->imesh_instance()->getVtxCoord(List_i[i], coords[0], coords[1], coords[2]);
			IBERRCHK(g_err, "Trouble get the xyz coordinates for node 3.");
			g_err = mk_core()->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), coords[0], coords[1], coords[2], uv_3[0], uv_3[1]);
			IBERRCHK(g_err, "Trouble get the uv parametric coordinates from xyz coordinates for node 3.");

			g_err = mk_core()->imesh_instance()->getVtxCoord(List_ii[i], coords[0], coords[1], coords[2]);
			IBERRCHK(g_err, "Trouble get the xyz coordinates for node 2.");
			g_err = mk_core()->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), coords[0], coords[1], coords[2], uv_2[0], uv_2[1]);
			IBERRCHK(g_err, "Trouble get the uv parametric coordinates from xyz coordinates for node 2.");

			double uv[2], r, s;
			r = ((uv_2[0]+uv_3[0])/2 - umin)/(umax - umin);
			s = ((uv_0[1]+uv_1[1])/2 - vmin)/(vmax - vmin);
			uv[0] = parametricTFI2D(r, s, uv_0[0], uv_1[0], uv_3[0], uv_2[0]);
			uv[1] = parametricTFI2D(r, s, uv_0[1], uv_1[1], uv_3[1], uv_2[1]);

			g_err = mk_core()->igeom_instance()->getEntUVtoXYZ(ent->geom_handle(), uv[0], uv[1], coords[0], coords[1], coords[2]);
			IBERRCHK(g_err, "Trouble get the xyz coordinates from uv parametric coordinates for interior node.");	

			iMesh::Error m_err = mk_core()->imesh_instance()->createVtx(coords[0], coords[1], coords[2], InteriorNodes[(j-1)*NodeList_i.size() + i -1]);
			IBERRCHK(m_err, "Trouble get create the interior node.");			
		}
	}

	//finish creating the interior nodes
	iBase_EntitySetHandle entityset;
	r_err = mk_core()->irel_pair()->getEntSetRelation(ent->geom_handle(), 0, entityset);
	if (r_err)
	{
		//create the entityset
		iMesh::Error m_err = mk_core()->imesh_instance()->createEntSet(true, entityset);
		IBERRCHK(m_err, "Trouble create the entity set.");

		r_err = mk_core()->irel_pair()->setEntSetRelation(ent->geom_handle(), entityset);
		IBERRCHK(r_err, "Trouble create the association between the geometry and mesh entity set.");
	}
	
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&InteriorNodes[0],InteriorNodes.size(), entityset);
	IBERRCHK(m_err, "Trouble add an array of entities to the mesh entity set.");
	
	//create the quads
	std::vector<iBase_EntityHandle> Quads((NodeList_j.size()+1)*(NodeList_i.size()+1));		
	for (unsigned int j = 0; j < (NodeList_j.size()+1); j++)
	{
		std::vector<iBase_EntityHandle> qNodes(4);
		if (j == 0)
		{//starting row boundary
			//first quad on first colume and first row
			qNodes[0] = List_i[0];
			qNodes[1] = List_j[1];
			qNodes[2] = InteriorNodes[0];
			qNodes[3] = List_i[1];
				
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[0]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
				
			for (unsigned int i = 1; i < (NodeList_i.size()); i++)
			{
				qNodes[0] =  List_i[i];
				qNodes[1] =  InteriorNodes[i-1];
				qNodes[2] =  InteriorNodes[i];
				qNodes[3] =  List_i[i+1];
				
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+i]);
				IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			}
			
			qNodes[0] = List_i[List_i.size()-2];
			qNodes[1] = InteriorNodes[List_i.size()-3];
			qNodes[2] = List_jj[1];
			qNodes[3] = List_jj[0];

			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[NodeList_i.size()]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
		}
		
		else if (j == NodeList_j.size())
		{//ending row boundary
			qNodes[0] = List_j[j];
			qNodes[1] = List_j[j+1];
			qNodes[2] = List_ii[1];
			qNodes[3] = InteriorNodes[(j-1)*(NodeList_i.size())];
				
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
				
			for (unsigned int i = 1; i < (NodeList_i.size()); i++)
			{
				qNodes[0] = InteriorNodes[(j-1)*(NodeList_i.size())+ i - 1];
				qNodes[1] = List_ii[i];
				qNodes[2] = List_ii[i+1];
				qNodes[3] = InteriorNodes[(j-1)*(NodeList_i.size())+ i];
				
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+i]);
				IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			}
			
			qNodes[0] = InteriorNodes[(j-1)*(NodeList_i.size())+NodeList_i.size() - 1];
			qNodes[1] = List_ii[j];
			qNodes[2] = List_ii[j+1];
			qNodes[3] = List_jj[List_jj.size()-2];
				
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[(NodeList_j.size()+1)*(NodeList_i.size()+1) - 1]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			
		}
		else
		{
			//first column
			qNodes[0] = List_j[j];
			qNodes[1] = List_j[j+1];
			qNodes[2] = InteriorNodes[j*NodeList_i.size()];
			qNodes[3] = InteriorNodes[(j-1)*NodeList_i.size()];
			
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			for (unsigned int i = 1; i < (NodeList_i.size()); i++)
			{
				qNodes[0] = InteriorNodes[(j-1)*NodeList_i.size() + i-1];
				qNodes[1] = InteriorNodes[j*NodeList_i.size() + i-1];
				qNodes[2] = InteriorNodes[j*NodeList_i.size() + i];
				qNodes[3] = InteriorNodes[(j-1)*NodeList_i.size() + i];
				
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+i]);
				IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			}

			//end column
			qNodes[0] = InteriorNodes[(j-1)*NodeList_i.size() + NodeList_i.size() - 1];
			qNodes[1] = InteriorNodes[j*NodeList_i.size() + NodeList_i.size() - 1];
			qNodes[2] = List_jj[j+1];
			qNodes[3] = List_jj[j];
			
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+NodeList_i.size()]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
		}

	}

	//finish creating the quads
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&Quads[0], Quads.size(), entityset);
	IBERRCHK(m_err, "Trouble add an array of quads to the mesh entity set.");
	
	return 1;
}



//****************************************************************************//
// function   : parametricTFI2D
// Date       : Feb 15, 2011
// Description: do the transfinite interpolation in (pt_0s, pt_1s), (pt_r0, pt_r1)
//***************************************************************************//
double TFIMapping::parametricTFI2D(double r, double s, double pt_0s, double pt_1s, double pt_r0, double pt_r1)
{
			//                             pt_r1
			//		  		 |		 
			//		pt_0s---------Node------------pt_1s	
			//		  	 	 |		 
			//		               pt_r0

	assert(r>= 0 && r <= 1.0);
	assert(s>= 0 && s <= 1.0);
	double pt_rs;

	//interpolate the pt_rs based on pt_r0, pt_r1, pt_0s and pt_1s
	pt_rs = 0.5*( linear_interpolation( s, pt_r0, pt_r1 ) + linear_interpolation(r, pt_0s, pt_1s) );

	return pt_rs;
}


//****************************************************************************//
// function   : linear_interpolation 
// Date       : Feb 15, 2011
// Description: interpolate linearly between x0 and x1
//***************************************************************************//
double TFIMapping::linear_interpolation(double r, double x0, double x1)
{
	assert(r >=0 && r <= 1.0);
	double pt= (1-r)*x0 + r*x1;
	return pt;
}

void TFIMapping::buildAssociation()
{
    std::vector<iBase_EntitySetHandle> entitySets;
    iBase_EntitySetHandle geom_root_set, mesh_root_set;

    // Get the root sets of the geometry and mesh.
    mesh_root_set = mk_core()->igeom_instance()->getRootSet();
    geom_root_set = mk_core()->imesh_instance()->getRootSet();

    //get the global geometrical id tag, global mesh id tag, global dimension id tag
    iBase_TagHandle geom_id_tag, mesh_id_tag, geom_dim_tag;    
    iGeom::Error g_err = mk_core()->igeom_instance()->getTagHandle("GLOBAL_ID", geom_id_tag);
    IBERRCHK(g_err, "Trouble get global geometry dimension id tag.");
    iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("GLOBAL_ID", mesh_id_tag);
    IBERRCHK(m_err, "Trouble get global mesh dimension id tag.");
    m_err = mk_core()->imesh_instance()->getTagHandle("GEOM_DIMENSION", geom_dim_tag);
    IBERRCHK(m_err, "Trouble get the geometric dimension id tag.");



    // Get all the entitySet in the mesh
    m_err = mk_core()->imesh_instance()->getEntSets(mesh_root_set, 0, entitySets);
    IBERRCHK(m_err, "Trouble get the geometric dimension id tag.");


    	//int ncount;
    	//iBase_EntityHandle gEntity;



    // Map all the geometric nodes
    //get size of nodes
    std::vector<iBase_EntityHandle> gNodes;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_VERTEX, gNodes);
    IBERRCHK(g_err, "Trouble get the geometric node entities.");
    
    int geom_id;
    std::map<int, iBase_EntityHandle> mapNodes;
    for (unsigned int i=0; i < gNodes.size(); i++)
    {
    	g_err = mk_core()->igeom_instance()->getIntData(gNodes[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of node entities.");
        mapNodes[geom_id] = gNodes[i];
    }
    

    // Map all the geometric edges.resize
    std::vector<iBase_EntityHandle> gEdges;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_EDGE, gEdges);
    IBERRCHK(g_err, "Trouble get the geometric edge entities.");

    
    std::map<int, iBase_EntityHandle> mapEdges;
    for (unsigned int i = 0; i < gEdges.size(); i++)
    {
        g_err = mk_core()->igeom_instance()->getIntData(gEdges[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of edge entities.");
        mapEdges[geom_id] = gEdges[i];
    }

    // Map all the geometric faces ...
    std::vector<iBase_EntityHandle> gFaces;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_FACE, gFaces);
    IBERRCHK(g_err, "Trouble get the geometric face entities.");

    std::map<int, iBase_EntityHandle> mapFaces;
    for (unsigned int i = 0; i < gFaces.size(); i++)
    {
        g_err = mk_core()->igeom_instance()->getIntData(gFaces[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of face entities.");
        mapFaces[geom_id] = gFaces[i];
    }

    // Map all the geometric cells ...
    std::vector<iBase_EntityHandle> gCells;
    g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_REGION, gCells);
    IBERRCHK(g_err, "Trouble get the geometric cell entities.");

    std::map<int, iBase_EntityHandle> mapCells;
    for (unsigned int i = 0; i < gCells.size(); i++)
    {
        g_err = mk_core()->igeom_instance()->getIntData(gCells[i], geom_id_tag, geom_id);
        IBERRCHK(g_err, "Trouble get the int data of cell entities.");
        mapCells[geom_id] = gCells[i];
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create Vertex Assocations:
    ///////////////////////////////////////////////////////////////////////////////
    cout << " Building Vertex Associations " << endl;
    int numNodes, ncount = 0, geom_dim;
    iBase_EntityHandle gEntity;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_VERTEX, numNodes);
    IBERRCHK(m_err, "Trouble get the number of vertices.");
    std::vector<iBase_EntityHandle> mNodes;
    
    int numAssociations = 0;
    
    for (unsigned int i = 0; i < entitySets.size(); i++)
    {
        mNodes.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_VERTEX, iMesh_ALL_TOPOLOGIES, mNodes);
        IBERRCHK(m_err, "Trouble get the number of vertices.");	
        if (mNodes.size() && (int(mNodes.size()) != numNodes))
        {
            ncount += mNodes.size();

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the mesh id tag with int value for entity set.");
	
            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
            IBERRCHK(m_err, "Trouble set the geom dim tag with int value for mesh entity set.");

            gEntity = 0;
            switch (geom_dim)
            {
            case 0:
            	if (mapNodes.find(geom_id) != mapNodes.end())
                {
                    gEntity = mapNodes[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            	
            case 1:

                if (mapEdges.find(geom_id) != mapEdges.end())
                {
                    gEntity = mapEdges[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                    gEntity = mapFaces[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Face not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Cell not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            default:
                cout << "Error: Invalid geometric dimension " << geom_dim << endl;
                exit(0);
            }

            if (gEntity)
            {
                iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
            }   
        }
    }

    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create Edge Assocations:
    ///////////////////////////////////////////////////////////////////////////////
    cout << " Building Edge Associations " << endl;

    int numEdges;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_EDGE, numEdges);
    IBERRCHK(m_err, "Trouble get the number of edges.");	

    std::vector<iBase_EntityHandle> mEdges;

    numAssociations = 0;
    ncount = 0;
    for (unsigned int i = 0; i < entitySets.size(); i++)
    {
        mEdges.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES, mEdges);
        IBERRCHK(m_err, "Trouble get the edge entities from the entity set.");

        if (mEdges.size() && (int(mEdges.size()) != numEdges))
        {
            ncount += mEdges.size();

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the geom_id for entity set in the edge association.");

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
            IBERRCHK(m_err, "Trouble set the geom_dim for entity set in the edge association.");
	    
            gEntity = 0;
            switch (geom_dim)
            {
            
            case 1:

                if (mapEdges.find(geom_id) != mapEdges.end())
                {
                    gEntity = mapEdges[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                    gEntity = mapFaces[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Face not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Cell not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            default:
                cout << "Error: Invalid geometric dimension " << geom_dim << endl;
                exit(0);
            }

            if (gEntity)
            {
                iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
            }
        }
    }

    if (numAssociations != int(mapEdges.size()))
        cout << "Warning: There are more edge entitySet than geometric edges " << endl;
    
    //////////////////////////////////////////////////////////////////////////////
    // Face Association
    //////////////////////////////////////////////////////////////////////////////
    cout << " Building Face Associations " << endl;

    std::vector<iBase_EntityHandle> mFaces;

    int numFaces;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_FACE, numFaces);
    IBERRCHK(m_err, "Trouble get the number of faces.");

    ncount = 0;
    numAssociations = 0;
    for (unsigned int i = 0; i < entitySets.size(); i++)
    {
        mFaces.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES, mFaces);
        IBERRCHK(m_err, "Trouble get the mesh face entities.");

        if (mFaces.size() && (int(mFaces.size()) != numFaces))
        {
            ncount += mFaces.size();
            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the geom_id for entity set in the face association.");

            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
            IBERRCHK(m_err, "Trouble set the geom_dim for entity set in the face association.");

            gEntity = 0;
            switch (geom_dim)
            {
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                {
                    gEntity = mapFaces[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                { IBERRCHK(m_err, "Trouble set the geom_id for entity set in the face association.");
                    cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                    exit(0);
                }
                break;
            }
            if (gEntity)
	    {
                iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
	    }        
	}
    }
    
    if (numAssociations != int(mapFaces.size()))
    {
	cout << "Warning: There are more face entitySet than geometric faces " << endl;
    }
	
    //////////////////////////////////////////////////////////////////////////////
    // Cell Association
    //////////////////////////////////////////////////////////////////////////////

    std::vector<iBase_EntityHandle> mCells;

    int numCells;
    m_err = mk_core()->imesh_instance()->getNumOfType(mesh_root_set, iBase_REGION, numCells);
    IBERRCHK(m_err, "Trouble get the number of cells.");

    ncount = 0; IBERRCHK(m_err, "Trouble set the geom_id for entity set in the face association.");
    for (unsigned int i = 0; i < entitySets.size(); i++)
    {
        mCells.clear();
        m_err = mk_core()->imesh_instance()->getEntities(entitySets[i], iBase_REGION, iMesh_ALL_TOPOLOGIES, mCells);
        IBERRCHK(m_err, "Trouble get the mesh cell entities.");

        if (mCells.size() && (int(mCells.size()) != numCells))
        {
            ncount += mCells.size();
            m_err = mk_core()->imesh_instance()->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
            IBERRCHK(m_err, "Trouble set the geom_id for entity set in the cell association.");

            if (mapCells.find(geom_id) != mapCells.end())
            {
                if (mapCells.find(geom_id) != mapCells.end())
                {
                    gEntity = mapCells[geom_id];
                    iRel::Error r_err = mk_core()->irel_pair()->setEntSetRelation(gEntity, entitySets[i]);
                    IBERRCHK(r_err, "Trouble set the association between the entity handle and mesh entity set.");	
                }
            }
        }
    }
}



}

