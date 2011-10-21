#include "meshkit/TFIMapping.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "moab/ReadUtilIface.hpp"
#ifdef HAVE_MESQUITE
#include "meshkit/MeshImprove.hpp"
#endif

#include <vector>
#include <iostream>
#include <math.h>
#include <map>

const double EPS = 1.0e-6;


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
		edges.clear();
		me->boundary(1, edges, &senses, &edge_sizes);
		//assert(4==(int)edges.size());

		

		SurfMapping(me);
		
      		//ok, we are done, commit to ME
    		me->commit_mesh(mit->second, COMPLETE_MESH);
  	}

}



/***********************************************************************************/
/*function   : SurfMapping                                                         */
/*Date       : Mar 3, 2011                                                         */
/*Description: Generate the mesh on the linking surface by using TFI               */
/*  prepare to generate the surface by using TFI mapping interpolation             */
/*  1. Get the mesh(edge mesh) from the bounding geometric edges                   */
/*  2. Find the corresponding relationship between edges, vertices                 */
/*  3. Check the nodes' corresponding relationship on the 4 bounding edges         */
/*  4. Do the TFI interpolation for interior nodes' location                       */
/***********************************************************************************/
int TFIMapping::SurfMapping(ModelEnt *ent)
{
	std::vector<iBase_EntityHandle> edges, gNode;	

	//get the four geometrical edges
	iGeom::Error g_err = ent->igeom_instance()->getEntAdj(ent->geom_handle(), iBase_EDGE, edges);
	IBERRCHK(g_err, "Trouble get the adjacent geometric edges on a surface.");

	//get the four corners on a surface
	g_err = ent->igeom_instance()->getEntAdj(ent->geom_handle(), iBase_VERTEX, gNode);
	IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a surface.");

	//create a tag for linking surface
	iBase_TagHandle taghandle;
	g_err = ent->igeom_instance()->getTagHandle("TFIMapping", taghandle);
	if (g_err){
		g_err = ent->igeom_instance()->createTag("TFIMapping", 1, iBase_INTEGER, taghandle);
		IBERRCHK(g_err, "Trouble create the taghandle for the surface.");
	}
	
	//it has already been checked that there are 4 edges bounding a surface in the execute function
	std::vector<iBase_EntitySetHandle> EdgeMeshSets(4);
	std::vector<int>  num(edges.size());
	//get the edge mesh entity sets for 4 edges
	for (unsigned int i = 0; i < edges.size(); i++){
		iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(edges[i], 0, EdgeMeshSets[i]);
		IBERRCHK(r_err, "Trouble get@ the mesh entity set of geometrical edge.");
		//get the number of line segments on each individual geometric edge
		iMesh::Error m_err = mk_core()->imesh_instance()->getNumOfType(EdgeMeshSets[i], iBase_EDGE, num[i]);
		IBERRCHK(m_err, "Trouble get the line segments in the edge mesh.");		
		//set the new int value for the geometric edges
		g_err = ent->igeom_instance()->setIntData(edges[i], taghandle, i);
		IBERRCHK(g_err, "Trouble create the taghandle for the surface.");
		//set the int data for geometric nodes
		g_err = ent->igeom_instance()->setIntData(gNode[i], taghandle, i);
		IBERRCHK(g_err, "Trouble create the taghandle for the surface.");
	}

	//find the corresponding edge and node relationship: based on the internal constraints, there are two edges with equal number of line segments
	int num_i=0, num_j=0;
	std::map<int, int> edgePair, nodePair, adjEdgesOfNode, oppNodeOfNode;
	num_i = num[0];


	std::vector<iBase_EntityHandle> nNodes, nEdges;//geometrical edges and nodes
	int index_a, index_b;//index_a---index for edges[0]
	//node_index_a denotes one node in edge 0, node_index_b denotes the corresponding node on the opposite edge
	int node_index_a, node_index_b;	
	//pick up the geometric edge 0, find its two end nodes	
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
	
	assert(nEdges.size()>=2);

	//find the index for one of the adjacent edges, which connects to node node_index_a
	int tmpIndex = 0; 
	for (unsigned int i = 0; i < nEdges.size(); i++){//loop over adjacent edges wrt node_index_a
		g_err = ent->igeom_instance()->getIntData(nEdges[i], taghandle, index_b);
		//finding an adjacent edge with respect to node_index_a, which definitely belongs to 4 bounding geometric edges
		if (!g_err){//maybe there are some geometric edges which don't belong to the four bounding geometric edges
			tmpIndex = i;
			break;
		}
		else{
			if (i ==(nEdges.size() - 1))
				std::cout << "test fail1\n";
		}	
	}
	if (index_a != index_b){//already found the adjacent edge with respect to node_index_a
		//try to choose nEdge[0] to find the opposite node with respect to nNodes[0]
		adjEdgesOfNode[0] = index_b;
		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;
		tnodes.clear();	
		g_err = ent->igeom_instance()->getEntAdj(nEdges[tmpIndex], iBase_VERTEX, tnodes);
		IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
		assert(tnodes.size()==2);
				
		g_err = ent->igeom_instance()->getIntData(tnodes[0], taghandle, node_index_b);
		IBERRCHK(g_err, "Trouble get the int data of node.");
		if (node_index_b != node_index_a)
			oppNodeOfNode[node_index_a] = node_index_b;			
		else{				
			g_err = ent->igeom_instance()->getIntData(tnodes[1], taghandle, node_index_b);
			IBERRCHK(g_err, "Trouble get the int data of node.");
			oppNodeOfNode[node_index_a] = node_index_b;
		}	
	}
	else
	{//choose nEdges[1]
		for (unsigned int i = tmpIndex + 1; i < nEdges.size(); i++){
			g_err = ent->igeom_instance()->getIntData(nEdges[i], taghandle, index_b);
			//IBERRCHK(g_err, "Trouble get the int data on edge 0.");
			if ((!g_err)&&(index_a != index_b)){
				tmpIndex = i;
				break;
			}
			else{
				if (i ==(nEdges.size() - 1))
					std::cout << "test fail2\n";
			}
		}		

		adjEdgesOfNode[0] = index_b;
		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;		
		g_err = ent->igeom_instance()->getEntAdj(nEdges[tmpIndex], iBase_VERTEX, tnodes);
		IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
		
		g_err = ent->igeom_instance()->getIntData(tnodes[0], taghandle, node_index_b);
		IBERRCHK(g_err, "Trouble get the int data of node.");
		if (node_index_b != node_index_a)
			//oppNodeOfNode.insert( pair<int, int>(node_index_a, node_index_b));
			oppNodeOfNode[node_index_a] = node_index_b;
		else{				
			g_err = ent->igeom_instance()->getIntData(tnodes[1], taghandle, node_index_b);
			IBERRCHK(g_err, "Trouble get the int data of node.");
			//oppNodeOfNode.insert( pair<int, int>(node_index_a, node_index_b));
			oppNodeOfNode[node_index_a] = node_index_b;
		}	
	}
	//do the other side, nodes[1] on the geometrical edges edges[0]
	int node_index_c, node_index_d;

	//pick up nodes[1] find its adjacent edges
	nEdges.clear();	
	g_err = ent->igeom_instance()->getEntAdj(nNodes[1], iBase_EDGE, nEdges);
	IBERRCHK(g_err, "Trouble get the adjacent geometric edges on a node.");
	
	assert(nEdges.size()>=2);

	g_err = ent->igeom_instance()->getIntData(nNodes[1], taghandle, node_index_c);
	IBERRCHK(g_err, "Trouble get the int data fro node 0.");	

	//find the index for one of the adjacent edges 
	for (unsigned int i = 0; i < nEdges.size(); i++){
		g_err = ent->igeom_instance()->getIntData(nEdges[i], taghandle, index_b);
		//IBERRCHK(g_err, "Trouble get the int data on edge 0.");	
		if (!g_err){
			tmpIndex = i;
			break;
		}
		else{
			if (i ==(nEdges.size() - 1))
				std::cout << "test fail3\n";
		}	
	}	
	if (index_a != index_b){
		//choose nEdge[0] to find the opposite node to nNodes[0]
		adjEdgesOfNode[1] = index_b;
		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;		
		g_err = ent->igeom_instance()->getEntAdj(nEdges[tmpIndex], iBase_VERTEX, tnodes);
		IBERRCHK(g_err, "Trouble get the adjacent geometric nodes on a edge.");
		
		g_err = ent->igeom_instance()->getIntData(tnodes[0], taghandle, node_index_d);
		IBERRCHK(g_err, "Trouble get the int data of node.");
		if (node_index_d != node_index_c)
			oppNodeOfNode[node_index_c] = node_index_d;
		else{				
			g_err = ent->igeom_instance()->getIntData(tnodes[1], taghandle, node_index_d);
			IBERRCHK(g_err, "Trouble get the int data of node.");
			oppNodeOfNode[node_index_c] = node_index_d;
		}	
	}
	else{//choose nEdges[1]
		for (unsigned int i = tmpIndex + 1; i < nEdges.size(); i++){
			g_err = ent->igeom_instance()->getIntData(nEdges[i], taghandle, index_b);
			//IBERRCHK(g_err, "Trouble get the int data on edge 0.");
			if ((!g_err)&&(index_a != index_b)){
				tmpIndex = i;
				break;
			}
			else
			{
				if (i ==(nEdges.size() - 1))
					std::cout << "test fail4\n";
			}
		}
		adjEdgesOfNode[1] = index_b;

		//extract the opposite node w.r.t node 0
		std::vector<iBase_EntityHandle> tnodes;		
		g_err = ent->igeom_instance()->getEntAdj(nEdges[tmpIndex], iBase_VERTEX, tnodes);
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

	//calculate the interior nodes based on TFI
	std::vector<iBase_EntityHandle> InteriorNodes((NodeList_j.size())*(NodeList_i.size()));
	for (unsigned int j = 1; j < (List_j.size()-1); j++){
		//                             Node 2 (List_ii)
		//		  		 |		 
		//	Node 0 (List_j)---------Node------------Node 1 (List_jj)	
		//		  	 	 |		 
		//		               Node 3 (List_i)
		double coords_i[3], coords_ii[3], coords_j[3], coords_jj[3], coords[3];
		
		g_err = mk_core()->imesh_instance()->getVtxCoord(List_j[j], coords_j[0], coords_j[1], coords_j[2]);
		IBERRCHK(g_err, "Trouble get the xyz coordinates for node 0.");

		g_err = mk_core()->imesh_instance()->getVtxCoord(List_jj[j], coords_jj[0], coords_jj[1], coords_jj[2]);
		IBERRCHK(g_err, "Trouble get the xyz coordinates for node 1.");
		
		for (unsigned int i = 1; i < (List_i.size()-1); i++){
			g_err = mk_core()->imesh_instance()->getVtxCoord(List_i[i], coords_i[0], coords_i[1], coords_i[2]);
			IBERRCHK(g_err, "Trouble get the xyz coordinates for node 3.");

			g_err = mk_core()->imesh_instance()->getVtxCoord(List_ii[i], coords_ii[0], coords_ii[1], coords_ii[2]);
			IBERRCHK(g_err, "Trouble get the xyz coordinates for node 2.");
			
			//calculate the parameter based on the distance
			double r, s, pts[3];
			ParameterCalculate(r, s, coords_j, coords_jj, coords_i, coords_ii, pts);

			parametricTFI3D(&pts[0], s, r, coords_j, coords_jj, coords_i, coords_ii);

			g_err = mk_core()->igeom_instance()->getEntClosestPtTrimmed(ent->geom_handle(), pts[0], pts[1], pts[2], coords[0], coords[1], coords[2]);
			if (g_err){
			    g_err = mk_core()->igeom_instance()->getEntClosestPt(ent->geom_handle(), pts[0], pts[1], pts[2], coords[0], coords[1], coords[2]);
			}
			IBERRCHK(g_err, "Trouble get the closest xyz coordinates on the linking surface.");	

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
			/*			
			qNodes[0] = List_i[0];
			qNodes[1] = List_j[1];
			qNodes[2] = InteriorNodes[0];
			qNodes[3] = List_i[1];
			*/
			qNodes[0] = List_i[0];
			qNodes[1] = List_i[1];
			qNodes[2] = InteriorNodes[0];
			qNodes[3] = List_j[1];			
			
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[0]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
				
			for (unsigned int i = 1; i < (NodeList_i.size()); i++)
			{
				/*
				qNodes[0] =  List_i[i];
				qNodes[1] =  InteriorNodes[i-1];
				qNodes[2] =  InteriorNodes[i];
				qNodes[3] =  List_i[i+1];
				*/
				qNodes[0] =  List_i[i];
				qNodes[1] =  List_i[i+1];
				qNodes[2] =  InteriorNodes[i];
				qNodes[3] =  InteriorNodes[i-1];				
				
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+i]);
				IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			}
			/*
			qNodes[0] = List_i[List_i.size()-2];
			qNodes[1] = InteriorNodes[List_i.size()-3];
			qNodes[2] = List_jj[1];
			qNodes[3] = List_jj[0];
			*/
			qNodes[0] = List_i[List_i.size()-2];
			qNodes[1] = List_jj[0];
			qNodes[2] = List_jj[1];
			qNodes[3] = InteriorNodes[List_i.size()-3];
			
			

			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[NodeList_i.size()]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
		}
		
		else if (j == NodeList_j.size())
		{//ending row boundary
			/*
			qNodes[0] = List_j[j];
			qNodes[1] = List_j[j+1];
			qNodes[2] = List_ii[1];
			qNodes[3] = InteriorNodes[(j-1)*(NodeList_i.size())];
			*/
			qNodes[0] = List_j[j];
			qNodes[1] = InteriorNodes[(j-1)*(NodeList_i.size())];
			qNodes[2] = List_ii[1];
			qNodes[3] = List_j[j+1];			

			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
				
			for (unsigned int i = 1; i < (NodeList_i.size()); i++)
			{
				/*
				qNodes[0] = InteriorNodes[(j-1)*(NodeList_i.size())+ i - 1];
				qNodes[1] = List_ii[i];
				qNodes[2] = List_ii[i+1];
				qNodes[3] = InteriorNodes[(j-1)*(NodeList_i.size())+ i];
				*/
				qNodes[0] = InteriorNodes[(j-1)*(NodeList_i.size())+ i - 1];
				qNodes[1] = InteriorNodes[(j-1)*(NodeList_i.size())+ i];
				qNodes[2] = List_ii[i+1];
				qNodes[3] = List_ii[i];
	
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+i]);
				IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			}
			/*
			qNodes[0] = InteriorNodes[(j-1)*(NodeList_i.size())+NodeList_i.size() - 1];
			qNodes[1] = List_ii[List_ii.size()-2];
			qNodes[2] = List_jj[List_jj.size()-1];
			qNodes[3] = List_jj[List_jj.size()-2];
			*/
			qNodes[0] = InteriorNodes[(j-1)*(NodeList_i.size())+NodeList_i.size() - 1];
			qNodes[1] = List_jj[List_jj.size()-2];
			qNodes[2] = List_jj[List_jj.size()-1];
			qNodes[3] = List_ii[List_ii.size()-2];

			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[(NodeList_j.size()+1)*(NodeList_i.size()+1) - 1]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			
		}
		else
		{
			//first column
			/*			
			qNodes[0] = List_j[j];
			qNodes[1] = List_j[j+1];
			qNodes[2] = InteriorNodes[j*NodeList_i.size()];
			qNodes[3] = InteriorNodes[(j-1)*NodeList_i.size()];
			*/
			qNodes[0] = List_j[j];
			qNodes[1] = InteriorNodes[(j-1)*NodeList_i.size()];
			qNodes[2] = InteriorNodes[j*NodeList_i.size()];
			qNodes[3] = List_j[j+1];
			
			
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			for (unsigned int i = 1; i < (NodeList_i.size()); i++)
			{
				/*				
				qNodes[0] = InteriorNodes[(j-1)*NodeList_i.size() + i-1];
				qNodes[1] = InteriorNodes[j*NodeList_i.size() + i-1];
				qNodes[2] = InteriorNodes[j*NodeList_i.size() + i];
				qNodes[3] = InteriorNodes[(j-1)*NodeList_i.size() + i];
				*/
				qNodes[0] = InteriorNodes[(j-1)*NodeList_i.size() + i-1];
				qNodes[1] = InteriorNodes[(j-1)*NodeList_i.size() + i];
				qNodes[2] = InteriorNodes[j*NodeList_i.size() + i];
				qNodes[3] = InteriorNodes[j*NodeList_i.size() + i-1];
				
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+i]);
				IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
			}

			//end column
			/*
			qNodes[0] = InteriorNodes[(j-1)*NodeList_i.size() + NodeList_i.size() - 1];
			qNodes[1] = InteriorNodes[j*NodeList_i.size() + NodeList_i.size() - 1];
			qNodes[2] = List_jj[j+1];
			qNodes[3] = List_jj[j];
			*/
			qNodes[0] = InteriorNodes[(j-1)*NodeList_i.size() + NodeList_i.size() - 1];
			qNodes[1] = List_jj[j];
			qNodes[2] = List_jj[j+1];
			qNodes[3] = InteriorNodes[j*NodeList_i.size() + NodeList_i.size() - 1];
			
			m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(NodeList_i.size()+1)+NodeList_i.size()]);
			IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
		}

	}

	//finish creating the quads
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&Quads[0], Quads.size(), entityset);
	IBERRCHK(m_err, "Trouble add an array of quads to the mesh entity set.");

	//remove the tag 
	//g_err = mk_core()->igeom_instance()->destroyTag(taghandle, 1);
	//IBERRCHK(g_err, "Trouble remove the taghandle from mesh core.");
	
	g_err = ent->igeom_instance()->rmvArrTag(&edges[0], 4, taghandle);
	IBERRCHK(g_err, "Trouble remove the tag values from an array of entities.");

	//get the four corners on a surface
	g_err = ent->igeom_instance()->rmvArrTag(&gNode[0], 4, taghandle);
	IBERRCHK(g_err, "Trouble remove the tag values from an array of entities.");
	//rmvArrTag

	//Get the global id tag
	const char *tag = "GLOBAL_ID";
	iBase_TagHandle mesh_id_tag;
	m_err = mk_core()->imesh_instance()->getTagHandle(tag, mesh_id_tag);
	IBERRCHK(m_err, "Trouble get the mesh_id_tag for 'GLOBAL_ID'.");

	std::vector<iBase_EntityHandle> m_Nodes, m_Edges, m_Quads;

	//set the int data for Global ID tag
	iBase_EntitySetHandle root_set;
	int err;
	iMesh_getRootSet(mk_core()->imesh_instance()->instance(), &root_set, &err);
	assert(!err);
	m_err = mk_core()->imesh_instance()->getEntities(root_set, iBase_VERTEX, iMesh_POINT, m_Nodes);
	IBERRCHK(m_err, "Trouble get the node list from the mesh entity set.");	
	m_err = mk_core()->imesh_instance()->getEntities(root_set, iBase_EDGE, iMesh_LINE_SEGMENT, m_Edges);
	IBERRCHK(m_err, "Trouble get the edges from the mesh entity set.");	
	m_err = mk_core()->imesh_instance()->getEntities(root_set, iBase_FACE, iMesh_QUADRILATERAL, m_Quads);
	IBERRCHK(m_err, "Trouble get the faces from the mesh entity set.");
	
	for (unsigned int i = 0; i < m_Nodes.size(); i++){
		m_err = mk_core()->imesh_instance()->setIntData(m_Nodes[i], mesh_id_tag, i);
		IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");		
	}
	for (unsigned int i = 0; i < m_Edges.size(); i++){
		m_err = mk_core()->imesh_instance()->setIntData(m_Edges[i], mesh_id_tag, i);
		IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");		
	}
	for (unsigned int i = 0; i < m_Quads.size(); i++){
		m_err = mk_core()->imesh_instance()->setIntData(m_Quads[i], mesh_id_tag, i);
		IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");		
	}

	//SurfImprove(ent->geom_handle(), entityset, iBase_FACE);

#ifdef HAVE_MESQUITE
	MeshImprove meshopt(mk_core(), true, false, false, false);
    	meshopt.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);

#endif

	return 1;
}

//****************************************************************************//
// function   : ParameterCalculate
// Date       : Feb 15, 2011
// Description: calculate the parameters for TFI mapping
//***************************************************************************//
int TFIMapping::ParameterCalculate(double &r, double &s, double pt_0s[3], double pt_1s[3], double pt_r0[3], double pt_r1[3], double *pts)
{
    // equations P_0s + s*(P_1s - P_0s) = P_r0 + t*(P_r1 - P_r0)

    assert(((fabs(pt_0s[0]-pt_1s[0])>1.0e-5)||(fabs(pt_0s[1]-pt_1s[1]) > 1.0e-5)||(fabs(pt_0s[2]-pt_1s[2]) > 1.0e-5)));
    assert(((fabs(pt_r0[0]-pt_r1[0])>1.0e-5)||(fabs(pt_r0[1]-pt_r1[1]) > 1.0e-5)||(fabs(pt_r0[2]-pt_1s[2]) > 1.0e-5)));
    s = 0; r = 0;
    double pt_s[3], pt_r[3];
    if (!LineLineIntersect(pt_0s, pt_1s, pt_r0, pt_r1, &pt_s[0], &pt_r[0], s, r)){      
	throw Error(1, "2 3D lines don't intersect at a point.");

    }
    for (int i = 0; i < 3; i++)
	pts[i] = (pt_s[i] + pt_r[i])/2;

    return 1;

}



//****************************************************************************//
// function   : SurfMeshImprove
// Date       : Oct 20, 2011
// Desription :
//   Calculate the line segment PaPb that is the shortest route between
//   two lines P1P2 and P3P4. Calculate also the values of mua and mub where
//      Pa = P1 + mua (P2 - P1)
//      Pb = P3 + mub (P4 - P3)
//   Return FALSE if no solution exists.
//****************************************************************************//
bool TFIMapping::LineLineIntersect(double p1[3], double p2[3], double p3[3], double p4[3], double *pa, double *pb, double &mua, double &mub)
{
   double p13[3], p43[3], p21[3];
   double d1343, d4321, d1321, d4343, d2121;
   double numer,denom;

   for (int i = 0; i < 3; i++){
       p13[i] = p1[i] - p3[i];
       p43[i] = p4[i] - p3[i];
   }
   if ( (fabs(p43[0]) < EPS) && (fabs(p43[1]) < EPS) && (fabs(p43[2]) < EPS))
       return false;

   for (int i = 0; i < 3; i++)
       p21[i] = p2[i] - p1[i];
   if ((fabs(p21[0]) < EPS) && (fabs(p21[1]) < EPS) && (fabs(p21[2]) < EPS))
      return false;

   d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
   d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
   d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
   d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
   d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

   denom = d2121 * d4343 - d4321 * d4321;
   if (fabs(denom) < EPS)
      return false;
   numer = d1343 * d4321 - d1321 * d4343;

   mua = numer / denom;
   mub = (d1343 + d4321 * (mua)) / d4343;

   for (int i = 0; i < 3; i++){
       pa[i] = p1[i] + mua*p21[i];
       pb[i] = p3[i] + mub*p43[i];
   }

   return true;
}





//****************************************************************************//
// function   : SurfMeshImprove
// Date       : Oct 20, 2011
// Description: smooth the surface mesh on the linking surface by using Mesquite
// Because the linking surface is convex in general, the Laplace is used to smooth
// the surface mesh on the linking surface
//***************************************************************************//
void TFIMapping::SurfImprove(iBase_EntityHandle surface, iBase_EntitySetHandle surfMesh, iBase_EntityType entity_type)
{
#ifdef HAVE_MESQUITE

#endif

}


//****************************************************************************//
// function   : parametricTFI2D
// Date       : Feb 15, 2011
// Description: do the transfinite interpolation in (pt_0s, pt_1s), (pt_r0, pt_r1)
//***************************************************************************//
void TFIMapping::parametricTFI3D(double *pts, double r, double s, double pt_0s[3], double pt_1s[3], double pt_r0[3], double pt_r1[3])
{
			//                             pt_r1
			//		  		 |		 
			//		pt_0s---------Node------------pt_1s	
			//		  	 	 |		 
			//		               pt_r0
	//assert(r>= 0 && r <= 1.0);
	//assert(s>= 0 && s <= 1.0);
	
	for (int i = 0; i < 3; i++)//interpolate the pt_rs based on pt_r0, pt_r1, pt_0s and pt_1s
	    pts[i] = 0.5*( linear_interpolation( s, pt_r0[i], pt_r1[i]) + linear_interpolation(r, pt_0s[i], pt_1s[i]) );
}


//****************************************************************************//
// function   : linear_interpolation 
// Date       : Feb 15, 2011
// Description: interpolate linearly between x0 and x1
//***************************************************************************//
double TFIMapping::linear_interpolation(double r, double x0, double x1)
{
	//assert(r >=0 && r <= 1.0);
	double pt= (1-r)*x0 + r*x1;
	return pt;
}


}

