#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/SimpleArray.hpp"
#include "meshkit/MeshImprove.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>


#define Sign(u, v) ( (v)>=0.0 ? Abs(u) : -Abs(u) )
#define Max(u, v) ( (u)>(v)? (u) : (v) )
#define Abs(u) ((u)>0 ? (u) : (-u))
#define Min(u, v) ( (u)>(v)? (v) : (u) )


namespace MeshKit
{

//---------------------------------------------------------------------------//
//Entity Type initilization for OneToOneSwept meshing
moab::EntityType OneToOneSwept_tps[] = {moab::MBVERTEX, moab::MBQUAD, moab::MBHEX, moab::MBMAXTYPE};
const moab::EntityType* OneToOneSwept::output_types()
  { return OneToOneSwept_tps; }

//---------------------------------------------------------------------------//
// construction function for OneToOneSwept class
OneToOneSwept::OneToOneSwept(MKCore *mk_core, const MEntVector &me_vec) : MeshScheme(mk_core, me_vec)
{
	//buildAssociation();	

}

//---------------------------------------------------------------------------//
// PreprocessMesh function: find the corner node list, inner node list and edge node list for the mesh on the source surface
void OneToOneSwept::PreprocessMesh(ModelEnt *me)
{
	int entity_index=0;
	std::vector<iBase_EntityHandle> Nodes, Edges, Faces;
	iBase_EntitySetHandle SourceSets;

	const char *tag = "GLOBAL_ID";
	iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle(tag, mesh_id_tag);
	IBERRCHK(m_err, "Trouble get the mesh_id_tag for 'GLOBAL_ID'.");
	
	//get the vertex list on the source surface
	iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(sourceSurface, 0, SourceSets);
	IBERRCHK(r_err, "Trouble get the mesh entity set from the geometry entity handle.");	

	//get inner nodes not boundary nodes
	m_err = mk_core()->imesh_instance()->getEntities(SourceSets, iBase_VERTEX, iMesh_POINT, Nodes);
	IBERRCHK(m_err, "Trouble get the node list from the mesh entity set.");	
	
	iBase_TagHandle taghandle=0;
	m_err = mk_core()->imesh_instance()->createTag("source", 1, iBase_INTEGER, taghandle);
 	IBERRCHK(m_err, "Trouble create the tag handle in the mesh instance.");	
	
	iBase_TagHandle taghandle_tar=0;
	m_err = mk_core()->imesh_instance()->createTag("TargetMesh", 1, iBase_INTEGER, taghandle_tar);
	IBERRCHK(m_err, "Trouble create the taghandle for the target mesh.");

	
	
	int testnum;
	m_err = mk_core()->imesh_instance()->getNumOfType(SourceSets, iBase_FACE, testnum);
	IBERRCHK(m_err, "Trouble get the number of mesh faces in the mesh entity set.");	
	
	NodeList.resize(Nodes.size());

	for (unsigned int i=0; i < Nodes.size(); i++)
	{
		entity_index++;
		NodeList[entity_index-1].gVertexHandle = Nodes[i];
		NodeList[entity_index-1].index = entity_index - 1;
		
		m_err = mk_core()->imesh_instance()->getIntData(Nodes[i], mesh_id_tag, NodeList[entity_index-1].id);		
		IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the source surface.");

				
		m_err = mk_core()->imesh_instance()->getVtxCoord(Nodes[i], NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]);
		IBERRCHK(m_err, "Trouble get the mesh node coordinates on the source surface.");
		
		Point3D pts3={NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]};
		Point2D pts2;	
		getUVCoords(sourceSurface, pts3, pts2);
		NodeList[entity_index-1].uvCoords[0] = pts2.pu;
		NodeList[entity_index-1].uvCoords[1] = pts2.pv;

		NodeList[entity_index-1].onBoundary = false;
		NodeList[entity_index-1].onCorner = false;
		
		//set the int data to the entity
		m_err = mk_core()->imesh_instance()->setIntData(Nodes[i], taghandle, NodeList[entity_index-1].index);
		IBERRCHK(m_err, "Trouble set the int value for nodes in the mesh instance.");
	}

	//loop over the edges and find the boundary nodes	
	for (unsigned int i=0; i < gsEdgeList.size(); i++)
	{
		iBase_EntitySetHandle tmpSet;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[i].gEdgeHandle, 0, tmpSet);
		IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometry edge entity handle .");

		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(tmpSet, iBase_VERTEX, iMesh_POINT, Nodes);
		IBERRCHK(m_err, "Trouble get the nodes' list from the mesh entity set.");
		
		NodeList.resize(NodeList.size() + Nodes.size());
		
		for (unsigned int j=0; j < Nodes.size(); j++)
		{
			entity_index++;
			NodeList[entity_index-1].gVertexHandle = Nodes[j];
			NodeList[entity_index-1].index = entity_index-1;

			m_err = mk_core()->imesh_instance()->getIntData(Nodes[j], mesh_id_tag, NodeList[entity_index-1].id);		
			IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the source edge.");

			m_err = mk_core()->imesh_instance()->getVtxCoord(Nodes[j], NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]);
			IBERRCHK(m_err, "Trouble get the node coordinates from mesh node entity handle.");
		
			Point3D pts3={NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]};
			Point2D pts2;	
			getUVCoords(sourceSurface, pts3, pts2);
			NodeList[entity_index-1].uvCoords[0] = pts2.pu;
			NodeList[entity_index-1].uvCoords[1] = pts2.pv;

			NodeList[entity_index-1].onBoundary = true;
			NodeList[entity_index-1].onCorner = false;
			
			m_err = mk_core()->imesh_instance()->setIntData(Nodes[j], taghandle, NodeList[entity_index-1].index);
			IBERRCHK(m_err, "Trouble set the int value for mesh node entity handle.");		
		}	
	}	
	
	//loop over the corners and find the corner nodes
	for (unsigned int i=0; i < gsVertexList.size(); i++)
	{
		iBase_EntitySetHandle tmpSet;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gsVertexList[i].gVertexHandle, 0, tmpSet);
		IBERRCHK(r_err, "Trouble get the entity set from the geometry edge handle.");		
		
		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(tmpSet, iBase_VERTEX, iMesh_POINT, Nodes);
		IBERRCHK(m_err, "Trouble get the nodes' list from mesh entity set.");
		
		NodeList.resize(NodeList.size()+Nodes.size());
		
		for (unsigned int j=0; j < Nodes.size(); j++)
		{
			entity_index++;
			NodeList[entity_index-1].gVertexHandle = Nodes[j];
			NodeList[entity_index-1].index = entity_index-1;
			
			m_err = mk_core()->imesh_instance()->getIntData(Nodes[j], mesh_id_tag, NodeList[entity_index-1].id);		
			IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the source edge.");
			
			m_err = mk_core()->imesh_instance()->getVtxCoord(Nodes[j], NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]);
			IBERRCHK(m_err, "Trouble get the nodes' coordinates from node handle.");
		
			Point3D pts3={NodeList[entity_index-1].xyzCoords[0], NodeList[entity_index-1].xyzCoords[1], NodeList[entity_index-1].xyzCoords[2]};
			Point2D pts2;	
			getUVCoords(sourceSurface, pts3, pts2);
			NodeList[entity_index-1].uvCoords[0] = pts2.pu;
			NodeList[entity_index-1].uvCoords[1] = pts2.pv;

			NodeList[entity_index-1].onBoundary = false;
			NodeList[entity_index-1].onCorner = true;
			
			m_err = mk_core()->imesh_instance()->setIntData(Nodes[j], taghandle, NodeList[entity_index-1].index);
			IBERRCHK(m_err, "Trouble set the int value for nodes.");
		}		
	}

	//sort out the boundary nodes
	//nBoundaries.resize(1);
			


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//get the edge list for the source surf mesh, this is not the boundary mesh edge
	Edges.clear();
	m_err = mk_core()->imesh_instance()->getEntities(SourceSets, iBase_EDGE, iMesh_LINE_SEGMENT, Edges);
	IBERRCHK(m_err, "Trouble get the mesh edges from the mesh entity set.");
	
	EdgeList.resize(Edges.size());
	std::cout << "8888888888888888888888888888888888888888888888888888\n";
	for (unsigned int i=0; i < Edges.size(); i++)
	{
		EdgeList[i].gEdgeHandle = Edges[i];
		EdgeList[i].index = i;
		m_err = mk_core()->imesh_instance()->getIntData(Edges[i], mesh_id_tag, EdgeList[i].EdgeID);
		IBERRCHK(m_err, "Trouble get the int data from the mesh edges.");

		m_err = mk_core()->imesh_instance()->setIntData(Edges[i], taghandle, i);
		IBERRCHK(m_err, "Trouble set the int data from the mesh edges.");

		//test  888888888888888888888888888888888888888888888888888
		//int test;
		//m_err = mk_core()->imesh_instance()->getIntData(Edges[i], taghandle, test);
		//IBERRCHK(m_err, "Trouble get the int data from the mesh edges");
		//std::cout << "-----------\nindex = " << i << "\ttest index = " << test << std::endl;

		//test  888888888888888888888888888888888888888888888888888

		//get the nodes for the edge[i], use the function iMesh_isEntContained
		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntAdj(Edges[i], iBase_VERTEX, Nodes);
		IBERRCHK(m_err, "Trouble get the adjacent nodes from the mesh edges.");
		
		//loop over the nodes on the edge elements		
		for (unsigned int j=0; j < Nodes.size(); j++)
		{
			int tmpIndex=-1;
			m_err = mk_core()->imesh_instance()->getIntData(Nodes[j], taghandle, tmpIndex);
			IBERRCHK(m_err, "Trouble get the int data from the mesh nodes.");

			//find the corresponding nodes on the vertex list
			EdgeList[i].connect[j] = &NodeList[tmpIndex];
		}

		//determine whether the edge is a boundary edge element or an inner edge element
		/*		
		if (isEdgeBoundary(Edges[i]))
			EdgeList[i].onBoundary = true;
		else
			EdgeList[i].onBoundary = false;
		*/
		EdgeList[i].onBoundary = false;
	}
	std::cout << "edge size = " << EdgeList.size() << std::endl;
	
	//get the boundary mesh edges	
	//first get the geometrical boundary edges
	entity_index = EdgeList.size() - 1;
	std::vector<iBase_EntityHandle> gEdges;
	iGeom::Error g_err = mk_core()->igeom_instance()->getEntAdj(sourceSurface, iBase_EDGE, gEdges);
	IBERRCHK(g_err, "Trouble get the mesh face entities from mesh entity set.");

	//loop over the geometrical edges
	for (unsigned int i = 0; i < gEdges.size(); i++)
	{
		//get the mesh entity set with respect to the geometrical edge
		iBase_EntitySetHandle mEdgeSet;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gEdges[i], 0, mEdgeSet);
		IBERRCHK(r_err, "Trouble get the entity set from the geometry edge handle.");

		//get the mesh edge entities
		Edges.clear();
		m_err = mk_core()->imesh_instance()->getEntities(mEdgeSet, iBase_EDGE, iMesh_LINE_SEGMENT, Edges);
		IBERRCHK(m_err, "Trouble get the mesh edges from the mesh entity set.");
		
		EdgeList.resize(EdgeList.size() + Edges.size());

		//loop over the mesh edge entities and add them to the list 
		for (unsigned int j = 0; j < Edges.size(); j++)
		{
			entity_index++;
			
			EdgeList[entity_index].gEdgeHandle = Edges[j];
			EdgeList[entity_index].index = entity_index;
			m_err = mk_core()->imesh_instance()->getIntData(Edges[j], mesh_id_tag, EdgeList[entity_index].EdgeID);
			IBERRCHK(m_err, "Trouble get the int data from the mesh edges.");

			m_err = mk_core()->imesh_instance()->setIntData(Edges[j], taghandle, entity_index);
			IBERRCHK(m_err, "Trouble set the int data from the mesh edges.");

			//test  888888888888888888888888888888888888888888888888888
			//int test;
			//m_err = mk_core()->imesh_instance()->getIntData(Edges[i], taghandle, test);
			//IBERRCHK(m_err, "Trouble get the int data from the mesh edges");
			//std::cout << "-----------\nindex = " << i << "\ttest index = " << test << std::endl;

			//test  888888888888888888888888888888888888888888888888888

			//get the nodes for the edge[i], use the function iMesh_isEntContained
			Nodes.clear();
			m_err = mk_core()->imesh_instance()->getEntAdj(Edges[j], iBase_VERTEX, Nodes);
			IBERRCHK(m_err, "Trouble get the adjacent nodes from the mesh edges.");
			assert(Nodes.size()==2);		

			//loop over the nodes on the edge elements		
			for (unsigned int k=0; k < Nodes.size(); k++)
			{
				int tmpIndex=-1;
				m_err = mk_core()->imesh_instance()->getIntData(Nodes[k], taghandle, tmpIndex);
				IBERRCHK(m_err, "Trouble get the int data from the mesh nodes.");

				//find the corresponding nodes on the vertex list
				EdgeList[entity_index].connect[k] = &NodeList[tmpIndex];
			}

			//determine whether the edge is a boundary edge element or an inner edge element
			EdgeList[entity_index].onBoundary = true;			
			//if (isEdgeBoundary(Edges[j]))
			//	EdgeList[j].onBoundary = true;
			//else
			//	EdgeList[j].onBoundary = false;
		}

	}

	//std::cout << "8888888888888888888888888888888888888888888888888888\n";
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Faces.clear();
	m_err = mk_core()->imesh_instance()->getEntities(SourceSets, iBase_FACE, iMesh_ALL_TOPOLOGIES, Faces);
	IBERRCHK(m_err, "Trouble get the mesh face entities from mesh entity set.");
	
	FaceList.resize(Faces.size());
	for (unsigned int i=0; i < Faces.size(); i++)
	{
		FaceList[i].gFaceHandle = Faces[i];
		FaceList[i].index = i;
		m_err = mk_core()->imesh_instance()->getIntData(Faces[i], mesh_id_tag, FaceList[i].FaceID);
		IBERRCHK(m_err, "Trouble get the int data from mesh face entity.");

		//get the nodes on the face elements
		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntAdj(Faces[i], iBase_VERTEX, Nodes);
		IBERRCHK(m_err, "Trouble get the adjacent nodes from mesh face entity.");

		FaceList[i].connect.resize(Nodes.size());
		for (unsigned int j=0; j < Nodes.size(); j++)
		{			
			int tmpIndex=-1;
			m_err = mk_core()->imesh_instance()->getIntData(Nodes[j], taghandle, tmpIndex);
			IBERRCHK(m_err, "Trouble get the int data from node handle.");
			
			FaceList[i].connect[j] = &NodeList[tmpIndex];
		}

		m_err = mk_core()->imesh_instance()->setIntData(Faces[i], taghandle, i);
		IBERRCHK(m_err, "Trouble set the int data for quad mesh on the source surface.");
	}

	std::cout << "Face size = " << FaceList.size() << std::endl;

	//initialize the mesh size on the target surface
	TVertexList.resize(NodeList.size());
	TEdgeList.resize(EdgeList.size());
	TFaceList.resize(FaceList.size());

	std::cout << "--------------------test create mesh corners---------------------------\n";

	//Initialize the vertex meshing on the target surface, don't forget to add the mesh entities to the set
	for (unsigned int i = 0; i < gsVertexList.size(); i++)
	{
		//get the mesh entity node from the source surface
		std::cout << "geometrical node " << i << "\tx=" << gsVertexList[i].xyzCoords[0] << "\ty=" << gsVertexList[i].xyzCoords[1] << "\tz=" << gsVertexList[i].xyzCoords[2] << std::endl;
		
		iBase_EntitySetHandle tmpSet;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gsVertexList[i].gVertexHandle, 0, tmpSet);
		IBERRCHK(r_err, "Trouble get the mesh entity set from the geometry entity handle.");
		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(tmpSet, iBase_VERTEX, iMesh_POINT, Nodes);
		IBERRCHK(r_err, "Trouble get the mesh entity node from the geometry entity vertex handle.");

		std::cout << "---------------------------------------\n";
		for (unsigned int j = 0; j < Nodes.size(); j++){
			double xyzCoords[3];
			m_err = mk_core()->imesh_instance()->getVtxCoord(Nodes[j], xyzCoords[0], xyzCoords[1], xyzCoords[2]);
			IBERRCHK(m_err, "Trouble get the mesh node coordinates.");

			std::cout << "Node " << j << "\tx = " << xyzCoords[0] << "\ty = " << xyzCoords[1] << "\tz = " << xyzCoords[2] << std::endl; 
		}
		std::cout << "---------------------------------------\n";



		std::cout << "node size = " << Nodes.size() << std::endl;



		assert(Nodes.size()==1);
		m_err = mk_core()->imesh_instance()->getIntData(Nodes[0], taghandle, entity_index);
		IBERRCHK(r_err, "Trouble get the int data for mesh node on the source surface.");

		
		//get the mesh entity node from the target surface
		iBase_EntitySetHandle tmpSet_tar;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gtVertexList[cornerPairs[i]].gVertexHandle, 0, tmpSet_tar);
		IBERRCHK(r_err, "Trouble get the mesh entity set from the geometry entity handle.");
		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(tmpSet_tar, iBase_VERTEX, iMesh_POINT, Nodes);
		IBERRCHK(m_err, "Trouble get the mesh entity node from the geometry entity vertex handle.");
		//assert(Nodes.size()==1);
		if (Nodes.size()==0)
		{
			Nodes.resize(1);
			m_err = mk_core()->imesh_instance()->createVtx(gtVertexList[cornerPairs[i]].xyzCoords[0], gtVertexList[cornerPairs[i]].xyzCoords[1], gtVertexList[cornerPairs[i]].xyzCoords[2], Nodes[0]);
			IBERRCHK(m_err, "Trouble create the mesh node for the geometry entity vertex.");
			m_err = mk_core()->imesh_instance()->addEntToSet(Nodes[0], tmpSet_tar);
			IBERRCHK(m_err, "Trouble add the mesh node entity to the set.");
		}

		assert(Nodes.size()==1);
		TVertexList[entity_index].gVertexHandle = Nodes[0];
		TVertexList[entity_index].index = entity_index;
		
		m_err = mk_core()->imesh_instance()->setIntData(Nodes[0], taghandle_tar, entity_index);		
		IBERRCHK(m_err, "Trouble set the int data for mesh nodes on the target surface.");

				
		m_err = mk_core()->imesh_instance()->getVtxCoord(Nodes[0], TVertexList[entity_index].xyzCoords[0], TVertexList[entity_index].xyzCoords[1], TVertexList[entity_index].xyzCoords[2]);
		IBERRCHK(m_err, "Trouble get the mesh node coordinates on the target surface.");
		
		Point3D pts3={TVertexList[entity_index].xyzCoords[0], TVertexList[entity_index].xyzCoords[1], TVertexList[entity_index].xyzCoords[2]};
		Point2D pts2;	
		getUVCoords(targetSurface, pts3, pts2);
		TVertexList[entity_index].uvCoords[0] = pts2.pu;
		TVertexList[entity_index].uvCoords[1] = pts2.pv;

		TVertexList[entity_index].onBoundary = false;
		TVertexList[entity_index].onCorner = true;

		std::cout << "mesh corner " << i << "\tx=" << TVertexList[entity_index].xyzCoords[0] << "\ty=" << TVertexList[entity_index].xyzCoords[1] << "\tz=" << TVertexList[entity_index].xyzCoords[2] << std::endl;		
	}
	


}

//---------------------------------------------------------------------------//
//determine whether a mesh edge is on the boundary or not
int OneToOneSwept::isEdgeBoundary(iBase_EntityHandle gEdgeHandle)
{
	std::vector<iBase_EntityHandle> Faces;	
	Faces.clear();
	
	iMesh::Error m_err = mk_core()->imesh_instance()->getEntAdj(gEdgeHandle, iBase_FACE, Faces);
	IBERRCHK(m_err, "Trouble get the adjacent faces from edge handle.");

	if (Faces.size()==1)
		return 1;
	else
		return 0;
}

//---------------------------------------------------------------------------//
// deconstruction function for OneToOneSwept class
OneToOneSwept::~OneToOneSwept()
{

}

//---------------------------------------------------------------------------//
// setup function: define the size between the different layers
void OneToOneSwept::setup_this()
{
    //compute the number of intervals for the associated ModelEnts, from the size set on them
    //the sizing function they point to, or a default sizing function
  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit->first;

      //first check to see whether entity is meshed
    if (me->get_meshed_state() >= COMPLETE_MESH || me->mesh_intervals() > 0)
      continue;
    
    SizingFunction *sf = mk_core()->sizing_function(me->sizing_function_index());
    if (!sf && me -> mesh_intervals() < 0 && me -> interval_firmness() == DEFAULT &&
        mk_core()->sizing_function(0))
      sf = mk_core()->sizing_function(0);
    
    if (!sf && me -> mesh_intervals() < 0 && me -> interval_firmness() == DEFAULT)
    {
        //no sizing set, just assume default #intervals as 4
      me->mesh_intervals(4);
      me->interval_firmness(DEFAULT);
    }
    else
    {
        //check # intervals first, then size, and just choose for now
      if (sf->intervals() > 0)
      {
        if (me->constrain_even() && sf->intervals()%2)
          me -> mesh_intervals(sf->intervals()+1);
        else
          me -> mesh_intervals(sf->intervals());
        me -> interval_firmness(HARD);
      }
      else if (sf->size()>0)
      {
        int intervals = me->measure()/sf->size();
        if (!intervals) intervals++;
        if (me->constrain_even() && intervals%2) intervals++;
        me->mesh_intervals(intervals);
        me->interval_firmness(SOFT);
      }
      else
        throw Error(MK_INCOMPLETE_MESH_SPECIFICATION,  "Sizing function for edge had neither positive size nor positive intervals.");
    }
  }


}

//---------------------------------------------------------------------------//
// PreprocessGeom function: preprocess the geometry and prepare for sweeping
// e.g. specify the source surfaces, target surfaces, linking surfaces, mapping 
//      the relation between the source surface and target surface
void OneToOneSwept::PreprocessGeom(ModelEnt *me)
{
	//get geometry root set	
	geom_root_set = me->igeom_instance()->getRootSet();

	//get the geom_id_tag
	const char *tag = "GLOBAL_ID";
	iGeom::Error g_err = mk_core()->igeom_instance()->getTagHandle(tag, geom_id_tag);
	IBERRCHK(g_err, "Trouble get the geom_id_tag for 'GLOBAL_ID'.");
	
	std::vector<iBase_EntityHandle> gFaces;
	g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_FACE, gFaces);
	IBERRCHK(g_err, "Trouble get the geometrical faces from geom_root_set.");

	std::vector<iBase_EntityHandle> gVertices;	
	for (unsigned int i = 0; i < gFaces.size(); i++)
	{
		gVertices.clear();
		g_err = mk_core()->igeom_instance()->getEntAdj(gFaces[i], iBase_VERTEX, gVertices);
		IBERRCHK(g_err, "Trouble get the geometrical vertices with respect to a face.");
	
		std::cout << "-------------------------------------------------------------------------" << std::endl;
		std::cout << "Face index = " << i << std::endl;		
		double coords[3];
		for (unsigned int j = 0; j < gVertices.size(); j++)
		{
			g_err = mk_core()->igeom_instance()->getVtxCoord(gVertices[j], coords[0], coords[1], coords[2]);
			IBERRCHK(g_err, "Trouble get the coordinates of vertices.");
			
			std::cout << "vertex " << j << "\tx = " << coords[0] << "\ty = " << coords[1] << "\tz = " << coords[2] << std::endl;
		}
	}

	//select the source surface and target surface	
	sourceSurface = gFaces[index_src];
	targetSurface = gFaces[index_tar];

	//create a id tag handle for source surface and target surface
	iBase_TagHandle  src_id_tag, tar_id_tag, link_id_tag;
	g_err = mk_core()->igeom_instance()->createTag("SourceIdTag", 1, iBase_INTEGER, src_id_tag);
	IBERRCHK(g_err, "Trouble create id tag handle for source surface.");

	g_err = mk_core()->igeom_instance()->createTag("TargetIdTag", 1, iBase_INTEGER, tar_id_tag);
	IBERRCHK(g_err, "Trouble create id tag handle for target surface.");

	g_err = mk_core()->igeom_instance()->createTag("LinkIdTag", 1, iBase_INTEGER, link_id_tag);
	IBERRCHK(g_err, "Trouble create id tag handle for target surface.");	

	//initialize the vertices on the source surface
	gVertices.resize(0);
	g_err = mk_core()->igeom_instance()->getEntAdj(sourceSurface, iBase_VERTEX, gVertices);
	IBERRCHK(g_err, "Trouble get the adjacent vertices around the source surface.");
	gsVertexList.resize(gVertices.size());	
	for (unsigned int i = 0; i < gVertices.size(); i++)
	{
		//initialize the vertex on the source surface based on the vertex data structure defined on OneToOneSwept.hpp file		
		gsVertexList[i].index = i;
		gsVertexList[i].gVertexHandle = gVertices[i];
		g_err = mk_core()->igeom_instance()->getVtxCoord(gVertices[i], gsVertexList[i].xyzCoords[0], gsVertexList[i].xyzCoords[1], gsVertexList[i].xyzCoords[2]);
		IBERRCHK(g_err, "Trouble get the vertex coordinates on the source surface.");
		g_err = mk_core()->igeom_instance()->getIntData(gVertices[i], geom_id_tag, gsVertexList[i].id);
		IBERRCHK(g_err, "Trouble get the int data for vertices on the source surface.");
		g_err = mk_core()->igeom_instance()->getEntXYZtoUVHint(sourceSurface, gsVertexList[i].xyzCoords[0], gsVertexList[i].xyzCoords[1], gsVertexList[i].xyzCoords[2], gsVertexList[i].uvCoords[0], gsVertexList[i].uvCoords[1]);
		IBERRCHK(g_err, "Trouble get the parametric coordinates from xyz coordinates for vertices on the source surface.");
		//set vertex index for int data on the source surface
		g_err = mk_core()->igeom_instance()->setIntData(gVertices[i], src_id_tag, i);
		IBERRCHK(g_err, "Trouble set the int data for vertices on the source surface.");
				
	}

	//initialize the vertices on the target surface
	gVertices.clear();
	g_err = mk_core()->igeom_instance()->getEntAdj(targetSurface, iBase_VERTEX, gVertices);
	IBERRCHK(g_err, "Trouble get the adjacent vertices around the target surface.");
	gtVertexList.resize(gVertices.size());
	for (unsigned int i = 0; i < gVertices.size(); i++)
	{
		gtVertexList[i].index = i;
		gtVertexList[i].gVertexHandle = gVertices[i];
		g_err = mk_core()->igeom_instance()->getVtxCoord(gVertices[i], gtVertexList[i].xyzCoords[0], gtVertexList[i].xyzCoords[1], gtVertexList[i].xyzCoords[2]);
		IBERRCHK(g_err, "Trouble get the vertex coordinates on the target surface.");
		g_err = mk_core()->igeom_instance()->getIntData(gVertices[i], geom_id_tag, gtVertexList[i].id);
		IBERRCHK(g_err, "Trouble get the int data for vertices on the target surface.");
		g_err = mk_core()->igeom_instance()->getEntXYZtoUVHint(targetSurface, gtVertexList[i].xyzCoords[0], gtVertexList[i].xyzCoords[1], gtVertexList[i].xyzCoords[2], gtVertexList[i].uvCoords[0], gtVertexList[i].uvCoords[1]);
		IBERRCHK(g_err, "Trouble get the parametric coordinates from xyz coordinates for vertices on the target surface.");
		//set vertex index for int data on the source surface
		g_err = mk_core()->igeom_instance()->setIntData(gVertices[i], tar_id_tag, i);
		IBERRCHK(g_err, "Trouble set the int data for vertices on the target surface.");
	}
	

	//initialize the edges on the source surface
	std::vector<iBase_EntityHandle> gEdges;
	g_err = mk_core()->igeom_instance()->getEntAdj(sourceSurface, iBase_EDGE, gEdges);
	IBERRCHK(g_err, "Trouble get the adjacent edges around the source surface.");
	gsEdgeList.resize(gEdges.size());
	//create a set for storing the geometrical edges
	std::set<iBase_EntityHandle> src_edges;
	for (unsigned int i = 0; i < gEdges.size(); i++)
	{
		gsEdgeList[i].index = i;
		gsEdgeList[i].gEdgeHandle = gEdges[i];
		src_edges.insert(gEdges[i]);
		g_err = mk_core()->igeom_instance()->getIntData(gEdges[i], geom_id_tag, gsEdgeList[i].EdgeID);
		IBERRCHK(g_err, "Trouble get the int data for edges on the source surface.");
		gVertices.clear();
		g_err = mk_core()->igeom_instance()->getEntAdj(gEdges[i], iBase_VERTEX, gVertices);
		IBERRCHK(g_err, "Trouble get the adjacent vertices of edge on the source surface.");
		assert(gVertices.size()==2);
		for (unsigned int j = 0; j < gVertices.size(); j++)
		{
			int index_id;
			g_err = mk_core()->igeom_instance()->getIntData(gVertices[j], src_id_tag, index_id);
			IBERRCHK(g_err, "Trouble get the index id for vertices on the source surface.");
			gsEdgeList[i].connect[j] = &gsVertexList[index_id];		
		}
		//set the int data for edges on the source surface
		g_err = mk_core()->igeom_instance()->setIntData(gEdges[i], src_id_tag, i);
		IBERRCHK(g_err, "Trouble set the int data for edges on the source surface.");				
		
	}

	//initialize the edges on the target surface
	gEdges.clear();	
	g_err = mk_core()->igeom_instance()->getEntAdj(targetSurface, iBase_EDGE, gEdges);
	IBERRCHK(g_err, "Trouble get the adjacent edges around the target surface.");
	gtEdgeList.resize(gEdges.size());
	for (unsigned int i = 0; i < gEdges.size(); i++)
	{
		gtEdgeList[i].index = i;
		gtEdgeList[i].gEdgeHandle = gEdges[i];
		g_err = mk_core()->igeom_instance()->getIntData(gEdges[i], geom_id_tag, gtEdgeList[i].EdgeID);
		IBERRCHK(g_err, "Trouble get the int data for edges on the target surface.");
		gVertices.clear();
		g_err = mk_core()->igeom_instance()->getEntAdj(gEdges[i], iBase_VERTEX, gVertices);
		IBERRCHK(g_err, "Trouble get the adjacent vertices of edge on the target surface.");
		assert(gVertices.size()==2);
		for (unsigned int j = 0; j < gVertices.size(); j++)
		{
			int index_id;
			g_err = mk_core()->igeom_instance()->getIntData(gVertices[j], tar_id_tag, index_id);
			IBERRCHK(g_err, "Trouble get the index id for vertices on the target surface.");
			gtEdgeList[i].connect[j] = &gtVertexList[index_id];
		}
		//set the int data for edges on the source surface
		g_err = mk_core()->igeom_instance()->setIntData(gEdges[i], tar_id_tag, i);
		IBERRCHK(g_err, "Trouble set the int data for edges on the target surface.");		
	}

	std::cout << "-----------------------------\nthe vertices on the source surface\n";
	for (unsigned int i = 0; i < gsVertexList.size(); i++)
	    std::cout << "source nodes[" << i << "]\t= {" << gsVertexList[i].xyzCoords[0] << ", " << gsVertexList[i].xyzCoords[1] << ", " << gsVertexList[i].xyzCoords[2] << "}" << std::endl;
	std::cout << "-----------------------------\n";
	std::cout << "-----------------------------\nthe vertices on the target surface\n";
	for (unsigned int i = 0; i < gtVertexList.size(); i++)
	    std::cout << "source nodes[" << i << "]\t= {" << gtVertexList[i].xyzCoords[0] << ", " << gtVertexList[i].xyzCoords[1] << ", " << gtVertexList[i].xyzCoords[2] << "}" << std::endl;
	std::cout << "-----------------------------\n";


	//initialize the edges on the linking faces
	gEdges.clear();
	//get all the edges in the whole model
	g_err = mk_core()->igeom_instance()->getEntities(geom_root_set, iBase_EDGE, gEdges);
	IBERRCHK(g_err, "Trouble get the adjacent vertices of edge on the target surface.");
	int index = 0;	
	for (unsigned int i = 0; i < gEdges.size(); i++)
	{
		int EID;
		g_err = mk_core()->igeom_instance()->getIntData(gEdges[i], geom_id_tag, EID);
		IBERRCHK(g_err, "Trouble get the int data for edges on the target surface.");
		bool isInSrc = false, isInTar = false;
		//check whether an edge is in edge list on the source surface
		for (unsigned int j = 0; j < gsEdgeList.size(); j++)
		{
			if (EID == gsEdgeList[j].EdgeID)
			{
				isInSrc = true;
				break;
			}
		}
		//check whether an edge is in edge list on the target surface
		if (!isInSrc)
		{
			for (unsigned int j = 0; j < gtEdgeList.size(); j++)
			{
				if (EID == gtEdgeList[j].EdgeID)
				{
					isInTar = true;
					break;
				}
			}
		}
		if (!isInSrc && !isInTar)
		{
			index++;
			gLinkSides.resize(index);
			gLinkSides[index - 1].index = index - 1;
			gLinkSides[index - 1].EdgeID = EID;
			gLinkSides[index - 1].gEdgeHandle = gEdges[i];
			gVertices.clear();
			g_err = mk_core()->igeom_instance()->getEntAdj(gEdges[i], iBase_VERTEX, gVertices);
			IBERRCHK(g_err, "Trouble get the adjacent vertices of linking edge.");
			assert(gVertices.size()==2);
			//there should be one vertex on the source surface and the other vertex on the target surface
			//gVertices[0] = vertex on the source surface, gVertices[1] = vertex on the target surface

//////////////////////////////////////////////////
				std::cout << "--------------------\nedge = " << index << std::endl;
				double test[3];				
				g_err = mk_core()->igeom_instance()->getVtxCoord(gVertices[1], test[0], test[1], test[2]);
				IBERRCHK(g_err, "Trouble get coordinates.");
				std::cout << "test gVertices[1] x = " << test[0] << "\ty = " << test[1] << "\tz = " << test[2] << std::endl;
				
				g_err = mk_core()->igeom_instance()->getVtxCoord(gVertices[0], test[0], test[1], test[2]);
				IBERRCHK(g_err, "Trouble get coordinates.");
				std::cout << "test gVertices[0] x = " << test[0] << "\ty = " << test[1] << "\tz = " << test[2] << std::endl;

				
//////////////////////////////////////////////////


			int index_id, index_id_all;
			g_err = mk_core()->igeom_instance()->getIntData(gVertices[0], geom_id_tag, index_id_all);
			IBERRCHK(g_err, "Trouble get the int data for vertices on the linking sides.");
			g_err = mk_core()->igeom_instance()->getIntData(gVertices[0], src_id_tag, index_id);
			//IBERRCHK(g_err, "Trouble get the int data for vertices on the linking sides.");
			if ((!g_err)&&(gsVertexList[index_id].id == index_id_all))
			//if (gsVertexList[index_id].id == index_id_all)
			{//gVertices[0] = vertex on the source surface, gVertices[1] = vertex on the target surface
				gLinkSides[index - 1].connect[0] = &gsVertexList[index_id];
				g_err = mk_core()->igeom_instance()->getIntData(gVertices[1], tar_id_tag, index_id);
				IBERRCHK(g_err, "Trouble get the int data for vertices on the linking sides.");
				gLinkSides[index - 1].connect[1] = &gtVertexList[index_id];
			}
			else
			{//gVertices[1] = vertex on the source surface, gVertices[0] = vertex on the target surface				
				g_err = mk_core()->igeom_instance()->getIntData(gVertices[1], src_id_tag, index_id);				
				IBERRCHK(g_err, "Trouble get the int data for vertices on the linking sides.");
				gLinkSides[index - 1].connect[0] = &gsVertexList[index_id];

				g_err = mk_core()->igeom_instance()->getIntData(gVertices[0], tar_id_tag, index_id);				
				IBERRCHK(g_err, "Trouble get the int data for vertices on the linking sides.");
				gLinkSides[index - 1].connect[1] = &gtVertexList[index_id];
			}

			//create the mapping relationship for vertices between the source surface and target surface
			cornerPairs[gLinkSides[index - 1].connect[0]->index] = gLinkSides[index - 1].connect[1]->index;

			//set the int data for linking sides
			g_err = mk_core()->igeom_instance()->setIntData(gEdges[i], link_id_tag, index - 1);				
			IBERRCHK(g_err, "Trouble set the int data for vertices on the linking sides.");

			
		}
		
	}	

	//initialize the linking faces
	index = 0;
	for (unsigned int i = 0; i < gFaces.size(); i++)
	{
		if (((int)i != index_src) && ( (int)i != index_tar ))
		{
			index++;
			gLinkFaceList.resize(index);
			gLinkFaceList[index - 1].gFaceHandle = gFaces[i];
			gLinkFaceList[index - 1].index = index - 1;
			g_err = mk_core()->igeom_instance()->getIntData(gFaces[i], geom_id_tag, gLinkFaceList[index - 1].FaceID);
			IBERRCHK(g_err, "Trouble get the int data for linking surfaces.");
			gEdges.clear();
			g_err = mk_core()->igeom_instance()->getEntAdj(gFaces[i], iBase_EDGE, gEdges);
			IBERRCHK(g_err, "Trouble get the adjacent vertices of linking edge.");
			
			//data structure for linking surface
			//	connect[2]----------connEdges[3]----------connect[3]
			//	   |					       |
			//     connEdges[1]				 connEdges[2]
			//         |                                           |
			//	connect[0]----------connEdges[0]----------connect[1]

			//detecting the edges: which one is on the source, which one is on the target and which two are on the linking surface
			assert(gEdges.size()==4);
			gLinkFaceList[index - 1].connEdges.resize(gEdges.size());
			int index_id, index_id_all;
			int link_edge_index1 = -1, link_edge_index2 = -1;
			//detect the edge which is the intersection between source surface and linking surface			
			for (unsigned int j = 0; j < gEdges.size(); j++)
			{
				g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], geom_id_tag, index_id_all);
				IBERRCHK(g_err, "Trouble get the int data for vertices on the source surface.");
				g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], src_id_tag, index_id);
				//IBERRCHK(g_err, "Trouble get the int data for vertices on the source surface.");
				//if (gsEdgeList[index_id].EdgeID == index_id_all)
				if ((!g_err)&&(gsEdgeList[index_id].EdgeID == index_id_all))				
				{// the edge is the intersection between the source surface and linking surface
					gLinkFaceList[index - 1].connEdges[0] = &gsEdgeList[index_id];
				}
				else
				{
					//g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], geom_id_tag, index_id_all);
					//IBERRCHK(g_err, "Trouble get the int data for vertices on the target surface.");
					g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], tar_id_tag, index_id);
					//IBERRCHK(g_err, "Trouble get the int data for vertices on the target surface.");
					if ((!g_err)&&(gtEdgeList[index_id].EdgeID == index_id_all))
					{// the edge is the intersection between the target surface and linking surface
						gLinkFaceList[index - 1].connEdges[3] = &gtEdgeList[index_id];
					}
					else
					{//the edge is on the linking surface
						if ((link_edge_index1 == -1)&&(link_edge_index2 == -1))
						{
							g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], link_id_tag, index_id);
							IBERRCHK(g_err, "Trouble get the int data for edges on the linking surface.");
							link_edge_index1 = index_id;
						}
						else if ((link_edge_index1 > -1)&&(link_edge_index2 == -1))
						{
							g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], link_id_tag, index_id);
							IBERRCHK(g_err, "Trouble get the int data for edges on the linking surface.");
							link_edge_index2 = index_id;
						}					
					}
				}
				
			}

			//assigned 2 vertices on the edge of source surface to the list
			gLinkFaceList[index - 1].connect.resize(4);
			gVertices.clear();
			g_err = mk_core()->igeom_instance()->getEntAdj(gLinkFaceList[index - 1].connEdges[0]->gEdgeHandle, iBase_VERTEX, gVertices);
			IBERRCHK(g_err, "Trouble get the adjacent vertices of linking edge.");

			assert(gVertices.size()==2);
			g_err = mk_core()->igeom_instance()->getIntData(gVertices[0], src_id_tag, index_id);
			IBERRCHK(g_err, "Trouble get the int data for vertices on the source surface.");
			gLinkFaceList[index - 1].connect[0] = &gsVertexList[index_id];
			//assign one vertex on the target surface to the list
			gLinkFaceList[index - 1].connect[2] = &gtVertexList[cornerPairs[index_id]];
			
			g_err = mk_core()->igeom_instance()->getIntData(gVertices[1], src_id_tag, index_id);
			IBERRCHK(g_err, "Trouble get the int data for vertices on the source surface.");

			gLinkFaceList[index - 1].connect[1] = &gsVertexList[index_id];
			//assign the other vertex on the target surface to the list
			gLinkFaceList[index - 1].connect[3] = &gtVertexList[cornerPairs[index_id]];
	
			//ReOrder the linking edges
			if (gLinkSides[link_edge_index1].connect[0]->id == gLinkFaceList[index - 1].connect[0]->id)
			{
				gLinkFaceList[index - 1].connEdges[1] = &gLinkSides[link_edge_index1];
				gLinkFaceList[index - 1].connEdges[2] = &gLinkSides[link_edge_index2];
			}
			else
			{
				gLinkFaceList[index - 1].connEdges[1] = &gLinkSides[link_edge_index2];
				gLinkFaceList[index - 1].connEdges[2] = &gLinkSides[link_edge_index1];
			}

			//store the edge relationship 
			edgePairs[gLinkFaceList[index - 1].connEdges[0]->index] = gLinkFaceList[index - 1].connEdges[3]->index;		
		}
	}

	//detect the geometrical boundary loops
	//gBoundaries.resize(1);
	//gBoundaries[0].insert(gsEdgeList[0].index);
	unsigned int boundary_loop = 0;
	std::set<int> sumBoundaryEdges;
	std::cout << "total number of edges on the source surface is " << gsEdgeList.size() << std::endl;

	while(sumBoundaryEdges.size() < gsEdgeList.size())
	{
		std::cout << "------------Boundary Loop " << boundary_loop << "------------" << std::endl;		


		boundary_loop++;
		gBoundaries.resize(boundary_loop);
		//get the starting geometrical edge
		if (boundary_loop == 1)
		{
			int first_index = 0;
			//DetectFirstEdge(src_edges, first_index);
			gBoundaries[boundary_loop - 1].push_back(gsEdgeList[first_index].index);
			sumBoundaryEdges.insert(gsEdgeList[first_index].index);

			std::cout << "+++++++++++++++++++++test first edge on the outer boundary loop++++++++++++++++++++++++++++++++++\n";
			std::cout << "Node 1 x=" << gsEdgeList[first_index].connect[0]->xyzCoords[0] << "\ty=" << gsEdgeList[first_index].connect[0]->xyzCoords[1] << "\tz=" << gsEdgeList[first_index].connect[0]->xyzCoords[2] << std::endl;
			std::cout << "Node 2 x=" << gsEdgeList[first_index].connect[1]->xyzCoords[0] << "\ty=" << gsEdgeList[first_index].connect[1]->xyzCoords[1] << "\tz=" << gsEdgeList[first_index].connect[1]->xyzCoords[2] << std::endl;
			std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		}
		else
		{
			for (unsigned int k = 0; k < gsEdgeList.size(); k++)
			{
				if (sumBoundaryEdges.find(gsEdgeList[k].index) == sumBoundaryEdges.end())
				{
					gBoundaries[boundary_loop - 1].push_back(gsEdgeList[k].index);
					sumBoundaryEdges.insert(gsEdgeList[k].index);
					break;
				}
			}
		}

		//get the info from the first edge on the boundary loop
		int start_index = gsEdgeList[*(gBoundaries[boundary_loop - 1].begin())].connect[0]->index;
		int end_index = gsEdgeList[*(gBoundaries[boundary_loop - 1].begin())].connect[1]->index;
		int first_edge_index = *(gBoundaries[boundary_loop - 1].begin());
		int second_edge_index = *(gBoundaries[boundary_loop - 1].begin());

		//std::cout << "node x=" << gsVertexList[start_index].xyzCoords[0] << "\ty=" << gsVertexList[start_index].xyzCoords[1] << "\tz=" << gsVertexList[start_index].xyzCoords[2] << std::endl;
		//std::cout << "node x=" << gsVertexList[end_index].xyzCoords[0] << "\ty=" << gsVertexList[end_index].xyzCoords[1] << "\tz=" << gsVertexList[end_index].xyzCoords[2] << std::endl; 
		while(end_index != (gsEdgeList[*(gBoundaries[boundary_loop - 1].begin())].connect[0]->index))
		{
			//get the adjacent edges around a geometrical vertex
			std::vector<iBase_EntityHandle> g_edges;
			g_err = mk_core()->igeom_instance()->getEntAdj(gsVertexList[end_index].gVertexHandle, iBase_EDGE, g_edges);
			IBERRCHK(g_err, "Trouble get the geometrical edges around a vertex.");

			//check whether all the edges are on the source surface. If there is any edge which doesn't belong to source surface,
			//delete it
			std::set<iBase_EntityHandle>::iterator it;
			for (unsigned int i = 0; i < g_edges.size();)
			{
				it = src_edges.find(g_edges[i]);
				if (it == src_edges.end())
				{
					g_edges.erase(g_edges.begin()+i);
					i = 0;
					continue;
				}
				else
					i++;		
			}				

			assert(g_edges.size()==2);

			//get the next geometrical edge
			int index_id;
			g_err = mk_core()->igeom_instance()->getIntData(g_edges[0], src_id_tag, index_id);
			IBERRCHK(g_err, "Trouble get the int data for geometrical edge entity.");

			if (index_id == second_edge_index)
			{
				g_err = mk_core()->igeom_instance()->getIntData(g_edges[1], src_id_tag, index_id);
				IBERRCHK(g_err, "Trouble get the int data for geometrical edge entity.");

				//update the edge index
				first_edge_index = second_edge_index;
				second_edge_index = index_id;

				//add the new geometrical edge to the set
				gBoundaries[boundary_loop - 1].push_back(index_id);
				sumBoundaryEdges.insert(index_id);

				//update the vertex index
				start_index = end_index;

				index_id = gsEdgeList[index_id].connect[0]->index;

				if (index_id == end_index)
					end_index = gsEdgeList[index_id].connect[1]->index;
				else
					end_index = index_id;
			}
			else		
			{
				//update edge index
				first_edge_index = second_edge_index;
				second_edge_index = index_id;

				//add the new geometrical edge to the set
				gBoundaries[boundary_loop - 1].push_back(index_id);
				sumBoundaryEdges.insert(index_id);
			
				//update the vertex index
				start_index = end_index;

				index_id = gsEdgeList[index_id].connect[0]->index;

				if (index_id == end_index)
					end_index = gsEdgeList[index_id].connect[1]->index;
				else
					end_index = index_id;
			}

			//std::cout << "node x=" << gsVertexList[end_index].xyzCoords[0] << "\ty=" << gsVertexList[end_index].xyzCoords[1] << "\tz=" << gsVertexList[end_index].xyzCoords[2] << std::endl;
		}

		//update the number of edges added to the boundary loop
		//std::cout << "number of edges on the boundary loop [" << (boundary_loop - 1) << "] = " << gBoundaries[boundary_loop - 1].size() << std::endl;
	}
	DetectFirstEdge();	

}

//---------------------------------------------------------------------------//
//detect the first geometric edge on the outer boundary loop
void OneToOneSwept::DetectFirstEdge()
{
	//data structure for linking surface
	//	connect[2]----------connEdges[3]----------connect[3]
	//	   |					                     |
	//  connEdges[1]				              connEdges[2]
	//     |                                         |
	//	connect[0]----------connEdges[0]----------connect[1]
	if (int(gBoundaries.size())==1)
		return;
	std::vector< std::vector<double> > corner;
	corner.resize(gBoundaries.size());
	for (unsigned int i = 0; i < gBoundaries.size(); i++){
		corner[i].resize(6);
		for (int j = 0; j < 3; j++){
			corner[i][j] = -1.0e10;
			corner[i][j] = 1.0e10;
		}
		std::list<int>::iterator it = gBoundaries[i].begin();
		for (; it != gBoundaries[i].end(); it++){
			//get the bounding box for an individual edge
			double MinMax[6];
			iGeom::Error g_err = mk_core()->igeom_instance()->getEntBoundBox(gsEdgeList[*it].gEdgeHandle, MinMax[0], MinMax[1], MinMax[2], MinMax[3], MinMax[4], MinMax[5]);
			IBERRCHK(g_err, "Trouble get the bounding box for edge entity.");
			if (MinMax[0] < corner[i][0])
				corner[i][0] = MinMax[0];
			if (MinMax[1] < corner[i][1])
				corner[i][1] = MinMax[1];
			if (MinMax[2] < corner[i][2])
				corner[i][1] = MinMax[1];
			if (MinMax[3] < corner[i][3])
				corner[i][3] = MinMax[3];
			if (MinMax[4] < corner[i][4])
				corner[i][4] = MinMax[4];
			if (MinMax[5] < corner[i][5])
				corner[i][5] = MinMax[5];
		}
	}	
	int index = 0;
	for (unsigned int i = 1; i < gBoundaries.size(); i++){
		bool isTruexy = (((corner[i][0] < corner[index][0])&&(corner[i][1] < corner[index][1]))&&((corner[i][3] > corner[index][3])&&(corner[i][4] > corner[index][4])));
		bool isTruexz = (((corner[i][0] < corner[index][0])&&(corner[i][2] < corner[index][2]))&&((corner[i][3] > corner[index][3])&&(corner[i][5] > corner[index][5])));
		bool isTrueyz = (((corner[i][1] < corner[index][1])&&(corner[i][2] < corner[index][2]))&&((corner[i][4] > corner[index][4])&&(corner[i][5] > corner[index][5])));

		if (isTruexy || isTruexz || isTrueyz)
			index = i;
	}
	if (index == 0)
		return;
	else{
		//swap gBoundaries[0] with gBoundaries[index];
		std::list<int> tmp1, tmp2;
		for (std::list<int>::iterator it = gBoundaries[0].begin(); it != gBoundaries[0].end(); it++)
			tmp1.push_back(*it);
		for (std::list<int>::iterator it = gBoundaries[index].begin(); it != gBoundaries[index].end(); it++)
			tmp2.push_back(*it);
		gBoundaries[0].clear();
		gBoundaries[index].clear();
		for (std::list<int>::iterator it = tmp2.begin(); it != tmp2.end(); it++)
			gBoundaries[0].push_back(*it);
		for (std::list<int>::iterator it = tmp1.begin(); it != tmp1.end(); it++)
			gBoundaries[index].push_back(*it);
		

	}


}

//---------------------------------------------------------------------------//
// execute function: generate the all-hex mesh through sweeping from source 
// surface to target surface
void OneToOneSwept::execute_this()
{
	std::vector<double> coords;
	std::vector<moab::EntityHandle> nodes;

	//set up the adjTable
	iMesh::AdjTableType adjTable = mk_core()->imesh_instance()->getAdjTable();
	
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			std::cout << "index " << (i+1) << "\t" << (j+1) << "\t" << adjTable[i][j] << std::endl;		
		}
	}
	//adjTable[1][1] = (iBase::AdjacencyCost)1;
	
	
	//iMesh_setAdjTable(mk_core()->imesh_instance(), adjTable, 16, &err);
	
	//assert(!err);
	
	for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  	{
    		ModelEnt *me = mit -> first;
		if (me->get_meshed_state() >= COMPLETE_MESH)
			continue;

		PreprocessGeom(me);		

		me->boundary(0, nodes);
		
		PreprocessMesh(me);

    	//resize the coords based on the interval setting
    	numLayers = me->mesh_intervals();		

		//necessary steps for setting up the source surface and target surfaces    		
		TargetSurfProjection();

		mk_core()->save_mesh("test_targetsurf.vtk");
		
		//get the volume mesh entity set
		iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(me->geom_handle(), 0, volumeSet);
		IBERRCHK(r_err, "Trouble get the volume mesh entity set from the geometrical volume.");
		
		//do the linking surface meshing, create the hexs
		InnerLayerMeshing();		

      	//   ok, we are done, commit to ME
    	me->commit_mesh(mit->second, COMPLETE_MESH);	
  	}

}

//---------------------------------------------------------------------------//
// set the source surface function
void OneToOneSwept::SetSourceSurface(int index)
{
	index_src = index;
}

//---------------------------------------------------------------------------//
// set the target surface function
void OneToOneSwept::SetTargetSurface(int index)
{
	index_tar = index;
}


//****************************************************************************//
// function   : InnerLayerMeshing
// Author     : Shengyong Cai
// Date       : Feb 16, 2011
// Description: Generate all-hex mesh by general sweeping in the interior volume
//***************************************************************************//
int OneToOneSwept::InnerLayerMeshing()
{
	//first discretize the linking sides
	//temporarily discretize the linking sides based on the equally spacing
	iRel::Error r_err;
	iMesh::Error m_err;
	
	//create the vertices on the linking surface
	//determine whether there exists the mesh on the linking surface
	for (unsigned int i = 0; i < gLinkFaceList.size(); i++)
	{		
		ModelEnt link_surf(mk_core(), gLinkFaceList[i].gFaceHandle, 0, 0, 0);
		if (link_surf.get_meshed_state()>=COMPLETE_MESH)
		{
			iBase_EntitySetHandle mEdgeSet;
			r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].connEdges[1]->gEdgeHandle, 0, mEdgeSet);
			IBERRCHK(r_err, "Trouble get the mesh entity set from the geometrical linking edges.");

			int num_lines;
			m_err = mk_core()->imesh_instance()->getNumOfType(mEdgeSet, iBase_EDGE, num_lines);
			IBERRCHK(m_err, "Trouble get the number of line segments from mesh entity set.");

			if (num_lines != numLayers)
			{
				numLayers = num_lines;
			}
		}
	}
	

	//case 1: if numLayers = 1, then it is not necessary to create any vertices for linking surface, All the vertices have been created by source surface and target surface
	vector<vector <Vertex> > linkVertexList(numLayers-1, vector<Vertex>(NodeList.size()));
	
	LinkSurfMeshing(linkVertexList);	
	
	//create the inner vertex on the different layers
	InnerNodesProjection(linkVertexList);
	
	
	//create the quadrilateral face elements on the different layers and cell elements
	CreateElements(linkVertexList); 

	return 1;
}



//****************************************************************************//
// function   : CreateElements
// Author     : Shengyong Cai
// Date       : Feb 16, 2011
// Description: create hexahedral elements by connecting 8 nodes which have 
//              already been created by previous functions
//***************************************************************************//
int OneToOneSwept::CreateElements(vector<vector <Vertex> > &linkVertexList)
{
	//create the quadrilateral elements on the different layers
	//it is not necessary to create the quadrilateral elements on the different layers. Hex element can be created based on the eight nodes
	iMesh::Error m_err;
	
	vector<iBase_EntityHandle> mVolumeHandle(FaceList.size());
		
	for (int m=0; m < numLayers-1; m++)
	{
		if (m==0){
			for (unsigned int i=0; i < FaceList.size(); i++){
				vector<iBase_EntityHandle> connect(8);
				
				connect[0] = NodeList[(FaceList[i].getVertex(0))->index].gVertexHandle;
				connect[1] = NodeList[(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[2] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[3] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
				
				connect[4] = NodeList[(FaceList[i].getVertex(3))->index].gVertexHandle;
				connect[5] = NodeList[(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[6] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[7] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, mVolumeHandle[i]);
				IBERRCHK(m_err, "Trouble create the hexahedral elements.");
			}
			m_err = mk_core()->imesh_instance()->addEntArrToSet(&mVolumeHandle[0], FaceList.size(), volumeSet);
			IBERRCHK(m_err, "Trouble add an array of hexahedral elements to the entity set.");
			if (m == (numLayers-2))
			{
				for (unsigned int i=0; i < FaceList.size(); i++)
				{
					vector<iBase_EntityHandle> connect(8);
					
					connect[0] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
					connect[1] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[2] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[3] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
				
					connect[4] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
					connect[5] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[6] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[7] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
					
					m_err = mk_core()->imesh_instance()->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, mVolumeHandle[i]);
					IBERRCHK(m_err, "Trouble create the hexahedral elements.");
				}
				m_err = mk_core()->imesh_instance()->addEntArrToSet(&mVolumeHandle[0], FaceList.size(), volumeSet);
				IBERRCHK(m_err, "Trouble add an array of hexahedral elements to the entity set.");
			}
		}
		else{
			for (unsigned int i=0; i < FaceList.size(); i++){
				vector<iBase_EntityHandle> connect(8);
				
				connect[0] = linkVertexList[m-1][(FaceList[i].getVertex(0))->index].gVertexHandle;
				connect[1] = linkVertexList[m-1][(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[2] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
				connect[3] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
				
				connect[4] = linkVertexList[m-1][(FaceList[i].getVertex(3))->index].gVertexHandle;
				connect[5] = linkVertexList[m-1][(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[6] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
				connect[7] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
				m_err = mk_core()->imesh_instance()->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, mVolumeHandle[i]);
				IBERRCHK(m_err, "Trouble create the hexahedral elements.");
			}
			m_err = mk_core()->imesh_instance()->addEntArrToSet(&mVolumeHandle[0], FaceList.size(), volumeSet);
			IBERRCHK(m_err, "Trouble add an array of hexahedral elements to the entity set.");
			if (m == (numLayers-2))
			{
				for (unsigned int i=0; i < FaceList.size(); i++)
				{
					vector<iBase_EntityHandle> connect(8);
					
					connect[0] = linkVertexList[m][(FaceList[i].getVertex(0))->index].gVertexHandle;
					connect[1] = linkVertexList[m][(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[2] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
					connect[3] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
				
					connect[4] = linkVertexList[m][(FaceList[i].getVertex(3))->index].gVertexHandle;
					connect[5] = linkVertexList[m][(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[6] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
					connect[7] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
					m_err = mk_core()->imesh_instance()->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, mVolumeHandle[i]);
					IBERRCHK(m_err, "Trouble create the hexahedral elements.");
				}
				m_err = mk_core()->imesh_instance()->addEntArrToSet(&mVolumeHandle[0], FaceList.size(), volumeSet);
				IBERRCHK(m_err, "Trouble add an array of hexahedral elements to the entity set.");
			}
		}
	}
	
	return 1;
}


//****************************************************************************//
// function   : InnerNodesProjection
// Author     : Shengyong Cai
// Date       : Feb 16, 2011
// Description: generate the interior nodes between the Source surface and 
//              target surface.
//***************************************************************************//
int OneToOneSwept::InnerNodesProjection(vector<vector <Vertex> > &linkVertexList)
{
	iMesh::Error m_err;
	int numPts=0;
	Point3D sPtsCenter = {0, 0, 0}, tPtsCenter={0, 0, 0};
	std::vector<Point3D> PtsCenter(numLayers-1);
	std::vector<Point3D> sBoundaryNodes(0), tBoundaryNodes(0);
	std::vector<vector <Point3D> > iBoundaryNodes(numLayers-1, std::vector<Point3D>(0));
	
	//calculate the center coordinates
	for (int i=0; i< numLayers-1; i++)
	{
		PtsCenter[i].px = 0;
		PtsCenter[i].py = 0;
		PtsCenter[i].pz = 0;
	}
	for (unsigned int i=0; i < NodeList.size(); i++)
	{
		if (NodeList[i].onBoundary || NodeList[i].onCorner)
		{
			sPtsCenter.px = sPtsCenter.px + NodeList[i].xyzCoords[0];
			sPtsCenter.py = sPtsCenter.py + NodeList[i].xyzCoords[1];
			sPtsCenter.pz = sPtsCenter.pz + NodeList[i].xyzCoords[2];
			
			tPtsCenter.px = tPtsCenter.px + TVertexList[i].xyzCoords[0];
			tPtsCenter.py = tPtsCenter.py + TVertexList[i].xyzCoords[1];
			tPtsCenter.pz = tPtsCenter.pz + TVertexList[i].xyzCoords[2];
			
			numPts++;
			
			sBoundaryNodes.resize(numPts);
			tBoundaryNodes.resize(numPts);
			
			sBoundaryNodes[numPts-1].px = NodeList[i].xyzCoords[0];
			sBoundaryNodes[numPts-1].py = NodeList[i].xyzCoords[1];
			sBoundaryNodes[numPts-1].pz = NodeList[i].xyzCoords[2];	
			
			tBoundaryNodes[numPts-1].px = TVertexList[i].xyzCoords[0];
			tBoundaryNodes[numPts-1].py = TVertexList[i].xyzCoords[1];
			tBoundaryNodes[numPts-1].pz = TVertexList[i].xyzCoords[2];
			
			for (int j=0; j< numLayers-1; j++)
			{
				iBoundaryNodes[j].resize(numPts);
				PtsCenter[j].px = PtsCenter[j].px + linkVertexList[j][i].xyzCoords[0];
				PtsCenter[j].py = PtsCenter[j].py + linkVertexList[j][i].xyzCoords[1];
				PtsCenter[j].pz = PtsCenter[j].pz + linkVertexList[j][i].xyzCoords[2];
				
				iBoundaryNodes[j][numPts-1].px = linkVertexList[j][i].xyzCoords[0];
				iBoundaryNodes[j][numPts-1].py = linkVertexList[j][i].xyzCoords[1];
				iBoundaryNodes[j][numPts-1].pz = linkVertexList[j][i].xyzCoords[2];
			}
		}
	}
	sPtsCenter.px = sPtsCenter.px/double(numPts);
	sPtsCenter.py = sPtsCenter.py/double(numPts);
	sPtsCenter.pz = sPtsCenter.pz/double(numPts);
	
	tPtsCenter.px = tPtsCenter.px/double(numPts);
	tPtsCenter.py = tPtsCenter.py/double(numPts);
	tPtsCenter.pz = tPtsCenter.pz/double(numPts);
	
	//calculate the center coordinates for the ith layer 
	for (int i=0; i< numLayers-1; i++)
	{
		PtsCenter[i].px = PtsCenter[i].px/double(numPts);
		PtsCenter[i].py = PtsCenter[i].py/double(numPts);
		PtsCenter[i].pz = PtsCenter[i].pz/double(numPts);
	}
	
	//loop over different layers
	for (int i=0; i < numLayers-1; i++)
	{
		double sA[3][3], tA[3][3];
		double stransMatrix[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, ttransMatrix[3][3]={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
		double sInvMatrix[3][3], tInvMatrix[3][3];
		double sb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, tb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
		std::vector<Point3D> sBNodes(numPts);	//boundary nodes on the source surface
		std::vector<Point3D> tBNodes(numPts);	//boundary nodes on the target surface
		std::vector<std::vector <Point3D> > isBNodes(numLayers-1, std::vector<Point3D>(numPts)), itBNodes(numLayers-1, vector<Point3D>(numPts));  //boundary nodes on the different layer
		
		
		//transform the coordinates
		for (int j=0; j < numPts; j++)
		{
			//transform the coordinates on the source suface
			sBNodes[j].px = sBoundaryNodes[j].px;
			sBNodes[j].py = sBoundaryNodes[j].py;
			sBNodes[j].pz = sBoundaryNodes[j].pz;
			
			tBNodes[j].px = tBoundaryNodes[j].px;
			tBNodes[j].py = tBoundaryNodes[j].py;
			tBNodes[j].pz = tBoundaryNodes[j].pz;
			
			//from the source surface to layer j
			isBNodes[i][j].px = iBoundaryNodes[i][j].px;
			isBNodes[i][j].py = iBoundaryNodes[i][j].py;
			isBNodes[i][j].pz = iBoundaryNodes[i][j].pz;
			//from the target surface to layer j
			itBNodes[i][j].px = iBoundaryNodes[i][j].px;
			itBNodes[i][j].py = iBoundaryNodes[i][j].py;
			itBNodes[i][j].pz = iBoundaryNodes[i][j].pz;
			
			sBNodes[j].px = sBNodes[j].px - 2*sPtsCenter.px + PtsCenter[i].px;
			sBNodes[j].py = sBNodes[j].py - 2*sPtsCenter.py + PtsCenter[i].py;
			sBNodes[j].pz = sBNodes[j].pz - 2*sPtsCenter.pz + PtsCenter[i].pz;
			
			tBNodes[j].px = tBNodes[j].px - 2*tPtsCenter.px + PtsCenter[i].px;
			tBNodes[j].py = tBNodes[j].py - 2*tPtsCenter.py + PtsCenter[i].py;
			tBNodes[j].pz = tBNodes[j].pz - 2*tPtsCenter.pz + PtsCenter[i].pz;
					
			//transform the coordinates on the different layers
			isBNodes[i][j].px = isBNodes[i][j].px - sPtsCenter.px;
			isBNodes[i][j].py = isBNodes[i][j].py - sPtsCenter.py;
			isBNodes[i][j].pz = isBNodes[i][j].pz - sPtsCenter.pz;
			
			itBNodes[i][j].px = itBNodes[i][j].px - tPtsCenter.px;
			itBNodes[i][j].py = itBNodes[i][j].py - tPtsCenter.py;
			itBNodes[i][j].pz = itBNodes[i][j].pz - tPtsCenter.pz;
		}
		
		//calculate the temporary matrix
		for (int j=0; j < numPts; j++)
		{
			//transform matrix: from the source surface to layer j
			//first row entries in the temporary matrix
			stransMatrix[0][0] = stransMatrix[0][0] + sBNodes[j].px*sBNodes[j].px;
			stransMatrix[0][1] = stransMatrix[0][1] + sBNodes[j].px*sBNodes[j].py;
			stransMatrix[0][2] = stransMatrix[0][2] + sBNodes[j].px*sBNodes[j].pz;
			//second row entries in the temporary matrix
			stransMatrix[1][0] = stransMatrix[1][0] + sBNodes[j].py*sBNodes[j].px;
			stransMatrix[1][1] = stransMatrix[1][1] + sBNodes[j].py*sBNodes[j].py;
			stransMatrix[1][2] = stransMatrix[1][2] + sBNodes[j].py*sBNodes[j].pz;
			//third row entries in the temporary matrix
			stransMatrix[2][0] = stransMatrix[2][0] + sBNodes[j].pz*sBNodes[j].px;
			stransMatrix[2][1] = stransMatrix[2][1] + sBNodes[j].pz*sBNodes[j].py;
			stransMatrix[2][2] = stransMatrix[2][2] + sBNodes[j].pz*sBNodes[j].pz;
			//transform matrix: from the target surface to layer j
			//first row entries in the temporary matrix
			ttransMatrix[0][0] = ttransMatrix[0][0] + tBNodes[j].px*tBNodes[j].px;
			ttransMatrix[0][1] = ttransMatrix[0][1] + tBNodes[j].px*tBNodes[j].py;
			ttransMatrix[0][2] = ttransMatrix[0][2] + tBNodes[j].px*tBNodes[j].pz;
			//second row entries in the temporary matrix
			ttransMatrix[1][0] = ttransMatrix[1][0] + tBNodes[j].py*tBNodes[j].px;
			ttransMatrix[1][1] = ttransMatrix[1][1] + tBNodes[j].py*tBNodes[j].py;
			ttransMatrix[1][2] = ttransMatrix[1][2] + tBNodes[j].py*tBNodes[j].pz;
			//third row entries in the temporary matrix
			ttransMatrix[2][0] = ttransMatrix[2][0] + tBNodes[j].pz*tBNodes[j].px;
			ttransMatrix[2][1] = ttransMatrix[2][1] + tBNodes[j].pz*tBNodes[j].py;
			ttransMatrix[2][2] = ttransMatrix[2][2] + tBNodes[j].pz*tBNodes[j].pz;
			
			//transform matrix: from the source surface to layer j
			//first row entries in the temporary matrix
			sb[0][0] = sb[0][0] + isBNodes[i][j].px*sBNodes[j].px;
			sb[0][1] = sb[0][1] + isBNodes[i][j].px*sBNodes[j].py;
			sb[0][2] = sb[0][2] + isBNodes[i][j].px*sBNodes[j].pz;
			//second row entries in the temporary matrix
			sb[1][0] = sb[1][0] + isBNodes[i][j].py*sBNodes[j].px;
			sb[1][1] = sb[1][1] + isBNodes[i][j].py*sBNodes[j].py;
			sb[1][2] = sb[1][2] + isBNodes[i][j].py*sBNodes[j].pz;
			//third row entries in the temporary matrix
			sb[2][0] = sb[2][0] + isBNodes[i][j].pz*sBNodes[j].px;
			sb[2][1] = sb[2][1] + isBNodes[i][j].pz*sBNodes[j].py;
			sb[2][2] = sb[2][2] + isBNodes[i][j].pz*sBNodes[j].pz;
			//transform matrix: from the target surface to layer j
			//first row entries in the temporary matrix
			tb[0][0] = tb[0][0] + itBNodes[i][j].px*tBNodes[j].px;
			tb[0][1] = tb[0][1] + itBNodes[i][j].px*tBNodes[j].py;
			tb[0][2] = tb[0][2] + itBNodes[i][j].px*tBNodes[j].pz;
			//second row entries in the temporary matrix
			tb[1][0] = tb[1][0] + itBNodes[i][j].py*tBNodes[j].px;
			tb[1][1] = tb[1][1] + itBNodes[i][j].py*tBNodes[j].py;
			tb[1][2] = tb[1][2] + itBNodes[i][j].py*tBNodes[j].pz;
			//third row entries in the temporary matrix
			tb[2][0] = tb[2][0] + itBNodes[i][j].pz*tBNodes[j].px;
			tb[2][1] = tb[2][1] + itBNodes[i][j].pz*tBNodes[j].py;
			tb[2][2] = tb[2][2] + itBNodes[i][j].pz*tBNodes[j].pz;
			
		}		
		
		//first determine the affine mapping matrix is singular or not
		double sdetMatrix = stransMatrix[0][2]*stransMatrix[1][1]*stransMatrix[2][0] - stransMatrix[0][1]*stransMatrix[1][2]*stransMatrix[2][0] - stransMatrix[0][2]*stransMatrix[1][0]*stransMatrix[2][1] + stransMatrix[0][0]*stransMatrix[1][2]*stransMatrix[2][1] + stransMatrix[0][1]*stransMatrix[1][0]*stransMatrix[2][2] - stransMatrix[0][0]*stransMatrix[1][1]*stransMatrix[2][2];
		double tdetMatrix = ttransMatrix[0][2]*ttransMatrix[1][1]*ttransMatrix[2][0] - ttransMatrix[0][1]*ttransMatrix[1][2]*ttransMatrix[2][0] - ttransMatrix[0][2]*ttransMatrix[1][0]*ttransMatrix[2][1] + ttransMatrix[0][0]*ttransMatrix[1][2]*ttransMatrix[2][1] + ttransMatrix[0][1]*ttransMatrix[1][0]*ttransMatrix[2][2] - ttransMatrix[0][0]*ttransMatrix[1][1]*ttransMatrix[2][2];
		//transMatrix[0][0]*(transMatrix[1][1]*transMatrix[2][2]-transMatrix[2][1]*transMatrix[1][2]) - transMatrix[0][1]*(transMatrix[1][0]*transMatrix[2][2] - transMatrix[2][0]*transMatrix[1][2]) + transMatrix[0][2]*(transMatrix[1][0]*transMatrix[2][1]-transMatrix[1][1]*transMatrix[2][0]);
		assert(pow(sdetMatrix, 2)>1.0e-20);
		assert(pow(tdetMatrix, 2)>1.0e-20);
		
		////solve the affine mapping matrix, make use of inverse matrix to get affine mapping matrix
		sInvMatrix[0][0] = (stransMatrix[2][1]*stransMatrix[1][2] - stransMatrix[1][1]*stransMatrix[2][2])/sdetMatrix;
		
		sInvMatrix[0][1] = (stransMatrix[0][1]*stransMatrix[2][2] - stransMatrix[0][2]*stransMatrix[2][1])/sdetMatrix;
		
		sInvMatrix[0][2] = (stransMatrix[0][2]*stransMatrix[1][1] - stransMatrix[0][1]*stransMatrix[1][2])/sdetMatrix;
		
		sInvMatrix[1][0] = (stransMatrix[1][0]*stransMatrix[2][2] - stransMatrix[1][2]*stransMatrix[2][0])/sdetMatrix;
		
		sInvMatrix[1][1] = (stransMatrix[0][2]*stransMatrix[2][0] - stransMatrix[0][0]*stransMatrix[2][2])/sdetMatrix;
		
		sInvMatrix[1][2] = (stransMatrix[0][0]*stransMatrix[1][2] - stransMatrix[0][2]*stransMatrix[1][0])/sdetMatrix;
		
		sInvMatrix[2][0] = (stransMatrix[1][1]*stransMatrix[2][0] - stransMatrix[1][0]*stransMatrix[2][1])/sdetMatrix;
		
		sInvMatrix[2][1] = (stransMatrix[0][0]*stransMatrix[2][1] - stransMatrix[0][1]*stransMatrix[2][0])/sdetMatrix;
		
		sInvMatrix[2][2] = (stransMatrix[0][1]*stransMatrix[1][0] - stransMatrix[0][0]*stransMatrix[1][1])/sdetMatrix;
		
		//projection from the target surface
		tInvMatrix[0][0] = (ttransMatrix[2][1]*ttransMatrix[1][2] - ttransMatrix[1][1]*ttransMatrix[2][2])/tdetMatrix;
		
		tInvMatrix[0][1] = (ttransMatrix[0][1]*ttransMatrix[2][2] - ttransMatrix[0][2]*ttransMatrix[2][1])/tdetMatrix;
		
		tInvMatrix[0][2] = (ttransMatrix[0][2]*ttransMatrix[1][1] - ttransMatrix[0][1]*ttransMatrix[1][2])/tdetMatrix;
		
		tInvMatrix[1][0] = (ttransMatrix[1][0]*ttransMatrix[2][2] - ttransMatrix[1][2]*ttransMatrix[2][0])/tdetMatrix;
		
		tInvMatrix[1][1] = (ttransMatrix[0][2]*ttransMatrix[2][0] - ttransMatrix[0][0]*ttransMatrix[2][2])/tdetMatrix;
		
		tInvMatrix[1][2] = (ttransMatrix[0][0]*ttransMatrix[1][2] - ttransMatrix[0][2]*ttransMatrix[1][0])/tdetMatrix;
		
		tInvMatrix[2][0] = (ttransMatrix[1][1]*ttransMatrix[2][0] - ttransMatrix[1][0]*ttransMatrix[2][1])/tdetMatrix;
		
		tInvMatrix[2][1] = (ttransMatrix[0][0]*ttransMatrix[2][1] - ttransMatrix[0][1]*ttransMatrix[2][0])/tdetMatrix;
		
		tInvMatrix[2][2] = (ttransMatrix[0][1]*ttransMatrix[1][0] - ttransMatrix[0][0]*ttransMatrix[1][1])/tdetMatrix;
		
		
		sA[0][0] = sInvMatrix[0][0]*sb[0][0] + sInvMatrix[0][1]*sb[0][1] + sInvMatrix[0][2]*sb[0][2];
		sA[0][1] = sInvMatrix[1][0]*sb[0][0] + sInvMatrix[1][1]*sb[0][1] + sInvMatrix[1][2]*sb[0][2];
		sA[0][2] = sInvMatrix[2][0]*sb[0][0] + sInvMatrix[2][1]*sb[0][1] + sInvMatrix[2][2]*sb[0][2];
		
		sA[1][0] = sInvMatrix[0][0]*sb[1][0] + sInvMatrix[0][1]*sb[1][1] + sInvMatrix[0][2]*sb[1][2];
		sA[1][1] = sInvMatrix[1][0]*sb[1][0] + sInvMatrix[1][1]*sb[1][1] + sInvMatrix[1][2]*sb[1][2];
		sA[1][2] = sInvMatrix[2][0]*sb[1][0] + sInvMatrix[2][1]*sb[1][1] + sInvMatrix[2][2]*sb[1][2];
		
		sA[2][0] = sInvMatrix[0][0]*sb[2][0] + sInvMatrix[0][1]*sb[2][1] + sInvMatrix[0][2]*sb[2][2];
		sA[2][1] = sInvMatrix[1][0]*sb[2][0] + sInvMatrix[1][1]*sb[2][1] + sInvMatrix[1][2]*sb[2][2];
		sA[2][2] = sInvMatrix[2][0]*sb[2][0] + sInvMatrix[2][1]*sb[2][1] + sInvMatrix[2][2]*sb[2][2];
		
		tA[0][0] = tInvMatrix[0][0]*tb[0][0] + tInvMatrix[0][1]*tb[0][1] + tInvMatrix[0][2]*tb[0][2];
		tA[0][1] = tInvMatrix[1][0]*tb[0][0] + tInvMatrix[1][1]*tb[0][1] + tInvMatrix[1][2]*tb[0][2];
		tA[0][2] = tInvMatrix[2][0]*tb[0][0] + tInvMatrix[2][1]*tb[0][1] + tInvMatrix[2][2]*tb[0][2];
		
		tA[1][0] = tInvMatrix[0][0]*tb[1][0] + tInvMatrix[0][1]*tb[1][1] + tInvMatrix[0][2]*tb[1][2];
		tA[1][1] = tInvMatrix[1][0]*tb[1][0] + tInvMatrix[1][1]*tb[1][1] + tInvMatrix[1][2]*tb[1][2];
		tA[1][2] = tInvMatrix[2][0]*tb[1][0] + tInvMatrix[2][1]*tb[1][1] + tInvMatrix[2][2]*tb[1][2];
		
		tA[2][0] = tInvMatrix[0][0]*tb[2][0] + tInvMatrix[0][1]*tb[2][1] + tInvMatrix[0][2]*tb[2][2];
		tA[2][1] = tInvMatrix[1][0]*tb[2][0] + tInvMatrix[1][1]*tb[2][1] + tInvMatrix[1][2]*tb[2][2];
		tA[2][2] = tInvMatrix[2][0]*tb[2][0] + tInvMatrix[2][1]*tb[2][1] + tInvMatrix[2][2]*tb[2][2];
		
		//calculate the inner nodes for different layers
		for (unsigned int j=0; j < NodeList.size(); j++)
		{
			if (!(NodeList[j].onBoundary || NodeList[j].onCorner))
			{
				Point3D spts, tpts, pts;
				double s;
				iBase_EntityHandle nodeHandle;
				spts.px = sA[0][0]*(NodeList[j].xyzCoords[0] - 2*sPtsCenter.px + PtsCenter[i].px) + sA[0][1]*(NodeList[j].xyzCoords[1] - 2*sPtsCenter.py + PtsCenter[i].py) + sA[0][2]*(NodeList[j].xyzCoords[2] - 2*sPtsCenter.pz + PtsCenter[i].pz) + sPtsCenter.px;
				spts.py = sA[1][0]*(NodeList[j].xyzCoords[0] - 2*sPtsCenter.px + PtsCenter[i].px) + sA[1][1]*(NodeList[j].xyzCoords[1] - 2*sPtsCenter.py + PtsCenter[i].py) + sA[1][2]*(NodeList[j].xyzCoords[2] - 2*sPtsCenter.pz + PtsCenter[i].pz) + sPtsCenter.py;
				spts.pz = sA[2][0]*(NodeList[j].xyzCoords[0] - 2*sPtsCenter.px + PtsCenter[i].px) + sA[2][1]*(NodeList[j].xyzCoords[1] - 2*sPtsCenter.py + PtsCenter[i].py) + sA[2][2]*(NodeList[j].xyzCoords[2] - 2*sPtsCenter.pz + PtsCenter[i].pz) + sPtsCenter.pz;
				
				tpts.px = tA[0][0]*(TVertexList[j].xyzCoords[0] - 2*tPtsCenter.px + PtsCenter[i].px) + tA[0][1]*(TVertexList[j].xyzCoords[1] - 2*tPtsCenter.py + PtsCenter[i].py) + tA[0][2]*(TVertexList[j].xyzCoords[2] - 2*tPtsCenter.pz + PtsCenter[i].pz) + tPtsCenter.px;
				tpts.py = tA[1][0]*(TVertexList[j].xyzCoords[0] - 2*tPtsCenter.px + PtsCenter[i].px) + tA[1][1]*(TVertexList[j].xyzCoords[1] - 2*tPtsCenter.py + PtsCenter[i].py) + tA[1][2]*(TVertexList[j].xyzCoords[2] - 2*tPtsCenter.pz + PtsCenter[i].pz) + tPtsCenter.py;
				tpts.pz = tA[2][0]*(TVertexList[j].xyzCoords[0] - 2*tPtsCenter.px + PtsCenter[i].px) + tA[2][1]*(TVertexList[j].xyzCoords[1] - 2*tPtsCenter.py + PtsCenter[i].py) + tA[2][2]*(TVertexList[j].xyzCoords[2] - 2*tPtsCenter.pz + PtsCenter[i].pz) + tPtsCenter.pz;
				
				s = (i+1)/double(numLayers);
				pts.px = linear_interpolation(s, spts.px, tpts.px);
				pts.py = linear_interpolation(s, spts.py, tpts.py);
				pts.pz = linear_interpolation(s, spts.pz, tpts.pz);
				
				linkVertexList[i][j].xyzCoords[0] = pts.px;
				linkVertexList[i][j].xyzCoords[1] = pts.py;
				linkVertexList[i][j].xyzCoords[2] = pts.pz;
				
				m_err = mk_core()->imesh_instance()->createVtx(pts.px, pts.py, pts.pz, nodeHandle);
				IBERRCHK(m_err, "Trouble create the vertex entity.");
				linkVertexList[i][j].gVertexHandle = nodeHandle;
				m_err = mk_core()->imesh_instance()->addEntToSet(nodeHandle, volumeSet);
				IBERRCHK(m_err, "Trouble add the node handle to the entity set.");
			}	
		}								
	}
	
	return 1;
}


//****************************************************************************//
// function   : LinkSurfMeshing
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: Generate the mesh on the linking surface by using TFI
//***************************************************************************//
int OneToOneSwept::LinkSurfMeshing(vector<vector <Vertex> > &linkVertexList)
{
	//discretize the linking sides
	iBase_TagHandle taghandle;
	iMesh::Error m_err;
	iGeom::Error g_err;
	iRel::Error r_err;	
	
	m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
	IBERRCHK(m_err, "Trouble get tag handle of the source surface.");
	
	

	//Prepare to do the TFIMapping for linking surface
	MEntVector surfs, link_surfs;
  	mk_core()->get_entities_by_dimension(2, surfs);
	int index_id_src, index_id_tar;
	g_err = mk_core()->igeom_instance()->getIntData(sourceSurface, geom_id_tag, index_id_src);
	IBERRCHK(g_err, "Trouble get the int data for source surface.");
	g_err = mk_core()->igeom_instance()->getIntData(targetSurface, geom_id_tag, index_id_tar);
	IBERRCHK(g_err, "Trouble get the int data for target surface.");
	int index = 0;
	for (unsigned int i = 0; i < surfs.size(); i++)
	{
		int index_id_link;
		g_err = mk_core()->igeom_instance()->getIntData(surfs[i]->geom_handle(), geom_id_tag, index_id_link);
		IBERRCHK(g_err, "Trouble get the int data for linking surface.");
		if ((index_id_link != index_id_src) && (index_id_link != index_id_tar))
		{
			index++;
			link_surfs.resize(index);
			link_surfs[index - 1] = surfs[i];
			continue;		
		}
	}

	//Loop over the linking surface, and check whether the two linking edges are meshed or not
	for (unsigned int i = 0; i < link_surfs.size(); i++)
	{	
		MEntVector curves;
		curves.clear();
		link_surfs[i]->get_adjacencies(1, curves);
		
		//initialize the size function to prepare for generating the edge mesh for linking edges
		SizingFunction esize(mk_core(), numLayers, -1);
	  	link_surfs[i]->sizing_function_index(esize.core_index());
		EdgeMesher *em = (EdgeMesher*) mk_core()->construct_meshop("EdgeMesher", curves);

		em->setup_this();
		em->execute_this();
	}
	//oK we are done with the edge mesh for linking edges

	
	//extract the edge mesh for linking edge
	for (unsigned int i = 0; i < gLinkSides.size(); i++)
	{
		iBase_EntitySetHandle set;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkSides[i].gEdgeHandle, 0, set);
		IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometrical linking edge.");
		std::vector<iBase_EntityHandle> nodes;
		
		nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(set, iBase_VERTEX, iMesh_POINT, nodes);
		IBERRCHK(r_err, "Trouble get the nodes from linking edge mesh entity set.");

		assert((int)nodes.size()==(numLayers-1));
		
		//detect the edge sense
		int sense = -2;		
		g_err = mk_core()->igeom_instance()->getEgVtxSense(gLinkSides[i].gEdgeHandle, gLinkSides[i].connect[0]->gVertexHandle, gLinkSides[i].connect[1]->gVertexHandle, sense);
		IBERRCHK(g_err, "Trouble get the sense of vertex pair with respect to linking edge.");
		
		//get index for node, make sure that node index for linking edge is the same as the node index on the source surface
		int index_id;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkSides[i].connect[0]->gVertexHandle, 0, set);
		IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometrical linking edge.");
		std::vector<iBase_EntityHandle> corner;
		m_err = mk_core()->imesh_instance()->getEntities(set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, corner);
		IBERRCHK(m_err, "Trouble get the mesh node from the geometrical vertex.");		
		assert(corner.size()==1);		
		m_err = mk_core()->imesh_instance()->getIntData(corner[0], taghandle, index_id);
		IBERRCHK(m_err, "Trouble get the int data for mesh node.");		
		if (sense == 1)
		{
			for (unsigned int i = 0; i < nodes.size(); i++)
			{
				linkVertexList[i][index_id].gVertexHandle = nodes[i];
				m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[i][index_id].gVertexHandle, linkVertexList[i][index_id].xyzCoords[0], linkVertexList[i][index_id].xyzCoords[1], linkVertexList[i][index_id].xyzCoords[2]);
				IBERRCHK(m_err, "Trouble get the coordinates for mesh node.");
				linkVertexList[i][index_id].index = index_id;
				//m_err = mk_core()->imesh_instance()->getIntData(linkVertexList[i][index_id].gVertexHandle, mesh_id_tag, linkVertexList[i][index_id].id);
				//IBERRCHK(m_err, "Trouble get the int data for the vertex.");					
			}
		}
		else
		{
			for (unsigned int i = 0; i < nodes.size(); i++)
			{
				linkVertexList[nodes.size() - i - 1][index_id].gVertexHandle = nodes[i];
				m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[nodes.size() - i - 1][index_id].gVertexHandle, linkVertexList[nodes.size() - i - 1][index_id].xyzCoords[0], linkVertexList[nodes.size() - i - 1][index_id].xyzCoords[1], linkVertexList[nodes.size() - i - 1][index_id].xyzCoords[2]);
				IBERRCHK(m_err, "Trouble get the coordinates for mesh node.");
				linkVertexList[nodes.size() - i - 1][index_id].index = index_id;
				//m_err = mk_core()->imesh_instance()->getIntData(linkVertexList[nodes.size() - i - 1][index_id].gVertexHandle, mesh_id_tag, linkVertexList[nodes.size() - i - 1][index_id].id);
				//IBERRCHK(m_err, "Trouble get the int data for the vertex.");					
			}
		}
	}

	TFIMapping *tm = (TFIMapping*)mk_core()->construct_meshop("TFIMapping", link_surfs);
	tm->setup_this();
	tm->execute_this();
	//finish the meshing for linking surface

	mk_core()->save_mesh("test_linkingsurf.vtk");

	//extract the surface mesh from the linking surface
	for (unsigned int i = 0; i < gLinkFaceList.size(); i++)
	{
		iBase_EntitySetHandle set;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].gFaceHandle, 0, set);
		IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometrical linking edge.");
		std::vector<iBase_EntityHandle> nodes;
		
		nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(set, iBase_VERTEX, iMesh_POINT, nodes);
		IBERRCHK(r_err, "Trouble get the nodes from linking edge mesh entity set.");

		//need a transformation matrix to organize the nodes on the linking surface
		//TFIMapping starts the surface from edge[0] on the linking surface and node[0] on this edge[0]. Based on this, we can get the transformation Matrix
		std::vector<iBase_EntityHandle> edges;
		g_err = mk_core()->igeom_instance()->getEntAdj(gLinkFaceList[i].gFaceHandle, iBase_EDGE, edges);
		IBERRCHK(g_err, "Trouble get the adjacent edges around the linking surface.");
		assert(edges.size() == 4);
			//data structure for linking surface
			//	connect[2]----------connEdges[3]----------connect[3]
			//	   |					       |
			//     connEdges[1]				 connEdges[2]
			//         |                                           |
			//	connect[0]----------connEdges[0]----------connect[1]		

		if (edges[0] == (gLinkFaceList[i].connEdges[0]->gEdgeHandle))
		{//the edge on the source surface
			r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].connEdges[0]->gEdgeHandle, 0, set);
			IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometrical edge.");

			std::vector<iBase_EntityHandle> m_nodes;
			m_err = mk_core()->imesh_instance()->getEntities(set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, m_nodes);
			IBERRCHK(m_err, "Trouble get the number of nodes in the mesh entity set.");
			int sense = -2;
			g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[0], gLinkFaceList[i].connEdges[0]->connect[0]->gVertexHandle, gLinkFaceList[i].connEdges[0]->connect[1]->gVertexHandle, sense);
			IBERRCHK(g_err, "Trouble get the edge sense with respect to two vertices.");

			for (unsigned int k = 0; k < m_nodes.size(); k++){
				int index_id;
				if (sense ==1)
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[k], taghandle, index_id);
				else
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[m_nodes.size() - k - 1], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh node.");
				for (unsigned int j = 0; (int)j < (numLayers - 1); j++)	{						
					linkVertexList[j][index_id].gVertexHandle = nodes[j*m_nodes.size() + k];
					m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[j][index_id].gVertexHandle, linkVertexList[j][index_id].xyzCoords[0], linkVertexList[j][index_id].xyzCoords[1], linkVertexList[j][index_id].xyzCoords[2]);
					IBERRCHK(m_err, "Trouble get the coordinates of mesh node.");						
				}
			}
		}
		else if (edges[0] == (gLinkFaceList[i].connEdges[1]->gEdgeHandle))
		{//the edge on the linking edge list
			r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].connEdges[0]->gEdgeHandle, 0, set);
			IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometrical edge.");

			std::vector<iBase_EntityHandle> m_nodes;
			m_err = mk_core()->imesh_instance()->getEntities(set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, m_nodes);
			IBERRCHK(m_err, "Trouble get the number of nodes in the mesh entity set.");
			
			//g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[0], gLinkFaceList[i].connEdges[1]->connect[0]->gVertexHandle, gLinkFaceList[i].connEdges[1]->connect[1]->gVertexHandle, sense);
			std::vector<iBase_EntityHandle> gNodes;

			g_err = mk_core()->igeom_instance()->getEntAdj(gLinkFaceList[i].connEdges[1]->gEdgeHandle, iBase_VERTEX, gNodes);
			IBERRCHK(g_err, "Trouble get the adjacent vertices of the geometrical edge.");
			assert(gNodes.size()==2);

			int sense = -2;
			g_err = mk_core()->igeom_instance()->getEgVtxSense(gLinkFaceList[i].connEdges[0]->gEdgeHandle, gLinkFaceList[i].connEdges[0]->connect[0]->gVertexHandle, gLinkFaceList[i].connEdges[0]->connect[1]->gVertexHandle, sense);
			IBERRCHK(g_err, "Trouble get the edge sense with respect to two vertices.");
			
			for (unsigned int j = 0; j < m_nodes.size(); j++){
				int index_id = -1;
				if (sense == 1)
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[j], taghandle, index_id);
				else
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[m_nodes.size()-1-j], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh node.");			
			
				if ((gLinkFaceList[i].connEdges[1]->connect[0]->gVertexHandle) == gNodes[0]){
						int index_id = -1;
						if (sense == 1)
							m_err = mk_core()->imesh_instance()->getIntData(m_nodes[j], taghandle, index_id);
						else
							m_err = mk_core()->imesh_instance()->getIntData(m_nodes[m_nodes.size()-1-j], taghandle, index_id);
						IBERRCHK(m_err, "Trouble get the int data for mesh node.");
						for (unsigned int k = 0; (int)k < (numLayers - 1); k++){
							linkVertexList[k][index_id].gVertexHandle = nodes[j*(numLayers - 1) + k];
							m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[k][index_id].gVertexHandle, linkVertexList[k][index_id].xyzCoords[0], linkVertexList[k][index_id].xyzCoords[1], linkVertexList[k][index_id].xyzCoords[2]);
							IBERRCHK(m_err, "Trouble get the coordinates for mesh node.");
						}
				}
				else{
					for (int k = (numLayers - 2); k >= 0; k--){
						linkVertexList[k][index_id].gVertexHandle = nodes[j*(numLayers - 1) - k + numLayers - 2];
						m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[k][index_id].gVertexHandle, linkVertexList[k][index_id].xyzCoords[0], linkVertexList[k][index_id].xyzCoords[1], linkVertexList[k][index_id].xyzCoords[2]);
						IBERRCHK(m_err, "Trouble get the coordinates for mesh node.");
					}
				}
			}
		}
		else if (edges[0] == (gLinkFaceList[i].connEdges[2]->gEdgeHandle))
		{//the edge on the linking edge list
			r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].connEdges[0]->gEdgeHandle, 0, set);
			IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometrical edge.");

			std::vector<iBase_EntityHandle> m_nodes;
			m_err = mk_core()->imesh_instance()->getEntities(set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, m_nodes);
			IBERRCHK(m_err, "Trouble get the number of nodes in the mesh entity set.");
			
			//g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[0], gLinkFaceList[i].connEdges[1]->connect[0]->gVertexHandle, gLinkFaceList[i].connEdges[1]->connect[1]->gVertexHandle, sense);
			std::vector<iBase_EntityHandle> gNodes;

			g_err = mk_core()->igeom_instance()->getEntAdj(gLinkFaceList[i].connEdges[2]->gEdgeHandle, iBase_VERTEX, gNodes);
			IBERRCHK(g_err, "Trouble get the adjacent vertices of the geometrical edge.");
			assert(gNodes.size()==2);

			int sense = -2;
			g_err = mk_core()->igeom_instance()->getEgVtxSense(gLinkFaceList[i].connEdges[0]->gEdgeHandle, gLinkFaceList[i].connEdges[0]->connect[0]->gVertexHandle, gLinkFaceList[i].connEdges[0]->connect[1]->gVertexHandle, sense);
			IBERRCHK(g_err, "Trouble get the edge sense with respect to two vertices.");			
			
			for (unsigned int j = 0; j < m_nodes.size(); j++){
				int index_id = -1;
				if (sense == 1)
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[m_nodes.size()-1-j], taghandle, index_id);
				else
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[j], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh node.");
				if ((gLinkFaceList[i].connEdges[2]->connect[0]->gVertexHandle) == gNodes[0]){				
					for (unsigned int k = 0; int(k) < (numLayers - 1); k++){
						linkVertexList[k][index_id].gVertexHandle = nodes[j*(numLayers - 1) + k];
						m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[k][index_id].gVertexHandle, linkVertexList[k][index_id].xyzCoords[0], linkVertexList[k][index_id].xyzCoords[1], linkVertexList[k][index_id].xyzCoords[2]);
						IBERRCHK(m_err, "Trouble get the coordinates for mesh node.");
					}
				}
				else{
					
					for (int k = (numLayers - 2); k >= 0; k--){
						linkVertexList[k][index_id].gVertexHandle = nodes[j*(numLayers - 1) - k + numLayers - 2];
						m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[k][index_id].gVertexHandle, linkVertexList[k][index_id].xyzCoords[0], linkVertexList[k][index_id].xyzCoords[1], linkVertexList[k][index_id].xyzCoords[2]);
						IBERRCHK(m_err, "Trouble get the coordinates for mesh node.");
					}					
				}
			}
		}
		else{//the edge on the target surface
			iBase_TagHandle t_taghandle;
			m_err = mk_core()->imesh_instance()->getTagHandle("TargetMesh", t_taghandle);//need change a bit? use the taghandle instead?
			IBERRCHK(m_err, "Trouble get tag handle of the target surface.");

			r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].connEdges[3]->gEdgeHandle, 0, set);
			IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometrical edge.");

			std::vector<iBase_EntityHandle> m_nodes;
			m_err = mk_core()->imesh_instance()->getEntities(set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, m_nodes);
			IBERRCHK(m_err, "Trouble get the number of nodes in the mesh entity set.");
			int sense = -2;
			g_err = mk_core()->igeom_instance()->getEgVtxSense(edges[0], gLinkFaceList[i].connEdges[3]->connect[0]->gVertexHandle, gLinkFaceList[i].connEdges[3]->connect[1]->gVertexHandle, sense);
			IBERRCHK(g_err, "Trouble get the edge sense with respect to two vertices.");
			for (unsigned int k = 0; k < m_nodes.size(); k++){
				int index_id;
				if (sense ==1)
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[k], t_taghandle, index_id);
				else
					m_err = mk_core()->imesh_instance()->getIntData(m_nodes[m_nodes.size() - k - 1], t_taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh node.");
				for (unsigned int j = 0; (int)j < (numLayers - 1); j++)
				{								
					linkVertexList[numLayers -2 - j][index_id].gVertexHandle = nodes[j*m_nodes.size() + k];
					m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[numLayers -2 - j][index_id].gVertexHandle, linkVertexList[numLayers -2 - j][index_id].xyzCoords[0], linkVertexList[numLayers -2 - j][index_id].xyzCoords[1], linkVertexList[numLayers -2 - j][index_id].xyzCoords[2]);
					IBERRCHK(m_err, "Trouble get the coordinates of mesh node.");
				}
			}
		}				
	}
	
	return 1;
}



//****************************************************************************//
// function   : parametricTFI2D
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: do the transfinite interpolation in (pt_0s, pt_1s), (pt_r0, pt_r1)
//***************************************************************************//
double OneToOneSwept::parametricTFI2D(double r, double s, double pt_0s, double pt_1s, double pt_r0, double pt_r1)
{
	assert(r>= 0 && r <= 1.0);
	assert(s>= 0 && s <= 1.0);
	double pt_rs;

	//interpolate the pt_rs based on pt_r0, pt_r1, pt_0s and pt_1s
	pt_rs = 0.5*((1-s)*pt_r0 + s*pt_r1 + (1-r)*pt_0s + r*pt_1s);

	return pt_rs;
}


//****************************************************************************//
// function   : linear_interpolation 
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: interpolate linearly between x0 and x1
//***************************************************************************//
double OneToOneSwept::linear_interpolation(double r, double x0, double x1)
{
	assert(r >=0 && r <= 1.0);
	double pt= (1-r)*x0 + r*x1;
	return pt;
}

//****************************************************************************//
// function   : linear_interpolation 
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: function for obtaining the parametric coordinates from x,y,z coordinates
//***************************************************************************//
int OneToOneSwept::getUVCoords(iBase_EntityHandle gFaceHandle, Point3D pts3, Point2D &pts2)
{
	double xmin, ymin, zmin, xmax, ymax, zmax;

	iGeom::Error g_err = mk_core()->igeom_instance()->getEntBoundBox(gFaceHandle, xmin, ymin, zmin, xmax, ymax, zmax);
	IBERRCHK(g_err, "Trouble get the bounding box for the face entity.");
	if (pts3.px < xmin || pts3.px > xmax)
	{
        	cout << "Warning: Query point outside X Range [" << xmin << "," << xmax << "], x=" << pts3.px << endl;
	}
    	if (pts3.py < ymin || pts3.py > ymax)
        {
		cout << "Warning: Query point outside Y Range [" << ymin << "," << ymax << "], y=" << pts3.py << endl;
	}
    	if (pts3.pz < zmin || pts3.pz > zmax)
	{        
		cout << "Warning: Query point outside Z Range [" << zmin << "," << zmax << "], z=" << pts3.pz << endl;
	}
	g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gFaceHandle, pts3.px, pts3.py, pts3.pz, pts2.pu, pts2.pv);
	IBERRCHK(g_err, "Trouble get the parametric coordinates from x,y,z coordinates.");

	return 1;
}

//****************************************************************************//
// function   : linear_interpolation 
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: function for obtaining the x,y,z coordinates from parametric coordinates
//***************************************************************************//

int OneToOneSwept::getXYZCoords(iBase_EntityHandle gFaceHandle, Point3D &pts3, double uv[2])
{
	double umin, umax, vmin, vmax;
	iGeom::Error g_err = mk_core()->igeom_instance()->getEntUVRange(gFaceHandle, umin, vmin, umax, vmax);
	IBERRCHK(g_err, "Trouble get the parametric coordinate range.");	
	if ((uv[0]<umin)||(uv[0]>umax))
		cout << "Warning: U exceeds the range" << endl;
	if ((uv[1]<vmin)||(uv[1]>vmax))
		cout << "Warning: V exceeds the range" << endl;
	g_err = mk_core()->igeom_instance()->getEntUVtoXYZ(gFaceHandle, uv[0], uv[1], pts3.px, pts3.py, pts3.pz);
	IBERRCHK(g_err, "Trouble get the x,y,z coordinates from parametric coordinates.");	

	return 1;
}

//****************************************************************************//
// function   : TargetSurfProjection
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: map the mesh on the source surface to the target surface
//***************************************************************************//
int OneToOneSwept::TargetSurfProjection()
{
	iMesh::Error m_err;
	iGeom::Error g_err;
	iRel::Error r_err;

	//first check whether the target surface is meshed or not
	/*
	if (me->get_meshed_state() >= COMPLETE_MESH)
	{
		//get the node list, edge list, face list
		//find the corresponding edges, nodes
		return 1;
	}
	*/
	
	//get the tag handle for source surface
	iBase_TagHandle taghandle;
	m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
	IBERRCHK(m_err, "Trouble get the tag handle 'source'.");
	
	iBase_TagHandle taghandle_tar;
	m_err = mk_core()->imesh_instance()->getTagHandle("TargetMesh", taghandle_tar);
	IBERRCHK(m_err, "Trouble get the tag handle 'TargetMesh'.");
	
	
	MEntVector surfs;
	mk_core()->get_entities_by_dimension(2, surfs);
	ModelEnt *target_surf;	
	for (unsigned int i = 0; i < surfs.size(); i++)
	{
		//int index_id = -1;
		//g_err = mk_core()->igeom_instance()->getIntData(surfs[i]->geom_handle(), geom_id_tag, index_id);
		//IBERRCHK(g_err, "Trouble get the int data for surfaces.");		
		if (surfs[i]->geom_handle() == targetSurface)
		{
			target_surf = surfs[i];
			break;
		}
	}

	MEntVector curves;
	target_surf->get_adjacencies(1, curves);	
	assert(curves.size()==gsEdgeList.size());
	//define the mesh set for various edges on the source surface
	//loop over the various edges
	for (unsigned int i=0; i < gsEdgeList.size(); i++)
	{
		//get the mesh entityset for edge[i]
		std::vector<iBase_EntityHandle> nodes_src;
		iBase_EntitySetHandle mEdgeSet;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[i].gEdgeHandle, 0, mEdgeSet);
		IBERRCHK(r_err, "Trouble get the entity set for edge from source surface.");		
		
		//get the edge nodes for edge[i] mesh
		nodes_src.clear();
		m_err = mk_core()->imesh_instance()->getEntities(mEdgeSet, iBase_VERTEX, iMesh_POINT, nodes_src);
		IBERRCHK(m_err, "Trouble get the mesh edge entities.");

		int num_lines;
		m_err = mk_core()->imesh_instance()->getNumOfType(mEdgeSet, iBase_EDGE, num_lines);
		IBERRCHK(m_err, "Trouble get the number of line segments from mesh entity sets.");

		
		//initial size functon for edges, get the number of edges and assign it to the edge

		//do the edge mesher 
		SizingFunction esize(mk_core(), num_lines, -1);
  		target_surf->sizing_function_index(esize.core_index());

		//detect the edge on the target surface which corresponds to gsEdgeList[i]
		MEntVector edge_curve;
		edge_curve.resize(1);
		for (unsigned int j = 0; j < curves.size(); j++)
		{
			int index_id;
			g_err = mk_core()->igeom_instance()->getIntData(curves[j]->geom_handle(), geom_id_tag, index_id);
			IBERRCHK(g_err, "Trouble get the int data for the edge on the target surface.");
			if (index_id == gtEdgeList[edgePairs[i]].EdgeID)
			{	
				edge_curve[0] = curves[j];
				break;
			}
		}
		EdgeMesher *em = (EdgeMesher*) mk_core()->construct_meshop("EdgeMesher", edge_curve);
		//mk_core()->setup_and_execute();
		em->setup_this();
		em->execute_this();
		//done with the meshing for edge i on the target surface


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//std::cout << "edge " << i << "\tnode 0 x=" << gsEdgeList[i].connect[0]->xyzCoords[0] << "\ty=" << gsEdgeList[i].connect[0]->xyzCoords[1] << "\tz=" << gsEdgeList[i].connect[0]->xyzCoords[2] << "\tnode 1 x=" << gsEdgeList[i].connect[1]->xyzCoords[0] << "\ty=" << gsEdgeList[i].connect[1]->xyzCoords[1] << "\tz=" << gsEdgeList[i].connect[1]->xyzCoords[2] << std::endl;
		

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		//assign the edge mesh to the list: TVertexList, TEdgeList
		int sense_src, sense_tar;
		g_err = mk_core()->igeom_instance()->getEgVtxSense(gsEdgeList[i].gEdgeHandle, gsEdgeList[i].connect[0]->gVertexHandle, gsEdgeList[i].connect[1]->gVertexHandle, sense_src);
		IBERRCHK(g_err, "Trouble get the sense of edge with respect to two vertices on the source surface.");
		g_err = mk_core()->igeom_instance()->getEgVtxSense(gtEdgeList[edgePairs[i]].gEdgeHandle, gtVertexList[cornerPairs[gsEdgeList[i].connect[0]->index]].gVertexHandle, gtVertexList[cornerPairs[gsEdgeList[i].connect[1]->index]].gVertexHandle, sense_tar);
		IBERRCHK(g_err, "Trouble get the sense of edge with respect to two vertices on the target surface.");
		
		r_err = mk_core()->irel_pair()->getEntSetRelation(gtEdgeList[edgePairs[i]].gEdgeHandle, 0, mEdgeSet);
		IBERRCHK(r_err, "Trouble get the entity set for edge from target surface.");

		std::vector<iBase_EntityHandle> nodes_tar;
		nodes_tar.clear();
		m_err = mk_core()->imesh_instance()->getEntities(mEdgeSet, iBase_VERTEX, iMesh_POINT, nodes_tar);
		IBERRCHK(m_err, "Trouble get the mesh node entities for edge on the target surface.");

		assert(nodes_src.size()==nodes_tar.size());
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//for (unsigned int j = 0; j < nodes_src.size(); j++)
		//{
			//int test_id;
			//m_err = mk_core()->imesh_instance()->getIntData(nodes_tar[j], mesh_id_tag, test_id);		
			//IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the target surface.");
			//std::cout << "int data for node " << j << "\ttest_id = " << test_id << std::endl;


		//}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		if (sense_src == sense_tar)
		{
			//get the edge node index on source surface
			for (unsigned int j = 0; j < nodes_src.size(); j++)
			{
				int index_id;
				m_err = mk_core()->imesh_instance()->getIntData(nodes_src[j], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh node on the edge of source surface.");

				TVertexList[index_id].gVertexHandle = nodes_tar[j];
				TVertexList[index_id].index = index_id;
				
				m_err = mk_core()->imesh_instance()->setIntData(TVertexList[index_id].gVertexHandle, taghandle_tar, index_id);
				IBERRCHK(m_err, "Trouble set the int data for mesh node on the edge of target surface.");	

				//m_err = mk_core()->imesh_instance()->getIntData(nodes_tar[j], mesh_id_tag, TVertexList[index_id].id);		
				//IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the target surface.");

				m_err = mk_core()->imesh_instance()->getVtxCoord(nodes_tar[j], TVertexList[index_id].xyzCoords[0], TVertexList[index_id].xyzCoords[1], TVertexList[index_id].xyzCoords[2]);
				IBERRCHK(m_err, "Trouble get the mesh node coordinates on the target surface.");
		
				Point3D pts3={TVertexList[index_id].xyzCoords[0], TVertexList[index_id].xyzCoords[1], TVertexList[index_id].xyzCoords[2]};
				Point2D pts2;	
				getUVCoords(targetSurface, pts3, pts2);
				TVertexList[index_id].uvCoords[0] = pts2.pu;
				TVertexList[index_id].uvCoords[1] = pts2.pv;

				TVertexList[index_id].onBoundary = true;
				TVertexList[index_id].onCorner = false;	

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//std::cout << "*****************************************************************************same sense\n";
				//std::cout << "node coordinates on the source surface x=" << NodeList[index_id].xyzCoords[0] << "\ty=" << NodeList[index_id].xyzCoords[1] << "\tz=" << NodeList[index_id].xyzCoords[2] << std::endl;
				//std::cout << "node coordinates on the target surface x=" << TVertexList[index_id].xyzCoords[0] << "\ty=" << TVertexList[index_id].xyzCoords[1] << "\tz=" << TVertexList[index_id].xyzCoords[2] << std::endl;
				

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
			}
			
		}
		else
		{
			for (unsigned int j = 0; j < nodes_src.size(); j++)
			{
				int index_id;
				m_err = mk_core()->imesh_instance()->getIntData(nodes_src[nodes_src.size()-j-1], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh node on the edge of source surface.");

				TVertexList[index_id].gVertexHandle = nodes_tar[j];
				TVertexList[index_id].index = index_id;

				m_err = mk_core()->imesh_instance()->setIntData(TVertexList[index_id].gVertexHandle, taghandle_tar, index_id);
				IBERRCHK(m_err, "Trouble set the int data for mesh node on the edge of target surface.");
		
				//m_err = mk_core()->imesh_instance()->getIntData(nodes_tar[j], mesh_id_tag, TVertexList[index_id].id);		
				//IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the target surface.");

				
				m_err = mk_core()->imesh_instance()->getVtxCoord(nodes_tar[j], TVertexList[index_id].xyzCoords[0], TVertexList[index_id].xyzCoords[1], TVertexList[index_id].xyzCoords[2]);
				IBERRCHK(m_err, "Trouble get the mesh node coordinates on the target surface.");
		
				Point3D pts3={TVertexList[index_id].xyzCoords[0], TVertexList[index_id].xyzCoords[1], TVertexList[index_id].xyzCoords[2]};
				Point2D pts2;	
				getUVCoords(targetSurface, pts3, pts2);
				TVertexList[index_id].uvCoords[0] = pts2.pu;
				TVertexList[index_id].uvCoords[1] = pts2.pv;

				TVertexList[index_id].onBoundary = true;
				TVertexList[index_id].onCorner = false;	

			}
		}									
	}

	//set up int data for mesh nodes around the corner on the target surface
	std::vector<iBase_EntityHandle> g_Corners, m_Corners;
	g_err = mk_core()->igeom_instance()->getEntAdj(targetSurface, iBase_VERTEX, g_Corners);
	IBERRCHK(g_err, "Trouble get the geometrical vertices on the target surface.");
	for (unsigned int i = 0; i < g_Corners.size(); i++)
	{
		iBase_EntitySetHandle m_sets;
		r_err = mk_core()->irel_pair()->getEntSetRelation(g_Corners[i], 0, m_sets);
		IBERRCHK(r_err, "Trouble get the entity set for corners from target surface.");

		m_Corners.clear();
		m_err = mk_core()->imesh_instance()->getEntities(m_sets, iBase_VERTEX, iMesh_POINT, m_Corners);
		IBERRCHK(m_err, "Trouble get the mesh edge entities.");

		assert(m_Corners.size()==1);
	}

	for (unsigned int i = 0; i < TVertexList.size(); i++)
	{
		if (TVertexList[i].onCorner)
		{
			//m_err = mk_core()->imesh_instance()->setIntData(TVertexList[i].gVertexHandle, taghandle_tar, i);
			//IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
		}
	}
	//Until now, all the nodes have been created on the boundary edge.
	
	A_Matrix.resize(gBoundaries.size());
	for (unsigned int i = 0; i < gBoundaries.size(); i++)
		A_Matrix[i].resize(4);
	
	src_center.resize(gBoundaries.size());
	tar_center.resize(gBoundaries.size());	
	CalculateMatrices(A_Matrix, src_center, tar_center);	
	//we have done the transformation calculation for all the boundary loops

	iBase_EntitySetHandle entityset;  //this entityset is for storing the inner nodes on the target surface
	vector<iBase_EntityHandle>  newNodehandle(0), newEdgeHandle(0);
	

	r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, entityset);
	if (r_err) //there is no entityset associated with targetSurface
	{
		m_err = mk_core()->imesh_instance()->createEntSet(1, entityset);
		IBERRCHK(m_err, "Trouble create the entity set");
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Based on the transformation in the physical space, project the mesh nodes on the source surface to the target surface
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sPtsCenter[3] = {0.0, 0.0, 0.0}, tPtsCenter[3] = {0.0, 0.0, 0.0};
	std::vector< std::vector<double> > sBndNodes, tBndNodes;
	int num_pts = 0;
	//calculate the barycenter for the outmost boundary	
	for (unsigned int i=0; i < NodeList.size(); i++){
	    if (NodeList[i].onBoundary || NodeList[i].onCorner){//one more constraint is needed to make sure that it is on the outmost boundary
	       num_pts++;
	       for (int j = 0; j < 3; j++){
		   sPtsCenter[j] = sPtsCenter[j] + NodeList[i].xyzCoords[j];
		   tPtsCenter[j] = tPtsCenter[j] + TVertexList[i].xyzCoords[j];
	       }
	       sBndNodes.resize(num_pts);
	       tBndNodes.resize(num_pts);
	       sBndNodes[num_pts-1].resize(3);
	       tBndNodes[num_pts-1].resize(3);
	       for (int j = 0; j < 3; j++){
		   sBndNodes[num_pts-1][j] = NodeList[i].xyzCoords[j];
		   tBndNodes[num_pts-1][j] = TVertexList[i].xyzCoords[j];
	       }
	    }
	}
	for (int i = 0; i < 3; i++){
	    sPtsCenter[i] = sPtsCenter[i]/num_pts;
	    tPtsCenter[i] = tPtsCenter[i]/num_pts;
	}
	//done with the barycenter calculation

	double tmpMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
	double bMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
	//transform the coordinates
	std::vector< std::vector<double> > sBNodes(num_pts, std::vector<double>(3)), tBNodes(num_pts, std::vector<double>(3));
	for (int i = 0; i < num_pts; i++){
	    //temporarily store the boundary nodes' coordinates on the source surface and target surface
	    for (int j = 0; j < 3; j++){
		sBNodes[i][j] = sBndNodes[i][j];
		tBNodes[i][j] = tBndNodes[i][j];		
	    }

	    //transform the boundary nodes
	    for (int j = 0; j < 3; j++){
		sBNodes[i][j] = sBNodes[i][j] - 2*sPtsCenter[j] + tPtsCenter[j];
		tBNodes[i][j] = tBNodes[i][j] - sPtsCenter[j];
	    }
	}
	
	//calculate the transformation matrix: transform the nodes on the source surface to the target surface in the physical space
	for (int i = 0; i < num_pts; i++){
	    //3 row entries in the temporary matrix
	    for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
		    tmpMatrix[j][k] = tmpMatrix[j][k] + sBNodes[i][j]*sBNodes[i][k];

	    for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
		    bMatrix[j][k] = bMatrix[j][k] + tBNodes[i][j]*sBNodes[i][k];
	}

	//first determine the affine mapping matrix is singular or not
	double detValue = tmpMatrix[0][2]*tmpMatrix[1][1]*tmpMatrix[2][0] - tmpMatrix[0][1]*tmpMatrix[1][2]*tmpMatrix[2][0] - tmpMatrix[0][2]*tmpMatrix[1][0]*tmpMatrix[2][1] + tmpMatrix[0][0]*tmpMatrix[1][2]*tmpMatrix[2][1] + tmpMatrix[0][1]*tmpMatrix[1][0]*tmpMatrix[2][2] - tmpMatrix[0][0]*tmpMatrix[1][1]*tmpMatrix[2][2];
	assert(pow(detValue, 2)>1.0e-20);

	//solve the affine mapping matrix, make use of inverse matrix to get affine mapping matrix
	double InvMatrix[3][3];
	InvMatrix[0][0] = (tmpMatrix[2][1]*tmpMatrix[1][2] - tmpMatrix[1][1]*tmpMatrix[2][2])/detValue;
	InvMatrix[0][1] = (tmpMatrix[0][1]*tmpMatrix[2][2] - tmpMatrix[0][2]*tmpMatrix[2][1])/detValue;
	InvMatrix[0][2] = (tmpMatrix[0][2]*tmpMatrix[1][1] - tmpMatrix[0][1]*tmpMatrix[1][2])/detValue;
	InvMatrix[1][0] = (tmpMatrix[1][0]*tmpMatrix[2][2] - tmpMatrix[1][2]*tmpMatrix[2][0])/detValue;
	InvMatrix[1][1] = (tmpMatrix[0][2]*tmpMatrix[2][0] - tmpMatrix[0][0]*tmpMatrix[2][2])/detValue;
	InvMatrix[1][2] = (tmpMatrix[0][0]*tmpMatrix[1][2] - tmpMatrix[0][2]*tmpMatrix[1][0])/detValue;
	InvMatrix[2][0] = (tmpMatrix[1][1]*tmpMatrix[2][0] - tmpMatrix[1][0]*tmpMatrix[2][1])/detValue;
	InvMatrix[2][1] = (tmpMatrix[0][0]*tmpMatrix[2][1] - tmpMatrix[0][1]*tmpMatrix[2][0])/detValue;
	InvMatrix[2][2] = (tmpMatrix[0][1]*tmpMatrix[1][0] - tmpMatrix[0][0]*tmpMatrix[1][1])/detValue;

	double transMatrix[3][3];
	for (int i = 0; i < 3; i++)
	    for (int j = 0; j < 3; j++)
		transMatrix[i][j] = InvMatrix[j][0]*bMatrix[i][0] + InvMatrix[j][1]*bMatrix[i][1] + InvMatrix[j][2]*bMatrix[i][2];
	//finish calculating the interior nodes' location


	//calculate the inner nodes on the target surface
	for (unsigned int i=0; i < NodeList.size(); i++){
	    if (!(NodeList[i].onBoundary || NodeList[i].onCorner)){
	       double xyz[3];
	       for (int j = 0; j < 3; j++){
		   		xyz[j] = transMatrix[j][0]*(NodeList[i].xyzCoords[0] - 2*sPtsCenter[0] + tPtsCenter[0]) + transMatrix[j][1]*(NodeList[i].xyzCoords[1] - 2*sPtsCenter[1] + tPtsCenter[1]) + transMatrix[j][2]*(NodeList[i].xyzCoords[2] - 2*sPtsCenter[2] + tPtsCenter[2]) + sPtsCenter[j];
	       //calculate the closest point on the target surface with respect to the projected point
			//std::cout << "index = " << i << "\tx = " << xyz[0] << "\ty = " << xyz[1] << "\tz = " << xyz[2] << std::endl;

	       iGeom::Error g_err = mk_core()->igeom_instance()->getEntClosestPt(targetSurface, xyz[0], xyz[1], xyz[2], TVertexList[i].xyzCoords[0], TVertexList[i].xyzCoords[1], TVertexList[i].xyzCoords[2]);
	       IBERRCHK(g_err, "Trouble get the closest point on the targets surface.");
	       

	       //create the node entities
	       iMesh::Error m_err = mk_core()->imesh_instance()->createVtx(TVertexList[i].xyzCoords[0], TVertexList[i].xyzCoords[1], TVertexList[i].xyzCoords[2], TVertexList[i].gVertexHandle);
	       IBERRCHK(m_err, "Trouble create the node entity.");
		
	       //get the uv coordinates for mesh nodes on the target surface
	       
	       TVertexList[i].onBoundary = false;
	       TVertexList[i].onCorner = false;

	       newNodehandle.resize(newNodehandle.size()+1);
	       newNodehandle[newNodehandle.size()-1] = TVertexList[i].gVertexHandle;
		
	       }	       
	    }
	}

	std::cout << "==========================================================\n";

	std::cout << "Transformation matrix \n";
	std::cout << "-----\n";
	for (int i = 0; i < 3; i++){
	    for (int j = 0; j < 3; j++)
		std::cout << transMatrix[i][j] << "\t";
	    std::cout << "\n";    

	}
	std::cout << "-----\n";

	//add the inner nodes to the entityset
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&newNodehandle[0], newNodehandle.size(), entityset);
	IBERRCHK(m_err, "Trouble add an array of nodes to the entityset.");

	//until now, all the nodes have been generated on the target surface

	//determine the numbering order for quadrilateral nodes
	int sense_out, sense_out1, sense_out2, sense_out3, sense_out4;
	g_err = mk_core()->igeom_instance()->getEgFcSense(gsEdgeList[0].gEdgeHandle, sourceSurface, sense_out1);
	IBERRCHK(g_err, "Trouble get the sense of edge with respect to the face.");
	g_err = mk_core()->igeom_instance()->getEgFcSense(gtEdgeList[edgePairs[gsEdgeList[0].index]].gEdgeHandle, targetSurface, sense_out2);
	IBERRCHK(g_err, "Trouble get the sense of edge with respect to the face.");
	g_err = mk_core()->igeom_instance()->getEgVtxSense(gsEdgeList[0].gEdgeHandle, gsEdgeList[0].connect[0]->gVertexHandle, gsEdgeList[0].connect[1]->gVertexHandle, sense_out3);
	IBERRCHK(g_err, "Trouble get the sense of vertex with respect to the edge.");	
	g_err = mk_core()->igeom_instance()->getEgVtxSense(gtEdgeList[edgePairs[gsEdgeList[0].index]].gEdgeHandle, gtVertexList[cornerPairs[gsEdgeList[0].connect[0]->index]].gVertexHandle, gtVertexList[cornerPairs[gsEdgeList[0].connect[1]->index]].gVertexHandle, sense_out4);
	IBERRCHK(g_err, "Trouble get the sense of vertex with respect to the edge.");	
	sense_out = sense_out1*sense_out2*sense_out3*sense_out4;


	//create the quadrilateral elements on the target surface
	vector<iBase_EntityHandle> mFaceHandle(FaceList.size());
	vector<iBase_EntityHandle> connect(FaceList.size()*4);
	for (unsigned int i=0; i < FaceList.size(); i++)
	{
		
		
		if (sense_out < 0)
		{
			connect[4*i+0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
			connect[4*i+1] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
			connect[4*i+2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
			connect[4*i+3] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
		}
		else
		{
			connect[4*i+0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
			connect[4*i+1] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
			connect[4*i+2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
			connect[4*i+3] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
		}
		
	}
	m_err = mk_core()->imesh_instance()->createEntArr(iMesh_QUADRILATERAL, &connect[0], connect.size(), &mFaceHandle[0]);
	IBERRCHK(g_err, "Trouble create the quadrilateral mesh.");	
		
		//add the face elements on the target surface to the list
	for (unsigned int i=0; i < FaceList.size(); i++)
	{
		TFaceList[i].index = FaceList[i].index;
		TFaceList[i].gFaceHandle = mFaceHandle[i];
		TFaceList[i].connect.resize(4);
		for (int j=0; j < FaceList[i].getNumNodes(); j++)
		{
			TFaceList[i].connect[j] = &TVertexList[FaceList[i].connect[j]->index];
		}	
	}






	//add the inner face elements to the entityset
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&mFaceHandle[0], FaceList.size(), entityset);
	IBERRCHK(m_err, "Trouble add an array of quadrilateral entities to the entity set.");	
	
	//build the association
	r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, entityset);
	if (r_err) //there is no entityset associated with region[0]
	{
		r_err = mk_core()->irel_pair()->setEntSetRelation(targetSurface, entityset);
		IBERRCHK(g_err, "Trouble set the association between the target surface entity and mesh entity set.");	
	}

#ifdef HAVE_MESQUITE

	SurfMeshOptimization();   

#endif

	return 1;
}



//-----------------------------------------------------------------------------------------
//calculate the transformation matrix for all the boundary loops
void OneToOneSwept::CalculateMatrices(std::vector< std::vector<double> > &matrix, std::vector<Point2D> &src_center, std::vector<Point2D> &tar_center)
{
	//std::vector< std::set<int> > index_nodes;
	//std::vector< std::set<int> > nBoundaries;
	assert(matrix.size()==gBoundaries.size());
		
	
	nBoundaries.resize(matrix.size());
	nBndEdges.resize(matrix.size());

	//get the tag handle for the mesh on the source surface
	iBase_TagHandle taghandle = 0;
	iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
	IBERRCHK(m_err, "Trouble get the tag handle 'source'.");

	
	//loop over the nodes on the different boundaries
	for (unsigned int i = 0; i < gBoundaries.size(); i++)
	{
		std::list<int>::iterator it;
		std::vector<iBase_EntityHandle> nodes;
		std::vector<iBase_EntityHandle> edges;
		std::set<int> corners;

		for (it = gBoundaries[i].begin(); it != gBoundaries[i].end(); it++)
		{
			iBase_EntitySetHandle mesh_set;
			int index_id;

			//get the nodes on the boundaries, this doesn't include the corners
			iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[*it].gEdgeHandle, 0, mesh_set);
			IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometry edge entity handle .");
			
			edges.clear();
			m_err = mk_core()->imesh_instance()->getEntities(mesh_set, iBase_EDGE, iMesh_LINE_SEGMENT, edges);
			IBERRCHK(m_err, "Trouble get the mesh edges from the entity set");
			
			for (unsigned int j = 0; j < edges.size(); j++){
				std::vector<iBase_EntityHandle> testnodes;
				double testcoords[3];				

				testnodes.clear();
				m_err = mk_core()->imesh_instance()->getEntAdj(edges[j], iBase_VERTEX, testnodes);
				IBERRCHK(m_err, "Trouble get the adjacent nodes with respect to edge mesh");

				assert(testnodes.size()==2);
				//std::cout << "----boundary mesh edge\n";				
				m_err = mk_core()->imesh_instance()->getVtxCoord(testnodes[0], testcoords[0], testcoords[1], testcoords[2]);
				IBERRCHK(m_err, "Trouble get the x y z coordinates for the mesh nodes");
				//std::cout << "node 0 x=" << testcoords[0] << "\ty=" << testcoords[1] << "\tz=" << testcoords[2];

				m_err = mk_core()->imesh_instance()->getVtxCoord(testnodes[1], testcoords[0], testcoords[1], testcoords[2]);
				IBERRCHK(m_err, "Trouble get the x y z coordinates for the mesh nodes");
				//std::cout << "\t\tnode 1 x=" << testcoords[0] << "\ty=" << testcoords[1] << "\tz=" << testcoords[2] << std::endl;
				
				m_err = mk_core()->imesh_instance()->getIntData(edges[j], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh edges");
				nBndEdges[i].insert(index_id);
			}

			nodes.clear();
			m_err = mk_core()->imesh_instance()->getEntities(mesh_set, iBase_VERTEX, iMesh_POINT, nodes);
			IBERRCHK(m_err, "Trouble get the nodes' list from the mesh entity set.");
			
			
			for (unsigned int j = 0; j < nodes.size(); j++)
			{
				m_err = mk_core()->imesh_instance()->getIntData(nodes[j], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for nodes on the boundaries.");
				nBoundaries[i].insert(index_id);				
			}

			//extract the mesh node on the corners
			for (unsigned int j = 0; j < 2; j++)
			{
				iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[*it].connect[j]->gVertexHandle, 0, mesh_set);
				IBERRCHK(r_err, "Trouble get the edge mesh entity set from the geometry edge entity handle .");

				nodes.clear();
				m_err = mk_core()->imesh_instance()->getEntities(mesh_set, iBase_VERTEX, iMesh_POINT, nodes);
				IBERRCHK(m_err, "Trouble get the nodes' list from the mesh entity set.");

				assert(nodes.size()==1);
				
				m_err = mk_core()->imesh_instance()->getIntData(nodes[0], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for nodes on the boundaries.");

				if (corners.find(index_id) == corners.end()){
					corners.insert(index_id);
					nBoundaries[i].insert(index_id);					
				}
				else
					continue;
			}	
		}
	}

	//ok we have done with all the boundary mesh nodes on all the boundary loops


	//calculate the affine center
	for (unsigned int i = 0; i < gBoundaries.size(); i++)
	{
		src_center[i].pu = 0.0;
		src_center[i].pv = 0.0;
		tar_center[i].pu = 0.0;
		tar_center[i].pv = 0.0;

		std::set<int>::iterator it;
		for (it = nBoundaries[i].begin(); it != nBoundaries[i].end(); it++)
		{
			src_center[i].pu += NodeList[*it].uvCoords[0];
			src_center[i].pv += NodeList[*it].uvCoords[1];
			tar_center[i].pu += TVertexList[*it].uvCoords[0];
			tar_center[i].pv += TVertexList[*it].uvCoords[1];
		}

		src_center[i].pu /= nBoundaries[i].size();
		src_center[i].pv /= nBoundaries[i].size();
		tar_center[i].pu /= nBoundaries[i].size();
		tar_center[i].pv /= nBoundaries[i].size();
		
	}
	//ok, we are done with affine centroid computation

	//compute the new parametric coordinates
	std::vector< std::list<Point2D> > srcUV, tarUV;
	Point2D tmpUV;
	srcUV.resize(gBoundaries.size());
	tarUV.resize(gBoundaries.size());

	for (unsigned int i = 0; i < gBoundaries.size(); i++)
	{
		std::set<int>::iterator it;
		for (it = nBoundaries[i].begin(); it != nBoundaries[i].end(); it++)
		{
			//update the nodes on the source surface			
			tmpUV.pu = NodeList[*it].uvCoords[0] - src_center[i].pu;
			tmpUV.pv = NodeList[*it].uvCoords[1] - src_center[i].pv;
			srcUV[i].push_back(tmpUV);

			//update the nodes on the target surface
			tmpUV.pu = TVertexList[*it].uvCoords[0] - tar_center[i].pu;
			tmpUV.pv = TVertexList[*it].uvCoords[1] - tar_center[i].pv;
			tarUV[i].push_back(tmpUV);
		}
	}
	
	//compute the transformation matrix
	std::vector< std::vector<double> > temp, b1;
	temp.resize(gBoundaries.size());
	b1.resize(gBoundaries.size());
	
	//get sum of boundary nodes' coordinate
	for (unsigned int i = 0; i < gBoundaries.size(); i++)
	{
		temp[i].resize(4);
		b1[i].resize(4);
	
		for (int j = 0; j < 4; j++){		
			temp[i][j] = 0.0;
			b1[i][j] = 0;
		}			
		
		assert(srcUV[i].size()==tarUV[i].size());
		std::list<Point2D>::iterator it_src, it_tar;
		for (it_src = srcUV[i].begin(), it_tar = tarUV[i].begin(); (it_src != srcUV[i].end())&&(it_tar != tarUV[i].end()); it_src++, it_tar++)
		{
			temp[i][0] += (*it_src).pu * (*it_src).pu;
			temp[i][1] += (*it_src).pu * (*it_src).pv;
			temp[i][2] += (*it_src).pu * (*it_src).pv;
			temp[i][3] += (*it_src).pv * (*it_src).pv;

			b1[i][0] += (*it_src).pu * (*it_tar).pu;
			b1[i][1] += (*it_src).pv * (*it_tar).pu;
			b1[i][2] += (*it_src).pu * (*it_tar).pv;
			b1[i][3] += (*it_src).pv * (*it_tar).pv;			
		}		
	}
	
	//Solve the equation to get affine mapping matrix A.resize
	//check whether the equation has the solution
	for (unsigned int i = 0; i < gBoundaries.size(); i++)
	{
		assert(temp[i][0]*temp[i][3]-temp[i][1]*temp[i][2]);
		assert(matrix[i].size()==4);
		matrix[i][0] = (temp[i][3]*b1[i][0] - temp[i][1]*b1[i][1])/(temp[i][0]*temp[i][3]-temp[i][1]*temp[i][2]);
		matrix[i][1] = (temp[i][0]*b1[i][1] - temp[i][2]*b1[i][0])/(temp[i][0]*temp[i][3]-temp[i][1]*temp[i][2]);
		matrix[i][2] = (temp[i][3]*b1[i][2] - temp[i][1]*b1[i][3])/(temp[i][0]*temp[i][3]-temp[i][1]*temp[i][2]);
		matrix[i][3] = (temp[i][0]*b1[i][3] - temp[i][2]*b1[i][2])/(temp[i][0]*temp[i][3]-temp[i][1]*temp[i][2]);

		for (int j = 0; j < 4; j++)
			if (fabs(matrix[i][j]) < 1.0e-8)
				matrix[i][j] = 0.0;
	}
		
}


void OneToOneSwept::CalculateCoeffs(std::vector<double> &ConstCoeffs, std::vector< std::vector<double> > coords, std::vector< std::vector<int> > list_edge, std::vector< std::vector<int> > triangles)
{
	for (unsigned int i = 0; i < list_edge.size(); i++){
		ConstCoeffs[i] = 0.0;
		if (!(list_edge[i][2] == -1)){
			//two adjacent vertices, list_edge[i][0], list_edge[i][1]
			int index = -1;

			for (int k = 0; k < 3; k++)
				if ((triangles[list_edge[i][2]][k] != list_edge[i][0])&&(triangles[list_edge[i][2]][k] != list_edge[i][1]))
					index = triangles[list_edge[i][2]][k];
			//the vertex: index
			double eLength[3] = {0.0, 0.0, 0.0};
			eLength[0] = sqrt(pow(coords[list_edge[i][0]][0] - coords[list_edge[i][1]][0], 2) + pow(coords[list_edge[i][0]][1] - coords[list_edge[i][1]][1], 2) + pow(coords[list_edge[i][0]][2] - coords[list_edge[i][1]][2], 2));
			eLength[1] = sqrt(pow(coords[list_edge[i][0]][0] - coords[index][0], 2) + pow(coords[list_edge[i][0]][1] - coords[index][1], 2) + pow(coords[list_edge[i][0]][2] - coords[index][2], 2));
			eLength[2] = sqrt(pow(coords[index][0] - coords[list_edge[i][1]][0], 2) + pow(coords[index][1] - coords[list_edge[i][1]][1], 2) + pow(coords[index][2] - coords[list_edge[i][1]][2], 2));

			double temp = (eLength[1]*eLength[1] + eLength[2]*eLength[2] - eLength[0]*eLength[0])/(2*eLength[1]*eLength[2]);

			double angle = acos(temp);
			ConstCoeffs[i] += fabs(1/tan(angle));											
		}
					
		//second triangle list_edge[*it][3]
		if (!(list_edge[i][3] == -1)){
			//two adjacent vertices, list_edge[*it][0], list_edge[*it][1]
			int index = -1;
			for (int k = 0; k < 3; k++)
				if ((triangles[list_edge[i][3]][k] != list_edge[i][0])&&(triangles[list_edge[i][3]][k] != list_edge[i][1]))
					index = triangles[list_edge[i][3]][k];
			//the vertex: index
			double eLength[3] = {0.0, 0.0, 0.0};
			eLength[0] = sqrt(pow(coords[list_edge[i][0]][0] - coords[list_edge[i][1]][0], 2) + pow(coords[list_edge[i][0]][1] - coords[list_edge[i][1]][1], 2) + pow(coords[list_edge[i][0]][2] - coords[list_edge[i][1]][2], 2));
			eLength[1] = sqrt(pow(coords[list_edge[i][0]][0] - coords[index][0], 2) + pow(coords[list_edge[i][0]][1] - coords[index][1], 2) + pow(coords[list_edge[i][0]][2] - coords[index][2], 2));
			eLength[2] = sqrt(pow(coords[index][0] - coords[list_edge[i][1]][0], 2) + pow(coords[index][1] - coords[list_edge[i][1]][1], 2) + pow(coords[index][2] - coords[list_edge[i][1]][2], 2));

			double temp = (eLength[1]*eLength[1] + eLength[2]*eLength[2] - eLength[0]*eLength[0])/(2*eLength[1]*eLength[2]);

			double angle = acos((eLength[1]*eLength[1] + eLength[2]*eLength[2] - eLength[0]*eLength[0])/(2*eLength[1]*eLength[2]));
			ConstCoeffs[i] += fabs(1/tan(angle));
		}
	}

}

//target surface mesh smoothing by Mesquite
void OneToOneSwept::SurfMeshOptimization()
{
    //create a tag to attach the coordinates to nodes
    iBase_TagHandle mapped_tag = 0;
    iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("COORDINATES_MAP", mapped_tag);
    if (m_err){
	m_err = mk_core()->imesh_instance()->createTag("COORDINATES_MAP", 3, iBase_DOUBLE, mapped_tag);
	IBERRCHK(m_err, "Trouble create a tag.");

    }
    //attach the coordinates to the nodes
    double tag_data[3*TVertexList.size()];
    std::vector<iBase_EntityHandle> vertexHandle(TVertexList.size());
    for (unsigned int i = 0; i < NodeList.size();i++){
	tag_data[3*i] = NodeList[i].xyzCoords[0];
	tag_data[3*i+1] = NodeList[i].xyzCoords[1];
	tag_data[3*i+2] = NodeList[i].xyzCoords[2];
	vertexHandle[i] = TVertexList[i].gVertexHandle;
    }
    m_err = mk_core()->imesh_instance()->setDblArrData(&vertexHandle[0], NodeList.size(), mapped_tag, &tag_data[0]);
    IBERRCHK(m_err, "Trouble set an array of int data to nodes.");

    //get the mesh entityset for target surface
    iBase_EntitySetHandle surfSets;
    iRel::Error r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, surfSets);
    IBERRCHK(r_err, "Trouble get the mesh entity set for the target surface.");
    //call the MeshImprove class to smooth the target surface mesh by using Mesquite
    MeshImprove meshopt(mk_core(), false, true, true, true);
    meshopt.SurfMeshImprove(targetSurface, surfSets, iBase_FACE);

    //update the new position for nodes on the target surfacce
    for (unsigned int i = 0; i < TVertexList.size(); i++){
	double coords[3];
	if (!((TVertexList[i].onBoundary)||(TVertexList[i].onCorner))){
	    m_err = mk_core()->imesh_instance()->getVtxCoord(TVertexList[i].gVertexHandle, coords[0], coords[1], coords[2]);
	    IBERRCHK(m_err, "Trouble get the node's coordinates on the target surface");
	    
	    for (int j = 0; j < 3; j++)
		TVertexList[i].xyzCoords[j] = coords[j];
	}
    }
    for (unsigned int i = 0; i < TVertexList.size(); i++){
	m_err = mk_core()->imesh_instance()->setVtxCoord(TVertexList[i].gVertexHandle, TVertexList[i].xyzCoords[0], TVertexList[i].xyzCoords[1], TVertexList[i].xyzCoords[2]);
	IBERRCHK(m_err, "Trouble set a new coordinates for nodes on the target surface");
    } 
}


}

