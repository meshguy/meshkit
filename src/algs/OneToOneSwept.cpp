#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/TFIMapping.hpp"
#include <iostream>
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
	buildAssociation();	

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//get the edge list for the source surf mesh
	Edges.clear();
	m_err = mk_core()->imesh_instance()->getEntities(SourceSets, iBase_EDGE, iMesh_LINE_SEGMENT, Edges);
	IBERRCHK(m_err, "Trouble get the mesh edges from the mesh entity set.");
	
	EdgeList.resize(Edges.size());
	for (unsigned int i=0; i < Edges.size(); i++)
	{
		EdgeList[i].gEdgeHandle = Edges[i];
		EdgeList[i].index = i;
		m_err = mk_core()->imesh_instance()->getIntData(Edges[i], mesh_id_tag, EdgeList[i].EdgeID);
		IBERRCHK(m_err, "Trouble get the int data from the mesh edges.");

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
		if (isEdgeBoundary(Edges[i]))
			EdgeList[i].onBoundary = true;
		else
			EdgeList[i].onBoundary = false;
	}
	
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
	}
	
	//changed the mesh state for source surface
	MEntVector surfs;
	mk_core()->get_entities_by_dimension(2, surfs);
	for (unsigned int i = 0; i < surfs.size(); i++)
	{
		//if (surfs[i]->geom_handle() == sourceSurface)
			
	}

	//initialize the mesh size on the target surface
	TVertexList.resize(NodeList.size());
	TEdgeList.resize(EdgeList.size());
	TFaceList.resize(FaceList.size());

	//Initialize the vertex meshing on the target surface
	for (unsigned int i = 0; i < gsVertexList.size(); i++)
	{
		//get the mesh entity node from the source surface
		iBase_EntitySetHandle tmpSet;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gsVertexList[i].gVertexHandle, 0, tmpSet);
		IBERRCHK(r_err, "Trouble get the mesh entity set from the geometry entity handle.");
		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(tmpSet, iBase_VERTEX, iMesh_POINT, Nodes);
		IBERRCHK(r_err, "Trouble get the mesh entity node from the geometry entity vertex handle.");
		assert(Nodes.size()==1);
		m_err = mk_core()->imesh_instance()->getIntData(Nodes[0], taghandle, entity_index);
		IBERRCHK(r_err, "Trouble get the int data for mesh node on the source surface.");
		
		//get the mesh entity node from the target surface
		r_err = mk_core()->irel_pair()->getEntSetRelation(gtVertexList[cornerPairs[i]].gVertexHandle, 0, tmpSet);
		IBERRCHK(r_err, "Trouble get the mesh entity set from the geometry entity handle.");
		Nodes.clear();
		m_err = mk_core()->imesh_instance()->getEntities(tmpSet, iBase_VERTEX, iMesh_POINT, Nodes);
		IBERRCHK(r_err, "Trouble get the mesh entity node from the geometry entity vertex handle.");
		//assert(Nodes.size()==1);
		if (Nodes.size()==0)
		{
			Nodes.resize(1);
			m_err = mk_core()->imesh_instance()->createVtx(gtVertexList[cornerPairs[i]].xyzCoords[0], gtVertexList[cornerPairs[i]].xyzCoords[1], gtVertexList[cornerPairs[i]].xyzCoords[2], Nodes[0]);
			IBERRCHK(m_err, "Trouble create the mesh node for the geometry entity vertex.");
			m_err = mk_core()->imesh_instance()->addEntToSet(Nodes[0], tmpSet);
			IBERRCHK(m_err, "Trouble add the mesh node entity to the set.");

		}

		TVertexList[entity_index].gVertexHandle = Nodes[0];
		TVertexList[entity_index].index = entity_index;
		
		m_err = mk_core()->imesh_instance()->getIntData(Nodes[0], mesh_id_tag, TVertexList[entity_index].id);		
		IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the target surface.");

				
		m_err = mk_core()->imesh_instance()->getVtxCoord(Nodes[0], TVertexList[entity_index].xyzCoords[0], TVertexList[entity_index].xyzCoords[1], TVertexList[entity_index].xyzCoords[2]);
		IBERRCHK(m_err, "Trouble get the mesh node coordinates on the target surface.");
		
		Point3D pts3={TVertexList[entity_index].xyzCoords[0], TVertexList[entity_index].xyzCoords[1], TVertexList[entity_index].xyzCoords[2]};
		Point2D pts2;	
		getUVCoords(targetSurface, pts3, pts2);
		TVertexList[entity_index].uvCoords[0] = pts2.pu;
		TVertexList[entity_index].uvCoords[1] = pts2.pv;

		TVertexList[entity_index].onBoundary = false;
		TVertexList[entity_index].onCorner = true;
		
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
	
	std::cout << "***********************************************************" << std::endl;
	for (unsigned int i = 0; i < gFaces.size(); i++)
	{
		std::vector<iBase_EntityHandle> gVertices;
		g_err = mk_core()->igeom_instance()->getEntAdj(gFaces[i], iBase_VERTEX, gVertices);
		IBERRCHK(g_err, "Trouble get the adjacent geometrical vertices to a geometrical face.");
		std::cout << "face index = " << i << std::endl;
		std::cout << "------------------------------------------------------" << std::endl;
		for (unsigned int j = 0; j < gVertices.size(); j++)
		{
			double testCoords[3];
			g_err = mk_core()->igeom_instance()->getVtxCoord(gVertices[j], testCoords[0], testCoords[1], testCoords[2]);
			IBERRCHK(g_err, "Trouble get the coordinates from the geometrical vertex.");
			std::cout << "vertex index = " << j << "\tx = " <<  testCoords[0] << "\ty = " << testCoords[1] << "\tz = " << testCoords[2] << std::endl; 
		}
		std::cout << "------------------------------------------------------" << std::endl;
	}
	std::cout << "***********************************************************" << std::endl;
	
	//select the source surface and target surface	
	//int index_src, index_tar;
	//std::cout << "Please select the source surface:";
	//std::cin >> index_src;
	sourceSurface = gFaces[index_src];
	//std::cout << "Please select the target surface:";
	//std::cin >> index_tar;
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
	std::vector<iBase_EntityHandle> gVertices;
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
	for (unsigned int i = 0; i < gEdges.size(); i++)
	{
		gsEdgeList[i].index = i;
		gsEdgeList[i].gEdgeHandle = gEdges[i];
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
}

//---------------------------------------------------------------------------//
// execute function: generate the all-hex mesh through sweeping from source 
// surface to target surface
void OneToOneSwept::execute_this()
{
	std::vector<double> coords;
	std::vector<moab::EntityHandle> nodes;
	
	for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  	{
    		ModelEnt *me = mit -> first;
		if (me->get_meshed_state() >= COMPLETE_MESH)
			continue;

		PreprocessGeom(me);
		
		
		/********************************************************************************************************/
		std::cout << "test geometry:\n"; 
		std::cout << "test vertex pair:\n"; 		
		for (unsigned int i = 0; i < gsVertexList.size(); i++)
			std::cout << "source vertex index = " << i  << "\t target vertex index = " << cornerPairs[i] << std::endl;
		std::cout << "test Edge pair:\n";
		for (unsigned int i = 0; i < gsEdgeList.size(); i++)
			std::cout << "source edge index = " << i  << "\t target edge index = " << edgePairs[i] << std::endl;

		/********************************************************************************************************/




		

		me->boundary(0, nodes);
		
		PreprocessMesh(me);

    		//resize the coords based on the interval setting
    		numLayers = me->mesh_intervals();

		

		//necessary steps for setting up the source surface and target surfaces
    		
		TargetSurfProjection();
		mk_core()->save_mesh("target.vtk");

		/***************************************************************************************************/
		//test
		std::cout << "test nodes on the edge of source surface and target surface\n";	
		for (unsigned int i = 0; i < NodeList.size(); i++)
		{
			std::cout << "source index = " << i << "\tx = " << NodeList[i].xyzCoords[0] << "\ty = " << NodeList[i].xyzCoords[1] << "\tz = " << NodeList[i].xyzCoords[2] << "\tOnBoundary = " << NodeList[i].onBoundary << "\tOnCorner = " << NodeList[i].onCorner << std::endl;
			std::cout << "target index = " << i << "\tx = " << TVertexList[i].xyzCoords[0] << "\ty = " << TVertexList[i].xyzCoords[1] << "\tz = " << TVertexList[i].xyzCoords[2] << "\tOnBoundary = " << TVertexList[i].onBoundary << "\tOnCorner = " << TVertexList[i].onCorner << std::endl;

		}
	
		/***************************************************************************************************/
		


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

	/*
	for (unsigned int i=0; i < gLinkSides.size(); i++)
	{
		iBase_EntitySetHandle testset;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkSides[i].gEdgeHandle, 0, testset);
		if (!r_err)
		{//there is the mesh for the linking sides
			std::vector<iBase_EntityHandle> lvertices;
			m_err = mk_core()->imesh_instance()->getEntities(testset, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, lvertices);
			IBERRCHK(m_err, "Trouble get the mesh entities from the entity set.");
			if (linking_num_pts==-1)
				linking_num_pts = lvertices.size();
			else
			{
				if (int(lvertices.size())!=linking_num_pts)
				{//there are the inconsistent mesh for the different linking surfs
					cout << "There are inconsistent meshes on the different linking surfs\nThe sweeeping algorithm will terminate" << endl;
					exit(1);
				}
			}
			numLayers = lvertices.size()+1;
		}
	}
	*/
	

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
		if (m==0)
		{
			for (unsigned int i=0; i < FaceList.size(); i++)
			{
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
		else
		{
			
			for (unsigned int i=0; i < FaceList.size(); i++)
			{
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
			
			isBNodes[i][j].px = iBoundaryNodes[i][j].px;
			isBNodes[i][j].py = iBoundaryNodes[i][j].py;
			isBNodes[i][j].pz = iBoundaryNodes[i][j].pz;
			
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
				
				s = (i+1)*1/double(numLayers);
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
	int LineIndex=0, faceIndex=0;
	std::vector<iBase_EntityHandle> lineHandle(0), faceHandle(0);
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

	for (unsigned int i = 0; i < gLinkSides.size(); i++)
	{
		iBase_EntitySetHandle TestSet;
		r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkSides[i].gEdgeHandle, 0, TestSet);
		IBERRCHK(r_err, "Trouble get the adjacent edges in a face.");
		std::vector<iBase_EntityHandle> TestPoints;
		TestPoints.clear();
		m_err = mk_core()->imesh_instance()->getEntities(TestSet, iBase_VERTEX, iMesh_POINT, TestPoints);
		IBERRCHK(m_err, "Trouble get the nodes.");

		double coords[3];
		std::cout << "Linking Side index = " << i << "\tnode size = " << TestPoints.size() << std::endl;
		for (unsigned int j = 0; j < TestPoints.size(); j++)
		{
			m_err = mk_core()->imesh_instance()->getVtxCoord(TestPoints[j], coords[0], coords[1], coords[2]);
			IBERRCHK(m_err, "Trouble get the nodes.");
			std::cout << j << "\tx=" << coords[0] << "\ty=" << coords[1] << "\tz=" << coords[2] << std::endl; 
		}
	}


	TFIMapping *tm = (TFIMapping*)mk_core()->construct_meshop("TFIMapping", link_surfs);
	tm->setup_this();
	tm->execute_this();
	//finish the meshing for linking surface

	mk_core()->save_mesh("temp.vtk");	

	//check whether this linking surf is meshed or not
	vector<bool> isMeshed(gLinkFaceList.size());
	for (unsigned int i=0; i < gLinkFaceList.size(); i++)
	{
		int t1=1, t2=1, t3=1;
		int sEdgeIndex, gLeftIndex, gRightIndex;
		iBase_EntitySetHandle tmpSet;
		std::vector<iBase_EntityHandle> aEdges, aNodes, aFaces, gEdges;
		std::vector<iBase_EntityHandle> bNodes1, bNodes2, bNodes3;
		iBase_EntitySetHandle aEdgeSet;
		iBase_EntityHandle preNode, nextNode;
		std::vector<iBase_EntityHandle> Connect;
		double xyzCoords[3];
		
		gEdges.clear();
		
		g_err = mk_core()->igeom_instance()->getEntAdj(gLinkFaceList[i].gFaceHandle, iBase_EDGE, gEdges);
		IBERRCHK(g_err, "Trouble get the adjacent edges in a face.");
		//find the corresponding edge on the linking surface, which is shared by source surface
		for (unsigned int j=0; j< gEdges.size(); j++)
		{
			int edgeId;
			g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], geom_id_tag, edgeId);
			IBERRCHK(g_err, "Trouble get the int data of edge entity handle.");
			//actually, the num of linking sides is equal to the num of edges on the source surface
			for (unsigned int k=0; ((k < gsEdgeList.size())&&(t1)); k++)
			{
				if (edgeId == gsEdgeList[k].EdgeID)
				{
					sEdgeIndex=k;
					t1=0;
					break;
				}
			}
			for (unsigned int k=0; (k < gLinkSides.size())&&(t2||t3); k++)
			{
				if ((edgeId == gLinkSides[k].EdgeID)&&(t2))
				{
					gLeftIndex=k;
					t2=0;
					break;
				}
				if ((edgeId == gLinkSides[k].EdgeID)&&(t3))
				{
					gRightIndex=k;
					t3=0;
					break;
				}
			}			
		}
		if (!((gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[0]->id)||(gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[1]->id)))
		{
			int temp=gRightIndex;
			gRightIndex = gLeftIndex;
			gLeftIndex = temp;
		}
		
		
		r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].gFaceHandle, 0, tmpSet);
		if(r_err)
		{
			isMeshed[i] = true;
			
			r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[sEdgeIndex].gEdgeHandle, 0, aEdgeSet);
			IBERRCHK(r_err, "Trouble set the association between the geometry edge and mesh entity set.");
			
			m_err = mk_core()->imesh_instance()->getEntities(aEdgeSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, aNodes);
			IBERRCHK(m_err, "Trouble get the entities from the entity set.");
			
			for (unsigned int mm=0; mm < aNodes.size(); mm++)
			{
				int VertexId;
				m_err = mk_core()->imesh_instance()->getIntData(aNodes[mm], taghandle, VertexId);
				IBERRCHK(m_err, "Trouble get teh int data from node entities.");
				
				nextNode = aNodes[mm];
				preNode = aNodes[mm];				
				
				m_err = mk_core()->imesh_instance()->getVtxCoord(nextNode, xyzCoords[0], xyzCoords[1], xyzCoords[2]);
				IBERRCHK(m_err, "Trouble get the x,y,z coordinates of node entity.");
				
				for (int m2=0; m2 < numLayers-1; m2++)
				{					
					int tmpId;
					aFaces.clear();
					m_err = mk_core()->imesh_instance()->getEntAdj(nextNode, iBase_FACE, aFaces);
					IBERRCHK(r_err, "Trouble get the adjacent faces of a nodes.");					
					m_err = mk_core()->imesh_instance()->getIntData(nextNode, mesh_id_tag, tmpId);
					IBERRCHK(m_err, "Trouble get the int data of node entity.");
					
					int tmpIndex=0;
					
					Connect.clear();
					if (m2==0)
					{
						tmpIndex = 0;
						for (unsigned int m1=0; m1 < aFaces.size(); m1++)
						{
							bool is_contained;
							m_err = mk_core()->imesh_instance()->isEntContained(tmpSet, aFaces[m1], is_contained);
							IBERRCHK(m_err, "Trouble determine whether the entity is contained in the entity set.");
							if (is_contained)
							{
								tmpIndex++;
								Connect.resize(tmpIndex);
								Connect[tmpIndex-1] = aFaces[m1];
							}						
						}
					}
					else
					{
						tmpIndex = 0;
						Connect.clear();
						int preIndex;
						
						m_err = mk_core()->igeom_instance()->getIntData(preNode, mesh_id_tag, preIndex);
						IBERRCHK(m_err, "Trouble get the int data of node entity.");
						
						for (unsigned int m1=0; m1 < aFaces.size(); m1++)
						{
							bool found = false;
							bNodes3.clear();
							m_err = mk_core()->imesh_instance()->getEntAdj(aFaces[m1], iBase_VERTEX, bNodes3);
							IBERRCHK(m_err, "Trouble get the int data of node entity.");
							
							for (unsigned int m3=0; m3 < bNodes3.size(); m3++)
							{
								int tmp;
								m_err = mk_core()->imesh_instance()->getIntData(bNodes3[m3], mesh_id_tag, tmp);
								IBERRCHK(m_err, "Trouble get the int data of node entity.");
								
								m_err = mk_core()->imesh_instance()->getVtxCoord(bNodes3[m3], xyzCoords[0], xyzCoords[1], xyzCoords[2]);
								IBERRCHK(m_err, "Trouble get the x,y,z coordinates of node entity.");
								if (preIndex == tmp)
								{
									found = true;
									break;
								}
								
							}
							if (!found)
							{
								tmpIndex++;
								Connect.resize(tmpIndex);
								Connect[tmpIndex-1] = aFaces[m1];
							}
						}
					}
					
					bNodes1.clear();
					m_err = mk_core()->imesh_instance()->getEntAdj(Connect[0], iBase_VERTEX, bNodes1);
					IBERRCHK(m_err, "Trouble get the adjacent nodes.");
					
					bNodes2.clear();
					
					m_err = mk_core()->imesh_instance()->getEntAdj(Connect[1], iBase_VERTEX, bNodes2);
					IBERRCHK(m_err, "Trouble get the adjacent nodes.");
					assert((bNodes1.size()==bNodes1.size())&&(bNodes1.size()==4));
					
					for (unsigned int m1=0; m1 < bNodes1.size(); m1++)
					{
						int tmp;
						m_err = mk_core()->imesh_instance()->getIntData(bNodes1[m1], mesh_id_tag, tmp);
						IBERRCHK(m_err, "Trouble get the int data of node entity.");
						if (tmpId==tmp)
						{
							tmpIndex = m1;
							Connect[0] = bNodes1[m1];
							break;
						}
					}
					
					int option1= (tmpIndex+1)%4, option2= (3+tmpIndex)%4, index1, index2;
					m_err = mk_core()->imesh_instance()->getIntData(bNodes1[option1], mesh_id_tag, index1);
					IBERRCHK(m_err, "Trouble get the int data of node entity.");
					m_err = mk_core()->imesh_instance()->getIntData(bNodes1[option2], mesh_id_tag, index2);
					IBERRCHK(m_err, "Trouble get the int data of node entity.");
	
					for (unsigned int m1 =0; m1 < bNodes2.size(); m1++)
					{
						int tmp;
						m_err = mk_core()->imesh_instance()->getIntData(bNodes2[m1], mesh_id_tag, tmp);
						IBERRCHK(m_err, "Trouble get the int data of node entity.");
						
						if (index1==tmp)
						{
							Connect[1] = bNodes2[m1];
							linkVertexList[m2][VertexId].gVertexHandle = bNodes2[m1];
							linkVertexList[m2][VertexId].index = tmp;
							
							m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[m2][VertexId].gVertexHandle, linkVertexList[m2][VertexId].xyzCoords[0], linkVertexList[m2][VertexId].xyzCoords[1], linkVertexList[m2][VertexId].xyzCoords[2]);
							IBERRCHK(m_err, "Trouble get the x,y,z coordinates of node entity.");
							break;
						}
						if (index2==tmp)
						{
							Connect[1] = bNodes2[m1];
							linkVertexList[m2][VertexId].gVertexHandle = bNodes2[m1];
							linkVertexList[m2][VertexId].index = tmp;
							
							m_err = mk_core()->imesh_instance()->getVtxCoord(linkVertexList[m2][VertexId].gVertexHandle, linkVertexList[m2][VertexId].xyzCoords[0], linkVertexList[m2][VertexId].xyzCoords[1], linkVertexList[m2][VertexId].xyzCoords[2]);
							IBERRCHK(m_err, "Trouble get the x,y,z coordinates of node entity.");
							break;
						}
					}					
					//preNode = Connect[1];
					nextNode = Connect[1];
					if (m2 > 0)
						preNode = linkVertexList[m2-1][VertexId].gVertexHandle;
				}
				
			}
		}
		else
		{
			isMeshed[i] = false;	
		}
		
	
	
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//discretize the linking surface		
	for (unsigned int i=0; i < gLinkFaceList.size(); i++)
	{
		//find the linking sides
		if (!isMeshed[i])
		{
			int t1=1, t2=1, t3=1, leftIndex, RightIndex, nodehandleIndex=0;//temporary variable for controlling the loop
			//leftIndex and RightIndex are used to store the index for mesh corners
			double suLeft, suRight;
			//suLeft------parametrical coordinatiMesh_addEntArrToSet(mesh, &mNodeHandle[0], numLayers-1, entityset, &err);es on the source boundary edge
			Point2D pt00, pt01, pt10, pt11;  //parametric coordinates for four corners on the linking surface
			std::vector<iBase_EntityHandle> gEdges, mNodes;
			int sEdgeIndex, gLeftIndex, gRightIndex;
			iBase_EntitySetHandle edgeNodeSet;
			vector<iBase_EntityHandle>  nodeHandle(0);
		
			gEdges.clear();
			g_err = mk_core()->igeom_instance()->getEntAdj(gLinkFaceList[i].gFaceHandle, iBase_EDGE, gEdges);
			IBERRCHK(m_err, "Trouble get the adjacent edges of a surface.");

			//find the corresponding edge on the linking surface, which is shared by source surface
			for (unsigned int j=0; j< gEdges.size(); j++)
			{
				int edgeId;
				g_err = mk_core()->igeom_instance()->getIntData(gEdges[j], geom_id_tag, edgeId);
				IBERRCHK(g_err, "Trouble get the int data of edge entity.");
				//actually, the num of linking sides is equal to the num of edges on the source surface
				for (unsigned int k=0; ((k < gsEdgeList.size())&&(t1)); k++)
				{
					if (edgeId == gsEdgeList[k].EdgeID)
					{
						sEdgeIndex=k;
						t1=0;
						break;
					}
				}
				for (unsigned int k=0; (k < gLinkSides.size())&&(t2||t3); k++)
				{
					if ((edgeId == gLinkSides[k].EdgeID)&&(t2))
					{
						gLeftIndex=k;
						t2=0;
						break;
					}
					if ((edgeId == gLinkSides[k].EdgeID)&&(t3))
					{
						gRightIndex=k;
						t3=0;
						break;
					}
				}			
			}
			if (!((gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[0]->id)||(gsEdgeList[sEdgeIndex].connect[0]->id==gLinkSides[gLeftIndex].connect[1]->id)))
			{
				int temp=gRightIndex;
				gRightIndex = gLeftIndex;
				gLeftIndex = temp;
			}
		
			//parametrize the linking surface
			//loop over the layers
			r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[sEdgeIndex].gEdgeHandle, 0, edgeNodeSet);
			IBERRCHK(r_err, "Trouble get the mesh entity set of an edge entity.");
			mNodes.clear();
			m_err = mk_core()->imesh_instance()->getEntities(edgeNodeSet, iBase_VERTEX, iMesh_POINT, mNodes);
			IBERRCHK(m_err, "Trouble get the mesh entities.");
		
			//get the parametric coordinates for four corners on the linking surface
			g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[2], pt00.pu, pt00.pv);
			IBERRCHK(g_err, "Trouble get the parametric coordinates u,v from x,y,z coordinates.");
			g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[2], pt10.pu, pt10.pv);
			IBERRCHK(g_err, "Trouble get the parametric coordinates u,v from x,y,z coordinates.");

			g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[0]->index]].xyzCoords[0], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[0]->index]].xyzCoords[1], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[0]->index]].xyzCoords[2], pt01.pu, pt01.pv);
			IBERRCHK(g_err, "Trouble get the parametric coordinates u,v from x,y,z coordinates.");			

			g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[1]->index]].xyzCoords[0], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[1]->index]].xyzCoords[1], gVertexList[cornerPairs[gsEdgeList[sEdgeIndex].connect[1]->index]].xyzCoords[2], pt11.pu, pt11.pv);
			IBERRCHK(g_err, "Trouble get the parametric coordinates u,v from x,y,z coordinates.");
		
			g_err = mk_core()->igeom_instance()->getEntXYZtoU(gsEdgeList[sEdgeIndex].gEdgeHandle, gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[0]->xyzCoords[2], suLeft);
			IBERRCHK(g_err, "Trouble get the parametric coordinates u from x,y,z coordinates.");
			g_err = mk_core()->igeom_instance()->getEntXYZtoU(gsEdgeList[sEdgeIndex].gEdgeHandle, gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[0], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[1], gsEdgeList[sEdgeIndex].connect[1]->xyzCoords[2], suRight);
			IBERRCHK(g_err, "Trouble get the parametric coordinates u from x,y,z coordinates.");
		
			//find the corresponding the bottom mesh corner and top mesh corner
			iBase_EntitySetHandle tmpSet;
			std::vector<iBase_EntityHandle> tmpNodes;
			r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[sEdgeIndex].connect[0]->gVertexHandle, 0, tmpSet);
			IBERRCHK(r_err, "Trouble get the entity set from the geometry.");
			m_err = mk_core()->imesh_instance()->getEntities(tmpSet, iBase_VERTEX, iMesh_POINT, tmpNodes);
			IBERRCHK(r_err, "Trouble get the entities from the entity set.");
			assert(tmpNodes.size()==1);
			m_err = mk_core()->imesh_instance()->getIntData(tmpNodes[0], taghandle, leftIndex);
			IBERRCHK(m_err, "Trouble get the int data from the node entity.");
		
			//get the index for mesh vertices on the source edge
			r_err = mk_core()->irel_pair()->getEntSetRelation(gsEdgeList[sEdgeIndex].connect[1]->gVertexHandle, 0, tmpSet);
			IBERRCHK(r_err, "Trouble get the entity set from the geometry.");			
			m_err = mk_core()->imesh_instance()->getEntities(tmpSet, iBase_VERTEX, iMesh_POINT, tmpNodes);
			IBERRCHK(m_err, "Trouble get the entities from the entity set.");
			assert(tmpNodes.size()==1);
			m_err = mk_core()->imesh_instance()->getIntData(tmpNodes[0], taghandle, RightIndex);
			IBERRCHK(m_err, "Trouble get the int data from the node entity.");
		
			vector<int> sIndex(mNodes.size());  //sIndex is used to store the index for mesh nodes on the source surface boundary
			for (unsigned int j=0; j < mNodes.size(); j++)  //inner node on the edge
			{
				int VertexId;
				double u;//parametric coordinates on the source boundary edge
			
				Point2D ptuv[2];
				m_err = mk_core()->imesh_instance()->getIntData(mNodes[j], taghandle, VertexId);
				IBERRCHK(m_err, "Trouble get the int data from the node entity.");
						
				sIndex[j] = VertexId;
			
				
			
				//now the geometrical bottom mesh corner is NodeList[leftIndex] and top mesh corner is TVertexList[RightIndex] 
				g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, NodeList[sIndex[j]].xyzCoords[0], NodeList[sIndex[j]].xyzCoords[1], NodeList[sIndex[j]].xyzCoords[2], ptuv[0].pu, ptuv[0].pv);
				IBERRCHK(g_err, "Trouble get the parametric coordinate u,v from the x,y,z coordinates.");
		
				g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, TVertexList[sIndex[j]].xyzCoords[0], TVertexList[sIndex[j]].xyzCoords[1], TVertexList[sIndex[j]].xyzCoords[2], ptuv[1].pu, ptuv[1].pv);
				IBERRCHK(g_err, "Trouble get the parametric coordinate u,v from the x,y,z coordinates.");
		
				//get the parametric coordinates mNodes[j] on the source boundary edge gsEdgeList[sEdgeIndex]
				g_err = mk_core()->igeom_instance()->getEntXYZtoU(gsEdgeList[sEdgeIndex].gEdgeHandle, NodeList[sIndex[j]].xyzCoords[0], NodeList[sIndex[j]].xyzCoords[1], NodeList[sIndex[j]].xyzCoords[2], u);
				IBERRCHK(g_err, "Trouble get the parametric coordinate u from the x,y,z coordinates.");
		
				//calculate the inner node in the linking surface
				for (int k=0; k < (numLayers-1); k++)
				{
					Point2D pt_0s, pt_1s, pt_r0, pt_r1;
					double r, s;
					s=(k+1)*1/double(numLayers);
			
					g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, linkVertexList[k][leftIndex].xyzCoords[0], linkVertexList[k][leftIndex].xyzCoords[1], linkVertexList[k][leftIndex].xyzCoords[2], pt_0s.pu, pt_0s.pv);
					IBERRCHK(g_err, "Trouble get the parametric coordinate u,v from the x,y,z coordinates.");
					g_err = mk_core()->igeom_instance()->getEntXYZtoUV(gLinkFaceList[i].gFaceHandle, linkVertexList[k][RightIndex].xyzCoords[0], linkVertexList[k][RightIndex].xyzCoords[1], linkVertexList[k][RightIndex].xyzCoords[2], pt_1s.pu, pt_1s.pv);
					IBERRCHK(g_err, "Trouble get the parametric coordinate u,v from the x,y,z coordinates.");
					r = (u - suLeft)/(suRight - suLeft);
			
					pt_r0.pu = ptuv[0].pu;
					pt_r0.pv = ptuv[0].pv;
					pt_r1.pu = ptuv[1].pu;
					pt_r1.pv = ptuv[1].pv;
			
					linkVertexList[k][sIndex[j]].uvCoords[0] = linear_interpolation(s, pt_r0.pu, pt_r1.pu);
					linkVertexList[k][sIndex[j]].uvCoords[1] = linear_interpolation(r, pt_0s.pv, pt_1s.pv);
			
					//get the Physical coordinates for inner nodes on the linking surface
					g_err = mk_core()->igeom_instance()->getEntUVtoXYZ(gLinkFaceList[i].gFaceHandle, linkVertexList[k][sIndex[j]].uvCoords[0], linkVertexList[k][sIndex[j]].uvCoords[1], linkVertexList[k][sIndex[j]].xyzCoords[0], linkVertexList[k][sIndex[j]].xyzCoords[1], linkVertexList[k][sIndex[j]].xyzCoords[2]);
					IBERRCHK(g_err, "Trouble get the x,y,z coordinates from parametric coordinate u,v.");
			
					//create the vertex in the mesh
					nodehandleIndex++;
					nodeHandle.resize(nodehandleIndex);
					m_err = mk_core()->imesh_instance()->createVtx(linkVertexList[k][sIndex[j]].xyzCoords[0], linkVertexList[k][sIndex[j]].xyzCoords[1], linkVertexList[k][sIndex[j]].xyzCoords[2], nodeHandle[nodehandleIndex-1]);
					IBERRCHK(m_err, "Trouble create the vertex entity.");
					//add new generated vertex to the list
					linkVertexList[k][sIndex[j]].gVertexHandle = nodeHandle[nodehandleIndex-1];
			
					//create the line segments on the linking surface in the vertical direction
					vector<iBase_EntityHandle> connect(2);
					if (k==0)
					{
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = NodeList[sIndex[j]].gVertexHandle;
						connect[1] = nodeHandle[nodehandleIndex-1];
						m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
						IBERRCHK(m_err, "Trouble create the line segments in the mesh.");
				
						if (k==(numLayers-2))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = TVertexList[sIndex[j]].gVertexHandle;
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
							IBERRCHK(g_err, "Trouble create the line segments in the mesh.");
						}
				
					}
					else
					{
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = nodeHandle[nodehandleIndex-2];
						connect[1] = nodeHandle[nodehandleIndex-1];
						m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
						IBERRCHK(m_err, "Trouble create the line segments in the mesh.");
						if (k==(numLayers-2))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = TVertexList[sIndex[j]].gVertexHandle;
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
							IBERRCHK(m_err, "Trouble create the line segments in the mesh.");
						}					
					}
					//create the inner line segments on the linking surface in the horizontal direction
					if (j==0)
					{
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
						connect[1] = nodeHandle[nodehandleIndex-1];
						m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
						IBERRCHK(g_err, "Trouble create the line segments in the mesh");
						if (j == (mNodes.size()-1))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
							IBERRCHK(m_err, "Trouble create the line segment entity.");
						}
					}
					else
					{
						//create the line segments between the adjacent vertical 
						LineIndex++;
						lineHandle.resize(LineIndex);
						connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
						connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
						m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
						IBERRCHK(m_err, "Trouble create the line segment entity.");
						if (j == (mNodes.size()-1))
						{
							LineIndex++;
							lineHandle.resize(LineIndex);
							connect[0] = nodeHandle[nodehandleIndex-1];
							connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT, &connect[0], 2, lineHandle[LineIndex-1]);
							IBERRCHK(m_err, "Trouble create the line segment entity.");
						}
					}
			
			
					//create the face elements for the linking surface
					connect.resize(4);
					int sense_FaceRegion;
					g_err = mk_core()->igeom_instance()->getEntNrmlSense(gLinkFaceList[i].gFaceHandle, volEntity, sense_FaceRegion);
					IBERRCHK(g_err, "Trouble create the normal sense of a face with respect to a volume.");
					int sense_EdgeFace;
					g_err = mk_core()->igeom_instance()->getEgFcSense(gsEdgeList[sEdgeIndex].gEdgeHandle, gLinkFaceList[i].gFaceHandle, sense_EdgeFace);
					IBERRCHK(g_err, "Trouble get the sense of an edge with respect to a face.");
					int sense_vtxEdge;
					g_err = mk_core()->igeom_instance()->getEgVtxSense(gsEdgeList[sEdgeIndex].gEdgeHandle, gsEdgeList[sEdgeIndex].connect[0]->gVertexHandle, gsEdgeList[sEdgeIndex].connect[1]->gVertexHandle, sense_vtxEdge);
					IBERRCHK(g_err, "Trouble create the sense of a vertex with respect to an edge.");
			
					if ((j==0))
					{
						if (k==0)
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{					
								connect[0] = NodeList[leftIndex].gVertexHandle;
								connect[1] = mNodes[j];
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][leftIndex].gVertexHandle;
							}
							else
							{
								connect[0] = NodeList[leftIndex].gVertexHandle;
								connect[1] = linkVertexList[k][leftIndex].gVertexHandle;
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = mNodes[j];
							}
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
							IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
					
							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{	
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[leftIndex].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = TVertexList[leftIndex].gVertexHandle;
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
							}
						}
						else
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{
								connect[0] = linkVertexList[k-1][leftIndex].gVertexHandle;
								connect[1] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][leftIndex].gVertexHandle;
							}
							else
							{
								connect[0] = linkVertexList[k-1][leftIndex].gVertexHandle;
								connect[1] = linkVertexList[k][leftIndex].gVertexHandle;				
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
							}
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
							IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;	
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[leftIndex].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][leftIndex].gVertexHandle;
									connect[1] = TVertexList[leftIndex].gVertexHandle;			
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
							}
						}
						if (j== (mNodes.size()-1))
						{
							if (k==0)
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = NodeList[RightIndex].gVertexHandle;	
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;		
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = NodeList[RightIndex].gVertexHandle;
								}
					
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
									IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
								}
							}
							else
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k-1][RightIndex].gVertexHandle;
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;			
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k-1][RightIndex].gVertexHandle;
								}
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;			
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
									IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
								}
					
							}
					
						}	
					}
					else
					{
						if (k==0)
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{
								connect[0] = NodeList[sIndex[j-1]].gVertexHandle;
								connect[1] = NodeList[sIndex[j]].gVertexHandle;						
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
							}
							else
							{
								connect[0] = NodeList[sIndex[j-1]].gVertexHandle;
								connect[1] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
												
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = NodeList[sIndex[j]].gVertexHandle;
							}
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
							IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
 							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
						
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[sIndex[j-1]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = TVertexList[sIndex[j-1]].gVertexHandle;				
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");	
							}
						}
						else
						{
							faceIndex++;
							faceHandle.resize(faceIndex);
							if (sense_EdgeFace*sense_vtxEdge > 0)
							{
								connect[0] = linkVertexList[k-1][sIndex[j-1]].gVertexHandle;
								connect[1] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
							}
							else
							{
								connect[0] = linkVertexList[k-1][sIndex[j-1]].gVertexHandle;
								connect[1] = linkVertexList[k][sIndex[j-1]].gVertexHandle;			
								connect[2] = linkVertexList[k][sIndex[j]].gVertexHandle;
								connect[3] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
							}
					
							m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
							IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
							if (k==(numLayers-2))
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;	
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = TVertexList[sIndex[j-1]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k][sIndex[j-1]].gVertexHandle;
									connect[1] = TVertexList[sIndex[j-1]].gVertexHandle;	
									connect[2] = TVertexList[sIndex[j]].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
							}
					
						}
						if (j== (mNodes.size()-1))
						{
							if (k==0)
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = NodeList[RightIndex].gVertexHandle;	
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = NodeList[sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
								
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = NodeList[RightIndex].gVertexHandle;
								}
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
						
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;		
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
									IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
								}
							}
							else
							{
								faceIndex++;
								faceHandle.resize(faceIndex);
								if (sense_EdgeFace*sense_vtxEdge > 0)
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k-1][RightIndex].gVertexHandle;
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k][sIndex[j]].gVertexHandle;
								}
								else
								{
									connect[0] = linkVertexList[k-1][sIndex[j]].gVertexHandle;
									connect[1] = linkVertexList[k][sIndex[j]].gVertexHandle;
							
									connect[2] = linkVertexList[k][RightIndex].gVertexHandle;
									connect[3] = linkVertexList[k-1][RightIndex].gVertexHandle;
								}
								m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
								IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
								if (k==(numLayers-2))
								{
									faceIndex++;
									faceHandle.resize(faceIndex);
									if (sense_EdgeFace*sense_vtxEdge > 0)
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = linkVertexList[k][RightIndex].gVertexHandle;
						
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = TVertexList[sIndex[j]].gVertexHandle;
									}
									else
									{
										connect[0] = linkVertexList[k][sIndex[j]].gVertexHandle;
										connect[1] = TVertexList[sIndex[j]].gVertexHandle;			
										connect[2] = TVertexList[RightIndex].gVertexHandle;
										connect[3] = linkVertexList[k][RightIndex].gVertexHandle;
									}
									m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], 4, faceHandle[faceIndex-1]);
									IBERRCHK(m_err, "Trouble create the quadrilateral element entity.");
								}					
							}					
						}
					}			
				}
			}
	
			//create the entityset for linking surface[i]
			iBase_EntitySetHandle entityset;
			r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].gFaceHandle, 0, entityset);
			if (r_err)	//there is no entityset associated with gLinkFaceList[i].gFaceHandle
			{
				m_err = mk_core()->imesh_instance()->createEntSet(1, entityset);
				IBERRCHK(m_err, "Trouble create the entity set.");
			}
			m_err = mk_core()->imesh_instance()->addEntArrToSet(&nodeHandle[0], nodehandleIndex, entityset);
			IBERRCHK(m_err, "Trouble add an array of nodes to entity set.");
			m_err = mk_core()->imesh_instance()->addEntArrToSet(&lineHandle[0], LineIndex, entityset);
			IBERRCHK(m_err, "Trouble add an array of line segments to the entity set.");
			m_err = mk_core()->imesh_instance()->addEntArrToSet(&faceHandle[0], faceIndex, entityset);
			IBERRCHK(m_err, "Trouble add an array of face elements to the entity set.");
	
			r_err = mk_core()->irel_pair()->getEntSetRelation(gLinkFaceList[i].gFaceHandle, 0, entityset);
			if (r_err)//there is no entityset associated with gLinkFaceList[i].gFaceHandle
			{
				//create the irel between the linking sides and entityset
				r_err = mk_core()->irel_pair()->setEntSetRelation(gLinkFaceList[i].gFaceHandle, entityset);
				IBERRCHK(m_err, "Trouble set the association between the geometry face and entity set.");
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
	
	//int target_id;
	//g_err = mk_core()->igeom_instance()->getIntData(targetSurface, geom_id_tag, target_id);
	//IBERRCHK(g_err, "Trouble get the int data for target surface.");
	
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
		
				m_err = mk_core()->imesh_instance()->getIntData(nodes_tar[j], mesh_id_tag, TVertexList[index_id].id);		
				IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the target surface.");

				
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
		else
		{
			for (unsigned int j = 0; j < nodes_src.size(); j++)
			{
				int index_id;
				m_err = mk_core()->imesh_instance()->getIntData(nodes_src[nodes_src.size()-j-1], taghandle, index_id);
				IBERRCHK(m_err, "Trouble get the int data for mesh node on the edge of source surface.");

				TVertexList[index_id].gVertexHandle = nodes_tar[j];
				TVertexList[index_id].index = index_id;
		
				m_err = mk_core()->imesh_instance()->getIntData(nodes_tar[j], mesh_id_tag, TVertexList[index_id].id);		
				IBERRCHK(m_err, "Trouble get the int data for mesh nodes on the target surface.");

				
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
	//Until now, all the nodes have been created on the boundary edge.	
	//get the parametric coordinates for nodes on the boundary edges and cornersfrom the source surface and target surface
	vector<Point2D> sPtsUV(0), tPtsUV(0);
	int index=0; 
	for (unsigned int i=0; i < NodeList.size(); i++)
	{
		if (NodeList[i].onBoundary || NodeList[i].onCorner)
		{
			index++;
			sPtsUV.resize(index);
			sPtsUV[index-1].pu = NodeList[i].uvCoords[0];
			sPtsUV[index-1].pv = NodeList[i].uvCoords[1];

			tPtsUV.resize(index);
			tPtsUV[index-1].pu = TVertexList[i].uvCoords[0];
			tPtsUV[index-1].pv = TVertexList[i].uvCoords[1];
		}
	}
	
	
	//create the A matrix for mapping (affine) the inner nodes.
	//this affine mapping is done in the parametric domain
	//first get the modified coefficents
	Point2D Sc={0.0,0.0}, Tc={0.0,0.0}; //modified coordinates in order not to be singular matrix
	for (int i=0; i < index; i++)
	{
		Sc.pu = Sc.pu + sPtsUV[i].pu;
		Sc.pv = Sc.pv + sPtsUV[i].pv;
		Tc.pu = Tc.pu + tPtsUV[i].pu;
		Tc.pv = Tc.pv + tPtsUV[i].pv;

	}

	Sc.pu = Sc.pu/double(index);
	Sc.pv = Sc.pv/double(index);
	Tc.pu = Tc.pu/double(index);
	Tc.pv = Tc.pv/double(index);

	//get the new coordinates for u&v in the parametric domain	
	for (int i=0; i < index; i++)
	{
		sPtsUV[i].pu = sPtsUV[i].pu - Sc.pu;
		sPtsUV[i].pv = sPtsUV[i].pv - Sc.pv;
		tPtsUV[i].pu = tPtsUV[i].pu - Tc.pu;
		tPtsUV[i].pv = tPtsUV[i].pv - Tc.pv;
	}
	
	double temp[2][2] = {{0.0, 0.0}, {0.0, 0.0}};//define an array variable for storing the sum of boundary nodes
	double b1[2] = {0.0,0.0}, b2[2] = {0.0,0.0};
	//get sum of boundary nodes' coordinate
	for (int i=0; i < index; i++)
	{
		temp[0][0] = temp[0][0] + sPtsUV[i].pu*sPtsUV[i].pu;
		temp[0][1] = temp[0][1] + sPtsUV[i].pu*sPtsUV[i].pv;
		temp[1][0] = temp[1][0] + sPtsUV[i].pu*sPtsUV[i].pv;
		temp[1][1] = temp[1][1] + sPtsUV[i].pv*sPtsUV[i].pv;

		b1[0] = b1[0] + sPtsUV[i].pu*tPtsUV[i].pu;
		b1[1] = b1[1] + sPtsUV[i].pv*tPtsUV[i].pu;
		b2[0] = b2[0] + sPtsUV[i].pu*tPtsUV[i].pv;
		b2[1] = b2[1] + sPtsUV[i].pv*tPtsUV[i].pv;
	}

	//Solve the equation to get affine mapping matrix A.resize
	//check whether the equation has the solution
	double A[2][2]= {{0, 0}, {0, 0}}; //affine map matrix
	assert((temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0])!=0);
	A[0][0] = (temp[1][1]*b1[0] - temp[0][1]*b1[1])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	A[0][1] = (temp[0][0]*b1[1]-temp[1][0]*b1[0])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	
	A[1][0] = (temp[1][1]*b2[0] - temp[0][1]*b2[1])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	A[1][1] = (temp[0][0]*b2[1]-temp[1][0]*b2[0])/(temp[0][0]*temp[1][1]-temp[0][1]*temp[1][0]);
	
	//the affine mapping matrix A is obtained

	//mapping the inner nodes on the source surface onto the target surface
	iBase_EntitySetHandle entityset;  //this entityset is for storing the inner nodes on the target surface
	vector<iBase_EntityHandle>  newNodehandle(0), newEdgeHandle(0);
	

	r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, entityset);
	if (r_err) //there is no entityset associated with targetSurface
	{
		m_err = mk_core()->imesh_instance()->createEntSet(1, entityset);
		IBERRCHK(m_err, "Trouble create the entity set");
	}

	index = 0;
	//create the inner nodes on the target surface
	for (unsigned int i=0; i < NodeList.size(); i++)
	{
		if ((!NodeList[i].onBoundary)&&(!NodeList[i].onCorner))//make sure that the node is the inner node
		{
			//first transform the vertex on the source surface: affine mapping
			double uv[2];
			Point3D pts; 
			uv[0] = A[0][0]*(NodeList[i].uvCoords[0] - Sc.pu) + A[0][1]*(NodeList[i].uvCoords[1]- Sc.pv) + Tc.pu;
			uv[1] = A[1][0]*(NodeList[i].uvCoords[0]- Sc.pu) + A[1][1]*(NodeList[i].uvCoords[1]- Sc.pv) + Tc.pv;
			
			getXYZCoords(targetSurface, pts, uv);//maybe there is a problem here,   notice
			//create the vertex on the target surface
			
			index++;
			newNodehandle.resize(index);
			m_err = mk_core()->imesh_instance()->createVtx(pts.px, pts.py, pts.pz, newNodehandle[index-1]);
			IBERRCHK(m_err, "Trouble create the vertex entity.");

			//add the new generated vertex on the target surface into the list
			TVertexList[i].xyzCoords[0] = pts.px;
			TVertexList[i].xyzCoords[1] = pts.py;
			TVertexList[i].xyzCoords[2] = pts.pz;
			TVertexList[i].uvCoords[0] = uv[0];
			TVertexList[i].uvCoords[1] = uv[1];
			TVertexList[i].gVertexHandle = newNodehandle[index-1];
			TVertexList[i].onCorner = false;
			TVertexList[i].onBoundary = false;
			TVertexList[i].index = NodeList[i].index;	
		}
	}
	//add the inner nodes to the entityset
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&newNodehandle[0], index, entityset);
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
	for (unsigned int i=0; i < FaceList.size(); i++)
	{
		vector<iBase_EntityHandle> connect(FaceList[i].getNumNodes());
		
		if (sense_out < 0)
		{
			connect[0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
			connect[1] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
			connect[2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
			connect[3] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
		}
		else
		{
			connect[0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
			connect[1] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
			connect[2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
			connect[3] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
		}
		m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &connect[0], FaceList[i].getNumNodes(), mFaceHandle[i]);
		IBERRCHK(g_err, "Trouble create the quadrilateral entity.");	
		
		//add the face elements on the target surface to the list
		TFaceList[i].index = FaceList[i].index;
		TFaceList[i].gFaceHandle = mFaceHandle[i];
		TFaceList[i].connect.resize(FaceList[i].getNumNodes());
		for (int j=0; j < FaceList[i].getNumNodes(); j++)
		{
			TFaceList[i].connect[j] = &TVertexList[FaceList[i].connect[j]->index];
		}	
	}
	//add the inner face elements to the entityset
	m_err = mk_core()->imesh_instance()->addEntArrToSet(&mFaceHandle[0], FaceList.size(), entityset);
	IBERRCHK(g_err, "Trouble add an array of quadrilateral entities to the entity set.");	
	
	//build the association
	r_err = mk_core()->irel_pair()->getEntSetRelation(targetSurface, 0, entityset);
	if (r_err) //there is no entityset associated with region[0]
	{
		r_err = mk_core()->irel_pair()->setEntSetRelation(targetSurface, entityset);
		IBERRCHK(g_err, "Trouble set the association between the target surface entity and mesh entity set.");	
	}

	return 1;
}

//****************************************************************************//
// function   : linear_interpolation 
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: build the association between the geometry and mesh
//***************************************************************************//
void OneToOneSwept::buildAssociation()
{

    int err;

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
    err = mk_core()->imesh_instance()->getTagHandle("GEOM_DIMENSION", geom_dim_tag);
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

