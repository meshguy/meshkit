#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/TFIMapping.hpp"
#ifdef HAVE_MESQUITE
#include "meshkit/MeshImprove.hpp"
#endif
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>
#ifdef HAVE_ARMADILLO
#include "SurfProHarmonicMap.hpp"
#endif

namespace MeshKit {

//---------------------------------------------------------------------------//
//Entity Type initialization for OneToOneSwept meshing
moab::EntityType OneToOneSwept_tps[] = { moab::MBVERTEX, moab::MBQUAD, moab::MBHEX, moab::MBMAXTYPE };
const moab::EntityType* OneToOneSwept::output_types()
{
  return OneToOneSwept_tps;
}

//---------------------------------------------------------------------------//
// construction function for OneToOneSwept class
OneToOneSwept::OneToOneSwept(MKCore *mk_core, const MEntVector &me_vec) :
  MeshScheme(mk_core, me_vec)
{
}

//---------------------------------------------------------------------------//
// deconstruction function for OneToOneSwept class
OneToOneSwept::~OneToOneSwept()
{
}

// these have to be called before setup, otherwise we don't know source and target
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

//---------------------------------------------------------------------------//
// setup function: define the size between the different layers
void OneToOneSwept::setup_this()
{
  //compute the number of intervals for the associated ModelEnts, from the size set on them
  //the sizing function they point to, or a default sizing function

  if (mentSelection.empty())
    return;
  ModelEnt * firstME = (mentSelection.begin())->first;
  igeom_inst = firstME->igeom_instance();
  /*moab::Interface **/
  mb = mk_core()->moab_instance();// we should get it also from
  // model entity, if we have multiple moab instances running around

  const char *tag = "GLOBAL_ID";
  iGeom::Error g_err = igeom_inst->getTagHandle(tag, geom_id_tag);
  IBERRCHK(g_err, "Trouble get the geom_id_tag for 'GLOBAL_ID'.");

  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit->first;
    //first check to see whether entity is meshed
    if (me->get_meshed_state() >= COMPLETE_MESH || me->mesh_intervals() > 0)
      continue;

    SizingFunction *sf = mk_core()->sizing_function(me->sizing_function_index());
    if (!sf && me -> mesh_intervals() < 0 && me -> interval_firmness() == DEFAULT && mk_core()->sizing_function(0))
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
        if (me->constrain_even() && sf->intervals() % 2)
          me -> mesh_intervals(sf->intervals() + 1);
        else
          me -> mesh_intervals(sf->intervals());
        me -> interval_firmness(HARD);
      }
      else if (sf->size() > 0)
      {
        // FIXME: This calculation is bogus if size means edge lengths.
        // This should divide by the measurement of the source surface.
        int intervals = me->measure() / sf->size();
        if (!intervals)
          intervals++;
        if (me->constrain_even() && intervals % 2)
          intervals++;
        me->mesh_intervals(intervals);
        me->interval_firmness(SOFT);
      }
      else
        throw Error(MK_INCOMPLETE_MESH_SPECIFICATION, "Sizing function for edge had neither positive size nor positive intervals.");
    }

    int numIntervals = me->mesh_intervals();
    if (!sf || sf->intervals() != numIntervals)
    {
      sf = new SizingFunction(mk_core(), numIntervals, -1.);
    }

    std::vector<iBase_EntityHandle> gFaces;
    g_err = igeom_inst->getEntAdj(me->geom_handle(), iBase_FACE, gFaces);
    IBERRCHK(g_err, "Trouble get the geometrical adjacent to the solid");
    //select the source surface and target surface
    // these were setup from the beginning, it could be different for each volume!!!!
    // we need a better way to decide/select target and source surfaces
    // this also assumes that the sourceSurface is fixed between setup and execute
    sourceSurface = gFaces[index_src];
    targetSurface = gFaces[index_tar];
    MEntVector linkingSurfaces;
    int index_id_src, index_id_tar;
    iGeom::Error g_err = igeom_inst->getIntData(sourceSurface, geom_id_tag, index_id_src);
    IBERRCHK(g_err, "Trouble get the int data for source surface.");
    g_err = igeom_inst->getIntData(targetSurface, geom_id_tag, index_id_tar);
    IBERRCHK(g_err, "Trouble get the int data for target surface.");
    std::cout << " source Global ID: " << index_id_src << " target Global ID: " << index_id_tar << "\n";

    for (unsigned int i = 0; i < gFaces.size(); i++)
    {
      int gid;
      g_err = igeom_inst->getIntData(gFaces[i], geom_id_tag, gid);
      IBERRCHK(g_err, "Trouble get the int data for source surface.");
      std::vector<iBase_EntityHandle> gEdges;
      g_err = igeom_inst->getEntAdj(gFaces[i], iBase_EDGE, gEdges);
      IBERRCHK(g_err, "Trouble get the geometrical edges adjacent to the surface");
      std::cout << " face index " << i << " with ID: " << gid << " has " << gEdges.size() << " edges\n";
    }
    MEntVector surfs;
    me->get_adjacencies(2, surfs);// all surfaces adjacent to the volume
    for (unsigned int i = 0; i < surfs.size(); i++)
    {
      int index_id_link;
      g_err = igeom_inst->getIntData(surfs[i]->geom_handle(), geom_id_tag, index_id_link);
      IBERRCHK(g_err, "Trouble get the int data for linking surface.");
      if ((index_id_link != index_id_src) && (index_id_link != index_id_tar))
      {
        MEntVector edges;
        surfs[i]->get_adjacencies(1, edges);// all surfaces adjacent to the volume
        std::cout << " linking surface with global ID: " << index_id_link << " has " << edges.size() << " edges\n";
        //assert((int)edges.size()==4);

        linkingSurfaces.push_back(surfs[i]);
      }
    }
    // now create a TFI mapping
    // this will trigger edge meshers for linking edges, and for target surface edges
    TFIMapping *tm = (TFIMapping*) mk_core()->construct_meshop("TFIMapping", linkingSurfaces);
    // For now, do not assign any sizing function at all
//    for (unsigned int i = 0; i < linkingSurfaces.size(); i++)
//    {
//      linkingSurfaces[i]->sizing_function_index(sf->core_index());
//    }
    mk_core()->insert_node(tm, (MeshOp*) this);
  }

  mk_core()->print_graph("AfterOneSetup.eps");

}

//---------------------------------------------------------------------------//
// execute function: generate the all-hex mesh through sweeping from source
// surface to target surface
void OneToOneSwept::execute_this()
{

  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit -> first;
    // maybe it will be decided later?
    //case 1: if numLayers = 1, then it is not necessary to create any vertices for linking surface, All the vertices have been created by source surface and target surface
    // bLayers will contain nodes in order, layer by layer
    std::vector<moab::EntityHandle> bLayers;
    BuildLateralBoundaryLayers(me, bLayers);
    //int sizeBLayer = bLayers.size()/(numLayers+1);
    volume = me->geom_handle();

	TargetSurfProjection(bLayers);// the top layer is the last list of nodes

    //get the volume mesh entity set
    iRel::Error r_err = mk_core()->irel_pair(me->iRelPairIndex())->getEntSetRelation(me->geom_handle(), 0, volumeSet);
    IBERRCHK(r_err,
        "Trouble get the volume mesh entity set from the geometrical volume.");

    vector<vector<Vertex> > linkVertexList(numLayers - 1, vector<Vertex> (NodeList.size()));
    // this will compute the interior points distribution, in each layer
    ProjectInteriorLayers(bLayers, linkVertexList);

    CreateElements(linkVertexList);

    mk_core()->save_mesh("BeforeVolumeImprove.exo");
#if HAVE_MESQUITE
    //MeshImprove meshImpr(mk_core());
    //meshImpr.VolumeMeshImprove(volumeSet, iBase_REGION);
#endif
    me->commit_mesh(mit->second, COMPLETE_MESH);
  }

}

// will also initialize NodeList, which comprises all <Vertex> data on the source sf
void OneToOneSwept::BuildLateralBoundaryLayers(ModelEnt * me, std::vector<moab::EntityHandle> & layers)
{
  // start with the source surface
  //ModelEnt * sourceME = ModelEnt()
  iGeom * igeom_inst = me->igeom_instance();
  iGeom::Error g_err;
  moab::ErrorCode rval = moab::MB_SUCCESS;
  int index_id_src, index_id_tar;

  g_err = igeom_inst->getIntData(sourceSurface, geom_id_tag, index_id_src);
  IBERRCHK(g_err, "Trouble get the int data for source surface.");
  g_err = igeom_inst->getIntData(targetSurface, geom_id_tag, index_id_tar);
  IBERRCHK(g_err, "Trouble get the int data for target surface.");

  iBase_TagHandle taghandle = 0;
  iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
  if (m_err)
  {
    m_err = mk_core()->imesh_instance()->createTag("source",
        1, iBase_INTEGER, taghandle);
    IBERRCHK(m_err,
        "Trouble creating the source tag handle in the mesh instance.");
  }
  // TODO: If the tag exists, we could verify that the tag is the proper one.
  // Probably the tag should be put in a namespace, and perhaps it should
  //   be deleted when we are finished with it.

  moab::Range quadsOnLateralSurfaces;
  ModelEnt * sME = NULL;

  MEntVector surfs;
  me->get_adjacencies(2, surfs); // all surfaces adjacent to the volume

  for (unsigned int i = 0; i < surfs.size(); i++)
  {
    int index_id_link;
    g_err = igeom_inst->getIntData(surfs[i]->geom_handle(), geom_id_tag, index_id_link);
    IBERRCHK(g_err, "Trouble get the int data for linking surface.");
    if ((index_id_link != index_id_src) && (index_id_link != index_id_tar))
    {
      moab::EntityHandle surfSet = surfs[i]->mesh_handle();
      // get all nodes from the surface, including the boundary nodes!
      rval = mb->get_entities_by_type(surfSet, moab::MBQUAD, quadsOnLateralSurfaces);
      MBERRCHK(rval, mb);
    }
    if (index_id_link == index_id_src)
      sME = surfs[i];
  }
  // the lateral layers will be build from these nodes
  moab::Range nodesOnLateralSurfaces;
  rval = mb->get_connectivity(quadsOnLateralSurfaces, nodesOnLateralSurfaces);
  MBERRCHK(rval, mb);

  // construct NodeList, nodes on source surface
  // list of node on boundary of source
  std::vector<moab::EntityHandle> bdyNodes;
  std::vector<int> group_sizes;
  sME-> boundary(0, bdyNodes, NULL, &group_sizes);// we do not need the senses, yet

  numLoops = (int) group_sizes.size();
  sizeBLayer = (int) bdyNodes.size();

  std::cout << "Execute one to one swept ....\n";
  std::cout << "Nodes on boundary: " << bdyNodes.size() << "\n";
  std::cout << "Number of loops: " << group_sizes.size() << "\n";
  std::cout << "Nodes on lateral surfaces: " << nodesOnLateralSurfaces.size() << "\n";
  std::cout << "Quads on lateral surfaces: " << quadsOnLateralSurfaces.size() << "\n";

  numLayers = 0;
  if (!bdyNodes.empty())
    numLayers = nodesOnLateralSurfaces.size() / bdyNodes.size() - 1;
  if (bdyNodes.size() * (numLayers + 1) != nodesOnLateralSurfaces.size())
  {
    std::cout << "Major problem: number of total nodes on boundary: " << nodesOnLateralSurfaces.size() << "\n";
    std::cout << " number of nodes in one layer on the boundary:" << bdyNodes.size() << "\n";
    std::cout << " we expect bdyNodes.size()*(numLayers+1) == nodesOnLateralSurfaces.size()" << bdyNodes.size() * (numLayers + 1)
        << " != " << nodesOnLateralSurfaces.size() << "\n";
  }
  // start from this list of boundary nodes (on the source)

  // get all nodes from the source surface
  moab::Range sourceNodes;
  moab::EntityHandle surfSet = sME->mesh_handle();
  // get all nodes from the surface, including the boundary nodes!
  moab::Range quads;
  rval = mb->get_entities_by_type(surfSet, moab::MBQUAD, quads);
  MBERRCHK(rval, mb);
  // we do not really care about corners, only about boundary
  rval = mb->get_connectivity(quads, sourceNodes);
  MBERRCHK(rval, mb);

  NodeList.resize(sourceNodes.size());
  int i = 0;
  for (moab::Range::iterator it = sourceNodes.begin(); it != sourceNodes.end(); it++, i++)
  {
    moab::EntityHandle node = *it;
    NodeList[i].gVertexHandle = (iBase_EntityHandle) sourceNodes[i];
    NodeList[i].index = i;

    rval = mb->get_coords(&node, 1, &(NodeList[i].xyz[0]));
    MBERRCHK(rval, mb);
    NodeList[i].onBoundary = false;

    //set the int data to the entity
    m_err = mk_core()->imesh_instance()->setIntData(NodeList[i].gVertexHandle, taghandle, NodeList[i].index);
    IBERRCHK(m_err, "Trouble set the int value for nodes in the mesh instance.");
  }
  // now loop over the boundary and mark the boundary nodes as such
  for (unsigned int j = 0; j < bdyNodes.size(); j++)
  {
    iBase_EntityHandle node = (iBase_EntityHandle) bdyNodes[j];
    int index;
    m_err = mk_core()->imesh_instance()->getIntData(node, taghandle, index);
    IBERRCHK(m_err, "Trouble set the int value for nodes in the mesh instance.");
    NodeList[index].onBoundary = true;// some could be corners, but it is not that important
  }

  // process faces on the source
  FaceList.resize(quads.size());// quads is actually a range....

  for (unsigned int i = 0; i < quads.size(); i++)
  {
    iBase_EntityHandle currentFace = (iBase_EntityHandle) quads[i];//
    FaceList[i].gFaceHandle = currentFace;
    //FaceList[i].index = i;

    //get the nodes on the face elements
    std::vector<iBase_EntityHandle> Nodes;
    m_err = mk_core()->imesh_instance()->getEntAdj(currentFace, iBase_VERTEX, Nodes);
    IBERRCHK(m_err, "Trouble get the adjacent nodes from mesh face entity.");

    FaceList[i].connect.resize(Nodes.size());
    for (unsigned int j = 0; j < Nodes.size(); j++)
    {
      int tmpIndex = -1;
      m_err = mk_core()->imesh_instance()->getIntData(Nodes[j], taghandle, tmpIndex);
      IBERRCHK(m_err, "Trouble get the int data from node handle.");

      FaceList[i].connect[j] = &NodeList[tmpIndex];// why not store tmp Index directly?
    }
	m_err = mk_core()->imesh_instance()->setIntData(currentFace, taghandle, i);
	IBERRCHK(m_err, "Trouble set the int data for quads on the source surface.");

   /* m_err = mk_core()->imesh_instance()->setIntData(currentFace, taghandle, i);
    IBERRCHK(m_err, "Trouble set the int data for quad mesh on the source surface.");*/
  }

  std::copy(bdyNodes.begin(), bdyNodes.end(), std::back_inserter(layers));

  // mark all of the boundary nodes
  markedMoabEnts.insert(bdyNodes.begin(), bdyNodes.end());

  for (int layer = 1; layer <= numLayers; layer++) // layer with index 0 is at the bottom/source
  {
    int start_index = 0;//this is the start index in group
    for (unsigned int k = 0; k < group_sizes.size(); k++)
    {
      // first group starts at index 0, the rest are accumulated
      if (k > 0)
        start_index += group_sizes[k - 1];
      moab::EntityHandle node1 = layers[(layer - 1) * sizeBLayer + start_index];
      moab::EntityHandle node2 = layers[(layer - 1) * sizeBLayer + start_index + 1];
      moab::EntityHandle node3, node4;
      /*
       *     node3 -- node4 -->  (next layer)
       *       |        |
       *     node1 -- node2 ---> (guide layer)
       */
      rval = NodeAbove(node1, node2, quadsOnLateralSurfaces, node3, node4);
      MBERRCHK(rval, mb);
      // check that node3 and node4 are not yet marked
      if (markedMoabEnts.find(node3) != markedMoabEnts.end())
        MBERRCHK(moab::MB_FAILURE, mb);
      if (markedMoabEnts.find(node4) != markedMoabEnts.end())
        MBERRCHK(moab::MB_FAILURE, mb);
      layers.push_back(node3);//start of a new layer
      layers.push_back(node4);//node above node 2
      // mark node3 and node4
      markedMoabEnts.insert(node3);
      markedMoabEnts.insert(node4);
      // we have at least 2 nodes in every loop, otherwise we are in deep trouble
      //start
      int currentIndex = start_index + 2;
      while (currentIndex < start_index + group_sizes[k])
      {
        node1 = node2;
        node2 = layers[(layer - 1) * sizeBLayer + currentIndex];// actually, the previous layer
        node3 = node4;
        rval = FourthNodeInQuad(node1, node2, node3, quadsOnLateralSurfaces, node4);
        MBERRCHK(rval, mb);
        if (markedMoabEnts.find(node4) != markedMoabEnts.end())
          MBERRCHK(moab::MB_FAILURE, mb);
        layers.push_back(node4);
        markedMoabEnts.insert(node4);
        currentIndex++;
      }
    }// end group
  }// end layer
  for (int lay = 0; lay <= numLayers; lay++)
  {
    std::cout << "layer : " << lay << "\n";
    for (int i = 0; i < sizeBLayer; i++)
    {
      std::cout << " " << layers[lay * sizeBLayer + i];
      if (i % 10 == 9)
        std::cout << "\n";
    }
    std::cout << "\n";
  }
}

moab::ErrorCode OneToOneSwept::NodeAbove(moab::EntityHandle node1, moab::EntityHandle node2, moab::Range & latQuads,
    moab::EntityHandle & node3, moab::EntityHandle & node4)
{
  // look for node above node 1 in an unused quad
  moab::EntityHandle nds[2] = { node1, node2 };
  // find all quads connected to these 2 nodes; find one in the range that is not used
  moab::Range adj_entities;
  moab::ErrorCode rval = mb->get_adjacencies(nds, 2, 2, false, adj_entities);
  MBERRCHK(rval, mb);
  // default is -1
  moab::EntityHandle nextQuad = 0;// null
  for (unsigned int i = 0; i < adj_entities.size(); i++)
  {
    if ((markedMoabEnts.find(adj_entities[i]) == markedMoabEnts.end()) &&
        (latQuads.find(adj_entities[i]) != latQuads.end()))
    {
      // the quadrilateral is not marked but is on the lateral/linking surface
      nextQuad = adj_entities[i];
      break;
    }
  }
  if (0 == nextQuad)
    MBERRCHK(moab::MB_FAILURE, mb);

  // decide on which side are nodes 1 and 2, then get the opposite side
  const moab::EntityHandle * conn4 = NULL;
  int nnodes;
  rval = mb->get_connectivity(nextQuad, conn4, nnodes);
  MBERRCHK(rval, mb);
  int index1 = -1;
  for (index1 = 0; index1 < 4; index1++)
    if (node1 == conn4[index1])
      break;
  if (4 == index1)
    MBERRCHK(moab::MB_FAILURE, mb);
  if (node2 == conn4[(index1 + 1) % 4]) // quad is oriented node1, node2, node4, node3
  {
    node3 = conn4[(index1 + 3) % 4];
  }
  else if (node2 == conn4[(index1 +3) % 4]) // quad is oriented node2, node1, node3, node4
  {
    node3 = conn4[(index1 + 1) % 4];
  }
  else
    MBERRCHK(moab::MB_FAILURE, mb);// something is really wrong
  node4 = conn4[(index1 + 2) % 4];

  // mark the quad to be sure not to use it again
  markedMoabEnts.insert(nextQuad);

  return moab::MB_SUCCESS;
}

moab::ErrorCode OneToOneSwept::FourthNodeInQuad(moab::EntityHandle node1, moab::EntityHandle node2, moab::EntityHandle node3,
    moab::Range & latQuads, moab::EntityHandle & node4)
{
  // look for the fourth node
  moab::EntityHandle nds[3] = { node1, node2, node3 };
  // find all quads connected to these 3 nodes; find one in the range that is not used
  moab::Range adj_entities;
  moab::ErrorCode rval = mb->get_adjacencies(nds, 3, 2, false, adj_entities);
  MBERRCHK(rval, mb);

  moab::EntityHandle nextQuad = 0;// null
  for (unsigned int i = 0; i < adj_entities.size(); i++)
  {
    if ((markedMoabEnts.find(adj_entities[i]) == markedMoabEnts.end()) &&
        (latQuads.find(adj_entities[i]) != latQuads.end()))
    {
      // the quadrilateral is not marked but is on the lateral/linking surface
      nextQuad = adj_entities[i];
      break;
    }
  }
  if (0 == nextQuad)
    MBERRCHK(moab::MB_FAILURE, mb);

  // decide on which side are nodes 1 and 2, then get the opposite side
  const moab::EntityHandle * conn4 = NULL;
  int nnodes;
  rval = mb->get_connectivity(nextQuad, conn4, nnodes);
  MBERRCHK(rval, mb);

  node4 = 0;
  for (int index = 0; index < 4; index++)
  {
    moab::EntityHandle c4 = conn4[index];
    if (node1 != c4 && node2 != c4 && node3 != c4)
    {
      node4 = c4;
      break;
    }
  }
  if (0 == node4)
    MBERRCHK(moab::MB_FAILURE, mb);

  // mark the quad to be sure not to use it again
  markedMoabEnts.insert(nextQuad);

  return moab::MB_SUCCESS;

}
#ifdef HAVE_ARMADILLO
void OneToOneSwept::SurfMeshHarmonic(iBase_EntityHandle vol, vector<iBase_EntityHandle> &newNodehandle)
{
	
	SurfProHarmonicMap surf_pro(mk_core(), sourceSurface, targetSurface, vol);
	surf_pro.match();
	surf_pro.setMeshData(NodeList, TVertexList, FaceList);
	surf_pro.projection();
	surf_pro.getMeshData(TVertexList);

	iMesh::Error m_err;
	for (unsigned int i = 0; i < TVertexList.size(); i++){
		if (TVertexList[i].onBoundary)
			continue;
		m_err = mk_core()->imesh_instance()->createVtx(TVertexList[i].xyz[0], TVertexList[i].xyz[1], TVertexList[i].xyz[2], TVertexList[i].gVertexHandle);
		IBERRCHK(m_err, "Trouble create a mesh nodes on the target surface!");
		newNodehandle.push_back(TVertexList[i].gVertexHandle);
	}	
}
#endif

bool OneToOneSwept::isConcave()
{
	//use the quad mesh on the source surface to determine whether it is concave or not
	iBase_TagHandle taghandle = 0;
  	iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
  	IBERRCHK(m_err, "Trouble get the tag handle in the mesh instance.");
	for (std::vector<Vertex>::iterator it = NodeList.begin(); it != NodeList.end(); it++){
		if (!it->onBoundary)
			continue;
		//get the adjacent quads around a vertex
		double angle = 0.0;
		std::vector<iBase_EntityHandle> adj_quads;
		m_err = mk_core()->imesh_instance()->getEntAdj(it->gVertexHandle, iBase_FACE, adj_quads);
		IBERRCHK(m_err, "Trouble get the adjacent quads around a vertex!");
		for (unsigned int i = 0; i < adj_quads.size(); i++){
			int face_index = -1;
			m_err = mk_core()->imesh_instance()->getIntData(adj_quads[i], taghandle, face_index);
			if (m_err) continue;//this is a quad entity from the linking surface
			int node_index = 0, pre_index, next_index;
			for (int j = 0; j < 4; j++)
				if (FaceList[face_index].connect[j]->index == it->index){
					node_index = j;
					break;
				}
			Vector3D v1(0.0), v2(0.0);
			pre_index = (FaceList[face_index].getVertex((node_index+4-1)%4))->index;
			next_index = (FaceList[face_index].getVertex((node_index+1)%4))->index;
			v1 = NodeList[pre_index].xyz - it->xyz;
			v2 = NodeList[next_index].xyz - it->xyz;
			double dotproduct = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
			angle += acos(dotproduct/sqrt((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])*(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])));
		}
		if (angle > M_PI){
			return true;	
		}
	}

	return false;	
}


//****************************************************************************//
// function   : TargetSurfProjection
// Author     : Shengyong Cai
// Date       : Feb 15, 2011
// Description: map the mesh on the source surface to the target surface
//***************************************************************************//
int OneToOneSwept::TargetSurfProjection(std::vector<moab::EntityHandle> & bLayers)
{
  iMesh::Error m_err;
  //iGeom::Error g_err;
  iRel::Error r_err;
  moab::ErrorCode rval;

  MEntVector surfs;
  mk_core()->get_entities_by_dimension(2, surfs);
  ModelEnt *target_surf = 0;
  for (unsigned int i = 0; i < surfs.size(); i++)
  {

    if (surfs[i]->geom_handle() == targetSurface)
    {
      target_surf = surfs[i];
      break;
    }
  }
  int irelIndx = target_surf->iRelPairIndex();
  iGeom * igeom_inst = mk_core()->igeom_instance(target_surf->iGeomIndex());
  // we could have got it from model tag, too
  // something along this:
  // err = igeom_instance(geom_index)->getData(*eit, iGeomModelTags[geom_index], &this_me);

  TVertexList.resize(NodeList.size());
  // make everything not on boundary, first; the true boundary flag will be set later
  for (unsigned int i = 0; i < TVertexList.size(); i++)
  {
    TVertexList[i].onBoundary = false;
  }

  iBase_TagHandle taghandle_tar = 0; 
  m_err = mk_core()->imesh_instance()->getTagHandle("TargetMesh",
      taghandle_tar);
  if (m_err)
  {
    m_err = mk_core()->imesh_instance()->createTag("TargetMesh",
        1, iBase_INTEGER, taghandle_tar);
    IBERRCHK(m_err,
        "Trouble create the target mesh tag handle in the mesh instance.");
  }
  // TODO: If the tag exists, we could verify that the tag is the proper one.
  // Probably the tag should be put in a namespace, and perhaps it should
  //   be deleted when we are finished with it.

  iBase_TagHandle taghandle = 0;
  m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
  IBERRCHK(m_err, "Trouble get tag handle of the source surface.");

  // we know the first layer of nodes (on the source)
  //
  // the first 0 , 1, ..., sizeBlayer - 1 are on bottom, were bdyNodes before
  //
  for (int i = 0; i < sizeBLayer; i++)
  {
    iBase_EntityHandle sBNode = (iBase_EntityHandle) bLayers[i];
    int index;
    m_err = mk_core()->imesh_instance()->getIntData(sBNode, taghandle, index);
    IBERRCHK(m_err, "Trouble get the int value for nodes in the mesh instance.");
    TVertexList[index].onBoundary = true;// some could be corners, but it is not that important
    moab::EntityHandle topNode = bLayers[numLayers * sizeBLayer + i];// node right above
    iBase_EntityHandle tBNode = (iBase_EntityHandle) topNode;// it is right above node i in loop
    TVertexList[index].gVertexHandle = tBNode;
    // now get the coordinates of the tBNode (node on boundary of target)
    rval = mb->get_coords(&topNode, 1, &(TVertexList[index].xyz[0]));
    MBERRCHK(rval, mb);

  }

  iBase_EntitySetHandle entityset; //this entityset is for storing the inner nodes on the target surface
  vector<iBase_EntityHandle> newNodehandle;

  r_err = mk_core()->irel_pair(irelIndx)->getEntSetRelation(targetSurface, 0, entityset);
  if (r_err) //there is no entityset associated with targetSurface; this should be an error
  {
    m_err = mk_core()->imesh_instance()->createEntSet(1, entityset);
    IBERRCHK(m_err, "Trouble create the entity set");
  }

 // bool condition_harmonic = isConcave();
  bool is_exe_harmonic = false;
#ifdef HAVE_ARMADILLO
  //target surface mesh projection based on harmonic mapping
  if (condition_harmonic){
	is_exe_harmonic = true;
	SurfMeshHarmonic(volume, newNodehandle);
  }
#endif

   if (!is_exe_harmonic){
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  // Based on the transformation in the physical space, project the mesh nodes on the source surface to the target surface
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  Vector3D sPtsCenter(0.0), tPtsCenter(0.0);
	  std::vector<Vector3D> sBndNodes(sizeBLayer), tBndNodes(sizeBLayer);
	  int num_pts = 0;
	  //calculate the barycenters for the source and target boundaries
	  for (unsigned int i = 0; i < NodeList.size(); i++)
	  {
		if (NodeList[i].onBoundary)
		{
		  sPtsCenter += NodeList[i].xyz;
		  tPtsCenter += TVertexList[i].xyz;
		  sBndNodes[num_pts] = NodeList[i].xyz;
		  tBndNodes[num_pts] = TVertexList[i].xyz;
		  num_pts++;
		}
	  }
	  if (sizeBLayer != num_pts)
		MBERRCHK(moab::MB_FAILURE, mb);
	  sPtsCenter = sPtsCenter / num_pts;
	  tPtsCenter = tPtsCenter / num_pts;
	  //done with the barycenter calculation

	  // compute transformation matrix from source boundary to target boundary
	  Matrix3D transMatrix;
	  computeTransformationFromSourceToTarget(sBndNodes, sPtsCenter, tBndNodes, tPtsCenter, transMatrix);

	  //calculate the inner nodes on the target surface
	  for (unsigned int i = 0; i < NodeList.size(); i++)
	  {
		if (!NodeList[i].onBoundary)
		{
		  Vector3D xyz;
		  xyz = transMatrix * (NodeList[i].xyz - 2 * sPtsCenter + tPtsCenter)+sPtsCenter;

		  iGeom::Error g_err = igeom_inst->getEntClosestPt(targetSurface, xyz[0], xyz[1], xyz[2], TVertexList[i].xyz[0],
		      TVertexList[i].xyz[1], TVertexList[i].xyz[2]);
		  IBERRCHK(g_err, "Trouble get the closest point on the targets surface.");

		  //create the node entities
		  iMesh::Error m_err = mk_core()->imesh_instance()->createVtx(TVertexList[i].xyz[0], TVertexList[i].xyz[1],
		      TVertexList[i].xyz[2], TVertexList[i].gVertexHandle);
		  IBERRCHK(m_err, "Trouble create the node entity.");
		  TVertexList[i].onBoundary = false;
		  newNodehandle.push_back(TVertexList[i].gVertexHandle);
		}
	  }
  }
  //add the inner nodes to the entityset
  m_err = mk_core()->imesh_instance()->addEntArrToSet(&newNodehandle[0], newNodehandle.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of nodes to the entityset.");

  //until now, all the nodes have been generated on the target surface

  // we need to decide the orientation with respect to the faces on the source (which are oriented correctly)
  // look at the first 2 nodes in the boundary, in target
  int sense_out = 1;
  std::vector<moab::EntityHandle> bdyTNodes;
  std::vector<int> group_sizes;
  target_surf-> boundary(0, bdyTNodes, NULL, &group_sizes);// we do not need the senses, yet. Or do we?
  // find the first node and second node in the bLayers[] list corresponding to the top
  moab::EntityHandle node1 = bdyTNodes[1], node2 = bdyTNodes[2];
  int index = -1;
  for (index = 0; index < sizeBLayer; index++)
  {
    if (bLayers[numLayers * sizeBLayer + index] == node1)
      break;
  }
  if (index == sizeBLayer)
    MBERRCHK(moab::MB_FAILURE, mb);

  int prevIndex = (index + sizeBLayer -1) % sizeBLayer;
  int nextIndex = (index + 1) % sizeBLayer;
  if (bLayers[numLayers * sizeBLayer + prevIndex] == node2)
    sense_out = -1;
  else if (bLayers[numLayers * sizeBLayer + nextIndex] == node2)
    sense_out = 1;
  else
    MBERRCHK(moab::MB_FAILURE, mb);// serious error in the logic, can't decide orientation


  //create the quadrilateral elements on the target surface
  vector<iBase_EntityHandle> mFaceHandle(FaceList.size());
  vector<iBase_EntityHandle> connect(FaceList.size() * 4);
  for (unsigned int i = 0; i < FaceList.size(); i++)
  {
    if (sense_out < 0)
    {
      connect[4 * i + 0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
      connect[4 * i + 1] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
      connect[4 * i + 2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
      connect[4 * i + 3] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
    }
    else
    {
      connect[4 * i + 0] = TVertexList[(FaceList[i].getVertex(0))->index].gVertexHandle;
      connect[4 * i + 1] = TVertexList[(FaceList[i].getVertex(1))->index].gVertexHandle;
      connect[4 * i + 2] = TVertexList[(FaceList[i].getVertex(2))->index].gVertexHandle;
      connect[4 * i + 3] = TVertexList[(FaceList[i].getVertex(3))->index].gVertexHandle;
    }

  }
  m_err = mk_core()->imesh_instance()->createEntArr(iMesh_QUADRILATERAL, &connect[0], connect.size(), &mFaceHandle[0]);
  IBERRCHK(m_err, "Trouble create the quadrilateral mesh.");

  //add the inner face elements to the entityset
  m_err = mk_core()->imesh_instance()->addEntArrToSet(&mFaceHandle[0], FaceList.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of quadrilateral entities to the entity set.");

  mk_core()->save_mesh("BeforeHex.vtk");
#ifdef HAVE_MESQUITE

  SurfMeshOptimization();

#endif

  return 1;
}
// input: list of nodes on source, boundary center, list of nodes on target, target center
// output: 3x3 matrix A such that
//  target= A * ( source - 2*sc + tc) + sc
void OneToOneSwept::computeTransformationFromSourceToTarget(std::vector<Vector3D> & sNodes, Vector3D & sc,
      std::vector<Vector3D> & tNodes, Vector3D & tc, Matrix3D & transMatrix)
{
  Matrix3D tmpMatrix(0.0);
  Matrix3D bMatrix(0.0);
  //transform the coordinates
  int size = (int)sNodes.size();
  assert(sNodes.size() == tNodes.size());
  for (int i = 0; i < size; i++)
  {
    //transform the boundary nodes
    sNodes[i] = sNodes[i] - 2 * sc + tc;
    tNodes[i] = tNodes[i] - sc;
  }

  //calculate the transformation matrix: transform the nodes on the source surface to the target surface in the physical space
  for (int i = 0; i < size; i++)
  {
    tmpMatrix = tmpMatrix + sNodes[i]*transpose(sNodes[i]);
    bMatrix = bMatrix + tNodes[i]*transpose(sNodes[i]);
  }

  //first determine the affine mapping matrix is singular or not
  double detValue = det(tmpMatrix);
  (void) detValue;
  assert(pow(detValue, 2)>1.0e-20);

  //solve the affine mapping matrix, make use of inverse matrix to get affine mapping matrix
  Matrix3D InvMatrix = inverse(tmpMatrix);
  transMatrix = bMatrix*InvMatrix;

}
// use similar code to TargetSurfProjection, but do not project on surface...
int OneToOneSwept::ProjectInteriorLayers(std::vector<moab::EntityHandle> & boundLayers, vector<vector<Vertex> > & linkVertexList)
{
  iMesh::Error m_err;
  Vector3D sPtsCenter(0.), tPtsCenter(0.);
  std::vector<Vector3D> PtsCenter(numLayers - 1);
  std::vector<Vector3D> sBoundaryNodes(0), tBoundaryNodes(0);

  iBase_TagHandle taghandle = 0;
  m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
  IBERRCHK(m_err, "Trouble getting the source tag handle");

  std::vector<Vector3D> latCoords;
  latCoords.resize(boundLayers.size() - 2 * sizeBLayer);// to store the coordinates of the lateral points
  // (exclude the source and target, they are already in NodeList and TVertexList)...

  moab::ErrorCode rval = mb->get_coords(&(boundLayers[sizeBLayer]), latCoords.size(), &(latCoords[0][0])); //these are the nodes for lateral sides
  MBERRCHK(rval, mb);

  //calculate the center coordinates
  for (int i = 0; i < numLayers - 1; i++)
    PtsCenter[i].set(0.);

  int numPts = sizeBLayer;// a lot of repetition ; do we really need
  sBoundaryNodes.resize(numPts);// another useless copy
  tBoundaryNodes.resize(numPts);// another useless copy
  for (int j = 0; j < sizeBLayer; j++)
  {
    iBase_EntityHandle sBNode = (iBase_EntityHandle) boundLayers[j];
    int i;
    // this is the index in NodeList...
    m_err = mk_core()->imesh_instance()->getIntData(sBNode, taghandle, i);
    IBERRCHK(m_err, "Trouble get the int value for nodes in the mesh instance.");

    sPtsCenter += NodeList[i].xyz;
    tPtsCenter += TVertexList[i].xyz;

    sBoundaryNodes[j] = NodeList[i].xyz;
    tBoundaryNodes[j] = TVertexList[i].xyz;

    for (int lay = 0; lay < numLayers - 1; lay++)
    {
      // interior layers
      int indexNodeOnSide = lay * sizeBLayer + j;// all p3d will be above index i
      Vector3D & p3d = latCoords[indexNodeOnSide];
      PtsCenter[lay] += p3d;
      linkVertexList[lay][i].xyz = p3d;
      linkVertexList[lay][i].gVertexHandle = (iBase_EntityHandle) boundLayers[indexNodeOnSide + sizeBLayer];
    }
  }

  sPtsCenter = sPtsCenter / numPts;
  tPtsCenter = tPtsCenter / numPts;

  //calculate the center coordinates for the ith layer
  for (int i = 0; i < numLayers - 1; i++)
    PtsCenter[i] = PtsCenter[i] / numPts;

  std::vector<Vector3D> sNodes(numPts); //boundary nodes on the source surface
  std::vector<Vector3D> tNodes(numPts); //boundary nodes on the target surface
  std::vector<Vector3D> isBNodes(numPts), itBNodes(numPts); //boundary nodes on the different layers

  //loop over different layers
  for (int i = 0; i < numLayers - 1; i++)
  {
    for (int j = 0; j < numPts; j++)
    {
      sNodes[j] = sBoundaryNodes[j];
      tNodes[j] = tBoundaryNodes[j];
      //coordinates on the different layers
      isBNodes[j] = itBNodes[j] = latCoords[i*sizeBLayer + j];
    }

    Matrix3D sA, tA;
    // from source to layer i
    computeTransformationFromSourceToTarget(sNodes, sPtsCenter, isBNodes, PtsCenter[i], sA);
    // from target to layer i
    computeTransformationFromSourceToTarget(tNodes, tPtsCenter, itBNodes, PtsCenter[i], tA);
    //calculate the inner nodes for different layers
    double s = (i + 1) / double(numLayers); // interpolation factor
    for (unsigned int j = 0; j < NodeList.size(); j++)
    {
      if (!(NodeList[j].onBoundary))
      {
        Vector3D spts, tpts, pts;
        iBase_EntityHandle nodeHandle;
        spts = sA * (NodeList[j].xyz - 2 * sPtsCenter + PtsCenter[i])+sPtsCenter;
        tpts = tA * (TVertexList[j].xyz - 2 * tPtsCenter + PtsCenter[i])+tPtsCenter;
        // interpolate the 2 results
        pts = (1-s) * spts + s * tpts;

        linkVertexList[i][j].xyz = pts;

        m_err = mk_core()->imesh_instance()->createVtx(pts[0], pts[1], pts[2], nodeHandle);
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
// function   : CreateElements
// Author     : Shengyong Cai
// Date       : Feb 16, 2011
// Description: create hexahedral elements by connecting 8 nodes which have
//              already been created by previous functions
//***************************************************************************//
int OneToOneSwept::CreateElements(vector<vector<Vertex> > &linkVertexList)
{
  //create the quadrilateral elements on the different layers
  //it is not necessary to create the quadrilateral elements on the different layers. Hex element can be created based on the eight nodes
  iMesh::Error m_err;
  vector<iBase_EntityHandle> mVolumeHandle(FaceList.size());
  // first decide orientation, based on the FaceList orientation, and first node above
  // take one face on source, and first node above in layer 1
  Vertex * v1 = FaceList[0].connect[0];
  Vertex * v2 = FaceList[0].connect[1];
  Vertex * v3 = FaceList[0].connect[2];
  // vertex 4 is on layer 1, above vertex 1
  int ind1 = v1->index;
  Vertex & v4 = linkVertexList[0][ind1];// this is the node above vertex 1 in layer 1

  Vector3D normal1 = (v2->xyz-v1->xyz)*(v3->xyz-v1->xyz);
  double vol6= normal1 % (v4.xyz-v1->xyz);
  std::cout << "Orientation is ";
  int skip = 0;
  if (vol6 < 0)
  {
    std::cout <<"negative\n";
    skip=4;
  }
  else
    std::cout <<"positive\n";

  vector<iBase_EntityHandle> connect(8);
  for (int m = 0; m < numLayers; m++)
  {
    if (m == 0) // first layer, above source
    {
      for (unsigned int i = 0; i < FaceList.size(); i++)
      {
        for (int k=0; k<4; k++)
        {
          connect[k+skip]=NodeList[(FaceList[i].connect[k])->index].gVertexHandle;
          connect[(k+skip+4)%8] = linkVertexList[m][(FaceList[i].connect[k])->index].gVertexHandle;
        }
        m_err = mk_core()->imesh_instance()->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, mVolumeHandle[i]);
        IBERRCHK(m_err, "Trouble create the hexahedral element.");
      }
      m_err = mk_core()->imesh_instance()->addEntArrToSet(&mVolumeHandle[0], FaceList.size(), volumeSet);
      IBERRCHK(m_err, "Trouble add an array of hexahedral elements to the entity set.");
    }
    else if (m == (numLayers - 1)) // the top layer, below target
    {
      for (unsigned int i = 0; i < FaceList.size(); i++)
      {
        for (int k=0; k<4; k++)
        {
          connect[k+skip]=linkVertexList[m-1][(FaceList[i].connect[k])->index].gVertexHandle;
          connect[(k+skip+4)%8] = TVertexList[(FaceList[i].connect[k])->index].gVertexHandle;
        }
        m_err = mk_core()->imesh_instance()->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, mVolumeHandle[i]);
        IBERRCHK(m_err, "Trouble create the hexahedral elements.");
      }
      m_err = mk_core()->imesh_instance()->addEntArrToSet(&mVolumeHandle[0], FaceList.size(), volumeSet);
      IBERRCHK(m_err, "Trouble add an array of hexahedral elements to the entity set.");
    }
    else // intermediate layers, m is between 0 and num_layers-2
    {
      for (unsigned int i = 0; i < FaceList.size(); i++)
      {
        for (int k=0; k<4; k++)
        {
          connect[k+skip]=linkVertexList[m-1][(FaceList[i].connect[k])->index].gVertexHandle;
          connect[(k+skip+4)%8] = linkVertexList[m][(FaceList[i].connect[k])->index].gVertexHandle;
        }
        m_err = mk_core()->imesh_instance()->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, mVolumeHandle[i]);
        IBERRCHK(m_err, "Trouble create the hexahedral elements.");
      }
      m_err = mk_core()->imesh_instance()->addEntArrToSet(&mVolumeHandle[0], FaceList.size(), volumeSet);
      IBERRCHK(m_err, "Trouble add an array of hexahedral elements to the entity set.");
    }
  }
  return 1;
}

#ifdef HAVE_MESQUITE
//target surface mesh smoothing by Mesquite
void OneToOneSwept::SurfMeshOptimization()
{
  MEntSelection::iterator mit = mentSelection.begin();
  ModelEnt *me = mit -> first;
  int irelPairIndex = me->iRelPairIndex();
  iGeom * igeom_inst = mk_core()->igeom_instance(me->iGeomIndex());
  //create a tag to attach the coordinates to nodes
  iBase_TagHandle mapped_tag = 0;
  iMesh::Error m_err = mk_core()->imesh_instance()->getTagHandle("MsqAltCoords", mapped_tag);
  if (m_err)
  {
    m_err = mk_core()->imesh_instance()->createTag("MsqAltCoords", 3, iBase_DOUBLE, mapped_tag);
    IBERRCHK(m_err, "Trouble create a tag.");

  }
  //attach the coordinates to the nodes
  double tag_data[3*TVertexList.size()];
  std::vector<iBase_EntityHandle> vertexHandle(TVertexList.size());
  for (unsigned int i = 0; i < NodeList.size();i++)
  {
    tag_data[3*i] = NodeList[i].xyz[0];
    tag_data[3*i+1] = NodeList[i].xyz[1];
    tag_data[3*i+2] = NodeList[i].xyz[2];
    vertexHandle[i] = TVertexList[i].gVertexHandle;
  }
  m_err = mk_core()->imesh_instance()->setDblArrData(&vertexHandle[0], NodeList.size(), mapped_tag, &tag_data[0]);
  IBERRCHK(m_err, "Trouble set an array of int data to nodes.");

  //get the mesh entityset for target surface
  iBase_EntitySetHandle surfSets;
  iRel::Error r_err = mk_core()->irel_pair(irelPairIndex)->getEntSetRelation(targetSurface, 0, surfSets);
  IBERRCHK(r_err, "Trouble get the mesh entity set for the target surface.");
  //call the MeshImprove class to smooth the target surface mesh by using Mesquite
  MeshImprove meshopt(mk_core(), false, true, true, true, igeom_inst);
  meshopt.SurfMeshImprove(targetSurface, surfSets, iBase_FACE);
}
#endif
}

