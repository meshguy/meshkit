#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/SimpleArray.hpp"
#ifdef HAVE_MESQUITE
#include "meshkit/MeshImprove.hpp"
#endif
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
    std::cout<<" source Global ID: " << index_id_src << " target Global ID: " << index_id_tar << "\n";

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
        std::cout<< " linking surface with global ID: " << index_id_link << " has " << edges.size() << " edges\n";
        assert((int)edges.size()==4);

        linkingSurfaces.push_back(surfs[i]);
      }
    }
    // now create a TFI mapping
    // this will trigger edge meshers for linking edges, and for target surface edges
    TFIMapping *tm = (TFIMapping*)mk_core()->construct_meshop("TFIMapping", linkingSurfaces);
    for (unsigned int i=0; i<linkingSurfaces.size(); i++)
    {
      linkingSurfaces[i]->sizing_function_index(sf->core_index());
    }
    mk_core()->insert_node(tm, (MeshOp*) this);
  }

  mk_core()->print_graph("AfterOneSetup.eps");

}

//---------------------------------------------------------------------------//
// execute function: generate the all-hex mesh through sweeping from source
// surface to target surface
void OneToOneSwept::execute_this()
{

  for (MEntSelection::iterator mit = mentSelection.begin(); mit
      != mentSelection.end(); mit++) {
    ModelEnt *me = mit -> first;
    if (me->get_meshed_state() >= COMPLETE_MESH)
      continue;
    numLayers = me->mesh_intervals();// maybe it will be decided later?
    //case 1: if numLayers = 1, then it is not necessary to create any vertices for linking surface, All the vertices have been created by source surface and target surface
    // bLayers will contain nodes in order, layer by layer
    std::vector<moab::EntityHandle> bLayers;
    BuildLateralBoundaryLayers(me, bLayers);
    //int sizeBLayer = bLayers.size()/(numLayers+1);

    TargetSurfProjection(bLayers);// the top layer is the last list of nodes

    //get the volume mesh entity set
    iRel::Error r_err =
        mk_core()->irel_pair(me->iRelPairIndex())->getEntSetRelation(
            me->geom_handle(), 0, volumeSet);
    IBERRCHK(r_err,
        "Trouble get the volume mesh entity set from the geometrical volume.");

    vector<vector<Vertex> > linkVertexList(numLayers - 1, vector<Vertex> (
        NodeList.size()));
    // this will compute the interior points distribution, in each layer
    ProjectInteriorLayers(bLayers, linkVertexList);

    CreateElements(linkVertexList);

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
  /*moab::Interface **/ mb = mk_core()->moab_instance();// we should get it also from
  // model entity, if we have multiple moab instances running around
  g_err = igeom_inst->getIntData(sourceSurface, geom_id_tag, index_id_src);
  IBERRCHK(g_err, "Trouble get the int data for source surface.");
  g_err = igeom_inst->getIntData(targetSurface, geom_id_tag, index_id_tar);
  IBERRCHK(g_err, "Trouble get the int data for target surface.");

  iBase_TagHandle taghandle=0;
  iMesh::Error m_err = mk_core()->imesh_instance()->createTag("source", 1, iBase_INTEGER, taghandle);
  IBERRCHK(m_err, "Trouble create the tag handle in the mesh instance.");

  MEntVector surfs;
  me->get_adjacencies(2, surfs);// all surfaces adjacent to the volume

  moab::Range quadsOnLateralSurfaces;

  ModelEnt * sME = NULL;
  me->get_adjacencies(2, surfs);// all surfaces adjacent to the volume
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

  numLoops = (int)group_sizes.size();

  std::cout << "Execute one to one swept ....\n";
  std::cout << "Nodes on boundary: " << bdyNodes.size()<< "\n";
  std::cout << "Number of loops: " << group_sizes.size() << "\n";
  std::cout << "Nodes on lateral surfaces: " << nodesOnLateralSurfaces.size() << "\n";
  std::cout << "Quads on lateral surfaces: " << quadsOnLateralSurfaces.size() << "\n";

  if (bdyNodes.size()*(numLayers+1)!=nodesOnLateralSurfaces.size())
  {
    std::cout << "Major problem: number of total nodes on boundary: "  <<nodesOnLateralSurfaces.size()<<"\n";
    std::cout << " number of nodes in one layer on the boundary:" << bdyNodes.size() << "\n";
    std::cout << " we expect bdyNodes.size()*(numLayers+1) == nodesOnLateralSurfaces.size()" << bdyNodes.size()*(numLayers+1) << " != " << nodesOnLateralSurfaces.size() << "\n";
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
  for (moab::Range::iterator it = sourceNodes.begin(); it!=sourceNodes.end(); it++, i++)
  {
    moab::EntityHandle node=*it;
    NodeList[i].gVertexHandle = (iBase_EntityHandle)sourceNodes[i];
    NodeList[i].index = i;

    rval = mb->get_coords(&node, 1, NodeList[i].xyzCoords);
    MBERRCHK(rval, mb);

    NodeList[i].onBoundary = false;
    NodeList[i].onCorner = false;

    //set the int data to the entity
    m_err = mk_core()->imesh_instance()->setIntData(NodeList[i].gVertexHandle, taghandle, NodeList[i].index);
    IBERRCHK(m_err, "Trouble set the int value for nodes in the mesh instance.");
  }
  // now loop over the boundary and mark the boundary nodes as such
  for (unsigned int j=0; j<bdyNodes.size(); j++)
  {
    iBase_EntityHandle node = (iBase_EntityHandle)bdyNodes[j];
    int index;
    m_err = mk_core()->imesh_instance()->getIntData(node, taghandle, index);
    IBERRCHK(m_err, "Trouble set the int value for nodes in the mesh instance.");
    NodeList[index].onBoundary = true;// some could be corners, but it is not that important
  }

  // process faces on the source
  FaceList.resize(quads.size());// quads is actually a range....

  for (unsigned int i=0; i < quads.size(); i++)
  {
    iBase_EntityHandle currentFace = (iBase_EntityHandle)quads[i];//
    FaceList[i].gFaceHandle = currentFace;
    FaceList[i].index = i;

    //get the nodes on the face elements
    std::vector<iBase_EntityHandle> Nodes;
    m_err = mk_core()->imesh_instance()->getEntAdj(currentFace, iBase_VERTEX, Nodes);
    IBERRCHK(m_err, "Trouble get the adjacent nodes from mesh face entity.");

    FaceList[i].connect.resize(Nodes.size());
    for (unsigned int j=0; j < Nodes.size(); j++)
    {
      int tmpIndex=-1;
      m_err = mk_core()->imesh_instance()->getIntData(Nodes[j], taghandle, tmpIndex);
      IBERRCHK(m_err, "Trouble get the int data from node handle.");

      FaceList[i].connect[j] = &NodeList[tmpIndex];// why not store tmp Index directly?
    }

    m_err = mk_core()->imesh_instance()->setIntData(currentFace, taghandle, i);
    IBERRCHK(m_err, "Trouble set the int data for quad mesh on the source surface.");
  }

  std::copy(bdyNodes.begin(), bdyNodes.end(), std::back_inserter(layers) );

  int deft = -1; // not set; 0 mean layer 1; maybe; we should maybe not use it?
  rval = mb->tag_get_handle( "mark", 1, moab::MB_TYPE_INTEGER, markTag, moab::MB_TAG_CREAT|moab::MB_TAG_DENSE, &deft);
  MBERRCHK(rval, mb);


  int sizeBLayer = (int)bdyNodes.size();
  int val=0;// layer level ; we will use this to verify the layers of boundary nodes
  for (i=0; i<sizeBLayer; i++)
  {
    rval = mb->tag_set_data(markTag, &(bdyNodes[i]), 1, &val);
    MBERRCHK(rval, mb);
  }

  for (int layer=1; layer <= numLayers; layer++) // layer with index 0 is at the bottom/source
  {
    int start_index = 0;//this is the start index in group
    for (unsigned int k=0; k<group_sizes.size(); k++)
    {
      // first group starts at index 0, the rest are accumulated
      if (k>0)
        start_index+=group_sizes[k-1];
      moab::EntityHandle node1 =layers[ (layer-1)*sizeBLayer +start_index];
      moab::EntityHandle node2 =layers[ (layer-1)*sizeBLayer +start_index+1];
      moab::EntityHandle node3, node4 ;
      /*
       *     node3 -- node4 -->  (next layer)
       *       |        |
       *     node1 -- node2 ---> (guide layer)
       */
      rval = NodeAbove(mb, node1, node2, quadsOnLateralSurfaces, node3, node4);
      MBERRCHK(rval, mb);
      // check the mark tag!
      rval = mb->tag_get_data(markTag, &node3, 1, &val);
      MBERRCHK(rval, mb);
      if (val!=-1)
        MBERRCHK(moab::MB_FAILURE, mb);
      rval = mb->tag_get_data(markTag, &node4, 1, &val);
      MBERRCHK(rval, mb);
      if (val!=-1)
        MBERRCHK(moab::MB_FAILURE, mb);
      layers.push_back( node3);//start of a new layer
      layers.push_back( node4);//node above node 2
      rval = mb->tag_set_data(markTag, &node3, 1, &layer);// layer level>0
      MBERRCHK(rval, mb);
      rval = mb->tag_set_data(markTag, &node4, 1, &layer);// layer level>0
      MBERRCHK(rval, mb);
      // we have at least 2 nodes in every loop, otherwise we are in deep trouble
      //start
      int currentIndex=start_index+2;
      while (currentIndex < start_index+group_sizes[k])
      {
        node1 = node2;
        node2 = layers[ (layer-1)*sizeBLayer + currentIndex];// actually, the previous layer
        node3 = node4;
        rval = FourthNodeInQuad(mb, node1, node2, node3, quadsOnLateralSurfaces, node4);
        MBERRCHK(rval, mb);
        rval = mb->tag_get_data(markTag, &node4, 1, &val);
        MBERRCHK(rval, mb);
        if (val!=-1)
          MBERRCHK(moab::MB_FAILURE, mb);
        layers.push_back( node4);
        rval = mb->tag_set_data(markTag, &node4, 1, &layer);// layer level>0
        MBERRCHK(rval, mb);
        currentIndex++;
      }
    }// end group
  }// end layer
  for (int lay=0; lay<=numLayers; lay++)
  {
    std::cout << "layer : " << lay <<"\n";
    for (int i =0; i<sizeBLayer; i++)
    {
      std::cout << " " << layers[lay*sizeBLayer + i];
      if (i%10 == 9)
        std::cout<<"\n";
    }
    std::cout<<"\n";
  }
}

moab::ErrorCode OneToOneSwept::NodeAbove(moab::Interface * mb, moab::EntityHandle node1,
    moab::EntityHandle node2, moab::Range & latQuads, moab::EntityHandle & node3,
    moab::EntityHandle & node4)
{
  // look for node above node 1 in an unused quad
  moab::EntityHandle nds[2] = {node1, node2};
  // find all quads connected to these 2 nodes; find one in the range that is not used
  moab::Range adj_entities;
  moab::ErrorCode rval = mb->get_adjacencies(nds, 2, 2, false,adj_entities);
  MBERRCHK(rval, mb);
  std::vector<int> tagVals(adj_entities.size());
  rval = mb->tag_get_data(markTag, adj_entities, &(tagVals[0]));
  MBERRCHK(rval, mb);
  // default is -1
  moab::EntityHandle nextQuad=0;// null
  for (unsigned int i=0; i<tagVals.size(); i++)
  {
    if ((tagVals[i]==-1 )&&  (latQuads.find(adj_entities[i])!=latQuads.end()) )
    {
      nextQuad = adj_entities[i];
      break;
    }
  }
  if (0==nextQuad)
    MBERRCHK(moab::MB_FAILURE, mb);

  // decide on which side are nodes 1 and 2, then get the opposite side
  const moab::EntityHandle * conn4 = NULL;
  int nnodes;
  rval = mb->get_connectivity(nextQuad, conn4, nnodes);
  MBERRCHK(rval, mb);
  int index1=-1;
  for (index1=0; index1<4; index1++)
    if (node1==conn4[index1])
      break;
  if (4==index1)
    MBERRCHK(moab::MB_FAILURE, mb);
  if (node2== conn4[(index1+1)%4]) // quad is oriented node1, node2, node4, node3
  {
    node3 = conn4[(index1-1)%4];
  }
  else if (node2 == conn4[(index1-1)%4]) // quad is oriented node2, node1, node3, node4
  {
    node3 = conn4[(index1+1)%4];
  }
  else
    MBERRCHK(moab::MB_FAILURE, mb);// something is really wrong
  node4 = conn4[(index1+2)%4];

  // mark the quad to not use it again
  int val=0; // maybe it should be the layer
  rval = mb->tag_set_data(markTag, &nextQuad, 1, &val);
  MBERRCHK(rval, mb);

  return moab::MB_SUCCESS;
}

moab::ErrorCode OneToOneSwept::FourthNodeInQuad(moab::Interface * mb, moab::EntityHandle node1,
      moab::EntityHandle node2, moab::EntityHandle node3,
      moab::Range & latQuads, moab::EntityHandle & node4)
{
  // look for the fourth node
  moab::EntityHandle nds[3] = {node1, node2, node3};
  // find all quads connected to these 3 nodes; find one in the range that is not used
  moab::Range adj_entities;
  moab::ErrorCode rval = mb->get_adjacencies(nds, 3, 2, false, adj_entities);
  MBERRCHK(rval, mb);

  std::vector<int> tagVals(adj_entities.size());
  rval = mb->tag_get_data(markTag, adj_entities, &(tagVals[0]));

  MBERRCHK(rval, mb);
  // default is -1
  moab::EntityHandle nextQuad=0;// null
  for (unsigned int i=0; i<tagVals.size(); i++)
  {
    if ((tagVals[i]==-1 )&&  (latQuads.find(adj_entities[i])!=latQuads.end()) )
    {
      nextQuad = adj_entities[i];
      break;
    }
  }
  if (0==nextQuad)
    MBERRCHK(moab::MB_FAILURE, mb);

  // decide on which side are nodes 1 and 2, then get the opposite side
  const moab::EntityHandle * conn4 = NULL;
  int nnodes;
  rval = mb->get_connectivity(nextQuad, conn4, nnodes);
  MBERRCHK(rval, mb);

  node4=0;
  for (int index=0; index<4; index++)
  {
    moab::EntityHandle c4 = conn4[index];
    if (node1!=c4 && node2 != c4 && node3!=c4)
    {
      node4 = c4;
      break;
    }
  }
  if (0==node4)
    MBERRCHK(moab::MB_FAILURE, mb);

  // mark the quad to not use it again
  int val=0; // maybe it should be the layer
  rval = mb->tag_set_data(markTag, &nextQuad, 1, &val);
  MBERRCHK(rval, mb);

  return moab::MB_SUCCESS;

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
  //moab::Interface *mb = mk_core()->moab_instance();
  int irelIndx = target_surf->iRelPairIndex();
  iGeom * igeom_inst = mk_core()->igeom_instance( target_surf->iGeomIndex());
  // we could have got it from model tag, too
  // something along this:
  // err = igeom_instance(geom_index)->getData(*eit, iGeomModelTags[geom_index], &this_me);

  TVertexList.resize(NodeList.size());
  // make everything not on boundary, first; the true boundary flag will be set later
  for (unsigned int i = 0; i < TVertexList.size(); i++)
  {
    TVertexList[i].onBoundary = false;
    TVertexList[i].onCorner = false; // it is not really needed!!! boundary is enough for our
    // purposes
  }

  iBase_TagHandle taghandle_tar=0;
  m_err = mk_core()->imesh_instance()->createTag("TargetMesh", 1, iBase_INTEGER, taghandle_tar);
  IBERRCHK(m_err, "Trouble create the taghandle for the target mesh.");

  iBase_TagHandle taghandle=0;
  m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
  IBERRCHK(m_err, "Trouble get tag handle of the source surface.");

  // we know that the first layer of nodes (on the source)
  int sizeBLayer = (int)bLayers.size()/(numLayers+1); // this could be a member data...
  // the first 0 , 1, ..., sizeBlayer - 1 are on bottom, were bdyNodes before
  //
  for (int i = 0; i<sizeBLayer; i++)
  {
    iBase_EntityHandle sBNode = (iBase_EntityHandle)bLayers[i];
    int index;
    m_err = mk_core()->imesh_instance()->getIntData(sBNode, taghandle, index);
    IBERRCHK(m_err, "Trouble get the int value for nodes in the mesh instance.");
    TVertexList[index].onBoundary = true;// some could be corners, but it is not that important
    moab::EntityHandle topNode = bLayers[numLayers*sizeBLayer + i];// node right above
    iBase_EntityHandle tBNode = (iBase_EntityHandle)topNode;// it is right above node i in loop
    TVertexList[index].gVertexHandle = tBNode;
    // now get the coordinates of the tBNode (node on boundary of target)
    rval = mb->get_coords(&topNode, 1, TVertexList[index].xyzCoords);
    MBERRCHK(rval, mb);

  }

  iBase_EntitySetHandle entityset;  //this entityset is for storing the inner nodes on the target surface
  vector<iBase_EntityHandle>  newNodehandle;


  r_err = mk_core()->irel_pair(irelIndx)->getEntSetRelation(targetSurface, 0, entityset);
  if (r_err) //there is no entityset associated with targetSurface; this should be an error
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
  for (unsigned int i = 0; i < NodeList.size(); i++) {
    if (NodeList[i].onBoundary /*|| NodeList[i].onCorner*/) {
      num_pts++;
      for (int j = 0; j < 3; j++) {
        sPtsCenter[j] = sPtsCenter[j] + NodeList[i].xyzCoords[j];
        tPtsCenter[j] = tPtsCenter[j] + TVertexList[i].xyzCoords[j];
      }
      sBndNodes.resize(num_pts);
      tBndNodes.resize(num_pts);
      sBndNodes[num_pts - 1].resize(3);
      tBndNodes[num_pts - 1].resize(3);
      for (int j = 0; j < 3; j++) {
        sBndNodes[num_pts - 1][j] = NodeList[i].xyzCoords[j];
        tBndNodes[num_pts - 1][j] = TVertexList[i].xyzCoords[j];
      }
    }
  }
  if (sizeBLayer!=num_pts)
    MBERRCHK(moab::MB_FAILURE, mb);

  for (int i = 0; i < 3; i++) {
    sPtsCenter[i] = sPtsCenter[i] / num_pts;
    tPtsCenter[i] = tPtsCenter[i] / num_pts;
  }
  //done with the barycenter calculation

  double tmpMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double bMatrix[3][3]   = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  //transform the coordinates
  std::vector<std::vector<double> > sBNodes(num_pts, std::vector<double>(3)),
      tBNodes(num_pts, std::vector<double>(3));
  for (int i = 0; i < num_pts; i++) {
    //temporarily store the boundary nodes' coordinates on the source surface and target surface
    for (int j = 0; j < 3; j++) {
      sBNodes[i][j] = sBndNodes[i][j];
      tBNodes[i][j] = tBndNodes[i][j];
    }

    //transform the boundary nodes
    for (int j = 0; j < 3; j++) {
      sBNodes[i][j] = sBNodes[i][j] - 2 * sPtsCenter[j] + tPtsCenter[j];
      tBNodes[i][j] = tBNodes[i][j] - sPtsCenter[j];
    }
  }

  //calculate the transformation matrix: transform the nodes on the source surface to the target surface in the physical space
  for (int i = 0; i < num_pts; i++) {
    //3 row entries in the temporary matrix
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        tmpMatrix[j][k] = tmpMatrix[j][k] + sBNodes[i][j] * sBNodes[i][k];

    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        bMatrix[j][k] = bMatrix[j][k] + tBNodes[i][j] * sBNodes[i][k];
  }

  //first determine the affine mapping matrix is singular or not
  double detValue =
        tmpMatrix[0][2] * tmpMatrix[1][1] * tmpMatrix[2][0]
      - tmpMatrix[0][1] * tmpMatrix[1][2] * tmpMatrix[2][0]
      - tmpMatrix[0][2] * tmpMatrix[1][0] * tmpMatrix[2][1]
      + tmpMatrix[0][0] * tmpMatrix[1][2] * tmpMatrix[2][1]
      + tmpMatrix[0][1] * tmpMatrix[1][0] * tmpMatrix[2][2]
      - tmpMatrix[0][0] * tmpMatrix[1][1] * tmpMatrix[2][2];
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
      transMatrix[i][j] = InvMatrix[j][0] * bMatrix[i][0]
       + InvMatrix[j][1] * bMatrix[i][1] + InvMatrix[j][2] * bMatrix[i][2];
  //finish calculating the interior nodes' location


  //calculate the inner nodes on the target surface
  for (unsigned int i = 0; i < NodeList.size(); i++) {
    if (!(NodeList[i].onBoundary /*|| NodeList[i].onCorner*/)) {
      double xyz[3];
      for (int j = 0; j < 3; j++)
        xyz[j] = transMatrix[j][0] * (NodeList[i].xyzCoords[0] - 2
            * sPtsCenter[0] + tPtsCenter[0]) + transMatrix[j][1]
            * (NodeList[i].xyzCoords[1] - 2 * sPtsCenter[1] + tPtsCenter[1])
            + transMatrix[j][2] * (NodeList[i].xyzCoords[2] - 2 * sPtsCenter[2]
                + tPtsCenter[2]) + sPtsCenter[j];
      //calculate the closest point on the target surface with respect to the projected point
      //std::cout << "index = " << i << "\tx = " << xyz[0] << "\ty = " << xyz[1] << "\tz = " << xyz[2] << std::endl;

      iGeom::Error g_err = igeom_inst->getEntClosestPt(targetSurface, xyz[0],
          xyz[1], xyz[2], TVertexList[i].xyzCoords[0],
          TVertexList[i].xyzCoords[1], TVertexList[i].xyzCoords[2]);
      IBERRCHK(g_err, "Trouble get the closest point on the targets surface.");

      //create the node entities
      iMesh::Error m_err = mk_core()->imesh_instance()->createVtx(
          TVertexList[i].xyzCoords[0], TVertexList[i].xyzCoords[1],
          TVertexList[i].xyzCoords[2], TVertexList[i].gVertexHandle);
      IBERRCHK(m_err, "Trouble create the node entity.");

      TVertexList[i].onBoundary = false;
      TVertexList[i].onCorner = false;

      newNodehandle.push_back(TVertexList[i].gVertexHandle);
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
  moab::EntityHandle node1 = bdyTNodes[1], node2=bdyTNodes[2];
  int index=-1;
  for (index=0; index<sizeBLayer; index++)
  {
    if (bLayers[numLayers*sizeBLayer+index]==node1)
      break;
  }
  if (index==sizeBLayer)
    MBERRCHK(moab::MB_FAILURE, mb);

  int prevIndex= (index-1) %sizeBLayer;
  int nextIndex= (index+1) %sizeBLayer;
  if (bLayers[numLayers*sizeBLayer+prevIndex] == node2)
    sense_out = -1;
  else if (bLayers[numLayers*sizeBLayer+nextIndex] == node2)
    sense_out = 1;
  else
    MBERRCHK(moab::MB_FAILURE, mb);// serious error in the logic, can't decide orientation


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
  IBERRCHK(m_err, "Trouble create the quadrilateral mesh.");

  //add the inner face elements to the entityset
  m_err = mk_core()->imesh_instance()->addEntArrToSet(&mFaceHandle[0], FaceList.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of quadrilateral entities to the entity set.");

  mk_core()->save_mesh("BeforeHex.vtk");
#ifdef HAVE_MESQUITE

  // this test is here because mesquite assumes that the first igeom instance is used
  // so, we will skip mesquite until we fix this
if (target_surf->iGeomIndex() == 0)
  SurfMeshOptimization();

#endif


  return 1;
}
// use similar code to TargetSurfProjection, but do not project on surface...
int OneToOneSwept::ProjectInteriorLayers(std::vector<moab::EntityHandle> & boundLayers,
    vector<vector <Vertex> > & linkVertexList)
{

  iMesh::Error m_err;
  Point3D sPtsCenter = {0, 0, 0}, tPtsCenter={0, 0, 0};
  std::vector<Point3D> PtsCenter(numLayers-1);
  std::vector<Point3D> sBoundaryNodes(0), tBoundaryNodes(0);


  // this could be saved, the tag handle for source...
  iBase_TagHandle taghandle=0;
  m_err = mk_core()->imesh_instance()->getTagHandle("source", taghandle);
  IBERRCHK(m_err, "Trouble getting the source tag handle");

  int sizeBLayer = boundLayers.size()/(numLayers+1);
  std::vector<vector <Point3D> > iBoundaryNodes(numLayers-1, std::vector<Point3D>(sizeBLayer));
  // we should use moab::CartVect, or MeshKit::Matrix<3,1> = MeshKit::Vector<3>
  // Point3D does not know how to do anything
  std::vector<Point3D> latCoords;
  latCoords.resize(boundLayers.size()-2*sizeBLayer);// to store the coordinates of the lateral points
  // (exclude the source and target, they are already in NodeList and TVertexList)...

  //moab::Interface *mb = mk_core()->moab_instance();
  moab::ErrorCode rval = mb->get_coords(&(boundLayers[sizeBLayer]), latCoords.size(),
      &(latCoords[0].px)) ; //these are the nodes for lateral sides
  MBERRCHK(rval, mb);

  //calculate the center coordinates
  for (int i=0; i< numLayers-1; i++)
  {
    PtsCenter[i].px = 0;
    PtsCenter[i].py = 0;
    PtsCenter[i].pz = 0;
  }
  int numPts = sizeBLayer;// a lot of repetition ; do we rfeally need
  sBoundaryNodes.resize(numPts);// another useless copy
  tBoundaryNodes.resize(numPts);// another useless copy
  for (int j=0; j < sizeBLayer; j++)
  {
    iBase_EntityHandle sBNode = (iBase_EntityHandle)boundLayers[j];
    int i;
    // this is the index in NodeList...
    m_err = mk_core()->imesh_instance()->getIntData(sBNode, taghandle, i);
    IBERRCHK(m_err, "Trouble get the int value for nodes in the mesh instance.");

    sPtsCenter.px = sPtsCenter.px + NodeList[i].xyzCoords[0];
    sPtsCenter.py = sPtsCenter.py + NodeList[i].xyzCoords[1];
    sPtsCenter.pz = sPtsCenter.pz + NodeList[i].xyzCoords[2];

    tPtsCenter.px = tPtsCenter.px + TVertexList[i].xyzCoords[0];
    tPtsCenter.py = tPtsCenter.py + TVertexList[i].xyzCoords[1];
    tPtsCenter.pz = tPtsCenter.pz + TVertexList[i].xyzCoords[2];

    sBoundaryNodes[j].px = NodeList[i].xyzCoords[0];
    sBoundaryNodes[j].py = NodeList[i].xyzCoords[1];
    sBoundaryNodes[j].pz = NodeList[i].xyzCoords[2];

    tBoundaryNodes[j].px = TVertexList[i].xyzCoords[0];
    tBoundaryNodes[j].py = TVertexList[i].xyzCoords[1];
    tBoundaryNodes[j].pz = TVertexList[i].xyzCoords[2];
    for (int lay=0; lay<numLayers-1; lay++)
    {
      // interior layers
      int indexNodeOnSide = lay*sizeBLayer+j;// all p3d will be above index i
      Point3D & p3d = latCoords[indexNodeOnSide];
      PtsCenter[lay].px = PtsCenter[lay].px + p3d.px;
      PtsCenter[lay].py = PtsCenter[lay].py + p3d.py;
      PtsCenter[lay].pz = PtsCenter[lay].pz + p3d.pz;

      iBoundaryNodes[lay][j].px = p3d.px;
      iBoundaryNodes[lay][j].py = p3d.py;
      iBoundaryNodes[lay][j].pz = p3d.pz;

      // another copy that is redundant
      linkVertexList[lay][i].xyzCoords[0] = p3d.px;
      linkVertexList[lay][i].xyzCoords[1] = p3d.py;
      linkVertexList[lay][i].xyzCoords[2] = p3d.pz;
      linkVertexList[lay][i].gVertexHandle = (iBase_EntityHandle)boundLayers[indexNodeOnSide+sizeBLayer];

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

  std::vector<Point3D> sBNodes(numPts); //boundary nodes on the source surface
  std::vector<Point3D> tBNodes(numPts); //boundary nodes on the target surface
  std::vector<std::vector <Point3D> > isBNodes(numLayers-1, std::vector<Point3D>(numPts)),
       itBNodes(numLayers-1, vector<Point3D>(numPts));  //boundary nodes on the different layer

  //loop over different layers
  for (int i=0; i < numLayers-1; i++)
  {
    double sA[3][3], tA[3][3];
    double stransMatrix[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, ttransMatrix[3][3]={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double sInvMatrix[3][3], tInvMatrix[3][3];
    double sb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, tb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};


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
// Description: function for obtaining the x,y,z coordinates from parametric coordinates
//***************************************************************************//
/*
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
}*/

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
    if (m_err){
	m_err = mk_core()->imesh_instance()->createTag("MsqAltCoords", 3, iBase_DOUBLE, mapped_tag);
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
    iRel::Error r_err = mk_core()->irel_pair(irelPairIndex)->getEntSetRelation(targetSurface, 0, surfSets);
    IBERRCHK(r_err, "Trouble get the mesh entity set for the target surface.");
    //call the MeshImprove class to smooth the target surface mesh by using Mesquite
    MeshImprove meshopt(mk_core(), false, true, true, true, igeom_inst);
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
#endif

}

