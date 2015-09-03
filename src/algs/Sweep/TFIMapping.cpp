#include "meshkit/TFIMapping.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp"
#include "EquipotentialSmooth.hpp"
#ifdef HAVE_MESQUITE
#include "meshkit/MeshImprove.hpp"
#endif

#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <algorithm>

const double EPS = 1.0e-6;

namespace MeshKit {

//---------------------------------------------------------------------------//
//Entity Type initialization for TFIMapping meshing
moab::EntityType TFIMapping_tps[] = { moab::MBVERTEX, moab::MBEDGE, moab::MBQUAD, moab::MBMAXTYPE };
const moab::EntityType* TFIMapping::output_types()
{
  return TFIMapping_tps;
}

//---------------------------------------------------------------------------//
// construction function for TFIMapping class
TFIMapping::TFIMapping(MKCore *mk_core, const MEntVector &me_vec) :
  MeshScheme(mk_core, me_vec)
{
  _shapeImprove=false;
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

  /* the only things we need to make sure :
   1) there are 4 edges, exactly: this is removed
   2) the opposite edges have the same meshcount
   - if some of the edges are meshed, we need to mesh the opposite edge with
   correct meshcount
   - if 2 opposite edges are meshed, verify the mesh count
   */
  // get iGeom instance from the first ment selection
  if (mentSelection.empty())
    return;

  //loop over the surfaces
  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit -> first;
    int dimME = me->dimension();
    if (dimME != 2)
      ECERRCHK(MK_FAILURE, "bad input for TFI Mapping, we can only mesh surfaces");
    //first check whether the surface is meshed or not
    if (me->get_meshed_state() >= COMPLETE_MESH)
      continue;
    int size_index = me->sizing_function_index();
    int mesh_count_surface = -1;
    if (size_index >= 0)
      mesh_count_surface = me->mk_core()->sizing_function(size_index)->intervals();

    // get the boundary loop of edges
    MEntVector boundEdges;
    std::vector<int> senses, group_sizes;
    me->ModelEnt::boundary(1, boundEdges, &senses, &group_sizes);
    //remove this constraint in case of the side-cylinder surface
    //if (boundEdges.size() != 4)
    //  ECERRCHK(MK_FAILURE, "bad input for TFI Mapping, we can only mesh surfaces with 4 edges");

    if (boundEdges.size()<4){
        ModelEnt *oppEdges[2] = {boundEdges[0], boundEdges[1]};
        MEntVector edgesToMesh;

        int mesh_count = mesh_count_surface;
        bool force = false;
        for (int j = 0; j < 2; j++){
            if (oppEdges[j]->get_meshed_state() >= COMPLETE_MESH){
                std::vector<moab::EntityHandle> medges;
                oppEdges[j]->get_mesh(1, medges, true);
                mesh_count = (int) medges.size();
                force = true;
            }
            else
            {
                int indexS = oppEdges[j]->sizing_function_index();
                if (indexS >= 0){
                    SizingFunction * sfe = mk_core()->sizing_function(indexS);
                    if (!force)
                    {
                        if (sfe->intervals() > 0)
                            mesh_count = sfe->intervals();
                        else if (sfe->size() > 0)
                            mesh_count = oppEdges[j]->measure() / sfe->size();
                        if (mesh_count % 2 && oppEdges[j]->constrain_even())
                            ++mesh_count;
                    }
                }
                edgesToMesh.push_back(oppEdges[j]);
            }
        }
		if (edgesToMesh.size() > 0)
      	{
		    EdgeMesher * em = (EdgeMesher*) me->mk_core()->construct_meshop("EdgeMesher", edgesToMesh);
		    if (mesh_count < 0)
		    {
		      std::cout << "mesh count not set properly on opposite edges, set it to 10\n";
		      mesh_count = 10; // 4 is a nice number, used in the default edge mesher;
		      // but I like 10 more
		    }

		    for (unsigned int j = 0; j < edgesToMesh.size(); j++)
		    {
                      int edgeMeshCount = 0;
                      int edgeSfIndex =
                          edgesToMesh[j]->sizing_function_index();
                      if (edgeSfIndex >= 0)
                      {
                        SizingFunction* edgeSf =
                            mk_core()->sizing_function(edgeSfIndex);
                        edgeMeshCount = edgeSf->intervals();
		      }
                      if (mesh_count != edgeMeshCount)
                      {
                        edgesToMesh[j]->mesh_intervals(mesh_count);
                      }
                      if (force)
                      {
                        // the opposite edge is already meshed, so the number
                        // of intervals is a hard constraint for this edge
                        edgesToMesh[j]->interval_firmness(HARD);
                      }
		      edgesToMesh[j]->add_meshop(em);
		    }
		    mk_core()->insert_node(em, (GraphNode*)this,
                        mk_core()->root_node());
      	}
    }
	else{

		// mesh edge 0 and 2 together, and 1 and 3 together (same mesh count)
		// look at all settings, to decide proper mesh count
		for (int k = 0; k <= 1; k++)
		{
		  // treat first edges 0 and 2, then 1 and 3
		  ModelEnt * oppEdges[2] = { boundEdges[k], boundEdges[k + 2] };
		  MEntVector edgesToMesh;// edges that are not meshed yet
		  //if one of them is meshed and the other not, use the same mesh count

		  // take the maximum of the proposed mesh counts, either from sizing function, or mesh intervals
		  int mesh_count = mesh_count_surface; // could be -1, still
		  bool force = false;
		  for (int j = 0; j < 2; j++)
		  {
		    if (oppEdges[j]->get_meshed_state() >= COMPLETE_MESH)
		    {
		      // in this case, force the other edge to have the same mesh count, do not take it from surface
		      std::vector<moab::EntityHandle> medges;
		      oppEdges[j]->get_mesh(1, medges, true);
		      mesh_count = (int) medges.size();
		      force = true;
		    }
		    else
		    {
		      int indexS = oppEdges[j]->sizing_function_index();
		      if (indexS >= 0)
		      {
		        SizingFunction * sfe = mk_core()->sizing_function(indexS);
                        if (!force)
                        {
                          // if a sizing function was set on an edge, use
                          // that rather than a mesh count from the surface
                          if (sfe->intervals() > 0)
                            mesh_count = sfe->intervals();
                          else if (sfe->size() > 0)
                            mesh_count = oppEdges[j]->measure() /  sfe->size();
                          if (mesh_count % 2 && oppEdges[j]->constrain_even())
                            ++mesh_count;
                        }
		      }
		      // push it to the list if it is not setup to another mesh op (edge mesher) already
		      //if (oppEdges[j]->is_meshops_list_empty())// it will create an EdgeMesher later
                      if ((j == 0 || (oppEdges[j] != oppEdges[0])) &&
                          oppEdges[j]->is_meshops_list_empty())
                      {
		        edgesToMesh.push_back(oppEdges[j]);
                      }
		    }
		  }
		  // decide on a mesh count now, if edgesToMesh.size()>0
		  if (edgesToMesh.size() > 0)
		  {

		    EdgeMesher * em = (EdgeMesher*) me->mk_core()->construct_meshop("EdgeMesher", edgesToMesh);
		    if (mesh_count < 0)
		    {
		      std::cout << "mesh count not set properly on opposite edges, set it to 10\n";
		      mesh_count = 10; // 4 is a nice number, used in the default edge mesher;
		      // but I like 10 more
		    }

		    for (unsigned int j = 0; j < edgesToMesh.size(); j++)
		    {
                      int edgeMeshCount = 0;
                      int edgeSfIndex =
                          edgesToMesh[j]->sizing_function_index();
                      if (edgeSfIndex >= 0)
                      {
                        SizingFunction* edgeSf =
                            mk_core()->sizing_function(edgeSfIndex);
                        edgeMeshCount = edgeSf->intervals();
		      }
                      if (mesh_count != edgeMeshCount)
                      {
                        edgesToMesh[j]->mesh_intervals(mesh_count);
                      }
                      if (force)
                      {
                        // the opposite edge is already meshed, so the number
                        // of intervals is a hard constraint for this edge
                        edgesToMesh[j]->interval_firmness(HARD);
                      }
		      edgesToMesh[j]->add_meshop(em);
		    }
		    mk_core()->insert_node(em, (GraphNode*)this,
                        mk_core()->root_node());
		  }
		} // end loop over pair of opposite edges
	}
  }// end loop over surfaces

  ensure_facet_dependencies(false);

  mk_core()->print_graph("AfterTFISetup.eps");
}

//---------------------------------------------------------------------------//
// execute function: generate the all-quad mesh through the TFI mapping
void TFIMapping::execute_this()
{

  //loop over the surfaces
  for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit -> first;
    //first check whether the surface is meshed or not
    if (me->get_meshed_state() >= COMPLETE_MESH){
		
#ifdef HAVE_MESQUITE
	iBase_EntitySetHandle entityset;
    iRel::Error r_err = mk_core()->irel_pair(me->iRelPairIndex())->getEntSetRelation(me->geom_handle(), 0, entityset);
	IBERRCHK(r_err, "Trouble get the entityset w.r.t a surface!");
    MeshImprove shapesmooth(mk_core(), false, false, true, false, mk_core()->igeom_instance(me->iGeomIndex()));
    shapesmooth.SurfMeshImprove(me->geom_handle(), entityset, iBase_FACE);
#endif     
	
	 continue;
	}
	
	MEntVector boundEdges;
    std::vector<int> senses, group_sizes;
    me->ModelEnt::boundary(1, boundEdges, &senses, &group_sizes);
    set<ModelEnt*> distinctBoundEdges;
    distinctBoundEdges.insert(boundEdges.begin(), boundEdges.end());
    if (distinctBoundEdges.size() == 4)
    	SurfMapping(me);
    else
	cylinderSurfMapping(me);

    //ok, we are done, commit to ME
    me->commit_mesh(mit->second, COMPLETE_MESH);
  }

}

int TFIMapping::cylinderSurfMapping(ModelEnt *ent)
{
  int irelPairIndex = ent->iRelPairIndex();

  // determine whether there is an edge along the linking surface
  MEntVector allBoundEdges;
  std::vector<int> allSenses, group_sizes;
  ent->ModelEnt::boundary(1, allBoundEdges, &allSenses, &group_sizes);
  map<ModelEnt*, int> boundEdgeCount;
  for (unsigned int i = 0; i < allBoundEdges.size(); ++i)
  {
    boundEdgeCount[allBoundEdges[i]] = 0;
  }
  for (unsigned int i = 0; i < allBoundEdges.size(); ++i)
  {
    ++boundEdgeCount[allBoundEdges[i]];
  }
  MEntVector boundEdges;
  std::vector<int> boundEdgeSenses;
  ModelEnt* linkSurfEdge = NULL;
  for (unsigned int i = 0; i < allBoundEdges.size(); ++i)
  {
    if (boundEdgeCount[allBoundEdges[i]] == 1)
    {
      boundEdges.push_back(allBoundEdges[i]);
      boundEdgeSenses.push_back(allSenses[i]);
    }
    else
    {
      linkSurfEdge = allBoundEdges[i];
    }
  }

  if (boundEdges.size() != 2)
  {
    ECERRCHK(MK_FAILURE, "Cylinder TFIMapping does not have exactly two distinct bounding edges");
  }

  std::cout << "TFIMapping on cylinder\n";

  std::vector<moab::EntityHandle> nList;
  std::vector<iBase_EntityHandle> List_i, List_ii;

  // get nodes of bounding edge 0 into List_i
  boundEdges[0]->get_mesh(0, nList, true);
  unsigned int ix = 0;
  int size_i = (int)nList.size()-1;
  for (ix = 0; ix < nList.size(); ix++)
  {
    List_i.push_back((iBase_EntityHandle) nList[ix]);
  }
  // we reverse the first boundary edge if it is "negative" in the loop
  if (boundEdgeSenses[0] == -1)
  {
    std::reverse(List_i.begin(), List_i.end());
  }
  nList.clear();

  // get nodes of bounding edge 1 into List_ii
  boundEdges[1]->get_mesh(0, nList, true);
  int size_ii = nList.size() - 1;
  for (ix = 0; ix < nList.size(); ix++)
  {
    List_ii.push_back((iBase_EntityHandle) nList[ix]);
  }
  // we reverse the second boundary edge if it is "positive" in the loop
  if (boundEdgeSenses[1] == 1)
  {
    std::reverse(List_ii.begin(), List_ii.end());
  }

  if (size_i != size_ii)
    ECERRCHK(MK_FAILURE, "Opposite edges have different mesh count, abort");

  // get nodes that are on the linking surface edge, if any,
  // and identify where they start on the source and target
  int linkingEdgeNodeI = -1;
  int linkingEdgeNodeII = -1;
  int offset = 0;
  std::vector<iBase_EntityHandle> linkEdgeNodeList;
  if (linkSurfEdge != NULL)
  {
    std::vector<moab::EntityHandle> lseNodes;
    linkSurfEdge->get_mesh(0, lseNodes, true);
    for (ix = 0; ix < lseNodes.size(); ++ix)
    {
      linkEdgeNodeList.push_back((iBase_EntityHandle)lseNodes[ix]);
    }
    iBase_EntityHandle nodeOnI = linkEdgeNodeList[0];
    iBase_EntityHandle nodeOnII = linkEdgeNodeList[linkEdgeNodeList.size() - 1];
    for (ix = 0; ix < List_i.size(); ix++)
    {
      if (nodeOnI == List_i[ix])
      {
        linkingEdgeNodeI = (int)ix;
      }
      else if (nodeOnI == List_ii[ix])
      {
        linkingEdgeNodeII = (int)ix;
      }
      if (nodeOnII == List_i[ix])
      {
        linkingEdgeNodeI = (int)ix;
      }
      else if (nodeOnII == List_ii[ix])
      {
        linkingEdgeNodeII = (int)ix;
      }
      if (linkingEdgeNodeI != -1 && linkingEdgeNodeII != -1)
      {
        offset = linkingEdgeNodeII - linkingEdgeNodeI;
        if (offset < 0)
        {
          offset += size_i;
        }
        break;
      }
    }
    if (linkingEdgeNodeI == -1 || linkingEdgeNodeII == -1)
    {
      ECERRCHK(MK_FAILURE, "Could not find vertices of linking surface edge on source and target.");
    }
    if (nodeOnI != List_i[linkingEdgeNodeI])
    {
      std::reverse(linkEdgeNodeList.begin(), linkEdgeNodeList.end());
      nodeOnI = linkEdgeNodeList[0];
      nodeOnII = linkEdgeNodeList[linkEdgeNodeList.size() - 1];
    }
  }
  // done with all the initalizations
	
  // get all the position vectors in 3D
  std::vector<Vector3D> pos_i(size_i), pos_ii(size_ii);
  iGeom::Error g_err =
      mk_core()->imesh_instance()->getVtxArrCoords(&(List_i[0]),
      size_i, iBase_INTERLEAVED, &(pos_i[0][0]));
  IBERRCHK(g_err, "Trouble get the xyz coordinates for nodes i");
  g_err = mk_core()->imesh_instance()->getVtxArrCoords(&(List_ii[0]),
      size_ii, iBase_INTERLEAVED, &(pos_ii[0][0]));
  IBERRCHK(g_err, "Trouble get the xyz coordinates for nodes i");

  // compute the interior nodes based on transforming the source and target
  // edges in position
  unsigned int numCreatedNodes = 0;
  unsigned int mesh_count = 0;
  if (!linkEdgeNodeList.empty())
  {
    mesh_count = linkEdgeNodeList.size() - 1;
    if (((int)mesh_count) != ent->mk_core()->sizing_function(
        ent->sizing_function_index())->intervals())
    {
      std::cout << "Warning: The number of nodes on the linking surface edge "
          << "does not match the number of intervals from the sizing "
          << "function on the surface.\n";
    }
    numCreatedNodes = (size_i - 1) * (mesh_count - 1);
  }
  else
  {
    mesh_count = ent->mk_core()->sizing_function(
        ent->sizing_function_index())->intervals();
    numCreatedNodes = size_i * (mesh_count - 1);
  }

  std::vector<iBase_EntityHandle> createdNodes(numCreatedNodes);
  std::vector<iBase_EntityHandle> interiorNodes(size_i * (mesh_count-1));

  Vector3D c0, c1;
  for (unsigned int k = 1; k < mesh_count; k++) //compute layer by layer
  {
    double interpolationFactor = 1.0-double(k)/double(mesh_count);
    for (int i = 0; i < size_i; i++){
      c0 = pos_i[i];
      c1 = pos_ii[(i + offset) % size_i];
      Vector3D pts = c0*interpolationFactor + c1*(1.0-interpolationFactor);
      Vector3D coords; 
      g_err = ent->igeom_instance()->getEntClosestPtTrimmed(ent->geom_handle(), pts[0], pts[1], pts[2], coords[0], coords[1], coords[2]);
      if (g_err)
      {
        g_err = ent->igeom_instance()->getEntClosestPt(ent->geom_handle(), pts[0], pts[1], pts[2], coords[0], coords[1], coords[2]);
      }
      IBERRCHK(g_err, "Trouble get the closest xyz coordinates on the linking surface.");
      int pastLinkEdgeOffset = 0;
      if (i == linkingEdgeNodeI)
      {
        // this node corresponds to one that should already exist on the edge
        interiorNodes[(k - 1)*size_i +i] = linkEdgeNodeList[k];
      }
      else if (i > linkingEdgeNodeI)
      {
        pastLinkEdgeOffset = -1;
      }
      iMesh::Error m_err =
          mk_core()->imesh_instance()->createVtx(coords[0], coords[1],
          coords[2], interiorNodes[(k-1)*size_i+i]);
      if (i != linkingEdgeNodeI)
      {
        createdNodes[(k - 1)*(size_i - 1) + i + pastLinkEdgeOffset] =
            interiorNodes[(k - 1)*size_i + i];
      }
      IBERRCHK(m_err, "Trouble create the interior node.");			
    }
  }

  //finish creating the interior nodes
  iBase_EntitySetHandle entityset;
  iRel::Error r_err = mk_core()->irel_pair(irelPairIndex)->getEntSetRelation(ent->geom_handle(), 0, entityset);
  if (r_err){
    //create the entityset
    iMesh::Error m_err = mk_core()->imesh_instance()->createEntSet(true, entityset);
    IBERRCHK(m_err, "Trouble create the entity set.");

    r_err = mk_core()->irel_pair(irelPairIndex)->setEntSetRelation(ent->geom_handle(), entityset);
    IBERRCHK(r_err, "Trouble create the association between the geometry and mesh entity set.");
  }

  iMesh::Error m_err = mk_core()->imesh_instance()->addEntArrToSet(&createdNodes[0], createdNodes.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of entities to the mesh entity set.");

  // copy nodes in a vector to create quads easier
  std::vector<iBase_EntityHandle> Nodes((mesh_count+1)*size_i);
  //create the int data for mesh nodes on the linking surface
  iBase_TagHandle mesh_tag;
  // TODO: Don't use a tag here, since multiple TFIMapping may occur at the
  // same time
  m_err = mk_core()->imesh_instance()->getTagHandle("MeshTFIMapping", mesh_tag);
  if (m_err)
  {
    m_err = mk_core()->imesh_instance()->createTag("MeshTFIMapping", 1, iBase_INTEGER, mesh_tag);
    IBERRCHK(m_err, "Trouble create the mesh_tag for the surface.");
  }
  int intdata = -1;
  for (int i = 0; i < size_i; i++){
    intdata++;
    m_err = mk_core()->imesh_instance()->setIntData(List_i[i], mesh_tag, intdata);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");

    Nodes[i] = List_i[i];

    m_err = mk_core()->imesh_instance()->setIntData(List_ii[i], mesh_tag, (intdata + mesh_count*size_i));
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    Nodes[size_i*mesh_count+i] = List_ii[i];
  }

  for (int ii = 0; ii < int(interiorNodes.size()); ii++){
    intdata++;
    m_err = mk_core()->imesh_instance()->setIntData(interiorNodes[ii], mesh_tag, intdata);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");

    Nodes[intdata] = interiorNodes[ii];
  }

  //we will always create them in the positive orientation, because we already reversed the Lists with nodes
  std::vector<iBase_EntityHandle> qNodes(4);//a generic quad
  std::vector<iBase_EntityHandle> Quads(size_i*mesh_count);
  for (unsigned int k = 0; k < mesh_count; k++){		
    for (int i = 0; i < size_i; i++){
      qNodes[0] = Nodes[ k*size_i + i ];
      qNodes[1] = Nodes[ k*size_i + (i + 1)%size_i ];
      qNodes[2] = Nodes[ (k+1)*size_i + (i + 1)%size_i ];
      qNodes[3] = Nodes[ (k+1)*size_i + i ];
      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[k*size_i+i]);
      IBERRCHK(m_err, "Trouble create the quadrilateral element.");
    }
  }

  //finish creating the quads
  m_err = mk_core()->imesh_instance()->addEntArrToSet(&Quads[0], Quads.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of quads to the mesh entity set.");
  //set int data for quads
  for (unsigned int i = 0; i < Quads.size(); i++)                               
  {                                                                             
    m_err = mk_core()->imesh_instance()->setIntData(Quads[i], mesh_tag, i);     
    IBERRCHK(m_err, "Trouble set the int data for quadrilateral elements.");    
  }

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

  for (unsigned int i = 0; i < m_Nodes.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->setIntData(m_Nodes[i], mesh_id_tag, i);
    IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");
  }
  for (unsigned int i = 0; i < m_Edges.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->setIntData(m_Edges[i], mesh_id_tag, i);
    IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");
  }
  for (unsigned int i = 0; i < m_Quads.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->setIntData(m_Quads[i], mesh_id_tag, i);
    IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");
  }

  //SurfImprove(ent->geom_handle(), entityset, iBase_FACE);
  mk_core()->save_mesh("InitialMapping.vtk");

  if (_shapeImprove)
  {
#ifdef HAVE_MESQUITE

    iGeom * ig_inst = mk_core()->igeom_instance(ent->iGeomIndex());

    MeshImprove meshopt(mk_core(), true, false, false, false, ig_inst);
    meshopt.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);

#endif

    mk_core()->save_mesh("AfterWinslow.vtk");
#ifdef HAVE_MESQUITE
    MeshImprove shapesmooth(mk_core(), false, false, true, false, ig_inst);
    shapesmooth.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);
#endif
  }

  //remove the mesh tag
  m_err = mk_core()->imesh_instance()->rmvArrTag(&Quads[0], Quads.size(), mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&List_i[0], List_i.size(), mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&List_ii[0], List_ii.size(), mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");

  if (interiorNodes.size() > 0)
  {
    m_err = mk_core()->imesh_instance()->rmvArrTag(&interiorNodes[0], interiorNodes.size(), mesh_tag);
    IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  }

  return 1; 
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

  int irelPairIndex = ent->iRelPairIndex();
  MEntVector boundEdges;
  std::vector<int> senses, group_sizes;
  ent->ModelEnt::boundary(1, boundEdges, &senses, &group_sizes);
  //remove this constraint in case of the side-cylinder case
  //if (boundEdges.size() != 4)
  //  ECERRCHK(MK_FAILURE, "bad input for TFI Mapping, we can only mesh surfaces with 4 edges");
  std::cout << "Surf mapping\n";

  std::vector<iBase_EntityHandle> List_i, List_j, List_ii, List_jj;

  std::vector<moab::EntityHandle> nList;
  /*
   corner[2]     NodeList_ii ->   corner[3]
   ^                                ^
   |                                |
   NodeList_j                    NodeList_jj
   ^                                ^
   |                                |
   corner[0]     NodeList_i  ->   corner[1]
   */
  //get the nodes on each edge
  // we want the start of node list i and j to be the same (corner 0)
  boundEdges[0]->get_mesh(0, nList, true); // include start and end vertices (corners)
  unsigned int ix = 0;
  int size_i=(int)nList.size();
  for (ix = 0; ix < nList.size(); ix++)
    List_i.push_back((iBase_EntityHandle) nList[ix]);
  // if sense is reverse for edge 0, reverse  list,
  if (senses[0] == -1)
    std::reverse(List_i.begin(), List_i.end());
  // so we know for sure corner 0 is at NodeList_i[0]!!
  nList.clear();
  boundEdges[1]->get_mesh(0, nList, true);
  int size_j=(int)nList.size();
  for (ix = 0; ix < nList.size(); ix++)
    List_jj.push_back((iBase_EntityHandle) nList[ix]);
  if (senses[1] == -1)
    std::reverse(List_jj.begin(), List_jj.end());

  nList.clear();
  boundEdges[2]->get_mesh(0, nList, true);
  for (ix = 0; ix < nList.size(); ix++)
    List_ii.push_back((iBase_EntityHandle) nList[ix]);
  if (senses[2] == 1) // we reverse it if this edge is "positive" in the loop
    std::reverse(List_ii.begin(), List_ii.end());

  nList.clear();
  boundEdges[3]->get_mesh(0, nList, true);
  for (ix = 0; ix < nList.size(); ix++)
    List_j.push_back((iBase_EntityHandle) nList[ix]);
  if (senses[3] == 1) // we reverse it if this edge is "positive" in the loop
    std::reverse(List_j.begin(), List_j.end());

  if (List_i.size() != List_ii.size())
    ECERRCHK(MK_FAILURE, "opposite edges have different mesh count, abort");
  if (List_j.size() != List_jj.size())
    ECERRCHK(MK_FAILURE, "opposite edges have different mesh count, abort");
  //ok, done with all the initializations

  // get all the vectors in 3d
  std::vector<Vector3D> pos_ii(List_ii.size());
  std::vector<Vector3D> pos_i(List_i.size());
  std::vector<Vector3D> pos_j(List_j.size());
  std::vector<Vector3D> pos_jj(List_jj.size());
  // iBase_INTERLEAVED
  /*getVtxArrCoords( const EntityHandle* vertex_handles,
                                   int vertex_handles_size,
                                   StorageOrder storage_order,
                                   double* coords_out ) const*/
  iGeom::Error g_err = mk_core()->imesh_instance()->getVtxArrCoords(&(List_ii[0]), List_ii.size(), iBase_INTERLEAVED, &(pos_ii[0][0]));
  IBERRCHK(g_err, "Trouble get the xyz coordinates for nodes ii.");
  g_err = mk_core()->imesh_instance()->getVtxArrCoords(&(List_i[0]), List_i.size(), iBase_INTERLEAVED, &(pos_i[0][0]));
  IBERRCHK(g_err, "Trouble get the xyz coordinates for nodes i.");
  g_err = mk_core()->imesh_instance()->getVtxArrCoords(&(List_j[0]), List_j.size(), iBase_INTERLEAVED, &(pos_j[0][0]));
  IBERRCHK(g_err, "Trouble get the xyz coordinates for nodes j.");
  g_err = mk_core()->imesh_instance()->getVtxArrCoords(&(List_jj[0]), List_jj.size(), iBase_INTERLEAVED, &(pos_jj[0][0]));
  IBERRCHK(g_err, "Trouble get the xyz coordinates for nodes jj.");
  //calculate the interior nodes based on transforming the top and bottom edges in position
  std::vector<iBase_EntityHandle> interiorNodes((List_j.size() - 2) * (List_i.size() - 2));
  // reminder
  /*
     corner[2]     NodeList_ii ->   corner[3]
     ^                                ^
     |                                |
     NodeList_j                    NodeList_jj
     ^                                ^
     |                                |
     corner[0]     NodeList_i  ->   corner[1]
  */
  Vector3D c0=pos_i[0], c1=pos_i[size_i-1], c2 = pos_ii[0], c3=pos_ii[size_i-1];

  Vector3D bc = 0.5*c0+0.5*c1;
  Vector3D tc = 0.5*c2+0.5*c3;
  for (int j = 1; j < size_j - 1; j++) // compute layer by layer
    // we will start from source (layer 0) to layer j (j>1, j< J-1)
    // also , we will look at the target, layer J-1 to layer j
  {

    // transformation from c0 and c1 to layer j
    Matrix3D tr1 , tr2;
    //target= A * ( source - 2*sc + tc) + sc
    Vector3D cj = 0.5*pos_j[j]+0.5*pos_jj[j]; // center of layer j
    computeTransformation(c0, c1, pos_j[j], pos_jj[j], tr1);
    // transformation from top, c2 and c3 to layer j
    computeTransformation(c2, c3, pos_j[j], pos_jj[j], tr2);

    double interpolationFactor = j/(size_j-1.);

    for (int i = 1; i < (size_i - 1); i++)
    {
      // transformation from bottom to layer j; source is b, target is j
      Vector3D res1= tr1*(pos_i[i] -2*bc+cj)+bc;
      // transformation from top to layer j; source is t, target is j
      Vector3D res2= tr2*(pos_ii[i] -2*tc+cj)+tc;
      // interpolate this result
      Vector3D pts = res1*(1-interpolationFactor) + res2*interpolationFactor;
      Vector3D coords;
      g_err = ent->igeom_instance()->getEntClosestPtTrimmed(ent->geom_handle(), pts[0], pts[1], pts[2], coords[0], coords[1],
          coords[2]);
      if (g_err)
      {
        g_err = ent->igeom_instance()->getEntClosestPt(ent->geom_handle(), pts[0], pts[1], pts[2], coords[0], coords[1], coords[2]);
      }
      IBERRCHK(g_err, "Trouble get the closest xyz coordinates on the linking surface.");

      iMesh::Error m_err = mk_core()->imesh_instance()->createVtx(coords[0], coords[1], coords[2], interiorNodes[(j - 1)
          * (size_i - 2) + i - 1]);
      IBERRCHK(m_err, "Trouble create the interior node.");
    }
  }

  //finish creating the interior nodes
  iBase_EntitySetHandle entityset;
  iRel::Error r_err = mk_core()->irel_pair(irelPairIndex)->getEntSetRelation(ent->geom_handle(), 0, entityset);
  if (r_err)
  {
    //create the entityset
    iMesh::Error m_err = mk_core()->imesh_instance()->createEntSet(true, entityset);
    IBERRCHK(m_err, "Trouble create the entity set.");

    r_err = mk_core()->irel_pair(irelPairIndex)->setEntSetRelation(ent->geom_handle(), entityset);
    IBERRCHK(r_err, "Trouble create the association between the geometry and mesh entity set.");
  }

  iMesh::Error m_err = mk_core()->imesh_instance()->addEntArrToSet(&interiorNodes[0], interiorNodes.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of entities to the mesh entity set.");

  // copy nodes in a vector to create the quads easier
  // they will be arranged in layers, from bottom (j=0) towards top (j=size_j-1)
  std::vector<iBase_EntityHandle> Nodes(size_j * size_i);

  //create the int data for mesh nodes on the linking surface
  iBase_TagHandle mesh_tag;
  m_err = mk_core()->imesh_instance()->getTagHandle("MeshTFIMapping", mesh_tag);
  if (m_err)
  {
    m_err = mk_core()->imesh_instance()->createTag("MeshTFIMapping", 1, iBase_INTEGER, mesh_tag);
    IBERRCHK(m_err, "Trouble create the mesh_tag for the surface.");
  }
  int intdata = -1;
  for (int i = 0; i < size_i; i++)
  {
    intdata++;
    m_err = mk_core()->imesh_instance()->setIntData(List_i[i], mesh_tag, intdata);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    // bottom row, j=0
    Nodes[i]=List_i[i];

    m_err = mk_core()->imesh_instance()->setIntData(List_ii[i], mesh_tag, intdata + size_i);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    // top row, j = size_j-1
    Nodes[ (size_i)*(size_j-1)+i] = List_ii[i];
  }
  intdata = 2 * size_i - 1;
  for (int j = 1; j < size_j - 1; j++)
  {
    intdata++;
    m_err = mk_core()->imesh_instance()->setIntData(List_j[j], mesh_tag, intdata);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    // right column, i=0
    Nodes[size_i*j] = List_j[j];

    m_err = mk_core()->imesh_instance()->setIntData(List_jj[j], mesh_tag, intdata + size_j - 2);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    // left column, i = size_i -1
    Nodes[size_i*j + size_i-1] = List_jj[j];
  }

  intdata = 2 * size_i + 2 * (size_j - 2) - 1;
  // it is clear that (size_i-2 > 0 and size_j-2 > 0) iff (interiorNodes.size()>0)
  for (unsigned int ii = 0; ii < interiorNodes.size(); ii++)
  {
    intdata++;
    m_err = mk_core()->imesh_instance()->setIntData(interiorNodes[ii], mesh_tag, intdata);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    int j = ii/(size_i-2) + 1;
    int i = (ii -(j-1)*(size_i-2)) + 1;
    // copy
    Nodes[ j*size_i + i] = interiorNodes[ii];
    // compute the row and column
  }

  // we will always create them in the positive orientation, because we already reversed the Lists
  // with nodes

  std::vector<iBase_EntityHandle> qNodes(4);// a generic quad
  std::vector<iBase_EntityHandle> Quads((size_j - 1) * (size_i - 1));
  for (int j=0; j <size_j-1; j++)
  {
    for (int i=0; i< size_i-1; i++)
    {
      qNodes[0] = Nodes[   j  *size_i+i  ];
      qNodes[1] = Nodes[   j  *size_i+i+1];
      qNodes[2] = Nodes[ (j+1)*size_i+i+1];
      qNodes[3] = Nodes[ (j+1)*size_i+i  ];
      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL, &qNodes[0], 4, Quads[j*(size_i-1)+i]);
      IBERRCHK(m_err, "Trouble create the quadrilateral element.");
    }
  }

  //finish creating the quads
  m_err = mk_core()->imesh_instance()->addEntArrToSet(&Quads[0], Quads.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of quads to the mesh entity set.");
  //set int data for quads
  for (unsigned int i = 0; i < Quads.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->setIntData(Quads[i], mesh_tag, i);
    IBERRCHK(m_err, "Trouble set the int data for quadrilateral elements.");
  }

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

  for (unsigned int i = 0; i < m_Nodes.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->setIntData(m_Nodes[i], mesh_id_tag, i);
    IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");
  }
  for (unsigned int i = 0; i < m_Edges.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->setIntData(m_Edges[i], mesh_id_tag, i);
    IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");
  }
  for (unsigned int i = 0; i < m_Quads.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->setIntData(m_Quads[i], mesh_id_tag, i);
    IBERRCHK(m_err, "Trouble set the int data for 'GLOBAL_ID'.");
  }

  //SurfImprove(ent->geom_handle(), entityset, iBase_FACE);
  mk_core()->save_mesh("InitialMapping.vtk");

  if (_shapeImprove)
  {
#ifdef HAVE_MESQUITE

    iGeom * ig_inst = mk_core()->igeom_instance(ent->iGeomIndex());

    MeshImprove meshopt(mk_core(), true, false, false, false, ig_inst);
    meshopt.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);

#endif
    //mk_core()->save_mesh("AfterLaplace.vtk");

    //if there is the parametric space, let Winslow smooth inside the parametric space
    SmoothWinslow(List_i, List_ii, List_j, List_jj, interiorNodes, Quads, mesh_tag, ent);

    mk_core()->save_mesh("AfterWinslow.vtk");
#ifdef HAVE_MESQUITE
    MeshImprove shapesmooth(mk_core(), false, false, true, false, ig_inst);
    shapesmooth.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);
#endif
  }
  //remove the mesh tag
  m_err = mk_core()->imesh_instance()->rmvArrTag(&Quads[0], Quads.size(), mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&List_i[0], List_i.size(), mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&List_ii[0], List_ii.size(), mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  if (List_j.size() > 2)
  {
    m_err = mk_core()->imesh_instance()->rmvArrTag(&List_j[0], List_j.size() - 2, mesh_tag);
    IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
    m_err = mk_core()->imesh_instance()->rmvArrTag(&List_jj[0], List_jj.size() - 2, mesh_tag);
    IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  }
  if (interiorNodes.size() > 0)
  {
    m_err = mk_core()->imesh_instance()->rmvArrTag(&interiorNodes[0], interiorNodes.size(), mesh_tag);
    IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  }

  return 1;
}
void TFIMapping::computeTransformation(Vector3D & A, Vector3D & B, Vector3D & C, Vector3D & D,
    Matrix3D & M)
{
  Matrix3D tmpMatrix;
  Matrix3D bMatrix;
  Vector3D c1=0.5*A+0.5*B;
  Vector3D c2=0.5*C+0.5*D;
  Vector3D normal=(A-D)*(B-C);// this should be not modified by the transformation
  // we add this just to increase the rank of tmpMatrix
  // the normal to the "plane" of interest should not be rotated at all
  // as if we say that we want A*normal = normal
  Vector3D s1=A-2*c1+c2;
  Vector3D s2=B-2*c1+c2;
  Vector3D t1=C-c1;
  Vector3D t2=D-c1;
  // so we are looking for M such that
  /*
   *  M*s1 = t1
   *  M*s2 = t2
   *  M*n  = n
   *
   *  In simple cases, M is identity
   */

  tmpMatrix.set_column(0, s1);
  tmpMatrix.set_column(1, s2);
  tmpMatrix.set_column(2, normal);
  bMatrix.set_column(0, t1);
  bMatrix.set_column(1, t2);
  bMatrix.set_column(2, normal);

  double detValue = det(tmpMatrix);
  (void) detValue;
  assert(detValue*detValue>1.e-20);

  //solve the affine mapping matrix, make use of inverse matrix to get affine mapping matrix
  Matrix3D InvMatrix = inverse(tmpMatrix);
  M = bMatrix*InvMatrix;

}
//smooth the quadrilateral mesh on the linking surface
void TFIMapping::SmoothWinslow(std::vector<iBase_EntityHandle> &List_i, std::vector<iBase_EntityHandle> &List_ii, std::vector<
    iBase_EntityHandle> &List_j, std::vector<iBase_EntityHandle> &List_jj, std::vector<iBase_EntityHandle> &interiorNodes,
    std::vector<iBase_EntityHandle> &quads, iBase_TagHandle &taghandle, ModelEnt *ent)
{
  std::vector<std::set<int> > AdjElements;
  std::vector<std::vector<int> > Quads;
  std::vector<std::vector<double> > coords;
  std::vector<bool> isBnd;
  std::vector<iBase_EntityHandle> nodes;
  std::vector<double> weight;

  bool isParameterized = false;

  iGeom::Error g_err = ent->igeom_instance()->isEntParametric(ent->geom_handle(), isParameterized);
  IBERRCHK(g_err, "Trouble check whether the surface is parameterized or not.");
  isParameterized = false;

  //resize the coords to store all the nodes's coordinates on the linking surface
  coords.resize(List_i.size() * List_j.size());
  isBnd.resize(coords.size());
  nodes.resize(coords.size());
  for (unsigned int i = 0; i < coords.size(); i++)
    coords[i].resize(3);

  iMesh::Error m_err;
  //input the boundary nodes
  for (unsigned int i = 0; i < List_i.size(); i++)
  {
    m_err = mk_core()->imesh_instance()->getVtxCoord(List_i[i], coords[i][0], coords[i][1], coords[i][2]);
    IBERRCHK(m_err, "Trouble get the vertex coordinates.");
    nodes[i] = List_i[i];
    isBnd[i] = true;

    m_err = mk_core()->imesh_instance()->getVtxCoord(List_ii[i], coords[List_i.size() + i][0], coords[List_i.size() + i][1],
        coords[List_i.size() + i][2]);
    IBERRCHK(m_err, "Trouble get the vertex coordinates.");
    if (isParameterized)
    {
      double uv[2];
      g_err = ent->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), coords[List_i.size() + i][0], coords[List_i.size() + i][1],
          coords[List_i.size() + i][2], uv[0], uv[1]);
      IBERRCHK(g_err, "Trouble get the uv from xyz.");
      coords[List_i.size() + i][0] = uv[0];
      coords[List_i.size() + i][1] = uv[1];
    }

    nodes[List_i.size() + i] = List_ii[i];
    isBnd[List_i.size() + i] = true;
  }
  if (int(List_j.size()) > 2)
  {
    for (unsigned int i = 1; i < (List_j.size() - 1); i++)
    {
      m_err = mk_core()->imesh_instance()->getVtxCoord(List_j[i], coords[2 * List_i.size() + i - 1][0], coords[2 * List_i.size()
          + i - 1][1], coords[2 * List_i.size() + i - 1][2]);
      IBERRCHK(m_err, "Trouble get the vertex coordinates.");
      nodes[2 * List_i.size() + i - 1] = List_j[i];
      isBnd[2 * List_i.size() + i - 1] = true;

      m_err = mk_core()->imesh_instance()->getVtxCoord(List_jj[i], coords[2 * List_i.size() + List_j.size() - 2 + i - 1][0],
          coords[2 * List_i.size() + List_j.size() - 2 + i - 1][1], coords[2 * List_i.size() + List_j.size() - 2 + i - 1][2]);
      IBERRCHK(m_err, "Trouble get the vertex coordinates.");
      nodes[2 * List_i.size() + List_j.size() - 2 + i - 1] = List_jj[i];
      isBnd[2 * List_i.size() + List_j.size() - 2 + i - 1] = true;
    }
  }
  //input the interior nodes
  if (interiorNodes.size() > 0)
  {
    for (unsigned int i = 0; i < interiorNodes.size(); i++)
    {
      m_err = mk_core()->imesh_instance()->getVtxCoord(interiorNodes[i],
          coords[2 * List_i.size() + 2 * (List_j.size() - 2) + i][0], coords[2 * List_i.size() + 2 * (List_j.size() - 2) + i][1],
          coords[2 * List_i.size() + 2 * (List_j.size() - 2) + i][2]);
      IBERRCHK(m_err, "Trouble get the vertex coordinates.");
      nodes[2 * List_i.size() + 2 * (List_j.size() - 2) + i] = interiorNodes[i];
      isBnd[2 * List_i.size() + 2 * (List_j.size() - 2) + i] = false;
    }
  }

  //update the AdjElements info
  //notice: during this process, adjacent quads will be returned around a node. The quads from source surface and target surface may be returned.
  AdjElements.resize(nodes.size());
  for (unsigned int i = 0; i < AdjElements.size(); i++)
  {
    if (!isBnd[i])
    {
      std::vector<iBase_EntityHandle> adjEnts;
      adjEnts.clear();
      m_err = mk_core()->imesh_instance()->getEntAdj(nodes[i], iBase_FACE, adjEnts);
      IBERRCHK(m_err, "Trouble get the adjacent quads wrt a node.");
      for (unsigned int j = 0; j < adjEnts.size(); j++)
      {
        int index_id = -1;
        m_err = mk_core()->imesh_instance()->getIntData(adjEnts[j], taghandle, index_id);
        IBERRCHK(m_err, "Trouble get int data for quads.");
        AdjElements[i].insert(index_id);
        //std::cout<< " i= " << i << " index_id:" << index_id << "\n";
      }
    }
  }

  //update the Quads' info
  Quads.resize(quads.size());
  for (unsigned int i = 0; i < Quads.size(); i++)
  {
    std::vector<iBase_EntityHandle> adjEnts;
    adjEnts.clear();
    m_err = mk_core()->imesh_instance()->getEntAdj(quads[i], iBase_VERTEX, adjEnts);
    IBERRCHK(m_err, "Trouble get the adjacent nodes wrt a quad.");

    assert(adjEnts.size()==4);
    Quads[i].resize(4);

    for (unsigned int j = 0; j < adjEnts.size(); j++)
    {
      int index_id = -1;
      m_err = mk_core()->imesh_instance()->getIntData(adjEnts[j], taghandle, index_id);
      IBERRCHK(m_err, "Trouble get int data for nodes.");
      Quads[i][j] = index_id;
    }
  }

  //detect the connectivity
  std::vector<std::vector<int> > connect(nodes.size(), std::vector<int>(9));
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    if (!isBnd[i])
    {
      //there are 4 adjacent quadrilateral elements around node i
      //std::cout << " element i:" << i << " AdjElements[i].size() " << AdjElements[i].size() << "\n";
      assert(AdjElements[i].size() == 4);
      std::set<int>::iterator it = AdjElements[i].begin();
      int st_index[4];
      //process 4 quad elements
      int j = -1;
      for (; it != AdjElements[i].end(); it++)
      {
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
      connect[i][2] = Quads[*it][(st_index[0] + 3) % 4];
      connect[i][8] = Quads[*it][(st_index[0] + 1) % 4];
      connect[i][1] = Quads[*it][(st_index[0] + 2) % 4];
      //finish processing the quad 1
      std::set<int>::iterator it1 = AdjElements[i].begin();
      it1++;
      for (j = 1; j < 4; j++, it1++)
      {
        if (connect[i][8] == Quads[*it1][(st_index[j] + 1) % 4])
        {
          connect[i][7] = Quads[*it1][(st_index[j] + 2) % 4];
          connect[i][6] = Quads[*it1][(st_index[j] + 3) % 4];
          break;
        }
        else if (connect[i][8] == Quads[*it1][(st_index[j] + 3) % 4])
        {
          connect[i][7] = Quads[*it1][(st_index[j] + 2) % 4];
          connect[i][6] = Quads[*it1][(st_index[j] + 1) % 4];
          break;
        }
        else
          continue;
      }
      //finish processing the quad 2
      std::set<int>::iterator it2 = AdjElements[i].begin();
      it2++;
      for (j = 1; it2 != AdjElements[i].end(); it2++, j++)
      {
        if (connect[i][2] == Quads[*it2][(st_index[j] + 1) % 4])
        {
          connect[i][3] = Quads[*it2][(st_index[j] + 2) % 4];
          connect[i][4] = Quads[*it2][(st_index[j] + 3) % 4];
          break;
        }
        else if (connect[i][2] == Quads[*it2][(st_index[j] + 3) % 4])
        {
          connect[i][3] = Quads[*it2][(st_index[j] + 2) % 4];
          connect[i][4] = Quads[*it2][(st_index[j] + 1) % 4];
          break;
        }
        else
          continue;
      }
      //finish processing the quad 4;
      std::set<int>::iterator it3 = AdjElements[i].begin();
      it3++;
      for (j = 1; it3 != AdjElements[i].end(); it3++, j++)
      {
        if ((it3 != it1) && (it3 != it2))
        {
          connect[i][5] = Quads[*it2][(st_index[j] + 2) % 4];
        }
        else
          continue;
      }
    }
  }
  //finish all the initialization

  EquipotentialSmooth smoother;

  //IsoLaplace smoother;

  smoother.SetupData(AdjElements, Quads, coords, isBnd, connect);
  smoother.Execute();

  //std::vector<std::vector<double> > coors;
  smoother.GetCoords(coords);

  //update the new position for nodes
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    if (!isBnd[i])
    {
      double tmp_coord[3] = { coords[i][0], coords[i][1], coords[i][2] };
      if (!isParameterized)
      {
        iGeom::Error g_err = ent->igeom_instance()->getEntClosestPt(ent->geom_handle(), coords[i][0], coords[i][1], coords[i][2],
            tmp_coord[0], tmp_coord[1], tmp_coord[2]);
        IBERRCHK(g_err, "Trouble get a closest point on the linking surface.");
      }
      else
      {
        iGeom::Error g_err = ent->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), coords[i][0], coords[i][1], tmp_coord[0],
            tmp_coord[1], tmp_coord[2]);
        IBERRCHK(g_err, "Trouble get the xyz from uv.");

      }
      m_err = mk_core()->imesh_instance()->setVtxCoord(nodes[i], tmp_coord[0], tmp_coord[1], tmp_coord[2]);
      IBERRCHK(m_err, "Trouble set the new coordinates for nodes.");
    }
  }

  //remove the unnecessary tag after smoothing
  m_err = mk_core()->imesh_instance()->rmvArrTag(&nodes[0], nodes.size(), taghandle);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&quads[0], quads.size(), taghandle);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  //m_err = mk_core()->imesh_instance()->destroyTag(taghandle, 1);
  //IBERRCHK(m_err, "Trouble destroy a tag.");
}

} // namespace MeshKit

