#include "meshkit/TFIMapping.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "moab/ReadUtilIface.hpp"
#include "EquipotentialSmooth.hpp"
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
//Entity Type initialization for TFIMapping meshing
moab::EntityType TFIMapping_tps[] = { moab::MBVERTEX, moab::MBEDGE,
    moab::MBQUAD, moab::MBMAXTYPE };
const moab::EntityType* TFIMapping::output_types()
{
  return TFIMapping_tps;
}

//---------------------------------------------------------------------------//
// construction function for TFIMapping class
TFIMapping::TFIMapping(MKCore *mk_core, const MEntVector &me_vec) :
  MeshScheme(mk_core, me_vec)
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

  /* the only things we need to make sure :
   1) there are 4 edges, exactly
   2) the opposite edges have the same meshcount
   - if some of the edges are meshed, we need to mesh the opposite edge with
   correct meshcount
   - if 2 opposite edges are meshed, verify the mesh count
   */
  // get iGeom instance from the first ment selection
  if (mentSelection.empty())
    return;

  //loop over the surfaces
  for (MEntSelection::iterator mit = mentSelection.begin(); mit
      != mentSelection.end(); mit++)
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
      mesh_count_surface
          = me->mk_core()->sizing_function(size_index)->intervals();

    // get the boundary loop of edges
    MEntVector boundEdges;
    std::vector<int> senses, group_sizes;
    me->ModelEnt::boundary(1, boundEdges, &senses, &group_sizes);
    if (boundEdges.size() != 4)
      ECERRCHK(MK_FAILURE, "bad input for TFI Mapping, we can only mesh surfaces with 4 edges");

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
            if (sfe->intervals() > 0 && !force)// if a sizing function was set on an edge, use it, do not
              // use the mesh count from surface
              // still, a mesh count from opposing edge is very powerful
              mesh_count = sfe->intervals();
          }
          // push it to the list if it is not setup to another mesh op (edge mesher) already
          //if (oppEdges[j]->is_meshops_list_empty())// it will create an EdgeMesher later
          edgesToMesh.push_back(oppEdges[j]);
        }
      }
      // decide on a mesh count now, if edgesToMesh.size()>0
      if (edgesToMesh.size() > 0)
      {

        EdgeMesher * em = (EdgeMesher*) me->mk_core()->construct_meshop(
            "EdgeMesher", edgesToMesh);
        if (mesh_count < 0)
        {
          std::cout << "mesh count not set properly on opposite edges, set it to 10\n";
          mesh_count = 10; // 4 is a nice number, used in the default edge mesher;
          // but I like 10 more
        }

        for (unsigned int j = 0; j < edgesToMesh.size(); j++)
        {
          edgesToMesh[j]->mesh_intervals(mesh_count);
          edgesToMesh[j]->interval_firmness(HARD);
          edgesToMesh[j]->add_meshop(em);
        }
        mk_core()->insert_node(em, (MeshOp*) this);

      }
    } // end loop over pair of opposite edges
  }// end loop over surfaces
  mk_core()->print_graph("AfterTFISetup.eps");
}

//---------------------------------------------------------------------------//
// execute function: generate the all-quad mesh through the TFI mapping
void TFIMapping::execute_this()
{

  //loop over the surfaces
  for (MEntSelection::iterator mit = mentSelection.begin(); mit
      != mentSelection.end(); mit++)
  {
    ModelEnt *me = mit -> first;
    //first check whether the surface is meshed or not
    if (me->get_meshed_state() >= COMPLETE_MESH)
      continue;

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

  int irelPairIndex = ent->iRelPairIndex();
  MEntVector boundEdges;
  std::vector<int> senses, group_sizes;
  ent->ModelEnt::boundary(1, boundEdges, &senses, &group_sizes);
  if (boundEdges.size() != 4)
    ECERRCHK(MK_FAILURE, "bad input for TFI Mapping, we can only mesh surfaces with 4 edges");

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
  for (ix = 0; ix < nList.size(); ix++)
    List_i.push_back((iBase_EntityHandle) nList[ix]);
  // if sense is reverse for edge 0, reverse  list,
  if (senses[0] == -1)
    std::reverse(List_i.begin(), List_i.end());
  // so we know for sure corner 0 is at NodeList_i[0]!!
  nList.clear();
  boundEdges[1]->get_mesh(0, nList, true);
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

  //calculate the interior nodes based on TFI
  iGeom::Error g_err;
  std::vector<iBase_EntityHandle> InteriorNodes((List_j.size() - 2)
      * (List_i.size() - 2));
  for (unsigned int j = 1; j < (List_j.size() - 1); j++)
  {
    //                     Node 2 (List_ii)
    //		  		 |
    //	Node 0 (List_j)---------Node------------Node 1 (List_jj)
    //		  	 	 |
    //		                  Node 3 (List_i)
    double coords_i[3], coords_ii[3], coords_j[3], coords_jj[3], coords[3];

    g_err = mk_core()->imesh_instance()->getVtxCoord(List_j[j], coords_j[0],
        coords_j[1], coords_j[2]);
    IBERRCHK(g_err, "Trouble get the xyz coordinates for node 0.");

    g_err = mk_core()->imesh_instance()->getVtxCoord(List_jj[j], coords_jj[0],
        coords_jj[1], coords_jj[2]);
    IBERRCHK(g_err, "Trouble get the xyz coordinates for node 1.");

    for (unsigned int i = 1; i < (List_i.size() - 1); i++)
    {
      g_err = mk_core()->imesh_instance()->getVtxCoord(List_i[i], coords_i[0],
          coords_i[1], coords_i[2]);
      IBERRCHK(g_err, "Trouble get the xyz coordinates for node 3.");

      g_err = mk_core()->imesh_instance()->getVtxCoord(List_ii[i],
          coords_ii[0], coords_ii[1], coords_ii[2]);
      IBERRCHK(g_err, "Trouble get the xyz coordinates for node 2.");

      //calculate the parameter based on the distance
      double r, s, pts[3];
      ParameterCalculate(r, s, coords_j, coords_jj, coords_i, coords_ii, pts,
          ent);

      std::cout << "r = " << r << "\ts = " << s << std::endl;

      parametricTFI3D(&pts[0], s, r, coords_j, coords_jj, coords_i, coords_ii);

      g_err = ent->igeom_instance()->getEntClosestPtTrimmed(ent->geom_handle(),
          pts[0], pts[1], pts[2], coords[0], coords[1], coords[2]);
      if (g_err)
      {
        g_err = ent->igeom_instance()->getEntClosestPt(ent->geom_handle(),
            pts[0], pts[1], pts[2], coords[0], coords[1], coords[2]);
      }
      IBERRCHK(g_err, "Trouble get the closest xyz coordinates on the linking surface.");

      iMesh::Error m_err = mk_core()->imesh_instance()->createVtx(coords[0],
          coords[1], coords[2], InteriorNodes[(j - 1) * (List_i.size() - 2) + i
              - 1]);
      IBERRCHK(m_err, "Trouble get create the interior node.");
    }
  }

  //finish creating the interior nodes
  iBase_EntitySetHandle entityset;
  iRel::Error r_err = mk_core()->irel_pair(irelPairIndex)->getEntSetRelation(
      ent->geom_handle(), 0, entityset);
  if (r_err)
  {
    //create the entityset
    iMesh::Error m_err = mk_core()->imesh_instance()->createEntSet(true,
        entityset);
    IBERRCHK(m_err, "Trouble create the entity set.");

    r_err = mk_core()->irel_pair(irelPairIndex)->setEntSetRelation(
        ent->geom_handle(), entityset);
    IBERRCHK(r_err, "Trouble create the association between the geometry and mesh entity set.");
  }

  iMesh::Error m_err = mk_core()->imesh_instance()->addEntArrToSet(
      &InteriorNodes[0], InteriorNodes.size(), entityset);
  IBERRCHK(m_err, "Trouble add an array of entities to the mesh entity set.");

  //create the int data for mesh nodes on the linking surface
  iBase_TagHandle mesh_tag;
  m_err = mk_core()->imesh_instance()->getTagHandle("MeshTFIMapping", mesh_tag);
  if (m_err)
  {
    m_err = mk_core()->imesh_instance()->createTag("MeshTFIMapping", 1,
        iBase_INTEGER, mesh_tag);
    IBERRCHK(m_err, "Trouble create the mesh_tag for the surface.");
  }
  int intdata = -1;
  for (unsigned int i = 0; i < List_i.size(); i++)
  {
    intdata++;
    m_err = mk_core()->imesh_instance()->setIntData(List_i[i], mesh_tag,
        intdata);
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");

    m_err = mk_core()->imesh_instance()->setIntData(List_ii[i], mesh_tag,
        intdata + List_i.size());
    IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
  }
  intdata = 2 * List_i.size() - 1;
  if (List_j.size() > 2)
  {
    for (unsigned int i = 1; i < List_j.size() - 1; i++)
    {
      intdata++;
      m_err = mk_core()->imesh_instance()->setIntData(List_j[i], mesh_tag,
          intdata);
      IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");

      m_err = mk_core()->imesh_instance()->setIntData(List_jj[i], mesh_tag,
          intdata + List_j.size() - 2);
      IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    }
  }
  intdata = 2 * List_i.size() + 2 * (List_j.size() - 2) - 1;
  if (InteriorNodes.size() > 0)
  {
    for (unsigned int i = 0; i < InteriorNodes.size(); i++)
    {
      intdata++;
      m_err = mk_core()->imesh_instance()->setIntData(InteriorNodes[i],
          mesh_tag, intdata);
      IBERRCHK(m_err, "Trouble set the int data for mesh nodes.");
    }

  }

  // we will always create them in the positive orientation, because we already reversed the Lists
  // with nodes

  std::vector<iBase_EntityHandle> Quads((List_j.size() - 1) * (List_i.size()
      - 1));
  for (unsigned int j = 0; j < (List_j.size() - 1); j++)
  {
    std::vector<iBase_EntityHandle> qNodes(4);
    if (j == 0)
    {//starting row boundary

      qNodes[0] = List_i[0];
      qNodes[1] = List_i[1];
      qNodes[2] = InteriorNodes[0];
      qNodes[3] = List_j[1];
      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
          &qNodes[0], 4, Quads[0]);
      IBERRCHK(m_err, "Trouble create the quadrilateral elements.");

      for (unsigned int i = 1; i < (List_i.size() - 2); i++)
      {

        qNodes[0] = List_i[i];
        qNodes[1] = List_i[i + 1];
        qNodes[2] = InteriorNodes[i];
        qNodes[3] = InteriorNodes[i - 1];
        m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
            &qNodes[0], 4, Quads[j * (List_i.size() - 1) + i]);
        IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
      }
      // the last one in first row
      qNodes[0] = List_i[List_i.size() - 2];
      qNodes[1] = List_jj[0];
      qNodes[2] = List_jj[1];
      qNodes[3] = InteriorNodes[List_i.size() - 3];

      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
          &qNodes[0], 4, Quads[List_i.size() - 2]);
      IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
    }

    else if (j == List_j.size() - 2)
    {//ending row boundary
      qNodes[0] = List_j[j];
      qNodes[1] = InteriorNodes[(j - 1) * (List_i.size() - 2)];
      qNodes[2] = List_ii[1];
      qNodes[3] = List_j[j + 1]; // == List_ii[0]!
      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
          &qNodes[0], 4, Quads[j * (List_i.size() - 1)]);
      IBERRCHK(m_err, "Trouble create the quadrilateral elements.");

      for (unsigned int i = 1; i < (List_i.size() - 2); i++)
      {

        qNodes[0] = InteriorNodes[(j - 1) * (List_i.size() - 2) + i - 1];
        qNodes[1] = InteriorNodes[(j - 1) * (List_i.size() - 2) + i];
        qNodes[2] = List_ii[i + 1];
        qNodes[3] = List_ii[i];

        m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
            &qNodes[0], 4, Quads[j * (List_i.size() - 1) + i]);
        IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
      }
      // last one , top row
      qNodes[0] = InteriorNodes[(j - 1) * (List_i.size() - 2) + List_i.size()
          - 3];// the last interior node
      qNodes[1] = List_jj[List_jj.size() - 2];
      qNodes[2] = List_jj[List_jj.size() - 1];
      qNodes[3] = List_ii[List_ii.size() - 2];

      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
          &qNodes[0], 4, Quads[(List_j.size() - 1) * (List_i.size() - 1) - 1]);
      IBERRCHK(m_err, "Trouble create the quadrilateral elements.");

    }
    else
    {

      qNodes[0] = List_j[j];
      qNodes[1] = InteriorNodes[(j - 1) * (List_i.size() - 2)];
      qNodes[2] = InteriorNodes[j * (List_i.size() - 2)];
      qNodes[3] = List_j[j + 1];

      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
          &qNodes[0], 4, Quads[j * (List_i.size() - 1)]);
      IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
      for (unsigned int i = 1; i < (List_i.size() - 2); i++)
      {

        qNodes[0] = InteriorNodes[(j - 1) * (List_i.size() - 2) + i - 1];
        qNodes[1] = InteriorNodes[(j - 1) * (List_i.size() - 2) + i];
        qNodes[2] = InteriorNodes[j * (List_i.size() - 2) + i];
        qNodes[3] = InteriorNodes[j * (List_i.size() - 2) + i - 1];

        m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
            &qNodes[0], 4, Quads[j * (List_i.size() - 1) + i]);
        IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
      }

      //end column

      qNodes[0] = InteriorNodes[(j - 1) * (List_i.size() - 2) + List_i.size()
          - 3];
      qNodes[1] = List_jj[j];
      qNodes[2] = List_jj[j + 1];
      qNodes[3] = InteriorNodes[j * (List_i.size() - 2) + List_i.size() - 3];

      m_err = mk_core()->imesh_instance()->createEnt(iMesh_QUADRILATERAL,
          &qNodes[0], 4, Quads[j * (List_i.size() - 1) + List_i.size() - 2]);
      IBERRCHK(m_err, "Trouble create the quadrilateral elements.");
    }

  }

  //finish creating the quads
  m_err = mk_core()->imesh_instance()->addEntArrToSet(&Quads[0], Quads.size(),
      entityset);
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
  m_err = mk_core()->imesh_instance()->getEntities(root_set, iBase_VERTEX,
      iMesh_POINT, m_Nodes);
  IBERRCHK(m_err, "Trouble get the node list from the mesh entity set.");
  m_err = mk_core()->imesh_instance()->getEntities(root_set, iBase_EDGE,
      iMesh_LINE_SEGMENT, m_Edges);
  IBERRCHK(m_err, "Trouble get the edges from the mesh entity set.");
  m_err = mk_core()->imesh_instance()->getEntities(root_set, iBase_FACE,
      iMesh_QUADRILATERAL, m_Quads);
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

#ifdef HAVE_MESQUITE

  iGeom * ig_inst = mk_core()->igeom_instance(ent->iGeomIndex());

  MeshImprove meshopt(mk_core(),  true, false, false, false, ig_inst);
  meshopt.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);

#endif
  //mk_core()->save_mesh("AfterLaplace.vtk");

  //if there is the parametric space, let Winslow smooth inside the parametric space
  SmoothWinslow(List_i, List_ii, List_j, List_jj, InteriorNodes, Quads,
      mesh_tag, ent);

  mk_core()->save_mesh("AfterWinslow.vtk");
#ifdef HAVE_MESQUITE

  MeshImprove shapesmooth(mk_core(), false, false, true, false, ig_inst);
  shapesmooth.SurfMeshImprove(ent->geom_handle(), entityset, iBase_FACE);

#endif

  //remove the mesh tag
  m_err = mk_core()->imesh_instance()->rmvArrTag(&Quads[0], Quads.size(),
      mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&List_i[0], List_i.size(),
      mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&List_ii[0], List_ii.size(),
      mesh_tag);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  if (List_j.size() > 2)
  {
    m_err = mk_core()->imesh_instance()->rmvArrTag(&List_j[0], List_j.size()
        - 2, mesh_tag);
    IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
    m_err = mk_core()->imesh_instance()->rmvArrTag(&List_jj[0], List_jj.size()
        - 2, mesh_tag);
    IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  }
  if (InteriorNodes.size() > 0)
  {
    m_err = mk_core()->imesh_instance()->rmvArrTag(&InteriorNodes[0],
        InteriorNodes.size(), mesh_tag);
    IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  }

  return 1;
}

//smooth the quadrilateral mesh on the linking surface
void TFIMapping::SmoothWinslow(std::vector<iBase_EntityHandle> &List_i,
    std::vector<iBase_EntityHandle> &List_ii,
    std::vector<iBase_EntityHandle> &List_j,
    std::vector<iBase_EntityHandle> &List_jj,
    std::vector<iBase_EntityHandle> &InteriorNodes, std::vector<
        iBase_EntityHandle> &quads, iBase_TagHandle &taghandle, ModelEnt *ent)
{
  std::vector<std::set<int> > AdjElements;
  std::vector<std::vector<int> > Quads;
  std::vector<std::vector<double> > coords;
  std::vector<bool> isBnd;
  std::vector<iBase_EntityHandle> nodes;
  std::vector<double> weight;

  bool isParameterized = false;

  iGeom::Error g_err = ent->igeom_instance()->isEntParametric(
      ent->geom_handle(), isParameterized);
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
    m_err = mk_core()->imesh_instance()->getVtxCoord(List_i[i], coords[i][0],
        coords[i][1], coords[i][2]);
    IBERRCHK(m_err, "Trouble get the vertex coordinates.");
    nodes[i] = List_i[i];
    isBnd[i] = true;

    m_err = mk_core()->imesh_instance()->getVtxCoord(List_ii[i],
        coords[List_i.size() + i][0], coords[List_i.size() + i][1],
        coords[List_i.size() + i][2]);
    IBERRCHK(m_err, "Trouble get the vertex coordinates.");
    if (isParameterized)
    {
      double uv[2];
      g_err = ent->igeom_instance()->getEntXYZtoUV(ent->geom_handle(),
          coords[List_i.size() + i][0], coords[List_i.size() + i][1],
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
      m_err = mk_core()->imesh_instance()->getVtxCoord(List_j[i], coords[2
          * List_i.size() + i - 1][0], coords[2 * List_i.size() + i - 1][1],
          coords[2 * List_i.size() + i - 1][2]);
      IBERRCHK(m_err, "Trouble get the vertex coordinates.");
      nodes[2 * List_i.size() + i - 1] = List_j[i];
      isBnd[2 * List_i.size() + i - 1] = true;

      m_err = mk_core()->imesh_instance()->getVtxCoord(List_jj[i], coords[2
          * List_i.size() + List_j.size() - 2 + i - 1][0], coords[2
          * List_i.size() + List_j.size() - 2 + i - 1][1], coords[2
          * List_i.size() + List_j.size() - 2 + i - 1][2]);
      IBERRCHK(m_err, "Trouble get the vertex coordinates.");
      nodes[2 * List_i.size() + List_j.size() - 2 + i - 1] = List_jj[i];
      isBnd[2 * List_i.size() + List_j.size() - 2 + i - 1] = true;
    }
  }
  //input the interior nodes
  if (InteriorNodes.size() > 0)
  {
    for (unsigned int i = 0; i < InteriorNodes.size(); i++)
    {
      m_err = mk_core()->imesh_instance()->getVtxCoord(InteriorNodes[i],
          coords[2 * List_i.size() + 2 * (List_j.size() - 2) + i][0], coords[2
              * List_i.size() + 2 * (List_j.size() - 2) + i][1], coords[2
              * List_i.size() + 2 * (List_j.size() - 2) + i][2]);
      IBERRCHK(m_err, "Trouble get the vertex coordinates.");
      nodes[2 * List_i.size() + 2 * (List_j.size() - 2) + i] = InteriorNodes[i];
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
      m_err = mk_core()->imesh_instance()->getEntAdj(nodes[i], iBase_FACE,
          adjEnts);
      IBERRCHK(m_err, "Trouble get the adjacent quads wrt a node.");
      for (unsigned int j = 0; j < adjEnts.size(); j++)
      {
        int index_id = -1;
        m_err = mk_core()->imesh_instance()->getIntData(adjEnts[j], taghandle,
            index_id);
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
    m_err = mk_core()->imesh_instance()->getEntAdj(quads[i], iBase_VERTEX,
        adjEnts);
    IBERRCHK(m_err, "Trouble get the adjacent nodes wrt a quad.");

    assert(adjEnts.size()==4);
    Quads[i].resize(4);

    for (unsigned int j = 0; j < adjEnts.size(); j++)
    {
      int index_id = -1;
      m_err = mk_core()->imesh_instance()->getIntData(adjEnts[j], taghandle,
          index_id);
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
        iGeom::Error g_err = ent->igeom_instance()->getEntClosestPt(
            ent->geom_handle(), coords[i][0], coords[i][1], coords[i][2],
            tmp_coord[0], tmp_coord[1], tmp_coord[2]);
        IBERRCHK(g_err, "Trouble get a closest point on the linking surface.");
      }
      else
      {
        iGeom::Error g_err = ent->igeom_instance()->getEntXYZtoUV(
            ent->geom_handle(), coords[i][0], coords[i][1], tmp_coord[0],
            tmp_coord[1], tmp_coord[2]);
        IBERRCHK(g_err, "Trouble get the xyz from uv.");

      }
      m_err = mk_core()->imesh_instance()->setVtxCoord(nodes[i], tmp_coord[0],
          tmp_coord[1], tmp_coord[2]);
      IBERRCHK(m_err, "Trouble set the new coordinates for nodes.");
    }
  }

  //remove the unnecessary tag after smoothing
  m_err = mk_core()->imesh_instance()->rmvArrTag(&nodes[0], nodes.size(),
      taghandle);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  m_err = mk_core()->imesh_instance()->rmvArrTag(&quads[0], quads.size(),
      taghandle);
  IBERRCHK(m_err, "Trouble remove the tag values from an array of entities.");
  //m_err = mk_core()->imesh_instance()->destroyTag(taghandle, 1);
  //IBERRCHK(m_err, "Trouble destroy a tag.");
}

//****************************************************************************//
// function   : ParameterCalculate
// Date       : Feb 15, 2011
// Description: calculate the parameters for TFI mapping
//***************************************************************************//
int TFIMapping::ParameterCalculate(double &r, double &s, double pt_0s[3],
    double pt_1s[3], double pt_r0[3], double pt_r1[3], double *pts,
    ModelEnt *ent)
{
  // equations P_0s + s*(P_1s - P_0s) = P_r0 + t*(P_r1 - P_r0)

  assert(((fabs(pt_0s[0]-pt_1s[0])>1.0e-5)||(fabs(pt_0s[1]-pt_1s[1]) > 1.0e-5)||(fabs(pt_0s[2]-pt_1s[2]) > 1.0e-5)));
  assert(((fabs(pt_r0[0]-pt_r1[0])>1.0e-5)||(fabs(pt_r0[1]-pt_r1[1]) > 1.0e-5)||(fabs(pt_r0[2]-pt_1s[2]) > 1.0e-5)));
  s = 0;
  r = 0;
  double pt_s[3], pt_r[3];
  bool isParameterized = false;
  iGeom::Error g_err = ent->igeom_instance()->isEntParametric(
      ent->geom_handle(), isParameterized);
  IBERRCHK(g_err, "Trouble check whether the surface is parameterized or not.");
  isParameterized = false;
  if (isParameterized)
  {
    double uv_0s[3] = { 0.0, 0.0, 0.0 }, uv_1s[3] = { 0.0, 0.0, 0.0 },
        uv_r0[3] = { 0.0, 0.0, 0.0 }, uv_r1[3] = { 0.0, 0.0, 0.0 };
    g_err = ent->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), pt_0s[0],
        pt_0s[1], pt_0s[2], uv_0s[0], uv_0s[1]);
    IBERRCHK(g_err, "Trouble get the UV coordinates from xyz.");
    g_err = ent->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), pt_1s[0],
        pt_1s[1], pt_1s[2], uv_1s[0], uv_1s[1]);
    IBERRCHK(g_err, "Trouble get the UV coordinates from xyz.");
    g_err = ent->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), pt_r0[0],
        pt_r0[1], pt_r0[2], uv_r0[0], uv_r0[1]);
    IBERRCHK(g_err, "Trouble get the UV coordinates from xyz.");
    g_err = ent->igeom_instance()->getEntXYZtoUV(ent->geom_handle(), pt_r1[0],
        pt_r1[1], pt_r1[2], uv_r1[0], uv_r1[1]);
    IBERRCHK(g_err, "Trouble get the UV coordinates from xyz.");

    if (!LineLineIntersect(uv_0s, uv_1s, uv_r0, uv_r1, &pt_s[0], &pt_r[0], s, r))
    {
      throw Error(1, "2 3D lines don't intersect at a point.");

    }
    double uv[3];
    for (int i = 0; i < 3; i++)
      uv[i] = (pt_s[i] + pt_r[i]) / 2;

    g_err = ent->igeom_instance()->getEntUVtoXYZ(ent->geom_handle(), uv[0],
        uv[1], pts[0], pts[1], pts[2]);
    IBERRCHK(g_err, "Trouble get the xyz from UV coordinates.");

    return 1;

  }
  else
  {
    if (!LineLineIntersect(pt_0s, pt_1s, pt_r0, pt_r1, &pt_s[0], &pt_r[0], s, r))
    {
      throw Error(1, "2 3D lines don't intersect at a point.");

    }
    for (int i = 0; i < 3; i++)
      pts[i] = (pt_s[i] + pt_r[i]) / 2;

    return 1;
  }

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
bool TFIMapping::LineLineIntersect(double p1[3], double p2[3], double p3[3],
    double p4[3], double *pa, double *pb, double &mua, double &mub)
{
  double p13[3], p43[3], p21[3];
  double d1343, d4321, d1321, d4343, d2121;
  double numer, denom;

  for (int i = 0; i < 3; i++)
  {
    p13[i] = p1[i] - p3[i];
    p43[i] = p4[i] - p3[i];
  }
  if ((fabs(p43[0]) < EPS) && (fabs(p43[1]) < EPS) && (fabs(p43[2]) < EPS))
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

  for (int i = 0; i < 3; i++)
  {
    pa[i] = p1[i] + mua * p21[i];
    pb[i] = p3[i] + mub * p43[i];
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
void TFIMapping::SurfImprove(iBase_EntityHandle surface,
    iBase_EntitySetHandle surfMesh, iBase_EntityType entity_type)
{
#ifdef HAVE_MESQUITE

#endif

}

//****************************************************************************//
// function   : parametricTFI2D
// Date       : Feb 15, 2011
// Description: do the transfinite interpolation in (pt_0s, pt_1s), (pt_r0, pt_r1)
//***************************************************************************//
void TFIMapping::parametricTFI3D(double *pts, double r, double s,
    double pt_0s[3], double pt_1s[3], double pt_r0[3], double pt_r1[3])
{
  //                             pt_r1
  //		  		 |
  //		pt_0s---------Node------------pt_1s
  //		  	 	 |
  //		               pt_r0
  //assert(r>= 0 && r <= 1.0);
  //assert(s>= 0 && s <= 1.0);

  for (int i = 0; i < 3; i++)//interpolate the pt_rs based on pt_r0, pt_r1, pt_0s and pt_1s
    pts[i] = 0.5 * (linear_interpolation(s, pt_r0[i], pt_r1[i])
        + linear_interpolation(r, pt_0s[i], pt_1s[i]));
}

//****************************************************************************//
// function   : linear_interpolation 
// Date       : Feb 15, 2011
// Description: interpolate linearly between x0 and x1
//***************************************************************************//
double TFIMapping::linear_interpolation(double r, double x0, double x1)
{
  //assert(r >=0 && r <= 1.0);
  double pt = (1 - r) * x0 + r * x1;
  return pt;
}

}

