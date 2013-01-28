#include "meshkit/EBMesher.hpp"
#include "meshkit/SCDMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "moab/ReadUtilIface.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/RegisterMeshOp.hpp"

#ifdef HAVE_MOAB
#include "moab/Core.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#include "moab/EntityType.hpp"
#include "moab/HomXform.hpp"
#include "moab/GeomTopoTool.hpp"
#endif

#ifdef HAVE_CGM
#include "CubitDefines.h"
#include "GeometryQueryTool.hpp"
#include "RefFace.hpp"
#endif

#define PI 3.14159265
#define ERRORRF(a) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return false;}}

double rayDir[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

const bool debug_ebmesh = false;
inline bool less_intersect(const IntersectDist d1, const IntersectDist d2) {
  return d1.distance < d2.distance;
}

inline bool equal_intersect(const IntersectDist d1, const IntersectDist d2) {
  return abs(d1.distance - d2.distance) < 10e-7;
}

namespace MeshKit 
{
// static registration of this  mesh scheme
moab::EntityType EBMesher_tps[] = {moab::MBVERTEX, moab::MBHEX, moab::MBMAXTYPE};
const moab::EntityType* EBMesher::output_types()
  { return EBMesher_tps; }

EBMesher::EBMesher(MKCore *mkcore, const MEntVector &me_vec,
                   double size, bool use_geom, int add_layer)
  : MeshScheme(mkcore, me_vec), m_nAddLayer(add_layer),  m_dInputSize(size), m_bUseGeom(use_geom)
{
  m_mesh = mkcore->imesh_instance()->instance();
  m_hRootSet = mkcore->imesh_instance()->getRootSet();
  m_nStatus = OUTSIDE;
  m_iStartHex = 0;
  m_hObbTree = NULL;
  m_hTreeRoot = 0;

#ifdef HAVE_MOAB
  // create tags with MOAB for dense tags
  int outside = 1;
  const void *out = &outside;
  m_elemStatusTag = get_tag("ELEM_STATUS_TAG", 1,
                            MB_TAG_DENSE, MB_TYPE_INTEGER, out);
  
  int length = 1;
  const void *leng = &length;
  m_edgeCutFracLengthTag = get_tag("EDGE_CUT_FRAC_LENGTH_TAG", // # of fractions
                                   3,
                                   //MB_TAG_SPARSE, // using dense for hdf5 exporting performance
                                   MB_TAG_DENSE,
                                   MB_TYPE_INTEGER, leng);

  m_edgeCutFracTag = get_various_length_tag("EDGE_CUT_FRAC_TAG",
                                            //MB_TAG_SPARSE,
                                            MB_TAG_DENSE,
                                            MB_TYPE_DOUBLE);

  int m_id = 1;
  const void *id = &m_id;
  m_volFracLengthTag = get_tag("VOL_FRAC_LENGTH_TAG", 1,
                               MB_TAG_SPARSE, MB_TYPE_INTEGER, id);
  
  m_volFracHandleTag = get_various_length_tag("VOL_FRAC_HANDLE_TAG",
                                          MB_TAG_SPARSE, MB_TYPE_INTEGER);

  m_volFracTag = get_various_length_tag("VOL_FRAC_TAG",
                                        MB_TAG_SPARSE, MB_TYPE_DOUBLE);
  
  // set bounding box size tag
  double bb_default[6] = { 0., 0., 0., 0., 0., 0. };
  m_bbTag = get_tag("BOUNDING_BOX_SIZE", 6,
                    MB_TAG_SPARSE, MB_TYPE_DOUBLE, bb_default);
  
  m_GeomTopoTool = new moab::GeomTopoTool( moab_instance() );
#endif
}

EBMesher::~EBMesher()
{
  delete m_GeomTopoTool;
}

bool EBMesher::can_mesh(ModelEnt *me) 
{
  if (me->dimension() == 3) return true;
  else return false;
}
    
bool EBMesher::add_modelent(ModelEnt *model_ent) 
{
    // make sure this represents geometric volumes
  if (3 != model_ent->dimension()) 
    throw Error(MK_WRONG_DIMENSION, "Wrong dimension entity added to EBMesher.");

  return MeshOp::add_modelent(model_ent);
}

MeshOp *EBMesher::get_scd_mesher()
{
  MeshOpProxy* proxy = NULL;
  proxy = MKCore::meshop_proxy("SCDMesh");
  if (proxy == NULL) throw Error(MK_FAILURE, "Couldn't find a MeshOp capable of producing the given mesh type.");
  return mk_core()->construct_meshop(proxy);
}

void EBMesher::setup_this()
{
  MeshOp *scd_mesher = NULL;
  std::vector<MeshOp*> meshops;
  bool inserted = false;

  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
    ModelEnt *me = (*sit).first;
    bool found = false;
    if (!scd_mesher) scd_mesher = get_scd_mesher();
    meshops.clear();
    me->get_meshops(meshops);
    int n_meshops = meshops.size();

    if (n_meshops > 0) {
      for (int i = 0; i < n_meshops; i++) {
        if (meshops[i] == scd_mesher) {
          found = true;
          break;
        }
      }
    }

    if (!found) { // if no scd meshop
      // add this entity to it, and if first, make sure it's added upstream
      scd_mesher->add_modelent(me);
      if (!inserted) {
        mk_core()->insert_node(scd_mesher, this);
        inserted = true;
      }
    }

    // no need to traverse all geometry for using whole geometry
    if (m_bUseWholeGeom) break;
  }

  SCDMesh *sm = (SCDMesh*) scd_mesher;
  sm->set_interface_scheme(SCDMesh::full);
  sm->set_grid_scheme(SCDMesh::cfMesh);
  
  sm->set_axis_scheme(SCDMesh::cartesian);
  sm->set_box_increase_ratio(m_boxIncrease); // add some extra layer to box

  if (m_bUseWholeGeom) sm->set_geometry_scheme(SCDMesh::all);
  else sm->set_geometry_scheme(SCDMesh::individual);

  // use mesh based geometry in SCDMesh and make tree to get box dimension
  if (m_bUseMeshGeom) {
    sm->use_mesh_geometry(true);
    sm->set_cart_box_min_max(m_minCoord, m_maxCoord, m_boxIncrease);
  }

  // set # of intervals for 3 directions
  std::vector<int> fine_i (m_nDiv[0], 1);
  sm->set_coarse_i_grid(m_nDiv[0]);
  sm->set_fine_i_grid(fine_i);
  std::vector<int> fine_j (m_nDiv[1], 1);
  sm->set_coarse_j_grid(m_nDiv[1]);
  sm->set_fine_j_grid(fine_j);
  std::vector<int> fine_k (m_nDiv[2], 1);
  sm->set_coarse_k_grid(m_nDiv[2]);
  sm->set_fine_k_grid(fine_k);
}

void EBMesher::execute_this() 
{
  iBase_ErrorType  err;
  GraphNode *scd_node = other_node(in_arcs());

#ifdef HAVE_MOAB
  double obb_time = .0;
  double intersection_time = .0;
  double set_info_time = .0;
  time_t time1, time2, time3, time4;
  unsigned long mem1, mem2, mem3, mem4;
  moab_instance()->estimated_memory_use(0, 0, 0, &mem1);
  moab::ErrorCode rval;

  if (debug_ebmesh) {
    rval = mk_core()->moab_instance()->write_mesh("input.vtk");
    MBERRCHK(rval, mk_core()->moab_instance());
  }

  time(&time1);
  if (m_bUseGeom && !m_bUseMeshGeom) { // construct obb tree for geometry surfaces and volumes by GeomTopoTool
    rval = m_GeomTopoTool->construct_obb_trees(m_bUseWholeGeom);
    MBERRCHK(rval, mk_core()->moab_instance());
  }
  time(&time2);
  obb_time += difftime(time2, time1);

  for (MEntSelection::iterator sit = mentSelection.begin(); sit != mentSelection.end(); sit++) {
    // 1. set division
    double box_min_max[6];
    if (m_bUseWholeGeom) {
      static_cast<SCDMesh*> (scd_node)->get_box_dimension(box_min_max, &box_min_max[3]);
    }
    else {
      err = mk_core()->imesh_instance()->getData(reinterpret_cast< iBase_EntityHandle >
                                                 (sit ->first->mesh_handle()), m_bbTag, &box_min_max[0]);
      IBERRCHK(err, "Failed to get hexes.\n");
    }
    set_division(box_min_max, &box_min_max[3]);

    // 2. set or construct obb tree
    time(&time1);
    set_tree_root(sit->first);
    time(&time2);
    obb_time += difftime(time2, time1);

    // 3. find intersected geometry surfaces by rays
    find_intersections();
    time(&time3);
    intersection_time += difftime(time3, time2);
    moab_instance()->estimated_memory_use(0, 0, 0, &mem3);
    
    // 4. set hex status and boundary hex cut fraction info
    set_tag_info();
    time(&time4);
    set_info_time += difftime(time4, time3);
    moab_instance()->estimated_memory_use(0, 0, 0, &mem4);

    if (m_bUseWholeGeom) break;
  }
#endif
  
  if (debug_ebmesh) {
    std::cout << "OBB_tree_construct_time: " << obb_time
              << ", intersection_time: " << intersection_time
              << ", set_info_time: " << set_info_time << std::endl;
    std::cout << "start_memory: " << mem1
              << ", OBB_tree_construct_moemory: " << mem2
              << ", intersection_memory: " << mem3
              << ", set_info_memory: " << mem4
              << std::endl;
  }
}

void EBMesher::set_num_interval(int* n_interval)
{
  for (int i = 0; i < 3; i++) {
    m_nDiv[i] = n_interval[i];
  }
}

#ifdef HAVE_MOAB
void EBMesher::set_tree_root(ModelEnt* me)
{
  moab::ErrorCode rval;
  if (m_bUseGeom) {
    if (m_hObbTree == NULL) {
      m_hObbTree = m_GeomTopoTool->obb_tree();
    }
    if (m_bUseWholeGeom) {
      if (!m_bUseMeshGeom) m_hTreeRoot = m_GeomTopoTool->get_one_vol_root();
    }
    else {
      rval = m_GeomTopoTool->get_root(me->mesh_handle(), m_hTreeRoot);
      MBERRCHK(rval, mk_core()->moab_instance());
    }
  }
  else { // facet data input case
    MBEntityHandle meshset;
    if (m_bUseWholeGeom) meshset = 0;
    else meshset = me->mesh_handle();

    // get triangles
    MBRange tris;
    rval = moab_instance()->get_entities_by_dimension(meshset, 2, tris);
    MBERRCHK(rval, mk_core()->moab_instance());
    
    // make tree
    if (m_hObbTree == NULL) {
      m_hObbTree = new MBOrientedBoxTreeTool( moab_instance() );
    }
    rval = m_hObbTree->build(tris, m_hTreeRoot);
    MBERRCHK(rval, mk_core()->moab_instance());
  }
}

void EBMesher::find_intersections()
{
  std::cout << "starting to find intersections." << std::endl;
  int i;
  //m_vnHexStatus.resize(m_nHex, 1); // initialize all hex as outside
  m_vnHexStatus.resize(m_nHex, 0); // initialize all hex as inside

  // fire rays to 3 directions
  for (i = 0; i < 3; i++) {
    //m_vnEdgeStatus[i].resize(m_nDiv[i]*m_nDiv[(i + 1)%3]*m_nDiv[(i + 2)%3], OUTSIDE);
    m_vnEdgeStatus[i].resize(m_nDiv[i]*m_nDiv[(i + 1)%3]*m_nDiv[(i + 2)%3], INSIDE);
    fire_rays(i);
  }
  std::cout << "finished to find intersections." << std::endl;
}

void EBMesher::export_mesh(const char* file_name, bool separate)
{ 
  // get all hexes
  time_t time1, time2, time3;
  time(&time1);
  int i;
  iBase_ErrorType err;
  std::vector<iBase_EntityHandle> hex_handles;
  err = mk_core()->imesh_instance()->getEntities(m_hRootSet, iBase_REGION,
                                                 iMesh_HEXAHEDRON, hex_handles);
  IBERRCHK(err, "Failed to get hexes.\n");
  /*
  int hex_size = hex_handles.size();
  int* hex_status = new int[hex_size];
  int* hex_status = NULL;
  std::vector<int> hex_status(hex_size);
  err = mk_core()->imesh_instance()->getIntArrData(&hex_handles[0], hex_size,
  m_elemStatusTag, &hex_status[0]);*/

  int error;
  int hex_size = hex_handles.size();
  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, &hex_handles[0],
                      hex_size, m_elemStatusTag,
                      &hex_status, &hex_status_alloc,
                      &hex_status_size, &error);
  IBERRCHK(error, "Failed to get hex status.\n");
  time(&time2);

  if (separate) {
    int n_inside_hex = 0;
    int n_outside_hex = 0;
    int n_boundary_hex = 0;
    int hex_stat;
    (void) hex_stat;
    std::vector<iBase_EntityHandle> insideHex, outsideHex, bndrHex;
    for (i = 0; i < hex_size; i++) {
      if (hex_status[i] == 0) {
        n_inside_hex++;
        hex_stat = 0;
        insideHex.push_back(hex_handles[i]);
      }
      else if (hex_status[i] == 1) {
        n_outside_hex++;
        hex_stat = 1;
        outsideHex.push_back(hex_handles[i]);
      }
      else if (hex_status[i] == 2) {
        n_boundary_hex++;
        hex_stat = 2;
        bndrHex.push_back(hex_handles[i]);
      }
      else {
        throw Error(MK_FAILURE, "Hex element status should be inside/outside/boundary.");
      }
    }
    
    std::cout << "# of exported inside hex:" << n_inside_hex
              << ", # of exported outside hex:" << n_outside_hex
              << ", # of exported boundary hex:" << n_boundary_hex
              << ", geom vol:"
              << n_inside_hex*m_dIntervalSize[0]*m_dIntervalSize[1]*m_dIntervalSize[2]
              << ", vox vol:" << hex_size*m_dIntervalSize[0]*m_dIntervalSize[1]*m_dIntervalSize[2]
              << std::endl;
    time(&time3);

    // save inside/outside/boundary elements separately
    if (n_inside_hex > 0) {
      write_mesh(file_name, 0, &insideHex[0], n_inside_hex);
    }
    if (n_boundary_hex > 0) {
      write_mesh(file_name, 2, &bndrHex[0], n_boundary_hex);
    }
    if (n_outside_hex > 0) {
      write_mesh(file_name, 1, &outsideHex[0], n_outside_hex);
    }

    if (debug_ebmesh) {
      std::cout << "hex_handle_get_time: "
                << difftime(time2, time1)
                << ", separate_write_time: "
                << difftime(time3, time2)
                << std::endl;
    }
  }
  else {
    write_mesh(file_name, 3, &hex_handles[0], hex_size);
    time3 = clock();
    if (debug_ebmesh) {
      std::cout << "hex_handle_get_time: "
                << difftime(time2, time1)
                << ", actual_write_time: "
                << difftime(time3, time2)
                << std::endl;
    }
  }
  delete [] hex_status;
}

void EBMesher::write_mesh(const char* file_name, int type,
                          iBase_EntityHandle* handles, int& n_elem)
{
  time_t time1, time2, time3;
  time(&time1);
  int is_list = 1;
  moab::ErrorCode rval;
  iBase_EntitySetHandle set;
  iBase_ErrorType err;

  err = mk_core()->imesh_instance()->createEntSet(is_list, set);
  IBERRCHK(err, "Couldn't create set.");

  err = mk_core()->imesh_instance()->addEntArrToSet(handles, n_elem, set);
  IBERRCHK(err, "Couldn't add hexes to set.");
  time(&time2);

  std::string out_name;
  std::stringstream ss;
  if (type == 0) ss << "inside_";
  else if (type == 1) ss << "outside_";
  else if (type == 2) ss << "boundary_";
  ss << file_name;
  ss >> out_name;
  
  rval = moab_instance()->write_mesh(out_name.c_str(),
                                     (const moab::EntityHandle*) &set, 1);
  MBERRCHK(rval, mk_core()->moab_instance());
  
  std::cout << "Elements are exported." << std::endl;
  time(&time3);

  if (debug_ebmesh) {
    std::cout << "set_creation_time: "
              << difftime(time2, time1)
              << ", write_time: "
              << difftime(time3, time2)
              << std::endl;
  }
}
#endif

EdgeStatus EBMesher::get_edge_status(const double dP, int& iSkip)
{
  if (m_nStatus == INSIDE) { // previous inside
    if (dP < m_dSecondP) {
      m_nMove = 0;
      iSkip = m_iSecondP;
      return INSIDE;
    }
    else {
      if (is_shared_overlapped_surf(m_iInter - 1)) m_nMove = 1;
      else m_nMove = 2;
        
      // skip shared and overlapped interesections
      while (m_vIntersection[m_iInter].distance -
             m_vIntersection[m_iInter - 1].distance < 1e-7 &&
             (unsigned int) (m_iInter + 1) < m_vIntersection.size()) m_iInter++;        
      
      return BOUNDARY;
    }
  }
  else if (m_nStatus == OUTSIDE) { // previous outside
    if (dP < m_dFirstP) {
      m_nMove = 0;
      return OUTSIDE;
    }
    else {
      if (dP < m_dSecondP) m_nMove = 0;
      else if (is_shared_overlapped_surf(m_iInter - 1)) m_nMove = 1;
      else m_nMove = 2;

      // skip shared and overlapped interesections
      while (m_vIntersection[m_iInter].distance -
             m_vIntersection[m_iInter - 1].distance < 1e-7 &&
             (unsigned int) (m_iInter + 1) < m_vIntersection.size()) m_iInter++;

      return BOUNDARY;
    }
  }
  else if (m_nStatus == BOUNDARY) { // previous boundary
    if (dP < m_dFirstP) {
      m_nMove = 0;
      iSkip = m_iFirstP;
      return OUTSIDE;
    }
    else {
      if (dP < m_dSecondP) {
        m_nMove = 0;
        if (m_prevPnt < m_dFirstP) return BOUNDARY;
        else {
          iSkip = m_iSecondP;
          return INSIDE;
        }
      }
      else {
        if (is_shared_overlapped_surf(m_iInter - 1)) m_nMove = 1;
        else m_nMove = 2;

        // skip shared and overlapped interesections
        while (m_vIntersection[m_iInter].distance -
               m_vIntersection[m_iInter - 1].distance < 1e-7 &&
               (unsigned int) (m_iInter + 1) < m_vIntersection.size()) m_iInter++;

        if (m_prevPnt < m_dSecondP) return BOUNDARY;
        else return OUTSIDE;
      }
    }
  }
  else {
    std::cerr << "Couldn't get edge status." << std::endl;
    return INSIDE;
  }
}

bool EBMesher::set_neighbor_hex_status(int dir, int i, int j, int k)
{
  int iElem;
  int otherDir1 = (dir + 1)%3;
  int otherDir2 = (dir + 2)%3;

  if (dir == 0) { // x coordinate ray
    m_iElem = j*m_nDiv[0]*m_nDiv[1] + i*m_nDiv[0] + k;
    if (!set_hex_status(m_iElem, m_nStatus, dir)) return false;
    iElem = m_iElem - m_nDiv[dir];
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
    iElem -= m_nDiv[dir]*m_nDiv[otherDir1];
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
    iElem += m_nDiv[dir];
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
  }
  else if (dir == 1) { // y coordinate ray
    m_iElem = i*m_nDiv[1]*m_nDiv[0] + k*m_nDiv[0] + j;
    if (!set_hex_status(m_iElem, m_nStatus, dir)) return false;
    iElem = m_iElem - 1;
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
    iElem -= m_nDiv[dir]*m_nDiv[otherDir2];
    if (!set_hex_status(iElem++, m_nStatus, dir)) return false;
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
  }
  else if (dir == 2) { // z coordinate ray
    m_iElem = k*m_nDiv[0]*m_nDiv[1] + j*m_nDiv[0] + i;
    if (!set_hex_status(m_iElem, m_nStatus, dir)) return false;
    iElem = m_iElem - 1;
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
    iElem -= m_nDiv[otherDir1];
    if (!set_hex_status(iElem++, m_nStatus, dir)) return false;
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
  }

  return true;
}

bool EBMesher::set_hex_status(int index, EdgeStatus value, int dir)
{
  if (index < 0 || index > m_nHex - 1) {
    return false;
  }

  if (m_vnHexStatus[index] != 2) {
    if (value == INSIDE) {
      m_vnHexStatus[index] = 0;
    }
    else if (value == OUTSIDE) {
      m_vnHexStatus[index] = 1;
    }
    else if (value == BOUNDARY) {
      m_vnHexStatus[index] = 2;
    }
  }

  return true;
}

bool EBMesher::set_edge_status(int dir)
{
  // set boundary cut information to edge
  std::vector<double> vdCutFraction;
  if (m_nMove == 0) vdCutFraction.push_back((m_curPnt - m_dFirstP)/m_dIntervalSize[dir] - 1.);
  else if (m_nMove == 1) vdCutFraction.push_back(1. - (m_curPnt - m_dSecondP)/m_dIntervalSize[dir]);
  else if (m_nMove == 2) {
    vdCutFraction.push_back(1. - (m_curPnt - m_dSecondP)/m_dIntervalSize[dir]);
    if (m_dFirstP < m_curPnt) {
      vdCutFraction.push_back((m_curPnt - m_dFirstP)/m_dIntervalSize[dir] - 1.);
    }
  }

  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.find(m_iElem);
  if (iter == m_mdCutFraction.end()) { // not exist
    CutFraction cFraction(dir, vdCutFraction);
    m_mdCutFraction[m_iElem] = cFraction;
  }
  else { // exist
    iter->second.add(dir, vdCutFraction);
  }

  //m_nFraction += vdCutFraction.size();

  return true;
}

#ifdef HAVE_MOAB
void EBMesher::set_division()
{
  int i, j;
  moab::CartVect box_center, box_axis1, box_axis2, box_axis3,
    min_cart_box(HUGE_VAL, HUGE_VAL, HUGE_VAL),
    max_cart_box(0., 0., 0.);
  
  moab::ErrorCode rval = m_hObbTree->box(m_hTreeRoot, box_center.array(),
                                     box_axis1.array(), box_axis2.array(),
                                     box_axis3.array());
  MBERRCHK(rval, mk_core()->moab_instance());

  // cartesian box corners
  moab::CartVect corners[8] = {box_center + box_axis1 + box_axis2 + box_axis3,
                               box_center + box_axis1 + box_axis2 - box_axis3,
                               box_center + box_axis1 - box_axis2 + box_axis3,
                               box_center + box_axis1 - box_axis2 - box_axis3,
                               box_center - box_axis1 + box_axis2 + box_axis3,
                               box_center - box_axis1 + box_axis2 - box_axis3,
                               box_center - box_axis1 - box_axis2 + box_axis3,
                               box_center - box_axis1 - box_axis2 - box_axis3};

  // get the max, min cartesian box corners
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 3; j++) {
      if (corners[i][j] < min_cart_box[j]) min_cart_box[j] = corners[i][j];
      if (corners[i][j] > max_cart_box[j]) max_cart_box[j] = corners[i][j];
    }
  }
  moab::CartVect length = max_cart_box - min_cart_box;
  
  // default value is adjusted to large geometry file
  // interval_size_estimate : 2*L*sqrt(2*PI*sqrt(2)/# of tris)
  if (m_dInputSize < 0.) {
    int n_tri;
    rval = moab_instance()->
      get_number_entities_by_dimension(reinterpret_cast<MBEntityHandle>(m_hRootSet),
                                       2, n_tri);
    MBERRCHK(rval, mk_core()->moab_instance());
    
    double box_length_ave = (length[0] + length[1] + length[2])/3.;
    m_dInputSize = 2.*box_length_ave*sqrt(8.886/n_tri);
  }

  for (i = 0; i < 3; i++) {
    m_nDiv[i] = length[i]/m_dInputSize;
    if (m_nDiv[i] < m_nAddLayer) m_nDiv[i] = m_nAddLayer; // make "m_nAddLayer" elements larger than bounding box
    m_dIntervalSize[i] = m_dInputSize;
    //if (m_nDiv[i]*.07 > m_nAddLayer) m_nDiv[i] += m_nDiv[i]*.07;
    //else m_nDiv[i] += m_nAddLayer;
    m_nDiv[i] += m_nAddLayer;
    m_nNode[i] = m_nDiv[i] + 1;
    m_origin_coords[i] = box_center[i] - .5*m_nDiv[i]*m_dIntervalSize[i];
  }

  m_nHex = m_nDiv[0]*m_nDiv[1]*m_nDiv[2];

  std::cout << "# of hex: " << m_nHex << ", interval size: "
            << m_dInputSize << std::endl;

  std::cout << "# of division: " << m_nDiv[0] << ","
            << m_nDiv[1] << "," << m_nDiv[2] << std::endl;
}

int EBMesher::set_division(double* min, double* max)
{
  for (int i = 0; i < 3; i++) {
    m_dIntervalSize[i] = (max[i] - min[i])/m_nDiv[i];
    m_nNode[i] = m_nDiv[i] + 1;
    m_origin_coords[i] = min[i];
  }

  m_nHex = m_nDiv[0]*m_nDiv[1]*m_nDiv[2];

  std::cout << "# of hex: " << m_nHex << ", interval_size: "
            << m_dIntervalSize[0] << ", " << m_dIntervalSize[1] << ", "
            << m_dIntervalSize[2] << std::endl;

  std::cout << "# of division: " << m_nDiv[0] << ","
            << m_nDiv[1] << "," << m_nDiv[2] << std::endl;

  return iBase_SUCCESS;
}

int EBMesher::make_scd_hexes()
{
  // create vertex and hex sequences
  int i;
  double max_coords[3];
  (void) max_coords;
  for (i = 0; i < 3; i++) {
    max_coords[i] = m_origin_coords[i] + m_nDiv[i]*m_dIntervalSize[i];
  }

  moab::HomCoord coord_min(0, 0, 0);
  moab::HomCoord coord_max(m_nDiv[0], m_nDiv[1], m_nDiv[2]);
  moab::EntitySequence* vertex_seq = NULL;
  moab::EntitySequence* cell_seq = NULL;
  moab::EntityHandle vs, cs;
  moab::Core *mbcore = dynamic_cast<moab::Core*>(moab_instance());

  moab::ErrorCode rval = mbcore->create_scd_sequence(coord_min, coord_max, MBVERTEX, 1, vs, vertex_seq);
  MBERRCHK(rval, mk_core()->moab_instance());
  
  mbcore->create_scd_sequence(coord_min, coord_max, MBHEX, 1, cs, cell_seq);
  MBERRCHK(rval, mk_core()->moab_instance());

  moab::HomCoord p2(coord_max.i(), coord_min.j(), coord_min.k());
  moab::HomCoord p3(coord_min.i(), coord_max.j(), coord_min.k()); 
  
  rval = mbcore->add_vsequence(vertex_seq, cell_seq, coord_min, coord_min,
                                        p2, p2, p3, p3);
  MBERRCHK(rval, mk_core()->moab_instance());

  int nTotNode = m_nNode[0]*m_nNode[1]*m_nNode[2];
  int ii, jj, kk;
  MBEntityHandle handle;
  double dumv[3];
  for (i = 0, handle = vs; i < nTotNode; i++, handle++) {
    ii = (i%(m_nNode[0]*m_nNode[1]))%m_nNode[0];
    jj = (i%(m_nNode[0]*m_nNode[1]))/m_nNode[0];
    kk = i/m_nNode[0]/m_nNode[1];
    dumv[0] = m_origin_coords[0] + ii*m_dIntervalSize[0];
    dumv[1] = m_origin_coords[1] + jj*m_dIntervalSize[1];
    dumv[2] = m_origin_coords[2] + kk*m_dIntervalSize[2];
    rval = moab_instance()->set_coords(&handle, 1, dumv);
    MBERRCHK(rval, mk_core()->moab_instance());
  }

  m_iStartHex = moab_instance()->id_from_handle(cs);

  return iBase_SUCCESS;
}

iBase_TagHandle EBMesher::get_tag(const char* name, int size,
                                     unsigned store, MBDataType type,
                                     const void* def_value,
                                     bool create_if_missing) 
{
  MBTag retval = 0;
  /*moab::ErrorCode result = moab_instance()->tag_create(name, size, store, type,
                                                   retval, def_value,
                                                   create_if_missing);*/
  store = store | moab::MB_TAG_CREAT;
  if(!create_if_missing)
    store = store | moab::MB_TAG_EXCL;
  moab::ErrorCode result = moab_instance()->tag_get_handle(name, size, type, retval, store,
                                                      def_value);
  if (create_if_missing && MB_SUCCESS != result) {
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;
  }
  
  return (iBase_TagHandle) retval;
}

iBase_TagHandle EBMesher::get_various_length_tag(const char* name,
                                                    unsigned store, MBDataType type)
{
  MBTag retval = 0;
  store = store | moab::MB_TAG_VARLEN | moab::MB_TAG_CREAT;
  moab::ErrorCode result = moab_instance()->
    tag_get_handle( name, 1, type, retval, store);
  if (MB_SUCCESS != result) {
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;
  }
  
  return (iBase_TagHandle) retval;
}

void EBMesher::set_tag_info()
{
  // get all hexes
  int i, j, k;
  iBase_ErrorType err;
  std::vector<iBase_EntityHandle> hex_handles;
  err = mk_core()->imesh_instance()->getEntities(m_hRootSet, iBase_REGION,
                                                 iMesh_HEXAHEDRON, hex_handles);
  IBERRCHK(err, "Failed to get hexes.\n");

  err = mk_core()->imesh_instance()->setIntArrData(&hex_handles[0], m_nHex,
                                                   m_elemStatusTag,
                                                   &m_vnHexStatus[0]);
  IBERRCHK(err, "Failed to set hex element status data.");

  // set cut fraction info to boundary hexes
  int nBndrHex = m_mdCutFraction.size();
  std::vector<iBase_EntityHandle> hvBndrHex(nBndrHex);
  int* frac_size = new int[nBndrHex];
  int* frac_leng = new int[3*nBndrHex];
  int nFracSize, nTempSize, iHex, nTolFrac = 0;
  int nDoubleSize = sizeof(double);
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.begin();
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  for (i = 0; iter != end_iter; iter++, i++) {
    iHex = iter->first;
    hvBndrHex[i] = hex_handles[iHex];

    nFracSize = 0;
    for (j = 0; j < 3; j++) {
      nTempSize = iter->second.fraction[j].size();
      nFracSize += nTempSize;
      frac_leng[3*i + j] = nTempSize;
    }
    frac_size[i] = nDoubleSize*nFracSize; // sizeof(double)*
    nTolFrac += nFracSize;
  }

  int iFrac = 0;
  double* fractions = new double[nTolFrac];
  std::vector<const void*> frac_data_pointer(nBndrHex);
  iter = m_mdCutFraction.begin();
  for (i = 0; iter != end_iter; iter++, i++) {
    for (j = 0; j < 3; j++) {
      nFracSize = iter->second.fraction[j].size();
      frac_data_pointer[i] = &fractions[iFrac];
      for (k = 0; k < nFracSize; k++) {
        fractions[iFrac++] = iter->second.fraction[j][k];
      }
    }
  }

  err = mk_core()->imesh_instance()->setIntArrData(&hvBndrHex[0], nBndrHex,
                                                   m_edgeCutFracLengthTag,
                                                   frac_leng);
  IBERRCHK(err, "Failed to set cut fraction sizes to hex.");

  moab::ErrorCode rval = moab_instance()->tag_set_by_ptr(reinterpret_cast<MBTag> (m_edgeCutFracTag),
                                                   reinterpret_cast<moab::EntityHandle*> (&hvBndrHex[0]),
                                                   nBndrHex, &frac_data_pointer[0], frac_size);
  MBERRCHK(rval, mk_core()->moab_instance());
  
  delete [] frac_size;
  delete [] frac_leng;
  delete [] fractions;
}

void EBMesher::fire_rays(int dir)
{
  // ray fire
  int i, j, k, l, index[3];
  double tolerance = 1e-12;
  double rayLength = m_nDiv[dir]*m_dIntervalSize[dir];
  int iNodeStart, iNodeEnd, nIntersect, nNodeSlice;
  double startPnt[3], endPnt[3];
  int otherDir1 = (dir + 1)%3;
  int otherDir2 = (dir + 2)%3;
 // harmless as this expression (does nothing) so that a compiler sees it is used.
  (void) iNodeEnd;
  (void) nNodeSlice;

  for (j = 1; j < m_nNode[otherDir2] - 1; j++) {
    for (i = 1; i < m_nNode[otherDir1] - 1; i++) {

      // get ray start and end points
      if (dir == 0) { // x coordinate ray
        iNodeStart = j*m_nNode[dir]*m_nNode[otherDir1] + i*m_nNode[dir];
        iNodeEnd = iNodeStart + m_nNode[dir] - 1;
        nNodeSlice = 1;
        index[0] = 0;
        index[1] = i;
        index[2] = j;
      }
      else if (dir == 1) { // y coordinate ray
        iNodeStart = i*m_nNode[otherDir2]*m_nNode[dir] + j;
        iNodeEnd = iNodeStart + m_nNode[otherDir2]*(m_nNode[dir] - 1);
        nNodeSlice = m_nNode[otherDir2];
        index[0] = j;
        index[1] = 0;
        index[2] = i;
      }
      else if (dir == 2) { // z coordinate ray
        iNodeStart = j*m_nNode[otherDir1] + i;
        iNodeEnd = iNodeStart + m_nNode[otherDir1]*m_nNode[otherDir2]*(m_nNode[dir] - 1);
        nNodeSlice = m_nNode[otherDir1]*m_nNode[otherDir2];
        index[0] = i;
        index[1] = j;
        index[2] = 0;
      }
      else IBERRCHK(iBase_FAILURE, "Ray direction should be 0 to 2.");
      
      for (l = 0; l < 3; l++) {
        if (l == dir) {
          startPnt[l] = m_origin_coords[l];
          endPnt[l] = m_origin_coords[l] + m_nDiv[dir]*m_dIntervalSize[l];
        }
        else {
          startPnt[l] = m_origin_coords[l] + index[l]*m_dIntervalSize[l];
          endPnt[l] = startPnt[l];
        }
      }
      
      k = 0;
      m_iInter = 0;
      if (!fire_ray(nIntersect, startPnt, endPnt, tolerance, dir, rayLength)) {
        IBERRCHK(iBase_FAILURE, "Couldn't fire ray.");
      }
      
      if (nIntersect > 0) {
        m_iFirstP = m_vIntersection[m_iInter].distance/m_dIntervalSize[dir];
        m_dFirstP = startPnt[dir] + m_vIntersection[m_iInter++].distance;
        m_iSecondP = m_vIntersection[m_iInter].distance/m_dIntervalSize[dir];
        m_dSecondP = startPnt[dir] + m_vIntersection[m_iInter++].distance;

        // set outside before the first hit
        for (k = 0; k < m_iFirstP; k++) {
          m_nStatus = OUTSIDE;
          
          if (!set_neighbor_hex_status(dir, i, j, k)) {
            IBERRCHK(iBase_FAILURE, "Couldn't set neighbor hex status.");
          }
        }
        
        for (; k < m_nNode[dir] - 1;) {
          int i_skip = 0;
          m_curPnt = startPnt[dir] + (k + 1)*m_dIntervalSize[dir];
          EdgeStatus preStat = m_nStatus;
          m_nStatus = get_edge_status(m_curPnt, i_skip);
          m_vnEdgeStatus[dir][i*m_nDiv[dir] + j*m_nDiv[dir]*m_nDiv[otherDir1] + k] = m_nStatus;

          // set status of all hexes sharing the edge
          if (!set_neighbor_hex_status(dir, i, j, k)) {
            IBERRCHK(iBase_FAILURE, "Couldn't set neighbor hex status.");
          }

          if (m_nMove > 0) {
            if (m_iInter < nIntersect) {
              if (!move_intersections(dir, nIntersect, startPnt)) {
                IBERRCHK(iBase_FAILURE, "Couldn't move intersections.");
              }
            }
            else {
              m_nMove = 1;
              if (m_nStatus == BOUNDARY && !set_edge_status(dir)) {
                IBERRCHK(iBase_FAILURE, "Couldn't set edge status.");
              }
              k++;
              break; // rest is all outside
            }
          }
          else if (m_nStatus == BOUNDARY && !set_edge_status(dir)) {
            IBERRCHK(iBase_FAILURE, "Couldn't set edge status.");
          }
          else if (m_nStatus == OUTSIDE && preStat == BOUNDARY) { // set outside until next hit
            k++;
            for (; k < i_skip; k++) {
              m_nStatus = OUTSIDE;
              if (!set_neighbor_hex_status(dir, i, j, k)) {
                IBERRCHK(iBase_FAILURE, "Couldn't set neighbor hex status.");
              }
            }
          }

          // set cut-cell edge status
          if (i_skip > 0) {
            m_prevPnt = startPnt[dir] + i_skip*m_dIntervalSize[dir];
            k = i_skip;
          }
          else {
            m_prevPnt = m_curPnt;
            k++;
          }
        }
      }
        
      // the rest are all outside
      for (; k < m_nNode[dir] - 1; k++) {
        m_nStatus = OUTSIDE;
        if (!set_neighbor_hex_status(dir, i, j, k)) {
          IBERRCHK(iBase_FAILURE, "Couldn't set neighbor hex status.");
        }
      }
    }
  }
}

bool EBMesher::fire_ray(int& nIntersect, double startPnt[3],
                      double endPnt[3], double tol, int dir,
                      double rayLength)
{
  m_vIntersection.clear();
  m_vhInterSurf.clear();
  m_vhInterFacet.clear();
  m_mhOverlappedSurf.clear();
  std::vector<double> temp_intersects;
  moab::ErrorCode rVal;
  if (m_bUseGeom) { // geometry input
    rVal = m_hObbTree->ray_intersect_sets(temp_intersects, m_vhInterSurf,
                                          m_vhInterFacet, m_hTreeRoot, tol,
                                          -1, startPnt, rayDir[dir], &rayLength);
  }
  else { // facet input
    std::vector<moab::EntityHandle> dum_facets_out;
    rVal = m_hObbTree->ray_intersect_triangles(temp_intersects, dum_facets_out,
                                               m_hTreeRoot, tol,
                                               startPnt, rayDir[dir], &rayLength);
  }
  
  nIntersect = temp_intersects.size();
  if (MB_SUCCESS != rVal) {
    std::cerr << "Failed : ray-triangle intersection." << std::endl;
    return false;
  }
  
  if (nIntersect > 0) {
    m_vIntersection.resize(nIntersect);
    
    for (int l = 0; l < nIntersect; l++) {
      IntersectDist temp_inter_dist(temp_intersects[l], l);
      m_vIntersection[l] = temp_inter_dist;
    }
    std::sort(m_vIntersection.begin(), m_vIntersection.end(), less_intersect);

    if (nIntersect > 1) {
      bool bMoveOnce;
      m_nIteration = 0;
      m_iOverlap = 0;
      
      if (m_bUseGeom) { // if there are geometry info
        remove_intersection_duplicates();
      }
      else if (is_ray_move_and_set_overlap_surf(bMoveOnce)) { // facet geom case
        if (!move_ray(nIntersect, startPnt, endPnt, tol, dir, bMoveOnce)) {
          std::cerr << "Number of Intersection between edges and ray should be even." << std::endl;
          return false;
        }
      }
      nIntersect = m_vIntersection.size();
    }
  }
  
  return true;
}

bool EBMesher::move_intersections(int n_dir, int n_inter, double start_pnt[3])
{
  if (m_nMove > 0) {
    if (m_iInter < n_inter) {
      if (m_nMove == 1) {
        do {
          if (m_nStatus == BOUNDARY && !set_edge_status(n_dir)) return false;
          m_iFirstP = m_iSecondP;
          m_dFirstP = m_dSecondP;
          m_iSecondP = m_vIntersection[m_iInter].distance/m_dIntervalSize[n_dir];
          m_dSecondP = start_pnt[n_dir] + m_vIntersection[m_iInter++].distance;
        }
        while (m_dSecondP < m_curPnt && m_iInter < n_inter);
      }
      else if (m_nMove == 2) {
        do {
          m_iFirstP = m_vIntersection[m_iInter].distance/m_dIntervalSize[n_dir];
          m_dFirstP = start_pnt[n_dir] + m_vIntersection[m_iInter++].distance;
          if (m_nStatus == BOUNDARY && !set_edge_status(n_dir)) return false;
          if (m_iInter < n_inter) {
            m_iSecondP = m_vIntersection[m_iInter].distance/m_dIntervalSize[n_dir];
            m_dSecondP = start_pnt[n_dir] + m_vIntersection[m_iInter++].distance;
          }
          else {
            m_iSecondP = m_iFirstP;
            m_dSecondP = m_dFirstP;
          }
        }
        while (m_dSecondP < m_curPnt && m_iInter < n_inter);
      }
    }
  }

  return true;
}

bool EBMesher::is_shared_overlapped_surf(int index)
{
  int nParent;
  iBase_ErrorType err;
  moab::EntityHandle hSurf;
  if (m_bUseGeom) {
    hSurf = m_vhInterSurf[m_vIntersection[index].index];
    err = mk_core()->imesh_instance()->getNumPrnt(reinterpret_cast<iBase_EntitySetHandle> (hSurf),
                                                  1, nParent);
    IBERRCHK(err, "Failed to get number of surface parents.\n");

    if (nParent > 1) return true;
  }
  else hSurf = m_vIntersection[index].index;

  return m_mhOverlappedSurf.count(hSurf) > 0;
}

void EBMesher::get_grid_and_edges_techX(double* boxMin, double* boxMax, int* nDiv,
                                      std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellSurfEdge,
                                      std::vector<int>& rvnInsideCell, bool isCornerExterior)
{
  int i, j, err, ii, jj, kk, iHex;
  for (i = 0; i < 3; i++) {
    boxMin[i] = m_origin_coords[i];
    boxMax[i] = m_origin_coords[i] + m_dIntervalSize[i]*m_nDiv[i];
    nDiv[i] = m_nDiv[i];
  }
  
  // get all hexes
  std::vector<iBase_EntityHandle> hex_handles;
  err = mk_core()->imesh_instance()->getEntities(m_hRootSet, iBase_REGION,
                                                 iMesh_HEXAHEDRON, hex_handles);
  IBERRCHK(err, "Failed to get hexes.\n");
  
  // get hex status
  int hex_size = hex_handles.size();
  /*
  //int* hex_status = new int[hex_size];
  //int* hex_status;
  //std::vector<int> hex_status(hex_size);
  int temp = 0;
  int* hex_status = &temp;
  err = mk_core()->imesh_instance()->getIntArrData(&hex_handles[0], hex_size,
                                                   m_elemStatusTag, hex_status);
  */
  int error;
  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, &hex_handles[0],
                      hex_size, m_elemStatusTag,
                      &hex_status, &hex_status_alloc,
                      &hex_status_size, &error);
  IBERRCHK(error, "Failed to get hex status.");

  // get inside and boundary hexes
  int nInside = 0;
  int nOutside = 0;
  int nBoundary = 0;
  rvnInsideCell.clear();
  for (i = 0; i < hex_size; i++) {
    if (hex_status[i] == 0) { // if inside
      iHex = moab_instance()->id_from_handle(reinterpret_cast<moab::EntityHandle>
                                             (hex_handles[i])) - m_iStartHex;
      rvnInsideCell.push_back((iHex%(m_nDiv[0]*m_nDiv[1]))%m_nDiv[0]);
      rvnInsideCell.push_back((iHex%(m_nDiv[0]*m_nDiv[1]))/m_nDiv[0]);
      rvnInsideCell.push_back(iHex/m_nDiv[0]/m_nDiv[1]);
      nInside++;
    }
    else if (hex_status[i] == 1) nOutside++;
    else if (hex_status[i] == 2) nBoundary++;
    else IBERRCHK(err, "Element status should be one of inside/outside/boundary."); 
  }
  std::cout << "# of inside, outside, boundary elem: " << nInside
            << ", " << nOutside << ", " << nBoundary << std::endl;
  delete [] hex_status;

  // get cut-cell fractions
  double dFrac[4];
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.begin();
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  
  for (; iter != end_iter; iter++) { // for each cut-cell
    // get i, j, k index from handle
    iHex = iter->first;
    ii = (iHex%(m_nDiv[0]*m_nDiv[1]))%m_nDiv[0];
    jj = (iHex%(m_nDiv[0]*m_nDiv[1]))/m_nDiv[0];
    kk = iHex/m_nDiv[0]/m_nDiv[1];

    // surface 1
    CutFraction cutFrac = iter->second;
    if (cutFrac.fraction[1].size() > 0) {
      dFrac[0] = cutFrac.fraction[1][0];
      if (dFrac[0] < 0.) dFrac[0] *= -1.;
    }
    else dFrac[0] = -1.;
    if (cutFrac.fraction[2].size() > 0) {
      dFrac[3] = cutFrac.fraction[2][0];
      if (dFrac[3] < 0.) dFrac[3] *= -1.;
    }
    else dFrac[3] = -1.;
    dFrac[1] = get_edge_fraction(iHex + m_nDiv[0], 2);
    dFrac[2] = get_edge_fraction(iHex + m_nDiv[0]*m_nDiv[1], 1);
    
    if (dFrac[0] > 0. || dFrac[1] > 0. || dFrac[2] > 0. || dFrac[3] > 0.) { // if surface is cut
      CutCellSurfEdgeKey surfKey(ii, jj, kk, 0);
      std::vector<double> edges;
      if (dFrac[0] < 0.) dFrac[0] = get_uncut_edge_fraction(ii, jj, kk, 1);
      if (dFrac[1] < 0.) dFrac[1] = get_uncut_edge_fraction(ii, jj + 1, kk, 2);
      if (dFrac[2] < 0.) dFrac[2] = get_uncut_edge_fraction(ii, jj, kk + 1, 1);
      if (dFrac[3] < 0.) dFrac[3] = get_uncut_edge_fraction(ii, jj, kk, 2);
      for (j = 0; j < 4; j++) edges.push_back(dFrac[j]);
      rmdCutCellSurfEdge[surfKey] = edges;
    }

    // surface 2
    cutFrac = iter->second;
    if (cutFrac.fraction[0].size() > 0) {
      dFrac[0] = cutFrac.fraction[0][0];
      if (dFrac[0] < 0.) dFrac[0] *= -1.;
    }
    else dFrac[0] = -1.;
    if (cutFrac.fraction[1].size() > 0) {
      dFrac[3] = cutFrac.fraction[1][0];
      if (dFrac[3] < 0.) dFrac[3] *= -1.;
    }
    else dFrac[3] = -1.;
    dFrac[1] = get_edge_fraction(iHex + 1, 2);
    dFrac[2] = get_edge_fraction(iHex + m_nDiv[0]*m_nDiv[1], 0);
    
    if (dFrac[0] > 0. || dFrac[1] > 0. || dFrac[2] > 0. || dFrac[3] > 0.) { // if surface is cut
      CutCellSurfEdgeKey surfKey(ii, jj, kk, 1);
      std::vector<double> edges;
      if (dFrac[0] < 0.) dFrac[0] = get_uncut_edge_fraction(ii, jj, kk, 0);
      if (dFrac[1] < 0.) dFrac[1] = get_uncut_edge_fraction(ii + 1, jj, kk, 2);
      if (dFrac[2] < 0.) dFrac[2] = get_uncut_edge_fraction(ii, jj, kk + 1, 0);
      if (dFrac[3] < 0.) dFrac[3] = get_uncut_edge_fraction(ii, jj, kk, 2);
      for (j = 0; j < 4; j++) edges.push_back(dFrac[j]);
      rmdCutCellSurfEdge[surfKey] = edges;
    }
    
    // surface 3
    cutFrac = iter->second;
    if (cutFrac.fraction[0].size() > 0) {
      dFrac[0] = cutFrac.fraction[0][0];
      if (dFrac[0] < 0.) dFrac[0] *= -1.;
    }
    else dFrac[0] = -1.;
    if (cutFrac.fraction[1].size() > 0) {
      dFrac[3] = cutFrac.fraction[1][0];
      if (dFrac[3] < 0.) dFrac[3] *= -1.;
    }
    else dFrac[3] = -1.;
    dFrac[1] = get_edge_fraction(iHex + 1, 1);
    dFrac[2] = get_edge_fraction(iHex + m_nDiv[0], 0);

    if (dFrac[0] > 0. || dFrac[1] > 0. || dFrac[2] > 0. || dFrac[3] > 0.) { // if surface is cut
      CutCellSurfEdgeKey surfKey(ii, jj, kk, 2);
      std::vector<double> edges;
      if (dFrac[0] < 0.) dFrac[0] = get_uncut_edge_fraction(ii, jj, kk, 0);
      if (dFrac[1] < 0.) dFrac[1] = get_uncut_edge_fraction(ii + 1, jj, kk, 1);
      if (dFrac[2] < 0.) dFrac[2] = get_uncut_edge_fraction(ii, jj + 1, kk, 0);
      if (dFrac[3] < 0.) dFrac[3] = get_uncut_edge_fraction(ii, jj, kk, 1);
      for (j = 0; j < 4; j++) edges.push_back(dFrac[j]);
      rmdCutCellSurfEdge[surfKey] = edges;
    }
  }
  
  if (debug_ebmesh) export_fraction_edges(rmdCutCellSurfEdge);
}

bool EBMesher::get_grid_and_edges(double* boxMin, double* boxMax, int* nDiv,
                                     std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellEdge,
                                     std::vector<int>& rvnInsideCell, bool isCornerExterior)
{
  int i, ii, jj, kk, iHex;
  for (i = 0; i < 3; i++) {
    boxMin[i] = m_origin_coords[i];
    boxMax[i] = m_origin_coords[i] + m_dIntervalSize[i]*m_nDiv[i];
    nDiv[i] = m_nDiv[i];
  }
  
  if (!get_inside_hex(rvnInsideCell)) return false;

  // get cut-cell fractions
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.begin();
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  
  for (; iter != end_iter; iter++) { // for each cut-cell
    // get i, j, k index from handle
    iHex = iter->first;
    ii = (iHex%(m_nDiv[0]*m_nDiv[1]))%m_nDiv[0];
    jj = (iHex%(m_nDiv[0]*m_nDiv[1]))/m_nDiv[0];
    kk = iHex/m_nDiv[0]/m_nDiv[1];

    for (i = 0; i < 3; i++) {
      std::vector<double> fractions(iter->second.fraction[i]);
      CutCellSurfEdgeKey edgeKey(ii, jj, kk, i);
      rmdCutCellEdge[edgeKey] = fractions;
    }
  }

  if (debug_ebmesh) return export_fraction_points(rmdCutCellEdge);
 
  return true;
}

bool EBMesher::get_inside_hex(std::vector<int>& rvnInsideCell)
{
  int i, err, iHex;
  
  // get all hexes
  std::vector<iBase_EntityHandle> hex_handles;
  err = mk_core()->imesh_instance()->getEntities(m_hRootSet, iBase_REGION,
                                                 iMesh_HEXAHEDRON, hex_handles);
  IBERRCHK(err, "Failed to get hexes.\n");
  
  // get hex status
  int hex_size = hex_handles.size();
  /*
  //int* hex_status = new int[hex_size];
  //int* hex_status;
  std::vector<int> hex_status(hex_size);
  err = mk_core()->imesh_instance()->getIntArrData(&hex_handles[0], hex_size,
                                                   m_elemStatusTag, &hex_status[0]);
  */
  int error;
  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, &hex_handles[0],
                      hex_size, m_elemStatusTag,
                      &hex_status, &hex_status_alloc,
                      &hex_status_size, &error);
  ERRORRF("Failed to get hex status.\n");

  // get inside and boundary hexes
  int nInside = 0;
  int nOutside = 0;
  int nBoundary = 0;
  rvnInsideCell.clear();
  for (i = 0; i < hex_size; i++) {
    if (hex_status[i] == 0) { // if inside
      iHex = moab_instance()->id_from_handle(reinterpret_cast<moab::EntityHandle>
                                             (hex_handles[i])) - m_iStartHex;
      rvnInsideCell.push_back((iHex%(m_nDiv[0]*m_nDiv[1]))%m_nDiv[0]);
      rvnInsideCell.push_back((iHex%(m_nDiv[0]*m_nDiv[1]))/m_nDiv[0]);
      rvnInsideCell.push_back(iHex/m_nDiv[0]/m_nDiv[1]);
      nInside++;
    }
    else if (hex_status[i] == 1) nOutside++;
    else if (hex_status[i] == 2) nBoundary++;
    else ERRORRF("Element status should be one of inside/outside/boundary.\n"); 
  }
  std::cout << "# of inside, outside, boundary elem: " << nInside
            << ", " << nOutside << ", " << nBoundary << std::endl;

  delete [] hex_status;
  return true;
}

bool EBMesher::export_fraction_edges(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellSurfEdge)
{
  // export fractions as edge
  double curPnt[3], ePnt[6];
  std::vector<iBase_EntityHandle> edge_handles;
  std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >::iterator itr = mdCutCellSurfEdge.begin();
  std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >::iterator e_itr = mdCutCellSurfEdge.end();
  for (; itr != e_itr; itr++) {
    curPnt[0] = m_origin_coords[0] + itr->first.i*m_dIntervalSize[0];
    curPnt[1] = m_origin_coords[1] + itr->first.j*m_dIntervalSize[1];
    curPnt[2] = m_origin_coords[2] + itr->first.k*m_dIntervalSize[2];
    std::vector<double> edges = itr->second;
    
    if (itr->first.l == 0) {
      if (edges[0] > 0. && edges[0] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1] + edges[0]*m_dIntervalSize[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[1] > 0. && edges[1] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1] + m_dIntervalSize[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2] + edges[1]*m_dIntervalSize[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[2] > 0. && edges[2] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2] + m_dIntervalSize[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1] + edges[2]*m_dIntervalSize[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[3] > 0. && edges[3] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2] + edges[3]*m_dIntervalSize[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
    }
    if (itr->first.l == 1) {
      if (edges[0] > 0. && edges[0] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0] + edges[0]*m_dIntervalSize[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[1] > 0. && edges[1] < 1.) {
        ePnt[0] = curPnt[0] + m_dIntervalSize[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2] + edges[1]*m_dIntervalSize[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[2] > 0. && edges[2] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2] + m_dIntervalSize[2];
        ePnt[3] = ePnt[0] + edges[2]*m_dIntervalSize[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[3] > 0. && edges[3] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2] + edges[3]*m_dIntervalSize[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
    }
    if (itr->first.l == 2) {
      if (edges[0] > 0. && edges[0] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0] + edges[0]*m_dIntervalSize[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[1] > 0. && edges[1] < 1.) {
        ePnt[0] = curPnt[0] + m_dIntervalSize[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1] + edges[1]*m_dIntervalSize[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[2] > 0. && edges[2] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1] + m_dIntervalSize[2];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0] + edges[2]*m_dIntervalSize[0];
        ePnt[4] = ePnt[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
      if (edges[3] > 0. && edges[3] < 1.) {
        ePnt[0] = curPnt[0];
        ePnt[1] = curPnt[1];
        ePnt[2] = curPnt[2];
        ePnt[3] = ePnt[0];
        ePnt[4] = ePnt[1] + edges[3]*m_dIntervalSize[1];
        ePnt[5] = ePnt[2];
        if (!make_edge(ePnt, edge_handles)) return false;
      }
    }
  }

  int is_list = 1;
  iBase_EntitySetHandle set;
  iBase_ErrorType err;
  err = mk_core()->imesh_instance()->createEntSet(is_list, set);
  IBERRCHK(err, "Couldn't create set.");
  
  err = mk_core()->imesh_instance()->addEntArrToSet(&edge_handles[0],
                                                    edge_handles.size(), set);
  IBERRCHK(err, "Couldn't add edges to set.");
  
  moab::ErrorCode rval = moab_instance()->write_mesh("edges.vtk",
                                                     (const moab::EntityHandle*) &set, 1);
  MBERRCHK(rval, mk_core()->moab_instance());

  return true;
}

bool EBMesher::export_fraction_points(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellEdge)
{
  // export fractions as edge
  int i, j, dir, nFrc;
  iBase_ErrorType err;
  double curPnt[3], fracPnt[3], frac;
  iBase_EntityHandle v_handle;
  std::vector<iBase_EntityHandle> vertex_handles;
  std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >::iterator itr = mdCutCellEdge.begin();
  std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >::iterator e_itr = mdCutCellEdge.end();
  for (; itr != e_itr; itr++) {
    curPnt[0] = m_origin_coords[0] + itr->first.i*m_dIntervalSize[0];
    curPnt[1] = m_origin_coords[1] + itr->first.j*m_dIntervalSize[1];
    curPnt[2] = m_origin_coords[2] + itr->first.k*m_dIntervalSize[2];
    dir = itr->first.l;
    nFrc = itr->second.size();
    
    for (i = 0; i < nFrc; i++) {
      frac = itr->second[i];
      if (frac > 0. && frac < 1.) {
        for (j = 0; j < 3; j++) {
          if (j == dir) fracPnt[j] = curPnt[j] + frac*m_dIntervalSize[j];
          else fracPnt[j] = curPnt[j];
        }
        err = mk_core()->imesh_instance()->createVtx(fracPnt[0], fracPnt[1],
                                                     fracPnt[2], v_handle);
        IBERRCHK(err, "Couldn't create vertex.");
        vertex_handles.push_back(v_handle);
      }
    }
  }
    
  int is_list = 1;
  moab::ErrorCode result;
  iBase_EntitySetHandle set;
  err = mk_core()->imesh_instance()->createEntSet(is_list, set);
  IBERRCHK(err, "Couldn't create set.");
  
  err = mk_core()->imesh_instance()->addEntArrToSet(&vertex_handles[0],
                                                    vertex_handles.size(), set);
  IBERRCHK(err, "Couldn't add vertices to set.");
  
  result = moab_instance()->write_mesh("frac_vertices.vtk",
                                       (const moab::EntityHandle*) &set, 1);
  if (MB_SUCCESS != result) {
    std::cerr << "Failed to write fraction vertices." << std::endl;
    return false;
  }

  return true;
}

bool EBMesher::make_edge(double ePnt[6], std::vector<iBase_EntityHandle>& edge_handles)
{
  iBase_ErrorType err;
  iBase_EntityHandle vertex_handle[2], edge_handle;
  iBase_EntityHandle* pVertexHandle = &vertex_handle[0];
  err = mk_core()->imesh_instance()->createVtxArr(2, iBase_INTERLEAVED,
                                                  ePnt, pVertexHandle);
  IBERRCHK(err, "Failed to create vertices.\n");

  err = mk_core()->imesh_instance()->createEnt(iMesh_LINE_SEGMENT,
                                               &vertex_handle[0], 2,
                                               edge_handle);
  IBERRCHK(err, "Failed to create edge.\n");

  edge_handles.push_back(edge_handle);

  return true;
}

double EBMesher::get_edge_fraction(int idHex, int dir)
{
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.find(idHex);
  if (iter != end_iter && iter->second.fraction[dir].size() > 0) {
    double frac = iter->second.fraction[dir][0];
    if (frac < 0.) frac *= -1.;
    return frac;
  }
  else return -1.;
}
#endif

double EBMesher::get_uncut_edge_fraction(int i, int j, int k, int dir)
{
  int iEdge;

  if (dir == 0) {
    if (j == m_nDiv[1] || k == m_nDiv[2]) return 0.; // return outside
    iEdge = k*m_nDiv[0]*m_nDiv[1] + j*m_nDiv[0] + i;
  }
  else if (dir == 1) {
    if (i == m_nDiv[0] || k == m_nDiv[2]) return 0.; // return outside
    iEdge = i*m_nDiv[1]*m_nDiv[2] + k*m_nDiv[1] + j;
  }
  else if (dir == 2) {
    if (i == m_nDiv[0] || j == m_nDiv[1]) return 0.; // return outside
    iEdge = j*m_nDiv[0]*m_nDiv[2] + i*m_nDiv[2] + k;
  }
  else return -1.;

  EdgeStatus edgeStatus = m_vnEdgeStatus[dir][iEdge];
  if (edgeStatus == INSIDE) return 1.;
  else if (edgeStatus == OUTSIDE) return 0.;
  else return -1.;
}

#ifdef HAVE_MOAB
bool EBMesher::move_ray(int& nIntersect, double* startPnt, double* endPnt,
                      double tol, int dir, bool bSpecialCase)
{
  //int i, nIteration;
  int i;
  bool bMove = true;

  if (bSpecialCase) m_iOverlap = 0;

  while (bMove) {
    if (m_nIteration > 10) {
      if (bSpecialCase) return true; // special case
      else if (m_bUseGeom) { // shared/overlapped, faceting case
        m_mhOverlappedSurf[m_iOverlap] = 0;
        m_mhOverlappedSurf[m_iOverlap + 1] = 0;
        return true;
      }
      else {
        return false;
      }
    }

    for (i = 0; i < 3; i++) {
      if (i != dir) {
        startPnt[i] += tol;
        endPnt[i] += tol;
      }
    }

    MBCartVect ray(endPnt[0] - startPnt[0], endPnt[1] - startPnt[1], endPnt[2] - startPnt[2]);
    double rayLength = ray.length();
    moab::ErrorCode rVal;
    ray.normalize();
    m_vIntersection.clear();
    m_vhInterSurf.clear();
    m_vhInterFacet.clear();
    
    std::vector<double> temp_intersects;
    if (m_bUseGeom) {
      rVal = m_hObbTree->ray_intersect_sets(temp_intersects, m_vhInterSurf,
                                            m_vhInterFacet, m_hTreeRoot, tol,
                                            -1, startPnt, ray.array(), &rayLength);
    }
    else { // facet input
      std::vector<moab::EntityHandle> dum_facets_out;
      rVal = m_hObbTree->ray_intersect_triangles(temp_intersects, dum_facets_out,
                                                 m_hTreeRoot, tol,
                                                 startPnt, ray.array(), &rayLength);
      m_vhInterSurf.resize(temp_intersects.size());
    }

    if (MB_SUCCESS != rVal) {
      std::cerr << "Failed : ray-triangle intersection." << std::endl;
      return false;
    }

    nIntersect = temp_intersects.size();
    m_vIntersection.resize(nIntersect);
    for (i = 0; i < nIntersect; i++) {
      IntersectDist temp_inter_dist(temp_intersects[i], i);
      m_vIntersection[i] = temp_inter_dist;
    }
    std::sort(m_vIntersection.begin(), m_vIntersection.end(), less_intersect);
    bMove = is_ray_move_and_set_overlap_surf(bSpecialCase);
    m_nIteration++;
  }

  std::cout << "ray is moved successfully." << std::endl;

  return true;
}
#endif

bool EBMesher::is_ray_move_and_set_overlap_surf(bool& bSpecialCase)
{
  int nInter = m_vIntersection.size() - 1;

  while (m_iOverlap < nInter) {
    if (m_vIntersection[m_iOverlap + 1].distance - m_vIntersection[m_iOverlap].distance < 1e-7) {
      if (m_bUseGeom) {
        moab::EntityHandle h1 = m_vhInterSurf[m_vIntersection[m_iOverlap].index];
        moab::EntityHandle h2 = m_vhInterSurf[m_vIntersection[m_iOverlap + 1].index];
        
        if (h1 == h2) { // remove too close case
          bSpecialCase = false;
          return true;
        }
        else if (m_nIteration < 10) { // when ray intesect shared edge by 2 triangles
          bSpecialCase = true;
          return true;
        }
        else { // overlapped surfaces
          m_mhOverlappedSurf[h1] = 0;
          m_mhOverlappedSurf[h2] = 0;
          m_nIteration = 0;
          m_iOverlap++;
        }
      }
      else {
        if (m_nIteration < 10) { // when ray intesect shared edge by 2 triangles
          bSpecialCase = true;
          return true;
        }
        else { // overlapped surfaces
          //bSpecialCase = false;
          //m_nIteration = 0;
          //m_iOverlap++;
          m_vIntersection.erase(m_vIntersection.begin() + m_iOverlap + 1);
          nInter = m_vIntersection.size();
          m_vhInterSurf.resize(nInter);
          //m_nIteration = 0;
        }
      }
    }
    else m_iOverlap++;
  }

  return false;
}

void EBMesher::remove_intersection_duplicates()
{
  int ind = 0;
  int nInter = m_vIntersection.size() - 1;

  while (ind < nInter) {
    if (m_vIntersection[ind + 1].distance - m_vIntersection[ind].distance < 1e-7) {
      moab::EntityHandle h1 = m_vhInterSurf[m_vIntersection[ind].index];
      moab::EntityHandle h2 = m_vhInterSurf[m_vIntersection[ind + 1].index];
      
      if (h1 != h2) { // overlapped/shared surfaces
        std::vector<iBase_EntitySetHandle> parents_out1, parents_out2;
        int err = mk_core()->imesh_instance()->getPrnts(reinterpret_cast<iBase_EntitySetHandle> (h1), 1, parents_out1);
        IBERRCHK(err, "Failed to get number of surface parents.\n");
        err = mk_core()->imesh_instance()->getPrnts(reinterpret_cast<iBase_EntitySetHandle> (h2), 1, parents_out2);
        IBERRCHK(err, "Failed to get number of surface parents.\n");
        if (parents_out1.size() == 1 && parents_out2.size() == 1) {
          if (parents_out1[0] != parents_out2[0]) { // parent volues are also different
            m_mhOverlappedSurf[h1] = 0;
            m_mhOverlappedSurf[h2] = 0;
          }
        }
      }
      m_vIntersection.erase(m_vIntersection.begin() + ind + 1);
      nInter--;
    }
    else ind++;
  }
}

bool EBMesher::get_volume_fraction(int vol_frac_div)
{
  std::cout << "starting get_volume_fraction." << std::endl;
  int i, j, k, l, n, iHex, dir, nIntersect, rayIndex[3], index[3],
    otherDir1, otherDir2, err, nParent;
  (void) otherDir1;
  (void) otherDir2;
  double tolerance = 1e-12, dDistance, elem_origin[3],
    elem_interval_size[3], startPnt[3], endPnt[3], rayMaxEnd[3];
  moab::ErrorCode rval;
  bool bOverlapSecond, bClosed;
  std::vector<iBase_EntityHandle> edge_handles;

  double dTotEdgeElem = 0.;
  for (i = 0; i < 3; i++) {
    rayMaxEnd[i] = m_origin_coords[i] + m_nDiv[i]*m_dIntervalSize[i];
    dTotEdgeElem += m_dIntervalSize[i]*(vol_frac_div + 1)*(vol_frac_div + 1);
  }

  // get all hexes
  std::vector<iBase_EntityHandle> hex_handles;
  err = mk_core()->imesh_instance()->getEntities(m_hRootSet, iBase_REGION,
                                                 iMesh_HEXAHEDRON, hex_handles);
  IBERRCHK(err, "Failed to get hexes.\n");

  int hex_size = hex_handles.size();
  /*
  std::vector<int> hex_status(hex_size);
  //int* hex_status = NULL;
  err = mk_core()->imesh_instance()->getIntArrData(&hex_handles[0], hex_size,
                                                   m_elemStatusTag, &hex_status[0]);
  */
  int error;
  std::vector<int> frac_length;
  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, &hex_handles[0],
                      hex_size, m_elemStatusTag,
                      &hex_status, &hex_status_alloc,
                      &hex_status_size, &error);
  ERRORRF("Failed to get hex status.\n");

  for (n = 0; n < hex_size; n++) {
    if (hex_status[n] == 2) { // just boundary
      std::map<moab::EntityHandle, VolFrac> vol_fraction;
      std::map<moab::EntityHandle, VolFrac>::iterator vf_iter;
      std::map<moab::EntityHandle, VolFrac>::iterator vf_end_iter;
      iHex = moab_instance()->id_from_handle(reinterpret_cast<moab::EntityHandle>
                                             (hex_handles[n])) - m_iStartHex;
      index[0] = (iHex%(m_nDiv[0]*m_nDiv[1]))%m_nDiv[0];
      index[1] = (iHex%(m_nDiv[0]*m_nDiv[1]))/m_nDiv[0];
      index[2] = iHex/m_nDiv[0]/m_nDiv[1];
      
      for (i = 0; i < 3; i++) {
        elem_origin[i] = m_origin_coords[i] + index[i]*m_dIntervalSize[i];
        elem_interval_size[i] = m_dIntervalSize[i]/vol_frac_div;
      }
      
      for (dir = 0; dir < 3; dir++) { // x, y, z directions
        otherDir1 = (dir + 1)%3;
        otherDir2 = (dir + 2)%3;
        
        for (j = 0; j < vol_frac_div + 1; j++) {
          for (i = 0; i < vol_frac_div + 1; i++) {
            // get ray start and end points
            if (dir == 0) { // x coordinate ray
              rayIndex[0] = 0;
              rayIndex[1] = i;
              rayIndex[2] = j;
            }
            else if (dir == 1) { // y coordinate ray
              rayIndex[0] = j;
              rayIndex[1] = 0;
              rayIndex[2] = i;
            }
            else if (dir == 2) { // z coordinate ray
              rayIndex[0] = i;
              rayIndex[1] = j;
              rayIndex[2] = 0;
            }
            else IBERRCHK(iBase_FAILURE, "Ray direction should be from 0 to 2.");
            
            for (k = 0; k < 3; k++) {
              if (k == dir) {
                startPnt[k] = elem_origin[k];
                endPnt[k] = startPnt[k] + m_dIntervalSize[k];
              }
              else {
                startPnt[k] = elem_origin[k] + rayIndex[k]*elem_interval_size[k];
                endPnt[k] = startPnt[k];
              }
            }
            
            // ray-tracing
            if (!fire_ray(nIntersect, startPnt, endPnt,
                          tolerance, dir, m_dIntervalSize[dir])) return iBase_FAILURE;
            
            if (nIntersect > 0) { // ray is intersected
              bOverlapSecond = false;
              bClosed = true;
              for (k = 0; k < nIntersect; k++) {
                std::vector<moab::EntityHandle> parents;
                //MBRange parents;
                rval = moab_instance()->get_parent_meshsets(m_vhInterSurf[m_vIntersection[k].index],
                                                            parents);
                if (rval != MB_SUCCESS) {
                  std::cerr << "Couldn't get parents." << std::endl;
                  return false;
                }

                nParent = parents.size();
                dDistance = m_vIntersection[k].distance;
                
                if (is_shared_overlapped_surf(k)) {
                  if (nParent == 2) { // shared surface
                    for (l = 0; l < 2; l++) {
                      if (l == 1) {
                        dDistance *= -1.;
                        bClosed = false;
                      }
                      else bClosed = true;

                      vf_iter = vol_fraction.find(parents[l]);
                      if (vf_iter == vol_fraction.end()) {
                        VolFrac temp_vf(dDistance, bClosed);
                        vol_fraction[parents[l]] = temp_vf;
                        //std::cout << "iHex=" << iHex << ", vh="
                        //        << parents[l] << ", dir=" << dir << ", dDistance1="
                        //        << dDistance << ", vol="
                        //        << temp_vf.vol << std::endl;
                      }
                      else {
                        vf_iter->second.vol += dDistance;
                        vf_iter->second.closed = bClosed;
                        //std::cout << "iHex=" << iHex << ", vh="
                        //        << vf_iter->first << ", dir=" << dir << ", dDistance1="
                        //                         << dDistance << ", vol="
                        //        << vf_iter->second.vol << std::endl;
                      }
                    }
                  }
                  else if (nParent == 1) { // overlapped surface
                    if (bOverlapSecond) { // second intersection
                      /*
                      for (int m = 0; m < 3; m++) {
                        ePnt[m] = startPnt[m];
                      }
                      ePnt[dir] += dDistance;
                      */
                      dDistance *= -1.;
                      bClosed = false;
                      bOverlapSecond = false;
                    }
                    else {
                      /*
                      // make edge
                      for (int m = 0; m < 3; m++) {
                        if (bClosed) ePnt[m] = startPnt[m];
                        ePnt[m + 3] = ePnt[m];
                      }
                      ePnt[dir + 3] += dDistance;
                      if (!make_edge(ePnt, edge_handles)) return iBase_FAILURE;
                      */
                      bOverlapSecond = true;// first intersection
                      bClosed = true;
                    }
                    
                    vf_iter = vol_fraction.find(parents[0]);
                    if (vf_iter == vol_fraction.end()) {
                      VolFrac temp_vf(dDistance, bClosed);
                      vol_fraction[parents[0]] = temp_vf;
                      //std::cout << "iHex=" << iHex << ", vh="
                      //        << parents[0] << ", dDistance2="
                      //        << dDistance << ", vol="
                      //        << temp_vf.vol << std::endl;
                    }
                    else {
                      vf_iter->second.vol += dDistance;
                      vf_iter->second.closed = bClosed;
                      //std::cout << "iHex=" << iHex << ", vh="
                      //        << vf_iter->first << ", dir=" << dir << ", dDistance2="
                      //        << dDistance << ", vol="
                      //        << vf_iter->second.vol << std::endl;
                    }
                  }
                  else return iBase_FAILURE;
                }
                else { // normal case
                  if (nParent != 1) return iBase_FAILURE;

                  if (!is_same_direct_to_ray(k, dir)) {
                    dDistance *= -1.; // outside
                    bClosed = false;
                  }
                  else bClosed = true;
                  
                  vf_iter = vol_fraction.find(parents[0]);
                  if (vf_iter == vol_fraction.end()) {
                    VolFrac temp_vf(dDistance, bClosed);
                    vol_fraction[parents[0]] = temp_vf;
                    //std::cout << "iHex=" << iHex << ", vh="
                    //        << parents[0] << ", dir=" << dir << ", dDistance3="
                    //        << dDistance << ", vol="
                    //        << temp_vf.vol << std::endl;
                  }
                  else {
                    vf_iter->second.vol += dDistance;
                    vf_iter->second.closed = bClosed;
                    //std::cout << "iHex=" << iHex << ", vh="
                    //        << vf_iter->first << ", dir=" << dir << ", dDistance3="
                    //        << dDistance << ", vol="
                    //        << vf_iter->second.vol << std::endl;
                  }
                }
              }

              // if fraction is not closed, add interval size to close
              vf_iter = vol_fraction.begin();
              vf_end_iter = vol_fraction.end();
              for (; vf_iter != vf_end_iter; vf_iter++) {
                if (!vf_iter->second.closed) {
                  vf_iter->second.vol += m_dIntervalSize[dir];
                  vf_iter->second.closed = true;
                  /*
                  for (int m = 0; m < 3; m++) {
                    ePnt[m + 3] = startPnt[m];
                  }
                  ePnt[dir + 3] += m_dIntervalSize[dir];
                  if (!make_edge(ePnt, edge_handles)) return iBase_FAILURE;
                  */
                  //std::cout << "iHex=" << iHex << ", vh="
                  //        << vf_iter->first << ", dir=" << dir << ", dDistance4="
                  //        << m_dIntervalSize[dir] << ", vol="
                  //        << vf_iter->second.vol << std::endl;
                }
              }
            }
            else { // decide if it is fully inside
              if (!fire_ray(nIntersect, startPnt, rayMaxEnd,
                            tolerance, dir, m_nDiv[dir]*m_dIntervalSize[dir])) return iBase_FAILURE;
              
              if (nIntersect > 0) { // fully inside
                if (is_shared_overlapped_surf(0) || // shared or overlapped
                    is_same_direct_to_ray(0, dir)) { // other inside case
                  MBRange parents;
                  rval = moab_instance()->get_parent_meshsets(m_vhInterSurf[m_vIntersection[0].index],
                                                              parents);
                  if (rval != MB_SUCCESS) {
                    std::cerr << "Couldn't get parents." << std::endl;
                    return false;
                  }

                  moab::EntityHandle parent_entity = parents.pop_front();
                  vf_iter = vol_fraction.find(parent_entity);
                  if (vf_iter == vf_end_iter) {
                    VolFrac temp_vf(m_dIntervalSize[dir], true);
                    vol_fraction[parent_entity] = temp_vf;
                    //std::cout << "iHex=" << iHex << ", vh="
                    //        << parents[0] << ", dir=" << dir << ", dDistance5="
                    //                             << dDistance << ", vol="
                    //        << temp_vf.vol << std::endl;
                  }
                  else {
                    vf_iter->second.vol += m_dIntervalSize[dir];
                    vf_iter->second.closed = bClosed;
                    //std::cout << "iHex=" << iHex << ", vh="
                    //        << vf_iter->first << ", dir=" << dir << ", dDistance5="
                    //                         << dDistance << ", vol="
                    //        << vf_iter->second.vol << std::endl;
                  }
                }
              }
            }
          }
        }
      }

      int vol_frac_length = vol_fraction.size();
      std::vector<int> vol_frac_id(vol_frac_length);
      std::vector<double> vol_frac(vol_frac_length);
      int vol_frac_id_size = vol_frac_length*sizeof(int);
      int vol_frac_size = vol_frac_length*sizeof(double);
      vf_iter = vol_fraction.begin();
      vf_end_iter = vol_fraction.end();
      for (int m = 0; vf_iter != vf_end_iter; vf_iter++, m++) {
        vol_frac_id[m] = vf_iter->first;
        vol_frac[m] = vf_iter->second.vol/dTotEdgeElem;
        //std::cout << "iHex=" << iHex << ", i=" << index[0]
        //        << ", j=" << index[1] << ", k=" << index[2]
        //        << ", vol_handle=" << vf_iter->first
        //        << ", vol=" << vf_iter->second.vol
        //        << ", vol_fraction=" << vf_iter->second.vol/dTotEdgeElem
        //        << std::endl;
      }
      const void* vol_frac_ids = &vol_frac_id[0];
      const void* vol_fracs = &vol_frac[0];

      // set volume fraction information as tag
      rval = moab_instance()->tag_set_data(reinterpret_cast<MBTag> (m_volFracLengthTag),
                                           reinterpret_cast<moab::EntityHandle*> (&hex_handles[n]),
                                           1, &vol_frac_length);
      MBERRCHK(rval, mk_core()->moab_instance());

      rval = moab_instance()->tag_set_by_ptr(reinterpret_cast<MBTag> (m_volFracHandleTag),
                                           reinterpret_cast<moab::EntityHandle*> (&hex_handles[n]),
                                           1, &vol_frac_ids, &vol_frac_id_size);
      MBERRCHK(rval, mk_core()->moab_instance());

      rval = moab_instance()->tag_set_by_ptr(reinterpret_cast<MBTag> (m_volFracTag),
                                           reinterpret_cast<moab::EntityHandle*> (&hex_handles[n]),
                                           1, &vol_fracs, &vol_frac_size);
      MBERRCHK(rval, mk_core()->moab_instance());

      if (debug_ebmesh) {
        for (int v = 0; v < vol_frac_length; v++) {
          std::cout << "#_boundary_hex_handles=" << hex_handles[n]
                    << ",vol_frac_handle[" << v << "]=" << vol_frac_id[v]
                    << ",vol_fracs[" << v << "]=" << vol_frac[v]
                    << std::endl;
        }
      }
    }
  }  

  std::cout << "end get_volume_fraction." << std::endl;
  delete [] hex_status;

  return true;
}

#ifdef HAVE_MOAB
bool EBMesher::is_same_direct_to_ray(int i, int dir)
{
  MBCartVect coords[3], normal(0.0), ray_dir(rayDir[dir]);
  const moab::EntityHandle *conn;
  int len;

  // get triangle normal vector
  moab::ErrorCode rval = moab_instance()->get_connectivity(m_vhInterFacet[m_vIntersection[i].index],
                                                           conn, len);
  assert(MB_SUCCESS == rval && 3 == len);
  
  rval = moab_instance()->get_coords(conn, 3, coords[0].array());
  assert(MB_SUCCESS == rval);

  coords[1] -= coords[0];
  coords[2] -= coords[0];
  normal += coords[1] * coords[2];
  normal.normalize();

  return angle(normal, ray_dir) < .5*PI;
}

void EBMesher::set_obb_tree_box_dimension()
{
  std::cout << "starting to construct obb tree." << std::endl;
  moab::ErrorCode rval = m_GeomTopoTool->construct_obb_trees(m_bUseWholeGeom);
  MBERRCHK(rval, mk_core()->moab_instance());

  m_hObbTree = m_GeomTopoTool->obb_tree();
  m_hTreeRoot = m_GeomTopoTool->get_one_vol_root();

  moab::CartVect box_center, box_axis1, box_axis2, box_axis3;
  for (int i = 0; i < 3; i++) {
    m_minCoord[i] = HUGE_VAL;
    m_maxCoord[i] = 0.;
  }
  rval = m_hObbTree->box(m_hTreeRoot, box_center.array(),
                         box_axis1.array(), box_axis2.array(),
                         box_axis3.array());
  MBERRCHK(rval, mk_core()->moab_instance());
  
  // cartesian box corners
  moab::CartVect corners[8] = {box_center + box_axis1 + box_axis2 + box_axis3,
                               box_center + box_axis1 + box_axis2 - box_axis3,
                               box_center + box_axis1 - box_axis2 + box_axis3,
                               box_center + box_axis1 - box_axis2 - box_axis3,
                               box_center - box_axis1 + box_axis2 + box_axis3,
                               box_center - box_axis1 + box_axis2 - box_axis3,
                               box_center - box_axis1 - box_axis2 + box_axis3,
                               box_center - box_axis1 - box_axis2 - box_axis3};
  
  // get the max, min cartesian box corners
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 3; j++) {
      if (corners[i][j] < m_minCoord[j]) m_minCoord[j] = corners[i][j];
      if (corners[i][j] > m_maxCoord[j]) m_maxCoord[j] = corners[i][j];
    }
  }
  std::cout << "finished to construct obb tree." << std::endl;
}
#endif
} // namespace MeshKit
