#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <sstream>
#include <algorithm>
#include <assert.h>

#include "MKDefines.h"
#include "EBMesh.hpp"

#ifdef MOAB
#include "moab/Core.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#include "moab/HomXform.hpp"
#endif

#ifdef CGM
#include "CubitDefines.h"
#include "GeometryQueryTool.hpp"
#include "RefFace.hpp"
#endif

#define PI 3.14159265

double rayDir[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

const bool debug = false;
inline bool less_intersect(const IntersectDist d1, const IntersectDist d2) {
  return d1.distance < d2.distance;
}

inline bool equal_intersect(const IntersectDist d1, const IntersectDist d2) {
  return abs(d1.distance - d2.distance) < 10e-7;
}

EBMesh::EBMesh(iMesh_Instance mesh, iBase_EntitySetHandle root_set,
			 double size, bool use_geom)
  : m_mesh(mesh), m_hRootSet(root_set), m_dInputSize(size)
{
  m_bUseGeom = use_geom;
  m_nStatus = OUTSIDE;
  m_iStartHex = 0;
  //m_nFraction = 0;

#ifdef MOAB
  // create tags with MOAB for dense tags
  int outside = 1;
  const void *out = &outside;
  m_elemStatusTag = get_tag("ELEM_STATUS_TAG", sizeof(int),
			    MB_TAG_DENSE, MB_TYPE_INTEGER, out);
  
  int length = 1;
  const void *leng = &length;
  m_edgeCutFracLengthTag = get_tag("EDGE_CUT_FRAC_LENGTH_TAG", // # of fractions
				   3*sizeof(int), MB_TAG_SPARSE,
				   MB_TYPE_INTEGER, leng);

  m_edgeCutFracTag = get_various_length_tag("EDGE_CUT_FRAC_TAG",
					    MB_TAG_SPARSE, MB_TYPE_DOUBLE);

  int m_id = 1;
  const void *id = &m_id;
  m_matFracIDTag = get_tag("MAT_FRAC_ID_TAG", sizeof(int),
			   MB_TAG_SPARSE, MB_TYPE_INTEGER, id);

  m_volFracTag = get_various_length_tag("VOL_FRAC_TAG",
					MB_TAG_SPARSE, MB_TYPE_DOUBLE);

  m_GeomTopoTool = new GeomTopoTool(reinterpret_cast<MBInterface*>(mesh));
#endif
}

EBMesh::~EBMesh()
{
  delete m_GeomTopoTool;
}

int EBMesh::do_mesh()
{
 #ifdef MOAB
  clock_t time1 = clock();
  unsigned long mem1, mem2, mem3, mem4, mem5;
  moab_instance()->estimated_memory_use(0, 0, 0, &mem1);
  MBErrorCode rval;

  if (debug) {
    rval = moab_instance()->write_mesh("input.vtk");
    MBERRORR("Couldn't write input mesh.", iBase_ERROR_MAP[rval]);
  }

  // 1. construct obb tree for all surfaces and volumes
  int err = construct_obb_tree();
  ERRORR("Couldn't construct obb tree.", err);
  clock_t time2 = clock();
  moab_instance()->estimated_memory_use(0, 0, 0, &mem2);

  // 2. set division
  err = set_division();
  ERRORR("Couldn't set division.", err);

  // 3. make hex vertices
  err = make_scd_hexes();
  ERRORR("Couldn't make scd hexes.", err);
  clock_t time3 = clock();
  moab_instance()->estimated_memory_use(0, 0, 0, &mem3);

  // 4. find intersected geometry surfaces by rays
  err = find_intersections();
  ERRORR("Couldn't find intersected surfaces.", err);
  clock_t time4 = clock();
  moab_instance()->estimated_memory_use(0, 0, 0, &mem4);
  
  // 5. set hex status and boundary hex cut fraction info
  err = set_tag_info();
  ERRORR("Couldn't set tag infor.", err);
  clock_t time5 = clock();
  moab_instance()->estimated_memory_use(0, 0, 0, &mem5);
#endif

  if (debug) {
    std::cout << "OBB_tree_construct_time: "
	      << (double) (time2 - time1)/CLOCKS_PER_SEC
	      << ", scd_hex_construct_time: "
	      << (double) (time3 - time2)/CLOCKS_PER_SEC
	      << ", intersection_time: "
	      << (double) (time4 - time3)/CLOCKS_PER_SEC
	      << ", set_info_time: "
	      << (double) (time5 - time4)/CLOCKS_PER_SEC
	      << std::endl;

    std::cout << "start_memory: " << mem1
	      << ", OBB_tree_construct_moemory: " << mem2
	      << ",hex_construct_memory : " << mem3
	      << ", intersection_memory: " << mem4
	      << ", set_info_memory: " << mem5
	      << std::endl;
  }

  return iBase_SUCCESS;
}

int EBMesh::construct_obb_tree()
{
  if (m_bUseGeom) {
    // construct obb tree for geometry surfaces and volumes by GeomTopoTool
    MBErrorCode rval = m_GeomTopoTool->construct_obb_trees(true);
    MBERRORR("Couldn't construct obb tree in GeomTopoTool.", iBase_ERROR_MAP[rval]);
    
    m_hObbTree = m_GeomTopoTool->obb_tree();
    m_hTreeRoot = m_GeomTopoTool->get_one_vol_root();
  }
  else { // facet data input case
    // get all triangles
    MBRange tris;
    MBErrorCode rval = moab_instance()->
      get_entities_by_dimension(reinterpret_cast<MBEntityHandle>(m_hRootSet), 2, tris);
    MBERRORR("Failed to get triangles.", iBase_ERROR_MAP[rval]);
    
    // make tree
    m_hObbTree = new MBOrientedBoxTreeTool(reinterpret_cast<MBInterface*> (m_mesh));
    rval = m_hObbTree->build(tris, m_hTreeRoot);
    MBERRORR("Failed to build tree.", iBase_ERROR_MAP[rval]);
  }

  return iBase_SUCCESS;
}

bool EBMesh::export_mesh(const char* file_name, bool separate)
{ 
  // get all hexes
  clock_t time1 = clock();
  int i, err;
  iBase_EntityHandle* hex_handles = NULL;
  int hex_allocated = 0;
  int hex_size = 0;
  iMesh_getEntities(m_mesh, m_hRootSet, iBase_REGION,
		    iMesh_HEXAHEDRON, &hex_handles,
		    &hex_allocated, &hex_size, &err);
  ERRORRF("Failed to get hexes.\n"); 

  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, hex_handles,
		      hex_size, m_elemStatusTag,
		      &hex_status, &hex_status_alloc,
		      &hex_status_size, &err);
  ERRORRF("Failed to get hex status.\n");
  clock_t time2 = clock();
  clock_t time3;

  if (separate) {
    int n_inside_hex = 0;
    int n_outside_hex = 0;
    int n_boundary_hex = 0;
    int hex_stat;
    std::vector<iBase_EntityHandle> insideHex, outsideHex, bndrHex;
    for (i = 0; i < hex_status_size; i++) {
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
    else ERRORRF("Hex element status should be inside/outside/boundary.\n");
    }
    
    std::cout << "# of exported inside hex:" << n_inside_hex
	      << ", # of exported outside hex:" << n_outside_hex
	      << ", # of exported boundary hex:" << n_boundary_hex
	      << ", geom vol:"
	      << n_inside_hex*m_dIntervalSize[0]*m_dIntervalSize[1]*m_dIntervalSize[2]
	      << ", vox vol:" << hex_size*m_dIntervalSize[0]*m_dIntervalSize[1]*m_dIntervalSize[2]
	      << std::endl;
    time3 = clock();

    // save inside/outside/boundary elements separately
    if (n_inside_hex > 0) {
      err = write_mesh(file_name, 0, &insideHex[0], n_inside_hex);
      ERRORRF("Couldn't write inside mesh.");
    }
    
    if (n_outside_hex > 0) {
      err = write_mesh(file_name, 1, &outsideHex[0], n_outside_hex);
      ERRORRF("Couldn't write outside mesh.");
    }
    
    if (n_boundary_hex > 0) {
      err = write_mesh(file_name, 2, &bndrHex[0], n_boundary_hex);
      ERRORRF("Couldn't write boundary mesh.");
    }
    
    if (debug) {
      std::cout << "hex_handle_get_time: "
		<< (double) (time2 - time1)/CLOCKS_PER_SEC
		<< ", separate_write_time: "
		<< (double) (time3 - time2)/CLOCKS_PER_SEC
		<< std::endl;
    }
  }
  else {
    err = write_mesh(file_name, 3, hex_handles, hex_size);
    ERRORRF("Couldn't write all hex mesh.");
  }

  free(hex_handles);
  free(hex_status);

  return true;
}

int EBMesh::write_mesh(const char* file_name, int type,
			    iBase_EntityHandle* handles, int& n_elem)
{
  clock_t time1 = clock();
  int is_list = 1, err;
  MBErrorCode rval;
  iBase_EntitySetHandle set;

  iMesh_createEntSet(m_mesh, is_list, &set, &err);
  ERRORR("Couldn't create set.\n", err);

  iMesh_addEntArrToSet(m_mesh, handles, n_elem, set, &err);
  ERRORR("Couldn't add hexes to set.\n", err);
  clock_t time2 = clock();

  std::string out_name;
  std::stringstream ss;
  if (type == 0) ss << "inside_";
  else if (type == 1) ss << "outside_";
  else if (type == 2) ss << "boundary_";
  ss << file_name;
  ss >> out_name;
  
  rval = moab_instance()->write_mesh(out_name.c_str(),
				     (const MBEntityHandle*) &set, 1);
  MBERRORR("Failed to write hex mesh.", iBase_ERROR_MAP[rval]);
  
  std::cout << "Elements are exported." << std::endl;
  clock_t time3 = clock();

  if (debug) {
    std::cout << "set_creation_time: "
	      << (double) (time2 - time1)/CLOCKS_PER_SEC
	      << ", write_time: "
	      << (double) (time3 - time2)/CLOCKS_PER_SEC
	      << std::endl;
  }

  return iBase_SUCCESS;
}

EdgeStatus EBMesh::get_edge_status(const double dP)
{
  if (m_nStatus == INSIDE) { // previous inside
    if (dP < m_dSecondP) {
      m_nMove = 0;
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
      return OUTSIDE;
    }
    else {
      if (dP < m_dSecondP) {
	m_nMove = 0;
	if (m_prevPnt < m_dFirstP) return BOUNDARY;
	else return INSIDE;
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
  }
  else {
    std::cerr << "Couldn't get edge status." << std::endl;
    return INSIDE;
  }
}

bool EBMesh::set_neighbor_hex_status(int dir, int i, int j, int k)
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

bool EBMesh::set_hex_status(int index, EdgeStatus value, int dir)
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

bool EBMesh::set_edge_status(int dir)
{
  // set boundary cut information to edge
  std::vector<double> vdCutFraction;
  if (m_nMove == 0) vdCutFraction.push_back(1. - (m_curPnt - m_dFirstP)/m_dIntervalSize[dir]);
  else if (m_nMove == 1) vdCutFraction.push_back(1. - (m_curPnt - m_dSecondP)/m_dIntervalSize[dir]);
  else if (m_nMove == 2) {
    vdCutFraction.push_back(1. - (m_curPnt - m_dSecondP)/m_dIntervalSize[dir]);
    if (m_dFirstP < m_curPnt) {
      vdCutFraction.push_back(1. - (m_curPnt - m_dFirstP)/m_dIntervalSize[dir]);
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

#ifdef MOAB
int EBMesh::make_scd_hexes()
{
  // create vertex and hex sequences
  int i;
  double max_coords[3];
  for (i = 0; i < 3; i++) {
    max_coords[i] = m_origin_coords[i] + m_nDiv[i]*m_dIntervalSize[i];
  }

  moab::HomCoord coord_min(0, 0, 0);
  moab::HomCoord coord_max(m_nDiv[0], m_nDiv[1], m_nDiv[2]);
  moab::EntitySequence* vertex_seq = NULL;
  moab::EntitySequence* cell_seq = NULL;
  MBEntityHandle vs, cs;
  moab::Core *mbcore = dynamic_cast<moab::Core*>(moab_instance());

  MBErrorCode rval = mbcore->create_scd_sequence(coord_min, coord_max, MBVERTEX, 1, vs, vertex_seq);
  MBERRORR("Failed to create structured vertices.", iBase_ERROR_MAP[rval]);
  
  mbcore->create_scd_sequence(coord_min, coord_max, MBHEX, 1, cs, cell_seq);
  MBERRORR("Failed to create structured hexes.", iBase_ERROR_MAP[rval]);

  moab::HomCoord p2(coord_max.i(), coord_min.j(), coord_min.k());
  moab::HomCoord p3(coord_min.i(), coord_max.j(), coord_min.k()); 
  
  rval = mbcore->add_vsequence(vertex_seq, cell_seq, coord_min, coord_min,
					p2, p2, p3, p3);
  MBERRORR("Failed to add vertex sequence.", iBase_ERROR_MAP[rval]);

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
    MBERRORR("Failed to set coords.", iBase_ERROR_MAP[rval]);
  }

  m_iStartHex = moab_instance()->id_from_handle(cs);

  return iBase_SUCCESS;
}

iBase_TagHandle EBMesh::get_tag(const char* name, int size,
				     MBTagType store, MBDataType type,
				     const void* def_value,
				     bool create_if_missing) 
{
  MBTag retval = 0;
  MBErrorCode result = moab_instance()->tag_create(name, size, store, type,
						   retval, def_value,
						   create_if_missing);
  if (create_if_missing && MB_SUCCESS != result) {
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;
  }
  
  return (iBase_TagHandle) retval;
}

iBase_TagHandle EBMesh::get_various_length_tag(const char* name,
						    MBTagType store, MBDataType type)
{
  MBTag retval = 0;
  MBErrorCode result = moab_instance()->
    tag_create_variable_length( name, store, type, retval);
  if (MB_SUCCESS != result) {
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;
  }
  
  return (iBase_TagHandle) retval;
}

int EBMesh::set_tag_info()
{
  // get all hexes
  int i, j, k, err;
  iBase_EntityHandle* hex_handles = NULL;
  int hex_allocated = 0;
  int hex_size = 0;
  iMesh_getEntities(m_mesh, m_hRootSet, iBase_REGION,
		    iMesh_HEXAHEDRON, &hex_handles,
		    &hex_allocated, &hex_size, &err);
  ERRORR("Failed to get hexes.\n", err);

  iMesh_setIntArrData(m_mesh, hex_handles, m_nHex, m_elemStatusTag,
		      &m_vnHexStatus[0], m_nHex, &err);
  ERRORR("Failed to set hex element status data.", err);

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

  iMesh_setIntArrData(m_mesh, &hvBndrHex[0],
		      nBndrHex, m_edgeCutFracLengthTag,
		      frac_leng, 3*nBndrHex, &err);
  ERRORR("Failed to set cut fraction sizes to hex.", err);

  MBErrorCode rval = moab_instance()->tag_set_data(reinterpret_cast<MBTag> (m_edgeCutFracTag),
						   reinterpret_cast<MBEntityHandle*> (&hvBndrHex[0]),
						   nBndrHex, &frac_data_pointer[0], frac_size);
  MBERRORR("Failed to set cut fraction infor to hex.", iBase_ERROR_MAP[rval]);
  
  delete [] frac_size;
  delete [] frac_leng;
  delete [] fractions;
  
  return iBase_SUCCESS;
}

static inline double axis_len( const double* axis )
{
  return std::sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
}

int EBMesh::set_division()
{
  int i, j;
  MBCartVect box_center, box_axis1, box_axis2, box_axis3,
    min_cart_box(HUGE_VAL, HUGE_VAL, HUGE_VAL),
    max_cart_box(0., 0., 0.);
  
  MBErrorCode rval = m_hObbTree->box(m_hTreeRoot, box_center.array(),
				     box_axis1.array(), box_axis2.array(),
				     box_axis3.array());
  MBERRORR("Error getting box information for tree root set.", iBase_ERROR_MAP[rval]);

  // cartesian box corners
  MBCartVect corners[8] = {box_center + box_axis1 + box_axis2 + box_axis3,
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
  MBCartVect length = max_cart_box - min_cart_box;
  
  // default value is adjusted to large geometry file
  // interval_size_estimate : 2*L*sqrt(2*PI*sqrt(2)/# of tris)
  if (m_dInputSize < 0.) {
    int n_tri;
    rval = moab_instance()->
      get_number_entities_by_dimension(reinterpret_cast<MBEntityHandle>(m_hRootSet),
				       2, n_tri);
    MBERRORR("Failed to get number of triangles.", iBase_ERROR_MAP[rval]);
    
    double box_length_ave = (length[0] + length[1] + length[2])/3.;
    m_dInputSize = 2.*box_length_ave*sqrt(8.886/n_tri);
  }

  for (i = 0; i < 3; i++) {
    m_nDiv[i] = length[i]/m_dInputSize;
    if (m_nDiv[i] < 5) m_nDiv[i] = 5; // make 5 elements larger than bounding box
    m_dIntervalSize[i] = m_dInputSize;
    if (m_nDiv[i]*.07 > 5) m_nDiv[i] += m_nDiv[i]*.07;
    else m_nDiv[i] += 5;
    m_nNode[i] = m_nDiv[i] + 1;
    m_origin_coords[i] = box_center[i] - .5*m_nDiv[i]*m_dIntervalSize[i];
  }

  m_nHex = m_nDiv[0]*m_nDiv[1]*m_nDiv[2];

  std::cout << "# of hex: " << m_nHex << ", interval size: "
	    << m_dInputSize << std::endl;

  return iBase_SUCCESS;
}

int EBMesh::find_intersections()
{
  int i, err;
  m_vnHexStatus.resize(m_nHex, 1); // initialize all hex as outside

  // fire rays to 3 directions
  for (i = 0; i < 3; i++) {
    m_vnEdgeStatus[i].resize(m_nDiv[i]*m_nDiv[(i + 1)%3]*m_nDiv[(i + 2)%3], OUTSIDE);
    err = fire_rays(i);
    if (err != iBase_SUCCESS) return err;
  }

  return iBase_SUCCESS;
}

int EBMesh::fire_rays(int dir)
{
  // ray fire
  int i, j, k, l, index[3];
  double tolerance = 1e-12;
  double rayLength = m_nDiv[dir]*m_dIntervalSize[dir];
  int err, iNodeStart, iNodeEnd, nIntersect, nNodeSlice;
  double startPnt[3], endPnt[3];
  int otherDir1 = (dir + 1)%3;
  int otherDir2 = (dir + 2)%3;

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
      else ERRORR("Ray direction should be 0 to 2.", iBase_FAILURE);

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
      if (!fire_ray(nIntersect, startPnt, endPnt,
		    tolerance, dir, rayLength)) return iBase_FAILURE;
      
      if (nIntersect > 0) {
	m_nStatus = OUTSIDE;
	m_dFirstP = startPnt[dir] + m_vIntersection[m_iInter++].distance;
	m_dSecondP = startPnt[dir] + m_vIntersection[m_iInter++].distance;
	m_prevPnt = startPnt[dir];

	for (; k < m_nNode[dir] - 1; k++) {
	  m_curPnt = startPnt[dir] + (k + 1)*m_dIntervalSize[dir];
	  m_nStatus = get_edge_status(m_curPnt);
	  m_vnEdgeStatus[dir][i*m_nDiv[dir] + j*m_nDiv[dir]*m_nDiv[otherDir1] + k] = m_nStatus;
	  
	  // set status of all hexes sharing the edge
	  if (!set_neighbor_hex_status(dir, i, j, k)) return iBase_FAILURE;
	  
	  if (m_nMove > 0) {
	    if (m_iInter < nIntersect) {
	      if (!move_intersections(dir, nIntersect, startPnt)) return iBase_FAILURE;
	    }
	    else {
	      m_nMove = 1;
	      if (m_nStatus == BOUNDARY && !set_edge_status(dir)) return iBase_FAILURE;
	      k++;
	      break; // rest is all outside
	    }
	  }
	  else if (m_nStatus == BOUNDARY && !set_edge_status(dir)) return iBase_FAILURE;
	  
	  // set cut-cell edge status
	  m_prevPnt = m_curPnt;
	}
      }
      
      // the rest are all outside
      for (; k < m_nNode[dir] - 1; k++) {
	m_nStatus = OUTSIDE;
	
	if (!set_neighbor_hex_status(dir, i, j, k)) {
	  return iBase_FAILURE;
	}
      }
    }
  }
  
  return iBase_SUCCESS;
}

bool EBMesh::fire_ray(int& nIntersect, double startPnt[3],
		      double endPnt[3], double tol, int dir,
		      double rayLength)
{
  m_vIntersection.clear();
  m_vhInterSurf.clear();
  m_vhInterFacet.clear();
  m_mhOverlappedSurf.clear();
  std::vector<double> temp_intersects;
  MBErrorCode rVal;
  
  if (m_bUseGeom) { // geometry input
    rVal = m_hObbTree->ray_intersect_sets(temp_intersects, m_vhInterSurf,
					  m_vhInterFacet, m_hTreeRoot, tol,
					  -1, startPnt, rayDir[dir], &rayLength);
  }
  else { // facet input
    rVal = m_hObbTree->ray_intersect_triangles(temp_intersects, m_hTreeRoot, tol,
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
    
    if (nIntersect > 1) { // when ray intersect shared edge of triangles
      bool bMoveOnce;
      m_nIteration = 0;
      m_iOverlap = 0;
      if (is_ray_move_and_set_overlap_surf(bMoveOnce)) {
	if (!move_ray(nIntersect, startPnt, endPnt, tol, dir, bMoveOnce)) {
	  std::cerr << "Number of Intersection between edges and ray should be even." << std::endl;
	  return false;
	}
      }
    }
  }

  return true;
}

bool EBMesh::move_intersections(int n_dir, int n_inter, double start_pnt[3])
{
  if (m_nMove > 0) {
    if (m_iInter < n_inter) {
      if (m_nMove == 1) {
	do {
	  if (m_nStatus == BOUNDARY && !set_edge_status(n_dir)) return false;
	  m_dFirstP = m_dSecondP;
	  m_dSecondP = start_pnt[n_dir] + m_vIntersection[m_iInter++].distance;
	}
	while (m_dSecondP < m_curPnt && m_iInter < n_inter);
      }
      else if (m_nMove == 2) {
	do {
	  m_dFirstP = start_pnt[n_dir] + m_vIntersection[m_iInter++].distance;
	  if (m_nStatus == BOUNDARY && !set_edge_status(n_dir)) return false;
	  if (m_iInter < n_inter) {
	    m_dSecondP = start_pnt[n_dir] + m_vIntersection[m_iInter++].distance;
	  }
	  else m_dSecondP = m_dFirstP;
	}
	while (m_dSecondP < m_curPnt && m_iInter < n_inter);
      }
    }
  }

  return true;
}

bool EBMesh::is_shared_overlapped_surf(int index)
{
  int nParent, err;
  MBEntityHandle hSurf;
  if (m_bUseGeom) {
    hSurf = m_vhInterSurf[m_vIntersection[index].index];
    iMesh_getNumPrnt(m_mesh,
		     reinterpret_cast<iBase_EntitySetHandle> (hSurf),
		     1, &nParent, &err);
    ERRORRF("Failed to get number of surface parents.\n");

    if (nParent > 1) return true;
  }
  else hSurf = m_vIntersection[index].index;

  return m_mhOverlappedSurf.count(hSurf) > 0;
}

bool EBMesh::get_grid_and_edges_techX(double* boxMin, double* boxMax, int* nDiv,
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
  iBase_EntityHandle* hex_handles = NULL;
  int hex_allocated = 0;
  int hex_size = 0;
  iMesh_getEntities(m_mesh, m_hRootSet, iBase_REGION,
		    iMesh_HEXAHEDRON, &hex_handles,
		    &hex_allocated, &hex_size, &err);
  ERRORRF("Failed to get hexes.\n");
  
  // get hex status
  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, hex_handles,
		      hex_size, m_elemStatusTag,
		      &hex_status, &hex_status_alloc,
		      &hex_status_size, &err);
  ERRORRF("Failed to get hex status.\n");

  // get inside and boundary hexes
  int nInside = 0;
  int nOutside = 0;
  int nBoundary = 0;
  rvnInsideCell.clear();
  for (i = 0; i < hex_size; i++) {
    if (hex_status[i] == 0) { // if inside
      iHex = moab_instance()->id_from_handle(reinterpret_cast<MBEntityHandle>
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
    if (cutFrac.fraction[1].size() > 0) dFrac[0] = cutFrac.fraction[1][0];
    else dFrac[0] = -1.;
    if (cutFrac.fraction[2].size() > 0) dFrac[3] = cutFrac.fraction[2][0];
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
    if (cutFrac.fraction[0].size() > 0) dFrac[0] = cutFrac.fraction[0][0];
    else dFrac[0] = -1.;
    if (cutFrac.fraction[1].size() > 0) dFrac[3] = cutFrac.fraction[1][0];
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
    if (cutFrac.fraction[0].size() > 0) dFrac[0] = cutFrac.fraction[0][0];
    else dFrac[0] = -1.;
    if (cutFrac.fraction[1].size() > 0) dFrac[3] = cutFrac.fraction[1][0];
    else dFrac[3] = -1.;
    dFrac[1] = get_edge_fraction(iHex + 1, 1);
    dFrac[2] = get_edge_fraction(iHex + m_nDiv[0], 0);

    if (dFrac[0] > 0. || dFrac[1] > 0. || dFrac[22] > 0. || dFrac[3] > 0.) { // if surface is cut
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

  free(hex_handles);
  free(hex_status);

  if (debug) return export_fraction_edges(rmdCutCellSurfEdge);
 
  return true;
}

bool EBMesh::get_grid_and_edges(double* boxMin, double* boxMax, int* nDiv,
				     std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellEdge,
				     std::vector<int>& rvnInsideCell, bool isCornerExterior)
{
  int i, ii, jj, kk, iHex;
  for (i = 0; i < 3; i++) {
    boxMin[i] = m_origin_coords[i];
    boxMax[i] = m_origin_coords[i] + m_dIntervalSize[i]*m_nDiv[i];
    nDiv[i] = m_nDiv[i];
  }
  
  if (!get_inside_boundary_hex(rvnInsideCell)) return false;

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

  if (debug) return export_fraction_points(rmdCutCellEdge);
 
  return true;
}

bool EBMesh::get_inside_boundary_hex(std::vector<int>& rvnInsideCell)
{
  int i, err, iHex;
  
  // get all hexes
  iBase_EntityHandle* hex_handles = NULL;
  int hex_allocated = 0;
  int hex_size = 0;
  iMesh_getEntities(m_mesh, m_hRootSet, iBase_REGION,
		    iMesh_HEXAHEDRON, &hex_handles,
		    &hex_allocated, &hex_size, &err);
  ERRORRF("Failed to get hexes.\n");
  
  // get hex status
  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, hex_handles,
		      hex_size, m_elemStatusTag,
		      &hex_status, &hex_status_alloc,
		      &hex_status_size, &err);
  ERRORRF("Failed to get hex status.\n");

  // get inside and boundary hexes
  int nInside = 0;
  int nOutside = 0;
  int nBoundary = 0;
  rvnInsideCell.clear();
  for (i = 0; i < hex_size; i++) {
    if (hex_status[i] == 0) { // if inside
      iHex = moab_instance()->id_from_handle(reinterpret_cast<MBEntityHandle>
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

  free(hex_handles);
  free(hex_status);

  return true;
}

bool EBMesh::export_fraction_edges(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellSurfEdge)
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

  int err, is_list = 1;
  iBase_EntitySetHandle set;
  
  iMesh_createEntSet(m_mesh, is_list, &set, &err);
  ERRORR("Couldn't create set.\n", err);
  
  iMesh_addEntArrToSet(m_mesh, &edge_handles[0], edge_handles.size(), set, &err);
  ERRORR("Couldn't add edges to set.\n", err);
  
  MBErrorCode rval = moab_instance()->write_mesh("edges.vtk",
						 (const MBEntityHandle*) &set, 1);
  MBERRORR("Couldn't write edges.\n", iBase_ERROR_MAP[rval]);

  return true;
}

bool EBMesh::export_fraction_points(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellEdge)
{
  // export fractions as edge
  int i, j, dir, nFrc, err;
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
	iMesh_createVtx(m_mesh, fracPnt[0], fracPnt[1], fracPnt[2],
			&v_handle, &err);
	ERRORR("Couldn't create vertex.\n", err);
	vertex_handles.push_back(v_handle);
      }
    }
  }
    
  int is_list = 1;
  MBErrorCode result;
  iBase_EntitySetHandle set;
  
  iMesh_createEntSet(m_mesh, is_list, &set, &err);
  ERRORR("Couldn't create set.\n", err);
  
  iMesh_addEntArrToSet(m_mesh, &vertex_handles[0], vertex_handles.size(), set, &err);
  ERRORR("Couldn't add vertices to set.\n", err);
  
  result = moab_instance()->write_mesh("frac_vertices.vtk",
				       (const MBEntityHandle*) &set, 1);
  if (MB_SUCCESS != result) {
    std::cerr << "Failed to write fraction vertices." << std::endl;
    return false;
  }

  return true;
}

bool EBMesh::make_edge(double ePnt[6], std::vector<iBase_EntityHandle>& edge_handles)
{
  int err, status;
  int vertex_alloc = sizeof(iBase_EntityHandle)*2;
  int vertex_size = 2;
  iBase_EntityHandle vertex_handle[2], edge_handle;
  iBase_EntityHandle* pVertexHandle = &vertex_handle[0];

  iMesh_createVtxArr(m_mesh, 2,
		     iBase_INTERLEAVED, ePnt, 6,
		     &pVertexHandle,
		     &vertex_alloc, &vertex_size, &err);
  ERRORRF("Failed to create vertices.\n");

  iMesh_createEnt(m_mesh, iMesh_LINE_SEGMENT, &vertex_handle[0], 2,
		  &edge_handle, &status, &err);
  ERRORRF("Failed to create edge.\n");

  edge_handles.push_back(edge_handle);

  return true;
}

double EBMesh::get_edge_fraction(int idHex, int dir)
{
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.find(idHex);
  if (iter != end_iter && iter->second.fraction[dir].size() > 0) return iter->second.fraction[dir][0];
  else return -1.;
}
#endif

double EBMesh::get_uncut_edge_fraction(int i, int j, int k, int dir)
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

bool EBMesh::move_ray(int& nIntersect, double* startPnt, double* endPnt,
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
      }
      else return false;             // error
    }

    for (i = 0; i < 3; i++) {
      if (i != dir) {
	startPnt[i] += tol;
	endPnt[i] += tol;
      }
    }

    MBCartVect ray(endPnt[0] - startPnt[0], endPnt[1] - startPnt[1], endPnt[2] - startPnt[2]);
    double rayLength = ray.length();
    MBErrorCode rVal;
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
      rVal = m_hObbTree->ray_intersect_triangles(temp_intersects, m_hTreeRoot, tol,
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

bool EBMesh::is_ray_move_and_set_overlap_surf(bool& bSpecialCase)
{
  int nInter = m_vIntersection.size() - 1;

  while (m_iOverlap < nInter) {
    if (m_vIntersection[m_iOverlap + 1].distance - m_vIntersection[m_iOverlap].distance < 1e-7) {
      if (m_bUseGeom) {
	MBEntityHandle h1 = m_vhInterSurf[m_vIntersection[m_iOverlap].index];
	MBEntityHandle h2 = m_vhInterSurf[m_vIntersection[m_iOverlap + 1].index];
	
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
	  m_nIteration = 0;
	  m_iOverlap++;
	}
      }
    }
    else m_iOverlap++;
  }

  return false;
}

bool EBMesh::get_volume_fraction(int vol_frac_div)
{
  int i, j, k, l, n, iHex, dir, nIntersect, rayIndex[3], index[3],
    otherDir1, otherDir2, err, nParent;
  double tolerance = 1e-12, dDistance, elem_origin[3],
    elem_interval_size[3], startPnt[3], endPnt[3], rayMaxEnd[3];
  MBErrorCode rval;
  bool bOverlapSecond, bClosed;
  std::vector<iBase_EntityHandle> edge_handles;

  double dTotEdgeElem = 0.;
  for (i = 0; i < 3; i++) {
    rayMaxEnd[i] = m_origin_coords[i] + m_nDiv[i]*m_dIntervalSize[i];
    dTotEdgeElem += m_dIntervalSize[i]*(m_nVolFracDiv + 1)*(m_nVolFracDiv + 1);
  }

  // get all hexes
  iBase_EntityHandle* hex_handles = NULL;
  int hex_allocated = 0;
  int hex_size = 0;
  iMesh_getEntities(m_mesh, m_hRootSet, iBase_REGION,
		    iMesh_HEXAHEDRON, &hex_handles,
		    &hex_allocated, &hex_size, &err);
  ERRORRF("Failed to get hexes.\n"); 

  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, hex_handles,
		      hex_size, m_elemStatusTag,
		      &hex_status, &hex_status_alloc,
		      &hex_status_size, &err);
  ERRORRF("Failed to get hex status.\n");

  for (n = 0; n < hex_status_size; n++) {
    if (hex_status[n] == 2) { // just boundary
      std::map<MBEntityHandle, VolFrac> vol_fraction;
      std::map<MBEntityHandle, VolFrac>::iterator vf_iter;
      std::map<MBEntityHandle, VolFrac>::iterator vf_end_iter;
      iHex = moab_instance()->id_from_handle(reinterpret_cast<MBEntityHandle>
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
	    else ERRORR("Ray direction should be from 0 to 2.", iBase_FAILURE);
	    
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
	      //std::cout << "startPnt=" << startPnt[0] << " "
	      //	<< startPnt[1] << " " << startPnt[2] << std::endl;
	      //std::cout << "endPnt=" << endPnt[0] << " "
	      //	<< endPnt[1] << " " << endPnt[2] << std::endl;
	      bOverlapSecond = false;
	      bClosed = true;
	      for (k = 0; k < nIntersect; k++) {
		std::vector<MBEntityHandle> parents;
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
			//  << parents[l] << ", dir=" << dir << ", dDistance1="
			//  << dDistance << ", vol="
			//  << temp_vf.vol << std::endl;
		      }
		      else {
			vf_iter->second.vol += dDistance;
			vf_iter->second.closed = bClosed;
			//std::cout << "iHex=" << iHex << ", vh="
			//  << vf_iter->first << ", dir=" << dir << ", dDistance1="
			//  << dDistance << ", vol="
			//  << vf_iter->second.vol << std::endl;
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
		      //	<< parents[0] << ", dDistance2="
		      //	<< dDistance << ", vol="
		      //	<< temp_vf.vol << std::endl;
		    }
		    else {
		      vf_iter->second.vol += dDistance;
		      vf_iter->second.closed = bClosed;
		      //std::cout << "iHex=" << iHex << ", vh="
		      //	<< vf_iter->first << ", dir=" << dir << ", dDistance2="
		      //	<< dDistance << ", vol="
		      //	<< vf_iter->second.vol << std::endl;
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
		    //      << parents[0] << ", dir=" << dir << ", dDistance3="
		    //      << dDistance << ", vol="
		    //      << temp_vf.vol << std::endl;
		  }
		  else {
		    vf_iter->second.vol += dDistance;
		    vf_iter->second.closed = bClosed;
		    //std::cout << "iHex=" << iHex << ", vh="
		    //      << vf_iter->first << ", dir=" << dir << ", dDistance3="
		    //      << dDistance << ", vol="
		    //      << vf_iter->second.vol << std::endl;
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
		  //    << vf_iter->first << ", dir=" << dir << ", dDistance4="
		  //    << m_dIntervalSize[dir] << ", vol="
		  //    << vf_iter->second.vol << std::endl;
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

		  MBEntityHandle parent_entity = parents.pop_front();
		  vf_iter = vol_fraction.find(parent_entity);
		  if (vf_iter == vf_end_iter) {
		    VolFrac temp_vf(m_dIntervalSize[dir], true);
		    vol_fraction[parent_entity] = temp_vf;
		    //std::cout << "iHex=" << iHex << ", vh="
		    //      << parents[0] << ", dir=" << dir << ", dDistance5="
		    //      << dDistance << ", vol="
		    //      << temp_vf.vol << std::endl;
		  }
		  else {
		    vf_iter->second.vol += m_dIntervalSize[dir];
		    vf_iter->second.closed = bClosed;
		    //std::cout << "iHex=" << iHex << ", vh="
		    //      << vf_iter->first << ", dir=" << dir << ", dDistance5="
		    //      << dDistance << ", vol="
		    //      << vf_iter->second.vol << std::endl;
		  }
		}
	      }
	    }
	  }
	}
      }
      /*
      vf_iter = vol_fraction.begin();
      vf_end_iter = vol_fraction.end();
      for (; vf_iter != vf_end_iter; vf_iter++) {
	std::cout << "iHex=" << iHex << ", i=" << index[0]
		  << ", j=" << index[1] << ", k=" << index[2]
		  << ", vol_handle=" << vf_iter->first
		  << ", vol=" << vf_iter->second.vol
		  << ", vol_fraction=" << vf_iter->second.vol/dTotEdgeElem
		  << std::endl;
		  }*/
    }
  }  

  free(hex_handles);
  delete [] hex_status;

  return true;
}
  
bool EBMesh::is_same_direct_to_ray(int i, int dir)
{
  MBCartVect coords[3], normal(0.0), ray_dir(rayDir[dir]);
  const MBEntityHandle *conn;
  int len;

  // get triangle normal vector
  MBErrorCode rval = moab_instance()->get_connectivity(m_vhInterFacet[m_vIntersection[i].index],
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
