#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <sstream>
#include <algorithm>
#include <assert.h>

#include "MKDefines.h"
#include "CutCellMesh.hpp"

#ifdef MOAB
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#endif

#ifdef CGM
#include "CubitDefines.h"
#include "GeometryQueryTool.hpp"
#include "RefFace.hpp"
#endif

#define PI 3.14159265

const bool debug = false;
inline bool less_intersect(const IntersectDist d1, const IntersectDist d2) {
  return d1.distance < d2.distance;
}

inline bool equal_intersect(const IntersectDist d1, const IntersectDist d2) {
  return abs(d1.distance - d2.distance) < 10e-7;
}

CutCellMesh::CutCellMesh(iMesh_Instance mesh, iBase_EntitySetHandle root_set, double size)
  : m_mesh(mesh), m_hRootSet(root_set), m_dInputSize(size)
{
  m_nStatus = OUTSIDE;
  m_iStartHex = 0;
  m_nFraction = 0;

#ifdef MOAB
  // create tags with MOAB to make them as dense tags
  int outside = 1;
  const void *out = &outside;
  m_elemStatusTag = get_tag("ELEM_STATUS_TAG", sizeof(int),
			    MB_TAG_DENSE, MB_TYPE_INTEGER, out);
  
  int length = 1;
  const void *leng = &length;
  m_edgeCutFracLengthTag = get_tag("EDGE_CUT_FRACTION_LENGTH_TAG",
				   3*sizeof(int), MB_TAG_DENSE,
				   MB_TYPE_INTEGER, leng);

  double fraction = 0.;
  const void *frac = &fraction;
  m_edgeCutFracTag = get_various_length_tag("EDGE_CUT_FRACTION_TAG",
					    MB_TAG_DENSE, MB_TYPE_DOUBLE);

  m_GeomTopoTool = new GeomTopoTool(reinterpret_cast<MBInterface*>(mesh));
#endif
}

CutCellMesh::~CutCellMesh()
{
}

int CutCellMesh::do_mesh(int exp, int file_type)
{
 #ifdef MOAB
  clock_t time1 = clock();
  struct rusage r_usage;
  long int rss1, rss2, rss3, rss4, rss5, rss6, rss7, rss8,
    msize1, msize2, msize3, msize4, msize5, msize6, msize7, msize8;
  static size_t pagesize = getpagesize();
  util_getrusage(r_usage);
  rss1 = r_usage.ru_maxrss*pagesize;
  msize1 = r_usage.ru_idrss*pagesize;
  MBErrorCode result;

  if (debug) {
    result = moab_instance()->write_mesh("input.vtk");
    if (MB_SUCCESS != result) {
      std::cerr << "Couldn't write input mesh.";
      return iBase_FAILURE;
    }
  }

  // 1. construct obb tree for all surfaces and volumes
  int err = construct_obb_tree();
  ERRORR("Couldn't construct obb tree.", err);
  clock_t time2 = clock();
  util_getrusage(r_usage);
  rss2 = r_usage.ru_maxrss*pagesize;
  msize2 = r_usage.ru_idrss*pagesize;

  // 2. set division
  err = set_division();
  ERRORR("Couldn't set division.", err);

  // 3. make hex vertices
  err = make_hex_vertices();
  ERRORR("Couldn't make hex vertices.", err);
  clock_t time3 = clock();
  util_getrusage(r_usage);
  rss3 = r_usage.ru_maxrss*pagesize;
  msize3 = r_usage.ru_idrss*pagesize; 
  
  // 4. make hexes
  err = make_hexes();
  ERRORR("Couldn't make hexes.", err);
  clock_t time4 = clock();
  util_getrusage(r_usage);
  rss4 = r_usage.ru_maxrss*pagesize;
  msize4 = r_usage.ru_idrss*pagesize;

  // 5. find intersected geometry surfaces by rays
  err = find_intersections();
  ERRORR("Couldn't find intersected surfaces.", err);
  clock_t time5 = clock();
  util_getrusage(r_usage);
  rss5 = r_usage.ru_maxrss*pagesize;
  msize5 = r_usage.ru_idrss*pagesize;
  
  // 6. set hex status and boundary hex cut fraction info
  err = set_tag_info();
  ERRORR("Couldn't set tag infor.", err);
  clock_t time6 = clock();
  util_getrusage(r_usage);
  rss6 = r_usage.ru_maxrss*pagesize;
  msize6 = r_usage.ru_idrss*pagesize;
#endif

  if (exp) {
    err = export_mesh(file_type);
    ERRORR("Couldn't print debug info.", err);
  }
  clock_t time7 = clock();
  util_getrusage(r_usage);
  rss7 = r_usage.ru_maxrss*pagesize;
  msize7 = r_usage.ru_idrss*pagesize;
  
  std::cout << "OBB_tree_construct_time: "
	    << (double) (time2 - time1)/CLOCKS_PER_SEC
	    << ", preparation time: "
	    << (double) (time3 - time2)/CLOCKS_PER_SEC
	    << ", hex_construct_time: "
	    << (double) (time4 - time3)/CLOCKS_PER_SEC
	    << ", intersection_time: "
	    << (double) (time5 - time4)/CLOCKS_PER_SEC
	    << ", set_info_time: "
	    << (double) (time6 - time5)/CLOCKS_PER_SEC
	    << ", export_time: "
	    << (double) (time7 - time6)/CLOCKS_PER_SEC
	    << std::endl;

  std::cout << "start_memory: " << rss1 << " " << msize1
	    << ", OBB_tree_construct_moemory: " << rss2 << " " << msize2
	    << ", preparation_moemory: " << rss3 << " " << msize3
	    << ",hex_construct_memory : " << rss4 << " " << msize4
	    << ", intersection_memory: " << rss5 << " " << msize5
	    << ", set_info_memory: " << rss6 << " " << msize6
	    << ", export_memory: " << rss7 << " " << msize7
	    << std::endl;
  
  return iBase_SUCCESS;
}

int CutCellMesh::construct_obb_tree()
{
  // construct obb tree for geometry surfaces and volumes by GeomTopoTool
  MBErrorCode rval = m_GeomTopoTool->construct_obb_trees(true);
  MBERRORR("Couldn't construct obb tree in GeomTopoTool.", iBase_ERROR_MAP[rval]);

  m_hObbTree = m_GeomTopoTool->obb_tree();
  m_hTreeRoot = m_GeomTopoTool->get_one_vol_root();
  
  return iBase_SUCCESS;
}

int CutCellMesh::export_mesh(int file_type)
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
  ERRORR("Failed to get hexes.\n", err); 

  int* hex_status = new int[hex_size];
  int hex_status_alloc = sizeof(int)*hex_size;
  int hex_status_size = 0;
  iMesh_getIntArrData(m_mesh, hex_handles,
		      hex_size, m_elemStatusTag,
		      &hex_status, &hex_status_alloc,
		      &hex_status_size, &err);
  ERRORR("Failed to get hex status.\n", err);
  clock_t time2 = clock();
  clock_t time3;

  if (debug) {
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
    else ERRORR("Hex element status should be inside/outside/boundary.\n", 1);
    }
    
    std::cout << "# of inside hex:" << n_inside_hex
	      << ", # of outside hex:" << n_outside_hex
	      << ", # of boundary hex:" << n_boundary_hex
	      << ", geom vol:"
	      << n_inside_hex*m_dIntervalSize[0]*m_dIntervalSize[1]*m_dIntervalSize[2]
	      << ", vox vol:" << hex_size*m_dIntervalSize[0]*m_dIntervalSize[1]*m_dIntervalSize[2]
	      << std::endl;
    time3 = clock();

    // save inside/outside/boundary elements separately
    if (n_inside_hex > 0) {
      err = write_mesh(0, file_type, &insideHex[0], n_inside_hex);
      ERRORR("Couldn't write inside mesh.", err);
    }
    
    if (n_outside_hex > 0) {
      err = write_mesh(1, file_type, &outsideHex[0], n_outside_hex);
      ERRORR("Couldn't write outside mesh.", err);
    }
    
    if (n_boundary_hex > 0) {
      err = write_mesh(2, file_type, &bndrHex[0], n_boundary_hex);
      ERRORR("Couldn't write boundary mesh.", err);
    }
    
    std::cout << "hex_handle_get_time: "
	      << (double) (time2 - time1)/CLOCKS_PER_SEC
	      << ", separate_write_time: "
	      << (double) (time3 - time2)/CLOCKS_PER_SEC
	      << std::endl;
  }
  else {
    err = write_mesh(3, file_type, hex_handles, hex_size);
    ERRORR("Couldn't write all hex mesh.", err);
  }

  free(hex_handles);
  free(hex_status);

  return iBase_SUCCESS;
}

int CutCellMesh::write_mesh(int type, int file_type,
			    iBase_EntityHandle* handles, int& n_elem)
{
  clock_t time1 = clock();
  int is_list = 1, err;
  MBErrorCode result;
  iBase_EntitySetHandle set;

  iMesh_createEntSet(m_mesh, is_list, &set, &err);
  ERRORR("Couldn't create set.\n", err);

  iMesh_addEntArrToSet(m_mesh, handles, n_elem, set, &err);
  ERRORR("Couldn't add hexes to set.\n", err);
  clock_t time2 = clock();

  std::string file_name;
  std::stringstream ss;
  if (type == 0) ss << "inside";
  else if (type == 1) ss << "outside";
  else if (type == 2) ss << "boundary";
  else if (type == 3) ss << "all_mesh";
  ss << m_dInputSize;
  ss >> file_name;
  if (file_type) file_name += std::string(".h5m");
  else file_name += std::string(".vtk");
  
  result = moab_instance()->write_mesh(file_name.c_str(),
				       (const MBEntityHandle*) &set, 1);
  if (MB_SUCCESS != result) {
    std::cerr << "Failed to write hex mesh." << std::endl;
    return iBase_FAILURE;
  }
  std::cout << "Elements are exported." << std::endl;
  clock_t time3 = clock();

  std::cout << "set_creation_time: "
	    << (double) (time2 - time1)/CLOCKS_PER_SEC
	    << ", write_time: "
	    << (double) (time3 - time2)/CLOCKS_PER_SEC
	    << std::endl;

  return iBase_SUCCESS;
}

int CutCellMesh::make_hex_vertices()
{
  int nNode = m_nNode[0]*m_nNode[1]*m_nNode[2];
  double* vertexCoords = new double[3*nNode];

  // make vertices of division
  int iNode;
  for (int k = 0; k < m_nNode[2]; k++) {
    for (int j = 0; j < m_nNode[1]; j++) {
      for (int i = 0; i < m_nNode[0]; i++) {
	iNode = 3*(k*m_nNode[0]*m_nNode[1] + j*m_nNode[0] + i);
	vertexCoords[iNode++] = m_origin_coords[0] + i*m_dIntervalSize[0];
	vertexCoords[iNode++] = m_origin_coords[1] + j*m_dIntervalSize[1];
	vertexCoords[iNode] = m_origin_coords[2] + k*m_dIntervalSize[2];
      }
    }
  }

  int err;
  int vertex_alloc = sizeof(iBase_EntityHandle)*nNode;
  int vertex_size = nNode;
  m_vhVertex.resize(nNode, NULL);
  iBase_EntityHandle* pVertexHandle = &m_vhVertex[0];

  iMesh_createVtxArr(m_mesh, nNode,
		     iBase_INTERLEAVED, vertexCoords, 3*nNode,
		     &pVertexHandle,
		     &vertex_alloc, &vertex_size, &err);
  delete [] vertexCoords;
  ERRORR("Failed to create vertices.\n", err);

  return iBase_SUCCESS;
}

int CutCellMesh::make_hexes()
{
  // get hex connectivity
  int err, iDiv, iNode;
  int nConn = 8*m_nHex;
  std::vector<iBase_EntityHandle> conn;
  conn.resize(nConn);

  for (int k = 0; k < m_nDiv[2]; k++) {
    for (int j = 0; j < m_nDiv[1]; j++) {
      iDiv = k*m_nDiv[0]*m_nDiv[1] + j*m_nDiv[0];
      for (int i = 0; i < m_nDiv[0]; i++) {
	conn[8*(iDiv + i)] = m_vhVertex[k*m_nNode[0]*m_nNode[1] + j*m_nNode[0] + i];
	conn[8*(iDiv + i) + 1] = m_vhVertex[k*m_nNode[0]*m_nNode[1] + j*m_nNode[0] + i + 1];
	conn[8*(iDiv + i) + 2] = m_vhVertex[k*m_nNode[0]*m_nNode[1] + (j + 1)*m_nNode[0] + i + 1];
	conn[8*(iDiv + i) + 3] = m_vhVertex[k*m_nNode[0]*m_nNode[1] + (j + 1)*m_nNode[0] + i];
	conn[8*(iDiv + i) + 4] = m_vhVertex[(k + 1)*m_nNode[0]*m_nNode[1] + j*m_nNode[0] + i];
	conn[8*(iDiv + i) + 5] = m_vhVertex[(k + 1)*m_nNode[0]*m_nNode[1] + j*m_nNode[0] + i + 1];
	conn[8*(iDiv + i) + 6] = m_vhVertex[(k + 1)*m_nNode[0]*m_nNode[1] + (j + 1)*m_nNode[0] + i + 1];
	conn[8*(iDiv + i) + 7] = m_vhVertex[(k + 1)*m_nNode[0]*m_nNode[1] + (j + 1)*m_nNode[0] + i];
      }
    }
  }

  // create hexes
  int entity_handles_allocated = sizeof(iBase_EntityHandle)*m_nHex,
    entity_handles_size = m_nHex,
    statusAllocated = m_nHex, statusSize = m_nHex;
  m_vhHex.resize(m_nHex, NULL);
  int* status = new int[m_nHex];
  iBase_EntityHandle* pEntityHandles = &m_vhHex[0];
  
  iMesh_createEntArr(m_mesh, iMesh_HEXAHEDRON, &conn[0], nConn,
		     &pEntityHandles,
		     &entity_handles_allocated,
		     &entity_handles_size,
		     &status,
		     &statusAllocated, &statusSize, &err);
  m_iStartHex = moab_instance()->id_from_handle(reinterpret_cast<MBEntityHandle>(m_vhHex[0]));

  delete [] status;
  ERRORR("Failed to create elements.", err);

  return iBase_SUCCESS;
}

EdgeStatus CutCellMesh::getEdgeStatus(const double dP)
{
  if (m_nStatus == INSIDE) { // previous inside
    if (dP < m_dSecondP) {
      m_nMove = 0;
      return INSIDE;
    }
    else {
      if (is_shared_overlapped_surf(m_iInter - 1)) {
	m_nMove = 1;
      }
      else m_nMove = 2;
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
      else if (is_shared_overlapped_surf(m_iInter - 1)) {
	m_nMove = 1;
      }
      else m_nMove = 2;
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
	if (is_shared_overlapped_surf(m_iInter - 1)) {
	  m_nMove = 1;
	}
	else m_nMove = 2;
	return BOUNDARY;
      }
    }
  }
  else {
    std::cerr << "Couldn't get edge status." << std::endl;
    return INSIDE;
  }
}

bool CutCellMesh::set_neighbor_hex_status(int dir, int i, int j, int k)
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

bool CutCellMesh::set_hex_status(int index, EdgeStatus value, int dir)
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

bool CutCellMesh::set_edge_status(int dir)
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

  m_nFraction += vdCutFraction.size();

  return true;
}

#ifdef MOAB
iBase_TagHandle CutCellMesh::get_tag(const char* name, int size,
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

iBase_TagHandle CutCellMesh::get_various_length_tag(const char* name,
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

int CutCellMesh::set_tag_info()
{
  int i, j, k, err;
  iMesh_setIntArrData(m_mesh, &m_vhHex[0], m_nHex, m_elemStatusTag,
		      &m_vnHexStatus[0], m_nHex, &err);
  ERRORR("Failed to set hex element status data.", err);

  // set cut fraction info to boundary hexes
  int nBndrHex = m_mdCutFraction.size();
  std::vector<iBase_EntityHandle> hvBndrHex(nBndrHex);
  int* frac_size = new int[nBndrHex];
  int* frac_leng = new int[3*nBndrHex];
  double dFrac;
  int nFracSize, nTempSize, nFracLeng, iHex, ii, jj, kk, nTolFrac = 0;
  int nDoubleSize = sizeof(double);
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.begin();
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  for (i = 0; iter != end_iter; iter++, i++) {
    iHex = iter->first;
    hvBndrHex[i] = m_vhHex[iHex];

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

int CutCellMesh::set_division()
{
  double box_center[3], box_axis1[3], box_axis2[3], box_axis3[3];
  MBErrorCode rval = m_hObbTree->box(m_hTreeRoot, box_center, box_axis1, box_axis2, box_axis3);
  MBERRORR("Error getting box information for tree root set.", iBase_ERROR_MAP[rval]);

  double lengths[3] = { axis_len(box_axis1), axis_len(box_axis2),
			axis_len(box_axis3) };
  
  // default value is adjusted to large geometry file
  // interval_size_estimate : 2*L*sqrt(2*PI*sqrt(2)/# of tris)
  if (m_dInputSize < 0.) {
    int n_tri;
    rval = moab_instance()->
      get_number_entities_by_dimension(reinterpret_cast<MBEntityHandle>(m_hRootSet),
				       2, n_tri);
    MBERRORR("Failed to get number of triangles.", iBase_ERROR_MAP[rval]);
    
    double box_length_ave = 2./3.*(lengths[0] + lengths[1] + lengths[2]);
    m_dInputSize = 2.*box_length_ave*sqrt(8.886/n_tri);
  }
    
  for (int i = 0; i < 3; i++) {
    m_nDiv[i] = 2.*lengths[i]/m_dInputSize;
    if (m_nDiv[i] < 5) m_nDiv[i] = 5;
    m_dIntervalSize[i] = m_dInputSize;
    if (m_nDiv[i]*.07 > 3) m_nDiv[i] += m_nDiv[i]*.07;
    else m_nDiv[i] += 3;
    m_nNode[i] = m_nDiv[i] + 1;
    m_origin_coords[i] = box_center[i] - .5*m_nDiv[i]*m_dIntervalSize[i];
  }

  m_nHex = m_nDiv[0]*m_nDiv[1]*m_nDiv[2];

  std::cout << "# of hex: " << m_nHex << ", interval size: " << m_dInputSize << std::endl;

  return iBase_SUCCESS;
}

int CutCellMesh::find_intersections()
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

int CutCellMesh::fire_rays(int dir)
{
  // ray fire
  int i, j, k, l, index[3];
  double tolerance = 1e-12;
  double rayLength = m_nDiv[dir]*m_dIntervalSize[dir];
  int err, iNodeStart, iNodeEnd, nIntersect, nNodeSlice;
  double startPnt[3], endPnt[3], rayDir[3];
  MBErrorCode rVal;
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
	rayDir[0] = 1.;
	rayDir[1] = 0.;
	rayDir[2] = 0.;
      }
      else if (dir == 1) { // y coordinate ray
	iNodeStart = i*m_nNode[otherDir2]*m_nNode[dir] + j;
	iNodeEnd = iNodeStart + m_nNode[otherDir2]*(m_nNode[dir] - 1);
	nNodeSlice = m_nNode[otherDir2];
	index[0] = j;
	index[1] = 0;
	index[2] = i;
	rayDir[0] = 0.;
	rayDir[1] = 1.;
	rayDir[2] = 0.;
      }
      else if (dir == 2) { // z coordinate ray
	iNodeStart = j*m_nNode[otherDir1] + i;
	iNodeEnd = iNodeStart + m_nNode[otherDir1]*m_nNode[otherDir2]*(m_nNode[dir] - 1);
	nNodeSlice = m_nNode[otherDir1]*m_nNode[otherDir2];
	index[0] = i;
	index[1] = j;
	index[2] = 0;
	rayDir[0] = 0.;
	rayDir[1] = 0.;
	rayDir[2] = 1.;
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

      if (dir == 2 && i == 8 && j == 10) {
	std::cout << "stop" << std::endl;
      }
      
      // ray-tracing
      m_vIntersection.clear();
      m_vhInterSurf.clear();
      m_vhInterFacet.clear();
      m_mhOverlappedSurf.clear();
      std::vector<double> temp_intersects;
      rVal = m_hObbTree->ray_intersect_sets(temp_intersects, m_vhInterSurf,
					    m_vhInterFacet, m_hTreeRoot, tolerance,
					    -1, startPnt, rayDir, &rayLength);
      nIntersect = temp_intersects.size();
      if (MB_SUCCESS != rVal) {
	std::cerr << "Failed : ray-triangle intersection." << std::endl;
	return iBase_FAILURE;
      }
      
      // put edge status info as tag
      k = 0;
      m_iInter = 0;

      if (nIntersect > 0) {
	bool bMoveRay = false;
	m_vIntersection.resize(nIntersect);

	for (l = 0; l < nIntersect; l++) {
	  IntersectDist temp_inter_dist(temp_intersects[l], l);
	  m_vIntersection[l] = temp_inter_dist;
	}
	std::sort(m_vIntersection.begin(), m_vIntersection.end(), less_intersect);

	if (nIntersect > 2) { // when ray intersect shared edge of triangles
	  bool bMoveOnce;
	  m_nIteration = 0;
	  m_iOverlap = 0;
	  if (is_ray_move_and_set_overlap_surf(startPnt, endPnt, bMoveOnce)) {
	    if (!move_ray(startPnt, endPnt, tolerance, dir, bMoveOnce)) {
	      std::cerr << "Number of Intersection between edges and ray should be even." << std::endl;
	      return iBase_FAILURE;
	    }
	  }

	  std::vector<IntersectDist>::iterator end = m_vIntersection.end();
	  std::vector<IntersectDist>::iterator new_end = std::unique(m_vIntersection.begin(),
								     end, equal_intersect);
	  m_vIntersection.erase(new_end, end);

	  nIntersect = m_vIntersection.size();
	}

	if (nIntersect > 0) {
	  m_nStatus = OUTSIDE;
	  m_dFirstP = startPnt[dir] + m_vIntersection[m_iInter++].distance;
	  m_dSecondP = startPnt[dir] + m_vIntersection[m_iInter++].distance;
	  m_prevPnt = startPnt[dir];

	  for (; k < m_nNode[dir] - 1; k++) {
	    m_curPnt = startPnt[dir] + (k + 1)*m_dIntervalSize[dir];
	    m_nStatus = getEdgeStatus(m_curPnt);
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

bool CutCellMesh::move_intersections(int n_dir, int n_inter, double start_pnt[3])
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

bool CutCellMesh::is_shared_overlapped_surf(int index)
{
  int nParent, err;
  MBEntityHandle hSurf = m_vhInterSurf[m_vIntersection[index].index];
  iMesh_getNumPrnt(m_mesh,
		   reinterpret_cast<iBase_EntitySetHandle> (hSurf),
		   1, &nParent, &err);
  ERRORRF("Failed to get number of surface parents.\n");

  if (nParent > 1) return true;

  return m_mhOverlappedSurf.count(hSurf) > 0;
}

bool CutCellMesh::get_grid_and_edges_techX(double* boxMin, double* boxMax, int* nDiv,
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
  std::cout << "# of inside, outside, boundary : " << nInside
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
    dFrac[3] = cutFrac.fraction[2][0];
    if (cutFrac.fraction[2].size() > 0) dFrac[3] = cutFrac.fraction[2][0];
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
    dFrac[3] = cutFrac.fraction[1][0];
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

bool CutCellMesh::get_grid_and_edges(double* boxMin, double* boxMax, int* nDiv,
				     std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellEdge,
				     std::vector<int>& rvnInsideCell, bool isCornerExterior)
{
  int i, err, ii, jj, kk, iHex;
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

bool CutCellMesh::get_inside_boundary_hex(std::vector<int>& rvnInsideCell)
{
  int i, j, err, ii, jj, kk, iHex;
  
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
  std::cout << "# of inside, outside, boundary : " << nInside
	    << ", " << nOutside << ", " << nBoundary << std::endl;

  return true;
}

bool CutCellMesh::export_fraction_edges(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellSurfEdge)
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

  int err;
  int is_list = 1;
  MBErrorCode result;
  iBase_EntitySetHandle set;
  
  iMesh_createEntSet(m_mesh, is_list, &set, &err);
  ERRORR("Couldn't create set.\n", err);
  
  iMesh_addEntArrToSet(m_mesh, &edge_handles[0], edge_handles.size(), set, &err);
  ERRORR("Couldn't add edges to set.\n", err);
  
  result = moab_instance()->write_mesh("edges.vtk",
				       (const MBEntityHandle*) &set, 1);
  if (MB_SUCCESS != result) {
    std::cerr << "Failed to write edges." << std::endl;
    return false;
  }

  return true;
}

bool CutCellMesh::export_fraction_points(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellEdge)
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

bool CutCellMesh::make_edge(double ePnt[6], std::vector<iBase_EntityHandle>& edge_handles)
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

double CutCellMesh::get_edge_fraction(int idHex, int dir)
{
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.find(idHex);
  if (iter != end_iter) return iter->second.fraction[dir][0];
  else return -1.;
}
#endif

double CutCellMesh::get_uncut_edge_fraction(int i, int j, int k, int dir)
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

bool CutCellMesh::move_ray(double* startPnt, double* endPnt,
			   double tol, int dir, bool bSpecialCase)
{
  //int i, nIteration;
  int i;
  bool bMove = true, temp_special;

  while (bMove) {
    if (m_nIteration > 10) {
      if (bSpecialCase) return true; // special case
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
    ray.normalize();
    m_vIntersection.clear();
    m_vhInterSurf.clear();
    m_vhInterFacet.clear();
    
    std::vector<double> temp_intersects;
    MBErrorCode rVal = m_hObbTree->ray_intersect_sets(temp_intersects, m_vhInterSurf,
						      m_vhInterFacet, m_hTreeRoot, tol,
						      -1, startPnt, ray.array(), &rayLength);
    if (MB_SUCCESS != rVal) {
      std::cerr << "Failed : ray-triangle intersection." << std::endl;
      return false;
    }

    int nInter = temp_intersects.size();
    m_vIntersection.resize(nInter);
    for (i = 0; i < nInter; i++) {
      IntersectDist temp_inter_dist(temp_intersects[i], i);
      m_vIntersection[i] = temp_inter_dist;
    }
    std::sort(m_vIntersection.begin(), m_vIntersection.end(), less_intersect);
    
    bMove = is_ray_move_and_set_overlap_surf(startPnt, endPnt, bSpecialCase);

    m_nIteration++;
  }

  std::cout << "ray is moved successfully." << std::endl;

  return true;
}

bool CutCellMesh::is_ray_move_and_set_overlap_surf(double* startPnt, double* endPnt, bool& bSpecialCase)
{
  int j, k, err, nInter = m_vIntersection.size() - 1;

  while (m_iOverlap < nInter) {
    if (m_vIntersection[m_iOverlap + 1].distance - m_vIntersection[m_iOverlap].distance < 1e-7) {
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
    else m_iOverlap++;
  }

  return false;
}

// get resource
void CutCellMesh::util_getrusage(struct rusage &r_usage)
{
  getrusage(RUSAGE_SELF, &r_usage);
  
  // this machine doesn't return rss - try going to /proc
  // print the file name to open
  if (r_usage.ru_maxrss == 0) {
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    close(file_ptr);
    file_str[file_len] = '\0';
    
    // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    static int page_size = getpagesize();
    int num_fields = sscanf(file_str, 
			    "%d " // pid
			    "%s " // comm
			    "%c " // state
			    "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
			    "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
			    "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
			    "%u %u " // timeout, itrealvalue
			    "%d " // starttime
			    "%u %u", // vsize, rss
			    &dum_int, 
			    dum_str, 
			    dum_str, 
			    &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
			    &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
			    &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
			    &dum_uint, &dum_uint, 
			    &dum_int,
			    &vm_size, &rss);
    if (num_fields == 24) {
      r_usage.ru_maxrss = rss/page_size;
      r_usage.ru_idrss = vm_size/page_size;
    }
  }
}
