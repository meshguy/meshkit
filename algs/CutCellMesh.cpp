#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <sstream>
#include <algorithm>

#include "CutCellMesh.hpp"

#ifdef MOAB
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBOrientedBox.hpp" // MBMatrix3.hpp is included
#endif

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
#define ERRORRF(a) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return false;}}

const bool debug = false;
static bool equal_to(double d1, double d2) { return abs(d1 - d2) < 10e-7; }

CutCellMesh::CutCellMesh(iMesh_Instance mesh, iBase_EntitySetHandle root_set, double size)
  : m_mesh(mesh), m_hRootSet(root_set), m_dInputSize(size)
{
  m_elem_status_tag_handle = NULL;
  m_cut_fraction_tag_handle = NULL;
  m_nStatus = OUTSIDE;
  m_iStartHex = 0;
}

CutCellMesh::~CutCellMesh()
{
}

int CutCellMesh::do_mesh(int exp, int file_type)
{
 #ifdef MOAB
  
  clock_t time1 = clock();
  struct rusage r_usage;
  long int rss1, rss2, rss3, rss4, rss5, rss7, rss8,
    msize1, msize2, msize3, msize4, msize5, msize7, msize8;
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
  
  // 1. get all triangles
  int err;
  MBRange tris;
  result = moab_instance()->
    get_entities_by_dimension(reinterpret_cast<MBEntityHandle>(m_hRootSet), 2, tris);
  if (MB_SUCCESS != result) {
    std::cerr << "Couldn't get triangles.";
    return iBase_FAILURE;
  }
  m_nTri = tris.size();
  clock_t time2 = clock();
  util_getrusage(r_usage);
  rss2 = r_usage.ru_maxrss*pagesize;
  msize2 = r_usage.ru_idrss*pagesize;

  // 2. build OBB trees for all triangles
  MBOrientedBoxTreeTool tool(reinterpret_cast<MBInterface*> (m_mesh));
  if (tool.build(tris, m_hTreeRoot) != MB_SUCCESS) {
    std::cerr << "Could'nt build tree." << std::endl;    
    return iBase_FAILURE;
  }
  clock_t time3 = clock();
  util_getrusage(r_usage);
  rss3 = r_usage.ru_maxrss*pagesize;
  msize3 = r_usage.ru_idrss*pagesize;
  
  // 3. build box
  MBOrientedBox box;
  result = tool.box(m_hTreeRoot, box);
  if (MB_SUCCESS != result) {
    std::cerr << "Error getting box for tree root set.";
    return iBase_FAILURE;
  }
  clock_t time4 = clock();
  util_getrusage(r_usage);
  rss4 = r_usage.ru_maxrss*pagesize;
  msize4 = r_usage.ru_idrss*pagesize;

  // 4. set division
  set_division(box);

  // 5. make hex vertices
  err = make_hex_vertices();
  ERRORR("Couldn't make hex vertices.", err);
  
  // 6. make hexes
  err = make_hexes();
  ERRORR("Couldn't make hexes.", err);
  clock_t time5 = clock();
  util_getrusage(r_usage);
  rss5 = r_usage.ru_maxrss*pagesize;
  msize5 = r_usage.ru_idrss*pagesize;
  
  // 7. find intersected geometry surfaces by rays
  err = find_intersections(tool);
  ERRORR("Couldn't find intersected surfaces.", err);
  clock_t time6 = clock();
  
  // 8. set hex status and boundary hex cut fraction info
  err = set_tag_info();
  ERRORR("Couldn't set tag infor.", err);
  clock_t time7 = clock();
  util_getrusage(r_usage);
  rss7 = r_usage.ru_maxrss*pagesize;
  msize7 = r_usage.ru_idrss*pagesize;
#endif

  if (exp) {
    err = export_mesh(file_type);
    ERRORR("Couldn't print debug info.", err);
  }
  clock_t time8 = clock();
  util_getrusage(r_usage);
  rss8 = r_usage.ru_maxrss*pagesize;
  msize8 = r_usage.ru_idrss*pagesize;

  std::cout << "triangle_construct_time: "
	    << (double) (time2 - time1)/CLOCKS_PER_SEC
	    << ", OBB_tree_construct_time: "
	    << (double) (time3 - time2)/CLOCKS_PER_SEC
	    << ", box_building_time: "
	    << (double) (time4 - time3)/CLOCKS_PER_SEC
	    << ", hex_construct_time: "
	    << (double) (time5 - time4)/CLOCKS_PER_SEC
	    << ", intersection_time: "
	    << (double) (time6 - time5)/CLOCKS_PER_SEC
	    << ", set_tag_info_time: "
	    << (double) (time7 - time6)/CLOCKS_PER_SEC
	    << ", mesh_set_info_time: "
	    << (double) (time7 - time3)/CLOCKS_PER_SEC
	    << ", export_time: "
	    << (double) (time8 - time7)/CLOCKS_PER_SEC
	    << std::endl;
  std::cout << "start_memory: " << rss1 << " " << msize1
	    << ", triangle_construct_memory: " << rss2 << " " << msize2
	    << ", OBB_tree_construct_moemory: " << rss3 << " " << msize3
	    << ", box_building_time: " << rss4 << " " << msize4
	    << ", hex_construct_time: " << rss5 << " " << msize5
	    << ", end_memory: " << rss7 << " " << msize7
	    << ", export_memory: " << rss8 << " " << msize8
	    << std::endl;

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
  
  char* hex_status = new char[hex_size];
  int hex_status_alloc = 0;
  int hex_status_size = 0;
  iMesh_getArrData(m_mesh, hex_handles,
		   hex_size, m_elem_status_tag_handle,
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
      iNode = 3*(k*m_nNode[0]*m_nNode[1] + j*m_nNode[0]);
      for (int i = 0; i < m_nNode[0]; i++) {
	vertexCoords[iNode + 3*i] = m_origin_coords[0] + i*m_dIntervalSize[0];
	vertexCoords[iNode + 3*i + 1] = m_origin_coords[1] + j*m_dIntervalSize[1];
	vertexCoords[iNode + 3*i + 2] = m_origin_coords[2] + k*m_dIntervalSize[2];
      }
    }
  }

  int err;
  int vertex_alloc = nNode;
  int vertex_size = nNode;
  m_vhVertex.resize(nNode, NULL);
  iBase_EntityHandle* pVertexHandle = &m_vhVertex[0];

  iMesh_createVtxArr(m_mesh, nNode,
		     iBase_BLOCKED, vertexCoords, 3*nNode,
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
      iNode = k*m_nNode[0]*m_nNode[1] + j*m_nNode[0];
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
  int entity_handles_allocated = m_nHex,
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

EdgeStatus CutCellMesh::getEdgeStatus(const double dP, bool& bMoveNext)
{
  if (m_nStatus == INSIDE) { // previous inside
    if (dP < m_dSecondP) {
      bMoveNext = false;
      return INSIDE;
    }
    else {
      bMoveNext = true;
      return BOUNDARY;
    }
  }
  else if (m_nStatus == OUTSIDE) { // previous outside
    if (dP < m_dFirstP) {
      bMoveNext = false;
      return OUTSIDE;
    }
    else {
      if (dP < m_dSecondP) bMoveNext = false;
      else {
	bMoveNext = true;
      }
      return BOUNDARY;
    }
  }
  else if (m_nStatus == BOUNDARY) { // previous boundary
    if (dP < m_dFirstP) {
      bMoveNext = false;
      return OUTSIDE;
    }
    else {
      if (dP < m_dSecondP) {
	bMoveNext = false;
	return INSIDE;
      }
      else {
	bMoveNext = true;
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
    //std::cout << "ijk: " << k << ", " << i << ", " << j
    //      << ", iElem: " << m_iElem << std::endl;
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
    //std::cout << "ijk: " << j << ", " << k << ", " << i
    //      << ", iElem: " << m_iElem << std::endl;
    if (!set_hex_status(m_iElem, m_nStatus, dir)) return false;
    iElem = m_iElem - 1;
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
    iElem -= m_nDiv[dir]*m_nDiv[otherDir2];
    if (!set_hex_status(iElem++, m_nStatus, dir)) return false;
    if (!set_hex_status(iElem, m_nStatus, dir)) return false;
  }
  else if (dir == 2) { // z coordinate ray
    m_iElem = k*m_nDiv[0]*m_nDiv[1] + j*m_nDiv[0] + i;
    //std::cout << "ijk: " << i << ", " << j << ", " << k
    //      << ", iElem: " << m_iElem << std::endl;
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
  double dCutFraction;
  if (m_bMove) dCutFraction = 1. - (m_curPnt - m_dSecondP)/m_dIntervalSize[dir];
  else dCutFraction = 1. - (m_curPnt - m_dFirstP)/m_dIntervalSize[dir];

  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.find(m_iElem);
  if (iter == m_mdCutFraction.end()) { // not exist
    CutFraction cFaction(dir, dCutFraction);
    m_mdCutFraction[m_iElem] = cFaction;
  }
  else { // exist
    iter->second.fraction[dir] = dCutFraction;
  }

  return true;
}

#ifdef MOAB
int CutCellMesh::set_tag_info()
{
  // set all hex status info as tag
  int i, err;
  MBTag new_tag = 0;
  char df_value = 1; // OUTSIDE
  MBErrorCode result = moab_instance()->tag_create("elem_status_tag",
						   sizeof(char), MB_TAG_DENSE,
						   new_tag, &df_value);
  if (MB_SUCCESS != result) {
    std::cerr << "Failed to create element status tag." << std::endl;
    return iBase_FAILURE;
  }
  m_elem_status_tag_handle = (iBase_TagHandle) new_tag;

  iMesh_setArrData(m_mesh, &m_vhHex[0], m_nHex, m_elem_status_tag_handle,
		   &m_vnHexStatus[0], m_nHex, &err);
  ERRORR("Failed to set hex element status data.", err);

  // set cut fraction info to boundary hexes
  std::vector<iBase_EntityHandle> hvBndrHex;
  int nBndrHex = m_mdCutFraction.size();
  hvBndrHex.resize(nBndrHex);
  double* daFraction = new double[3*nBndrHex];
  std::map<int, CutFraction>::iterator iter = m_mdCutFraction.begin();
  std::map<int, CutFraction>::iterator end_iter = m_mdCutFraction.end();
  for (i = 0; iter != end_iter; iter++, i++) {
    hvBndrHex[i] = m_vhHex[iter->first];
    for (int j = 0; j < 3; j++) {
      daFraction[3*i + j] = iter->second.fraction[j];
    }
  }

  iMesh_setDblArrData(m_mesh, &hvBndrHex[0], nBndrHex,
		      m_cut_fraction_tag_handle,
		      daFraction, 3*nBndrHex, &err);
  ERRORR("Failed to set cut fraction infor to hex.", err);

  delete [] daFraction;

  return iBase_SUCCESS;
}

void CutCellMesh::set_division(const MBOrientedBox& box)
{
  double box_length_ave = 2./3.*(box.length[0] + box.length[1] + box.length[2]);

  // default value is adjusted to large geometry file
  // interval_size_estimate : 2*L*sqrt(2*PI*sqrt(2)/# of tris)
  if (m_dInputSize < 0.) m_dInputSize = 2.*box_length_ave*sqrt(8.886/m_nTri);

  for (int i = 0; i < 3; i++) {
    m_nDiv[i] = 2.*box.length[i]/m_dInputSize;
    m_dIntervalSize[i] = 2.*box.length[i]/m_nDiv[i];
    if (m_nDiv[i]*.07 > 3) m_nDiv[i] += m_nDiv[i]*.07;
    else m_nDiv[i] += 3;
    m_nNode[i] = m_nDiv[i] + 1;
    m_origin_coords[i] = box.center[i] - .5*m_nDiv[i]*m_dIntervalSize[i];
  }

  m_nHex = m_nDiv[0]*m_nDiv[1]*m_nDiv[2];

  std::cout << "# of hex: " << m_nHex << ", interval size: " << m_dInputSize << std::endl;
}

int CutCellMesh::find_intersections(MBOrientedBoxTreeTool& tool)
{
  // create dense edge intersection tag
  int i, err;
  MBTag new_tag = 0;
  double df_value = 0.;
  MBErrorCode result = moab_instance()->tag_create("edge_cut_fraction_tag",
						   3*sizeof(double), MB_TAG_DENSE,
						   new_tag, &df_value);
  if (MB_SUCCESS != result) {
    std::cerr << "Failed to create edge cut fraction tag." << std::endl;
    return iBase_FAILURE;
  }
  m_cut_fraction_tag_handle = (iBase_TagHandle) new_tag;

  // initialize all hex as outside
  m_vnHexStatus.resize(m_nHex, 1);

  // fire rays to 3 directions
  for (i = 0; i < 3; i++) {
    m_vnEdgeStatus[i].resize(m_nDiv[i]*m_nDiv[(i + 1)%3]*m_nDiv[(i + 2)%3], OUTSIDE);
    err = fire_rays(tool, i);
    if (err != iBase_SUCCESS) return err;
  }

  return iBase_SUCCESS;
}

int CutCellMesh::fire_rays(MBOrientedBoxTreeTool& tool, int dir)
{
  // ray fire
  int i, j, k, n, index[3];
  double tolerance = 1e-20;
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

      for (n = 0; n < 3; n++) {
	if (n == dir) {
	  startPnt[n] = m_origin_coords[n];
	  endPnt[n] = m_origin_coords[n] + m_nDiv[dir]*m_dIntervalSize[n];
	}
	else {
	  startPnt[n] = m_origin_coords[n] + index[n]*m_dIntervalSize[n];
	  endPnt[n] = startPnt[n];
	}
      }
      
      // ray-tracing
      m_vdIntersection.clear();
      rVal = tool.ray_intersect_triangles(m_vdIntersection, m_hTreeRoot, tolerance,
					  startPnt, rayDir, &rayLength);
      if (MB_SUCCESS != rVal) {
	std::cerr << "Failed : ray-triangle intersection." << std::endl;
	return iBase_FAILURE;
      }
      
      // put edge status info as tag
      nIntersect = m_vdIntersection.size();
      k = 0;
      m_iInter = 0;

      if (nIntersect > 0) {
	std::sort(m_vdIntersection.begin(), m_vdIntersection.end());

	if (nIntersect % 2 != 0) { // when ray intersect shared edge of triangles
	  std::vector<double>::iterator end = m_vdIntersection.end();
	  std::vector<double>::iterator new_end = std::unique(m_vdIntersection.begin(), end, equal_to);
	  m_vdIntersection.erase(new_end, end);
	  
	  nIntersect = m_vdIntersection.size();

	  // when ray intersect region sharing many sliver triangles
	  if (nIntersect % 2 != 0) {
	    if (!move_ray(tool, startPnt, endPnt, tolerance, nIntersect, dir)) {
	      std::cerr << "Number of Intersection between edges and ray should be even." << std::endl;
	      return iBase_FAILURE;
	    }
	  }
	}

	m_nStatus = OUTSIDE;
	m_dFirstP = m_origin_coords[dir] + m_vdIntersection[m_iInter++];
	m_dSecondP = m_origin_coords[dir] + m_vdIntersection[m_iInter++];

	for (; k < m_nNode[dir] - 1; k++) {
	  m_curPnt =  m_origin_coords[dir] + (k + 1)*m_dIntervalSize[dir];
	  m_nStatus = getEdgeStatus(m_curPnt, m_bMove);
	  m_vnEdgeStatus[dir][i*m_nDiv[dir] + j*m_nDiv[dir]*m_nDiv[otherDir1] + k] = m_nStatus;

	  // set status of all hexes sharing the edge
	  if (!set_neighbor_hex_status(dir, i, j, k)) return iBase_FAILURE;

	  // set cut-cell edge status
	  if (m_nStatus == BOUNDARY && !set_edge_status(dir)) return iBase_FAILURE;

	  if (m_bMove) {
	    if (m_iInter < nIntersect) {
	      m_dFirstP = startPnt[dir] + m_vdIntersection[m_iInter++];
	      m_dSecondP = startPnt[dir] + m_vdIntersection[m_iInter++];
	    }
	    else {
	      k++;
	      break; // rest is all outside
	    }
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

bool CutCellMesh::get_grid_and_edges(double* boxMin, double* boxMax, int* nDiv,
				     std::map< CutCellSurfKey, std::vector<double>, LessThan >& rmdCutCellSurfEdge,
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
  char* hex_status = new char[hex_size];
  int hex_status_alloc = 0;
  int hex_status_size = 0;
  iMesh_getArrData(m_mesh, hex_handles,
		   hex_size, m_elem_status_tag_handle,
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
    dFrac[0] = cutFrac.fraction[1];
    dFrac[3] = cutFrac.fraction[2];
    dFrac[1] = get_edge_fraction(iHex + m_nDiv[0], 2);
    dFrac[2] = get_edge_fraction(iHex + m_nDiv[0]*m_nDiv[1], 1);
    
    if (dFrac[0] > 0. || dFrac[1] > 0. || dFrac[2] > 0. || dFrac[3] > 0.) { // if surface is cut
      CutCellSurfKey surfKey(ii, jj, kk, 0);
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
    dFrac[0] = cutFrac.fraction[0];
    dFrac[3] = cutFrac.fraction[2];
    dFrac[1] = get_edge_fraction(iHex + 1, 2);
    dFrac[2] = get_edge_fraction(iHex + m_nDiv[0]*m_nDiv[1], 0);
    
    if (dFrac[0] > 0. || dFrac[1] > 0. || dFrac[2] > 0. || dFrac[3] > 0.) { // if surface is cut
      CutCellSurfKey surfKey(ii, jj, kk, 1);
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
    dFrac[0] = cutFrac.fraction[0];
    dFrac[3] = cutFrac.fraction[1];
    dFrac[1] = get_edge_fraction(iHex + 1, 1);
    dFrac[2] = get_edge_fraction(iHex + m_nDiv[0], 0);

    if (dFrac[0] > 0. || dFrac[1] > 0. || dFrac[2] > 0. || dFrac[3] > 0.) { // if surface is cut
      CutCellSurfKey surfKey(ii, jj, kk, 2);
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

bool CutCellMesh::export_fraction_edges(std::map< CutCellSurfKey, std::vector<double>, LessThan >& mdCutCellSurfEdge)
{
  // export fractions as edge
  double curPnt[3], ePnt[6];
  std::vector<iBase_EntityHandle> edge_handles;
  std::map< CutCellSurfKey, std::vector<double>, LessThan >::iterator itr = mdCutCellSurfEdge.begin();
  std::map< CutCellSurfKey, std::vector<double>, LessThan >::iterator e_itr = mdCutCellSurfEdge.end();
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
}
 
bool CutCellMesh::make_edge(double ePnt[6], std::vector<iBase_EntityHandle>& edge_handles)
{
  int err, status;
  int vertex_alloc = 2;
  int vertex_size = 2;
  iBase_EntityHandle vertex_handle[2], edge_handle;
  //iBase_EntityHandle conn[2];
  iBase_EntityHandle* pVertexHandle = &vertex_handle[0];

  iMesh_createVtxArr(m_mesh, 2,
		     iBase_BLOCKED, ePnt, 6,
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
  if (iter != end_iter) return iter->second.fraction[dir];
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

// move ray when the first try is failed near extreme points
bool CutCellMesh::move_ray(MBOrientedBoxTreeTool& tool,
			   double* startPnt, double* endPnt,
			   double tol, int& nIntersect, int dir)
{
  int nIteration = 0;

  while (nIntersect % 2 == 0) {

    if (nIteration > 10) return false;

    m_vdIntersection.clear();
    for (int i = 0; i < 3; i++) {
      if (i != dir) {
	startPnt[i] += tol;
	endPnt[i] += tol;
      }
    }

    MBCartVect ray(endPnt[0] - startPnt[0], endPnt[1] - startPnt[1], endPnt[2] - startPnt[2]);
    double rayLength = ray.length();
    ray.normalize();
    m_vdIntersection.clear();
    
    MBErrorCode rVal = tool.ray_intersect_triangles(m_vdIntersection, m_hTreeRoot, tol,
					startPnt, ray.array(), &rayLength);
    if (MB_SUCCESS != rVal) {
      std::cerr << "Failed : ray-triangle intersection." << std::endl;
      return false;
    }

    nIntersect = m_vdIntersection.size();
    nIteration++;
  }

  std::cout << "ray is moved successfully." << std::endl;

  return true;
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
