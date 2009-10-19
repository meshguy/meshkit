#include <stdio.h>
#include <stdlib.h>
#include <iostream>

// TEST
#include <sstream>

#include <math.h>

#include "CutCellMesh.hpp"
#include "iMesh.h"

#ifdef DAGMC
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBOrientedBox.hpp" // MBMatrix3.hpp is included
#endif

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

const bool debug = true;

CutCellMesh::CutCellMesh(iMesh_Instance mesh, iBase_EntitySetHandle root_set)
  : m_mesh(mesh), m_hRootSet(root_set)
{
  /*
//#ifdef MOAB
  obbTag = get_tag( MB_OBB_TREE_TAG_NAME, sizeof(MBEntityHandle), MB_TAG_DENSE, MB_TYPE_HANDLE );
//  #endif*/
  m_nStatus = OUTSIDE;
}

CutCellMesh::~CutCellMesh()
{
}

int CutCellMesh::do_mesh()
{
 #ifdef DAGMC
  // get all triangles
  MBRange tris;
  MBErrorCode result = moab_instance()->
    get_entities_by_dimension(reinterpret_cast<MBEntityHandle>(m_hRootSet), 2, tris);
  if (MB_SUCCESS != result) return iBase_FAILURE;
  m_nTri = tris.size();

  // build OBB trees for all triangles
  MBOrientedBoxTreeTool tool(reinterpret_cast<MBInterface*> (m_mesh));
  if (tool.build(tris, m_hTreeRoot) != MB_SUCCESS) {
    std::cerr << "Could'nt build tree." << std::endl;    
    return iBase_FAILURE;
  }
  
  // build box
  MBOrientedBox box;
  result = tool.box(m_hTreeRoot, box);
  if (MB_SUCCESS != result) {
    std::cerr << "Error getting box for tree root set.";
    return iBase_FAILURE;
  }

  // set initial division
  set_initial_division(box);

  // make hex vertices
  int err = make_hex_vertices();
  ERRORR("Couldn't make hex vertices.", err);
  
  // make initial hexes
  err = make_initial_hexes();
  ERRORR("Couldn't make initial hexes.", err);

  // find intersected geometry surfaces by rays
  err = find_intersected_surfaces(tool);
  ERRORR("Couldn't find intersected surfaces.", err);
  
  // set hex status info as tag
  // ???? should be changed to dense char ????
  iBase_TagHandle elem_status_tag_handle = NULL;
  const char *tag_name = "elem_status_tag";
  int tag_name_size = 15;
  iMesh_createTag(m_mesh, tag_name, 1, iBase_INTEGER, 
                  &elem_status_tag_handle, &err, tag_name_size);
  ERRORR("Failed to create element status tag.", err);

  iMesh_setIntArrData(m_mesh, &m_vhHex[0], m_nHex, elem_status_tag_handle,
		      &m_vnHexStatus[0], m_nHex, &err);
  ERRORR("Failed to set hex element status data.", err);

  // TEST
  iBase_EntityHandle* hex_edges = NULL;
  int adj_entity_handles_allocated, adj_entity_handles_size;
  iMesh_getEntAdj(m_mesh, m_vhHex[0], iBase_EDGE, &hex_edges,
		  &adj_entity_handles_allocated, &adj_entity_handles_size,
                  &err);
  ERRORR("Failed to get edges of hex.", err);
  delete hex_edges;
#endif 
 
  return iBase_SUCCESS;
}

int CutCellMesh::make_hex_vertices()
{
  int nInitialNode = m_nNode[0]*m_nNode[1]*m_nNode[2];
  double* vertexCoords = new double[3*nInitialNode];

  // make vertices of initial division
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
  int vertex_alloc = nInitialNode;
  int vertex_size = nInitialNode;
  m_vhVertex.resize(nInitialNode, NULL);
  iBase_EntityHandle* pVertexHandle = &m_vhVertex[0];

  iMesh_createVtxArr(m_mesh, nInitialNode,
		     iBase_BLOCKED, vertexCoords, 3*nInitialNode,
		     &pVertexHandle,
		     &vertex_alloc, &vertex_size, &err);
  delete [] vertexCoords;
  ERRORR("Failed to create vertices.\n", err);

  return iBase_SUCCESS;
}

int CutCellMesh::make_initial_hexes()
{
  // make initial hexes
  int err, iDiv, iNode;
  m_nHex = m_nDiv[0]*m_nDiv[1]*m_nDiv[2];
  int nConn = 8*m_nHex;
  iBase_EntityHandle* conn = new iBase_EntityHandle[nConn];
  int ii = 0;
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
	
	for (int x = 0; x < 8; x++) {
	  double xx, yy, zz;
	  iBase_EntityHandle vertex_handle = conn[iDiv + 8*i + x];
	  
	  iMesh_getVtxCoord(m_mesh, vertex_handle, &xx, &yy, &zz, &err);
	  ERRORR("Failed to create m_vhVertex.", err);
	}
	ii++;
      }
    }
  }

  int entity_handles_allocated = m_nHex,
    entity_handles_size = m_nHex,
    statusAllocated = m_nHex, statusSize = m_nHex;
  m_vhHex.resize(m_nHex, NULL);
  int* status = new int[m_nHex];
  iBase_EntityHandle* pEntityHandles = &m_vhHex[0];
  
  iMesh_createEntArr(m_mesh, iMesh_HEXAHEDRON, conn, nConn,
		     &pEntityHandles,
		     &entity_handles_allocated,
		     &entity_handles_size,
		     &status,
		     &statusAllocated, &statusSize, &err);

  delete [] conn;
  delete [] status;
  ERRORR("Failed to create elements.", err);

  return iBase_SUCCESS;
}

EdgeStatus CutCellMesh::getEdgeStatus(const double dZ, bool bMoveNext)
{
  if (m_nStatus == INSIDE) { // previous inside
    if (dZ < m_dSecondZ) {
      bMoveNext = false;
      return INSIDE;
    }
    else {
      bMoveNext = true;
      return BOUNDARY;
    }
  }
  else if (m_nStatus == OUTSIDE) { // previous outside
    if (dZ < m_dFirstZ) {
      bMoveNext = false;
      return OUTSIDE;
    }
    else {
      if (dZ < m_dSecondZ) bMoveNext = false;
      else bMoveNext = true;
      return BOUNDARY;
    }
  }
  else if (m_nStatus == BOUNDARY) { // previous boundary
    if (dZ < m_dFirstZ) {
      bMoveNext = false;
      return OUTSIDE;
    }
    else {
      if (dZ < m_dSecondZ) {
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

bool CutCellMesh::set_hex_status(int index, int value)
{
  if (index < 0 || index > m_nHex - 1) {
    return false;
  }

  if (m_vnHexStatus[index] != 0) m_vnHexStatus[index] = value;
  return true;
}

#ifdef DAGMC
void CutCellMesh::set_initial_division(const MBOrientedBox& box)
{
  // get initial division
  // interval_size_estimate : L*sqrt(PI*2*sqrt(2)/# of tris)
  double box_length_ave = 2./3.*(box.length[0] + box.length[1] + box.length[2]);
  double interval_size_estimate = box_length_ave*sqrt(2*3.141592*sqrt(2)/m_nTri);

  for (int i = 0; i < 3; i++) {
    m_nDiv[i] = 2.*box.length[i]/interval_size_estimate; // twice of the # of division
    m_dIntervalSize[i] = 2.*box.length[i]/m_nDiv[i];
    m_nDiv[i] = m_nDiv[i] + 1;
    m_nNode[i] = m_nDiv[i] + 1;
    m_origin_coords[i] = box.center[i] - .5*m_nDiv[i]*m_dIntervalSize[i];
  }
}

int CutCellMesh::find_intersected_surfaces(MBOrientedBoxTreeTool& tool)
{
  // find intersected surfaces
  std::vector<double> distances;
  std::vector<MBEntityHandle> surfaces;
  double tolerance = 1e-6;
  double rayLength = m_nDiv[2]*m_dIntervalSize[2];
  int err, status, k, iNodeStart, iNodeEnd, nIntersect;
  double dIntersectZ;
  int m_nNodeSlice = m_nNode[0]*m_nNode[1];
  m_vnHexStatus.resize(m_nHex, 1); // initialize as outsize
  double startPnt[3], endPnt[3], curPnt[3];
  MBErrorCode rVal;
  bool bMove;
  int iElem;
  
  // create edge intersection tag
  iBase_TagHandle edge_intersection_tag_handle = NULL;
  const char *edge_intersection_tag_name = "edge_intersection_tag";
  int tag_name_size = 21;
  iMesh_createTag(m_mesh, edge_intersection_tag_name, 1, iBase_DOUBLE, 
                  &edge_intersection_tag_handle, &err, tag_name_size);
  ERRORR("Failed to create edge intersection tag.", err);

  for (int j = 1; j < m_nNode[1] - 1; j++) {
    for (int i = 1; i < m_nNode[0] - 1; i++) {
      
      // create edge entity
      iNodeStart = j*m_nNode[0] + i;
      iNodeEnd = m_nNode[0]*m_nNode[1]*(m_nNode[2] - 1) + j*m_nNode[0] + i;
      iBase_EntityHandle conn_edge[2] = {m_vhVertex[iNodeStart], m_vhVertex[iNodeEnd]};
      iBase_EntityHandle edgeHandle;
      iMesh_createEnt(m_mesh, iMesh_LINE_SEGMENT, conn_edge, 2,
		      &edgeHandle, &status, &err);
      ERRORR("Couldn't create edge element.", err);
		      
      // do ray tracing
      distances.clear();
      surfaces.clear();
      
      iMesh_getVtxCoord(m_mesh, m_vhVertex[iNodeStart], startPnt, startPnt + 1, startPnt + 2, &err);
      ERRORR("Couldn't get ray start vertex coordinates.", err);
      
      iMesh_getVtxCoord(m_mesh, m_vhVertex[iNodeEnd], endPnt, endPnt + 1, endPnt + 2, &err);
      ERRORR("Couldn't get ray end vertex coordinates.", err);
      
      MBCartVect ray(endPnt[0] - startPnt[0], endPnt[1] - startPnt[1], endPnt[2] - startPnt[2]);
      rayLength = ray.length();
      ray.normalize();
      m_vdIntersection.clear();


      rVal = tool.ray_intersect_triangles(m_vdIntersection, m_hTreeRoot, tolerance,
					       startPnt, ray.array(), &rayLength);
      if (MB_SUCCESS != rVal) {
	std::cerr << "Failed to build obb." << std::endl;
	return iBase_FAILURE;
      }
      
      // put edge status info as tag
      nIntersect = m_vdIntersection.size();
      k = 0;
      m_iInter = 0;
      if (nIntersect > 0) {
	/*
	if (nIntersect % 2 != 0) {
	  std::cerr << "Number of Intersection between edges and ray should be even." << std::endl;
	  return iBase_FAILURE;
	  }*/

	m_dFirstZ = startPnt[2] + m_vdIntersection[m_iInter++];
	m_dSecondZ = startPnt[2] + m_vdIntersection[m_iInter++];

	for (; k < m_nNode[2] - 1;k++) {
	  iMesh_getVtxCoord(m_mesh, m_vhVertex[iNodeStart + k*m_nNodeSlice + 2],
			    curPnt, curPnt + 1, curPnt + 2, &err);
	  ERRORR("Couldn't get cur vertex coordinates.", err);

	  m_nStatus = getEdgeStatus(curPnt[2], bMove);
	  iElem = k*m_nDiv[0]*m_nDiv[1] + j*m_nDiv[0] + i;

	  // set status of all hexes sharing the edge
	  if (m_nStatus == INSIDE) {
	    if (!set_hex_status(iElem--, -1)) return iBase_FAILURE;
	    if (!set_hex_status(iElem, -1)) return iBase_FAILURE;
	    iElem -= m_nDiv[0];
	    if (!set_hex_status(iElem++, -1)) return iBase_FAILURE;
	    if (!set_hex_status(iElem, -1)) return iBase_FAILURE;
	  }
	  else if (m_nStatus == OUTSIDE) {
	    if (!set_hex_status(iElem--, 1)) return iBase_FAILURE;
	    if (!set_hex_status(iElem, 1)) return iBase_FAILURE;
	    iElem -= m_nDiv[0];
	    if (!set_hex_status(iElem++, 1)) return iBase_FAILURE;
	    if (!set_hex_status(iElem, 1)) return iBase_FAILURE;
	  }
	  else { // boundary, set edge intersection info as tag

	    if (bMove) dIntersectZ = m_dSecondZ;
	    else dIntersectZ = m_dFirstZ;

	    iMesh_setDblData(m_mesh, edgeHandle,
			     edge_intersection_tag_handle,
			     dIntersectZ, &err);
	    ERRORR("Failed to set edge intersection position tag data.", err);

	    if (!set_hex_status(iElem--, 0)) return iBase_FAILURE;
	    if (!set_hex_status(iElem, 0)) return iBase_FAILURE;
	    iElem -= m_nDiv[0];
	    if (!set_hex_status(iElem++, 0)) return iBase_FAILURE;
	    if (!set_hex_status(iElem, 0)) return iBase_FAILURE;
	  }

	  if (bMove) {
	    if (m_iInter < nIntersect) {
	      m_dFirstZ = startPnt[2] + m_vdIntersection[m_iInter++]*rayLength;
	      m_dSecondZ = startPnt[2] + m_vdIntersection[m_iInter++]*rayLength;
	    }
	    else break; // rest is all outside
	  }
	}
      }
       
      // the rest are all outside
      //nEdgeStatus = 1;
      for (; k < m_nNode[2] - 1;k++) {
 	iElem = k*m_nDiv[0]*m_nDiv[1] + j*m_nDiv[0] + i;
	if (!set_hex_status(iElem--, 1)) return iBase_FAILURE;
	if (!set_hex_status(iElem, 1)) return iBase_FAILURE;
	iElem -= m_nDiv[0];
	if (!set_hex_status(iElem++, 1)) return iBase_FAILURE;
	if (!set_hex_status(iElem, 1)) return iBase_FAILURE;
      }
    }
  }
  
  return iBase_SUCCESS;
}
#endif
