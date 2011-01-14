// Author : Hong-Jun Kim
// Nov 16, 09
// It reads geometry file, makes facet triangles
// constructs tree structure for triangles and makes hexes bounding geometry
// With ray-tracing, find intersections and determine element inside/outside/boundary.
// Intersection fraction is stored to boundary elements.
// Element inside/outside/boundary status are stored as tag.

#ifndef EBMESH_HPP
#define EBMESH_HPP

#include <vector>
#include <map>
#include <sys/resource.h>

#include "iMesh.h"

#ifdef MOAB
#include "MBInterface.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBCartVect.hpp"
#include "MBOrientedBoxTreeTool.hpp"
#include "MBiMesh.hpp"
#endif

using moab::GeomTopoTool;

enum EdgeStatus {
  INSIDE,
  OUTSIDE,
  BOUNDARY
};

template <int size> struct EHARR
{ 
  iBase_EntityHandle h[size];
  iBase_EntityHandle& operator[](int i){ return h[i]; } 
  operator iBase_EntityHandle*() { return h; }
};

struct CutFraction {
  std::vector<double> fraction[3];

  CutFraction() {};

  CutFraction(int dir, const std::vector<double>& val) {
    add(dir, val);
  };
  
  void add(int dir, const std::vector<double>& val) {
    for (unsigned int i = 0; i < val.size(); i++) {
      fraction[dir].push_back(val[i]);
    }
  }
};

struct CutCellSurfEdgeKey {
  int i, j, k, l;

  CutCellSurfEdgeKey() {
    i = j = k = l = 0;
  };

  CutCellSurfEdgeKey(int ii, int jj, int kk, int ll) {
    i = ii;
    j = jj;
    k = kk;
    l = ll;
  };
};

struct IntersectDist {
  double distance;
  int index;

  IntersectDist() {};

  IntersectDist(double d, int i) {
    distance = d;
    index = i;
  };
};

struct VolFrac {
  double vol;
  bool closed;
  double ePnt[6];

  VolFrac() {};

  VolFrac(double f, int c) {
    vol = f;
    closed = c;
  };
};

struct LessThan
{
  bool operator() (const CutCellSurfEdgeKey key1, const CutCellSurfEdgeKey key2) const
  {
    if (key1.i < key2.i) return true;
    else if (key1.i > key2.i) return false;
    else {
      if (key1.j < key2.j) return true;
      else if (key1.j > key2.j) return false;
      else {
	if (key1.k < key2.k) return true;
	else if (key1.k > key2.k) return false;
	else {
	  if (key1.l < key2.l) return true;
	  else return false;
	}
      }
    }
  };
};

class EBMesh
{
public:

  EBMesh(iMesh_Instance mesh, iBase_EntitySetHandle root_set,
	 double size = -1., bool use_geom = true,
	 int add_layer = 3);
  
  ~EBMesh();

  // do EB mesh generation
  int do_mesh();

  // query function for techX
  bool get_grid_and_edges_techX(double* boxMin, double* boxMax, int* nDiv,
				std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellSurfEdge,
				std::vector<int>& rvnInsideCell, bool isCornerExterior = true);
  
  // query function to get multiple cut-cell edges
  bool get_grid_and_edges(double* boxMin, double* boxMax, int* nDiv,
			  std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellEdge,
			  std::vector<int>& rvnInsideCell, bool isCornerExterior = true);

  // get volume fraction for each material
  bool get_volume_fraction(int vol_frac_div);

  // export mesh to file
  bool export_mesh(const char* file_name, bool separate = false);
  
private:
  
  iBase_TagHandle m_elemStatusTag, m_edgeCutFracLengthTag,
    m_edgeCutFracTag, m_matFracIDTag, m_volFracTag;
  iMesh_Instance m_mesh;
  iBase_EntitySetHandle m_hRootSet;
  std::vector<IntersectDist> m_vIntersection;
  int m_nTri, m_nHex, m_iInter, m_nFraction, m_iStartHex, m_nMove, m_nAddLayer,
    m_nIteration, m_iOverlap, m_iElem, m_nNode[3], m_nDiv[3], m_nVolFracDiv,
    m_iFirstP, m_iSecondP;
  double m_dFirstP, m_dSecondP, m_curPnt, m_prevPnt,
    m_dIntervalSize[3], m_origin_coords[3], m_dInputSize;
  EdgeStatus m_nStatus;
  bool m_bMovedOnce, m_bUseGeom;
  std::vector<iBase_EntityHandle> m_vhVertex;
  std::vector<int> m_vnHexStatus;
  std::map<int, CutFraction> m_mdCutFraction;
  std::vector<EdgeStatus> m_vnEdgeStatus[3];
  //std::vector<double> m_vdCutFraction;
  
  EdgeStatus get_edge_status(const double dZ, int& iSkip);
  bool set_neighbor_hex_status(int dir, int i, int j, int k);
  bool set_hex_status(int index, EdgeStatus value, int dir);
  bool set_edge_status(int dir);
  int set_tag_info();
  int write_mesh(const char* file_name, int type,
		 iBase_EntityHandle* handles, int& n_elem);
  double get_edge_fraction(int idHex, int dir);
  double get_uncut_edge_fraction(int i, int j, int k, int dir);
  bool is_shared_overlapped_surf(int index);
  bool move_intersections(int n_dir, int n_inter, double start_pnt[3]);
  bool get_inside_boundary_hex(std::vector<int>& rvnInsideCell);

    // test functions
  bool export_fraction_edges(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellSurfEdge);
  bool export_fraction_points(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellEdge);
  bool make_edge(double ePnt[6], std::vector<iBase_EntityHandle>& edge_handles);

#ifdef MOAB
  MBInterface* moab_instance() {return reinterpret_cast<MBiMesh*> (m_mesh)->mbImpl;}
  iBase_TagHandle get_tag(const char* name, int size, MBTagType store, MBDataType type,
			  const void* def_value = NULL, bool create_if_missing = true);
  iBase_TagHandle get_various_length_tag(const char* name,
					 MBTagType store, MBDataType type);
  int make_scd_hexes(); // make structured hexes
  int make_uscd_hexes(); // make unstructured hexes
  int construct_obb_tree();
  int set_division();
  int find_intersections();
  int fire_rays(int dir);
  bool fire_ray(int& nIntersect, double startPnt[3],
		double endPnt[3], double tol, int dir,
		double rayLength);
  bool move_ray(int& nIntersect, double* startPnt, double* endPnt,
		double tol, int dir, bool bMoveOnce);
  bool is_ray_move_and_set_overlap_surf(bool& bMoveOnce);
  MBErrorCode surface_sense(MBEntityHandle volume, 
			    MBEntityHandle surface,
			    int& sense_out);
  bool is_same_direct_to_ray(int i, int dir);

  int m_nSurf;
  GeomTopoTool* m_GeomTopoTool;
  MBEntityHandle m_hTreeRoot;
  std::vector<MBEntityHandle> m_vhSurfSet;
  MBOrientedBoxTreeTool* m_hObbTree;
  std::vector<MBEntityHandle> triList, surfList;
  std::vector<double> distList;
  std::vector<MBEntityHandle> m_vhInterSurf;
  std::vector<MBEntityHandle> m_vhInterFacet;
  std::map<MBEntityHandle, int> m_mhOverlappedSurf;
#endif
};

#endif
