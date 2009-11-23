// Author : Hongjun Kim
// Nov 16, 09
// It reads geometry file, makes facet triangles
// constructs tree structure for triangles and makes hexes bounding geometry
// With ray-tracing, find intersections and determine element inside/outside/boundary.
// Intersection fraction is stored to boundary elements.
// Element inside/outside/boundary status are stored as tag.

#ifndef CUTCELLMESH_HPP
#define CUTCELLMESH_HPP

#include <vector>
#include <map>
#include <sys/resource.h>

#include "iMesh.h"

#ifdef MOAB
#include "MBInterface.hpp"
#include "MBOrientedBoxTreeTool.hpp"
#include "MBCartVect.hpp"
#endif

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
  double fraction[3];

  CutFraction() {
    for (int i = 0; i < 3; i ++) {
      fraction[i] = 0.;
    }
  }

  CutFraction(int dir, double val) {
    for (int i = 0; i < 3; i ++) {
      if (i == dir) fraction[i] = val;
      else fraction[i] = 0.;
    }
  };
};

class CutCellMesh
{
public:

  //explicit CutCellMesh(iMesh_Instance mesh,
  CutCellMesh(iMesh_Instance mesh,
	      iBase_EntitySetHandle root_set,
	      double size = -1.);
  
  //virtual ~CutCellMesh();
  ~CutCellMesh();

  iMesh_Instance mesh_impl() const
  {
    return m_mesh;
  }

  int do_mesh(int exp, int s_exp, int file_type);
  int set_division();
  int make_hex_vertices();
  int make_hexes();
  bool move_ray(MBOrientedBoxTreeTool& tool,
		double* startPnt, double* endPnt,
		double tol, int& nIntersect, int dir);
  
#ifdef MOAB
  MBInterface* moab_instance() {return reinterpret_cast<MBInterface*> (m_mesh);}
  void set_division(const MBOrientedBox& box);
  int find_intersections(MBOrientedBoxTreeTool& tool);
  int fire_rays(MBOrientedBoxTreeTool& tool, int dir);
#endif

private:

  iBase_TagHandle m_elem_status_tag_handle;
  iBase_TagHandle m_cut_fraction_tag_handle;
  iMesh_Instance m_mesh;
  iBase_EntitySetHandle m_hRootSet;
  std::vector<double> m_vdIntersection;

  int m_nNode[3], m_nDiv[3], m_iElem;
  double m_dIntervalSize[3], m_origin_coords[3], m_dInputSize;
  int m_nTri, m_nHex;
  double m_dFirstP, m_dSecondP;
  int m_iInter;
  EdgeStatus m_nStatus;

  std::vector<iBase_EntityHandle> m_vhVertex;
  std::vector<iBase_EntityHandle> m_vhHex;
  std::vector<char> m_vnHexStatus;
  std::map<iBase_EntityHandle, CutFraction> m_mdCutFraction;

  EdgeStatus getEdgeStatus(const double dZ, bool& bMoveNext);
  bool set_neighbor_hex_status(int dir, int i, int j, int k);
  bool set_hex_status(int index, EdgeStatus value);
  int set_tag_info();

  int print_debug();
  int export_mesh(int s_exp, int file_type);
  int write_mesh(int type, int file_type, iBase_EntityHandle* handles,
		 int& n_elem);
  void util_getrusage(struct rusage &r_usage);

#ifdef MOAB
  MBEntityHandle m_hTreeRoot;
  std::vector<MBEntityHandle> triList, surfList;
  std::vector<double> distList; 
#endif
};

#endif
