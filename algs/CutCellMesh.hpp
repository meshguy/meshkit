#ifndef CUTCELLMESH_HPP
#define CUTCELLMESH_HPP

#include "iMesh.h"

#ifdef DAGMC
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

class CutCellMesh
{
public:

  explicit CutCellMesh(iMesh_Instance mesh, iBase_EntitySetHandle root_set);
  virtual ~CutCellMesh();

  iMesh_Instance mesh_impl() const
  {
    return m_mesh;
  }

  int do_mesh();
  int set_initial_division();
  int make_hex_vertices();
  int make_initial_hexes();
  
#ifdef DAGMC
  MBInterface* moab_instance() {return reinterpret_cast<MBInterface*> (m_mesh);}
  void set_initial_division(const MBOrientedBox& box);
  int find_intersected_surfaces(MBOrientedBoxTreeTool& tool);
#endif

private:

  iMesh_Instance m_mesh;
  iBase_EntitySetHandle m_hRootSet;
  std::vector<double> m_vdIntersection;

  int m_nNode[3], m_nDiv[3];
  double m_dIntervalSize[3], m_origin_coords[3];
  int m_nTri, m_nHex;
  double m_dFirstZ, m_dSecondZ;
  int m_iInter;
  EdgeStatus m_nStatus;

  std::vector<iBase_EntityHandle> m_vhVertex;
  std::vector<iBase_EntityHandle> m_vhHex;
  std::vector<int> m_vnHexStatus;

  EdgeStatus getEdgeStatus(const double dZ, bool bMoveNext);
  bool set_hex_status(int index, int value);

#ifdef DAGMC
  MBEntityHandle m_hTreeRoot;
  std::vector<MBEntityHandle> triList, surfList;
  std::vector<double> distList; 

  //MBTag get_tag( const char* name, int size, MBTagType store, MBDataType type,
  //             bool create_if_missing = true);
#endif
};

#endif
