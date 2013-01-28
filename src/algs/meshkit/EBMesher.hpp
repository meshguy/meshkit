#ifndef MESHKIT_EB_MESHER_HPP
#define MESHKIT_EB_MESHER_HPP

#include <vector>
#include <map>
#include <sys/resource.h>

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/iMesh.hpp"

#ifdef HAVE_MOAB
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBOrientedBoxTreeTool.hpp"
#endif

//! edge status
enum EdgeStatus {
  INSIDE,
  OUTSIDE,
  BOUNDARY
};

//! stores boundary element edge cut fraction
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

//! boundary cut-cell surface edge key
//! made for TechX
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

//! stores intersection distances
struct IntersectDist {
  double distance;
  int index;

  IntersectDist() {};

  IntersectDist(double d, int i) {
    distance = d;
    index = i;
  };
};

//! stores volume fractions
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

namespace MeshKit {

class MKCore;
    

/** \class EBMesher EBMesher.hpp "meshkit/EBMesher.hpp"
 * \brief A meshing geometry as Cartesian structured mesh
 * It makes constructs tree structure with  triangles and
 * makes hexes bounding geometry
 * With ray-tracing, find intersections and determine element inside/outside/boundary.
 * Intersection fraction is stored to boundary elements.
 * Element inside/outside/boundary status are stored as tag.
 */
class EBMesher : public MeshScheme
{
public:

    //! Bare constructor
  EBMesher(MKCore *mkcore, const MEntVector &me_vec,
           double size = -1., bool use_geom = true,
           int add_layer = 3);
  
    //! Destructor
  virtual ~EBMesher();
  
  /**\brief Get class name */
  static const char* name() 
    { return "EBMesher"; }

  /**\brief Function returning whether this scheme can mesh entities of t
   *        the specified dimension.
   *\param dim entity dimension
   */
  static bool can_mesh(iBase_EntityType dim)
    { return iBase_REGION == dim; }
   
  /** \brief Function returning whether this scheme can mesh the specified entity
   * 
   * Used by MeshOpFactory to find scheme for an entity.
   * \param me ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
  static bool can_mesh(ModelEnt *me);
  
  /**\brief Get list of mesh entity types that can be generated.
   *\return array terminated with \c moab::MBMAXTYPE
   */
  static const moab::EntityType* output_types();
  
  /** \brief Return the mesh entity types operated on by this scheme
   * \return array terminated with \c moab::MBMAXTYPE
   */
  virtual const moab::EntityType* mesh_types_arr() const
  { return output_types(); }
  
  
  virtual bool add_modelent(ModelEnt *model_ent);
  
  /** \brief set # of divisions (x/y/z directions) of Cartesian box to use SCDMesh output
   * \param min Cartesian box min coordinates
   * \param max Cartesian box max coordinates
   * \return int if is working correctly
   */
  int set_division(double* min, double* max);
  
  /**\brief Setup is a no-op, but must be provided since it's pure virtual
   */  
  virtual void setup_this();
  
  /**\ The only setup/execute function we need, since meshing vertices is trivial
   */
  virtual void execute_this();
  
  /** \brief query function for techX
   * \param boxMin Cartesian box min coordinates returned
   * \param boxMax Cartesian box max coordinates returned
   * \param nDiv nDiv(x/yz directions) returned
   * \param rmdCutCellSurfEdge map of cut-cell surface index and edge cut information returned
   * \param rvnInsideCell Inside elements returned
   * \param isCornerExterior if box corner is exterior returned
   */
  void get_grid_and_edges_techX(double* boxMin, double* boxMax, int* nDiv,
                                std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellSurfEdge,
                                std::vector<int>& rvnInsideCell, bool isCornerExterior = true);
  
  /** \brief query function to get multiple cut-cell edges
   * \param boxMin Cartesian box min coordinates returned
   * \param boxMax Cartesian box max coordinates returned
   * \param nDiv nDiv(x/yz directions) returned
   * \param rmdCutCellEdge map of cut-cell surface index and edge cut information returned
   * \param rvnInsideCell Inside elements returned
   * \param isCornerExterior if box corner is exterior returned
   * \return if this function is working correctly
   */
  bool get_grid_and_edges(double* boxMin, double* boxMax, int* nDiv,
                          std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellEdge,
                          std::vector<int>& rvnInsideCell, bool isCornerExterior = true);
  
  /** \brief get volume fraction for each material
   * \param vol_frac_div resolution to get volume fraction
   * \return if this function is working correctly
   */
  bool get_volume_fraction(int vol_frac_div);
  
  /** \brief export mesh to file
   * \param file_name export file name
   * \param separate export file separately(inside/outside/boundary)
   */
  void export_mesh(const char* file_name, bool separate = false);

  /** \brief set number of intervals
   * \param n_interval number of interval array for x/y/z directions
   */
  void set_num_interval(int* n_interval);

  /** \brief set number of intervals
   * \param box_increase number of interval array for x/y/z directions
   */
  void increase_box(double box_increase);

  /** \brief set if mesh for whole geometry at once or individually
   * \param use if mesh for whole geometry at once
   */
  void use_whole_geom(int use);

  /** \brief set if mesh based geometry is used
   * \param use if mesh based geometry is used
   */
  void use_mesh_geometry(bool use);

  /** \brief construct obb tree with faceted geometry and set the box max/min coordinates
   */
  void set_obb_tree_box_dimension();
  
protected:
  
private:

  //! No copy constructor, since there's only meant to be one of these
  EBMesher(const EBMesher &);
  
  MeshOp* get_scd_mesher();

  //! No operator=, since there's only meant to be one of these
  EBMesher &operator=(const EBMesher &);

  iBase_TagHandle m_elemStatusTag, m_edgeCutFracLengthTag,
    m_edgeCutFracTag, m_volFracLengthTag, m_volFracHandleTag, m_volFracTag, m_bbTag;
  iMesh_Instance m_mesh;
  iBase_EntitySetHandle m_hRootSet;
  std::vector<IntersectDist> m_vIntersection;
  int m_nTri, m_nHex, m_iInter, m_nFraction, m_iStartHex, m_nMove, m_nAddLayer,
    m_nIteration, m_iOverlap, m_iElem, m_nNode[3], m_nDiv[3],
    m_iFirstP, m_iSecondP;
  double m_dFirstP, m_dSecondP, m_curPnt, m_prevPnt, m_boxIncrease,
    m_dIntervalSize[3], m_origin_coords[3], m_dInputSize,
    m_min[3], m_max[3];
  EdgeStatus m_nStatus;
  bool m_bMovedOnce, m_bUseGeom, m_bUseWholeGeom, m_bUseMeshGeom;
  std::vector<iBase_EntityHandle> m_vhVertex;
  std::vector<int> m_vnHexStatus;
  std::map<int, CutFraction> m_mdCutFraction;
  std::vector<EdgeStatus> m_vnEdgeStatus[3];
  
  /** \brief get hex edge status (inside/outside/boundary)
   * \param dZ edge end coordinate
   * \param iSkip how many index skipped for next intersection checking
   * \return EdgeStatus edge status
   */
  EdgeStatus get_edge_status(const double dZ, int& iSkip);

  /** \brief set hex status for neighboring elements
   * \param dir current fired ray direction
   * \param i index for i direction
   * \param j index for j direction
   * \param k index for k direction
   * \return bool if is working correctly
   */
  bool set_neighbor_hex_status(int dir, int i, int j, int k);

  /** \ brief set hex status by edge status
   * \param index hex index in m_vnHexStatus vector
   * \param value edge status
   * \param dir ray direction
   * \return bool if is working correctly
   */
  bool set_hex_status(int index, EdgeStatus value, int dir);

  /** \brief set edge status
   * \param dir ray direction
   * \return bool if is working correctly
   */
  bool set_edge_status(int dir);

  /** \brief set all produced mesh information as tag
   * \ hex status, edge-cut information.....
   */
  void set_tag_info();

  /** \brief wirte mesh
   * \param file_name
   * \param type element type (inside:0, outside:1, boundary:2)
   * \param handles element handles
   * \param n_elem # of elements
   * \return int if is working correctly
   */
  void write_mesh(const char* file_name, int type,
                  iBase_EntityHandle* handles, int& n_elem);

  /** \brief get edge fraction information
   * \param idHex index in m_mdCutFraction
   * \param dir ray direction
   * \return double edge fraction
   */
  double get_edge_fraction(int idHex, int dir);

  /** \brief get if the edge is fully inside(returns 1) or outside(returns 0)
   * \param i index for i direction
   * \param j index for j direction
   * \param k index for k direction
   * \param dir ray direction
   * \return double edge fraction
   */
  double get_uncut_edge_fraction(int i, int j, int k, int dir);

  /** \brief check if the intersected surface geometry is shared or overlapped
   * \param index intersection index in m_vIntersection vector
   * \return bool if is working correctly
   */
  bool is_shared_overlapped_surf(int index);

  /** \brief move intersection pairs to check element status
   * \param n_dir ray direction
   * \param n_inter index of intersection points
   * \param start_pnt ray starting point
   * \return bool if is working correctly
   */
  bool move_intersections(int n_dir, int n_inter, double start_pnt[3]);

  /** \brief get inside status elements
   * \param rvnInsideCell cell indices (i,j,k triple)
   * \return bool if is working correctly
   */
  bool get_inside_hex(std::vector<int>& rvnInsideCell);

  /** \brief check if ray is passing shared vertices or edges
   * \ And, check if the passing surface is overlapped one
   * \param bMoveOnce if ray is already moved before
   * \return bool if is working correctly
   */
  bool is_ray_move_and_set_overlap_surf(bool& bMoveOnce);

  void remove_intersection_duplicates();

  // test function 1 for debugging
  bool export_fraction_edges(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& rmdCutCellSurfEdge);

  // test functions 2 for debugging
  bool export_fraction_points(std::map< CutCellSurfEdgeKey, std::vector<double>, LessThan >& mdCutCellEdge);

  // test functions 3 for debugging
  bool make_edge(double ePnt[6], std::vector<iBase_EntityHandle>& edge_handles);

    //! Static variable, used in registration
  static int init;

#ifdef HAVE_MOAB
  /** \brief get MOAB instance
   * \return MOAB instance
   */
  MBInterface* moab_instance() {return mk_core()->moab_instance();}

  /** \brief get MOAB's Tag
   * \param name Tag name
   * \param size Tag size
   * \param store Tag type
   * \param type data type stored
   * \param def_value default value
   * \param create_if_missing create Tag if it is missed (flag)
   * \return Tag handle
   */
  iBase_TagHandle get_tag(const char* name, int size, unsigned flags, MBDataType type,
                          const void* def_value = NULL, bool create_if_missing = true);
  
  /** \brief get MOAB's various length Tag
   * \param name Tag name
   * \param store Tag type
   * \param type data type stored
   * \return Tag handle
   */
  iBase_TagHandle get_various_length_tag(const char* name,
                                         unsigned store, MBDataType type);

  /** \brief construct hexes in MOAB structured mesh format
   * \return int if is working correctly
   */
  int make_scd_hexes();

  /** \brief construct hexes in MOAB unstructured mesh format
   * \return int if is working correctly
   */
  int make_uscd_hexes();

  /** \brief construct OBB tree
   * \input geometry is faceted which are constructed to OBB tree
   * \return int if is working correctly
   */
  int construct_obb_tree();

  /** \brief construct and set OBB tree
   * \input Model entity pointer for tree construction
   */
  void set_tree_root(ModelEnt* me);

  /** \brief set # of divisions
   * \calcalate the best ones from the size of faceted triangles
   */
  void set_division();

  /** \brief find intersections by firing rays
   */
  void find_intersections();

  /** \brief fires rays to 3 directions
   * \param dir give the ray direction
   * \it calls fire_ray
   */
  void fire_rays(int dir);

  /** \brief fires ray
   * \param nIntersect # of intersections returned
   * \param startPnt ray starting point
   * \param endPnt ray ending point
   * \param tol tolerance to find intersection
   * \param dir ray direction
   * \param rayLength ray length
   * \return bool if is working correctly
   */
  bool fire_ray(int& nIntersect, double startPnt[3],
                double endPnt[3], double tol, int dir,
                double rayLength);

  /** \brief moves ray which passes any singluar point
   * \param nIntersect # of intersections returned
   * \param startPnt ray starting point
   * \param endPnt ray ending point
   * \param tol tolerance to find intersection
   * \param dir ray direction
   * \param bMoveOnce if ray is already moved before
   * \return bool if is working correctly
   */
  bool move_ray(int& nIntersect, double* startPnt, double* endPnt,
                double tol, int dir, bool bMoveOnce);

  /** \brief if the facet has the same direction to the ray
   * \param i index in m_vIntersection (intersection vector)
   * \param dir ray direction
   * \return bool if is working correctly
   */
  bool is_same_direct_to_ray(int i, int dir);

  // ! GeomTopoTool instance
  moab::GeomTopoTool* m_GeomTopoTool;

  // ! Tree root
  MBEntityHandle m_hTreeRoot;

  // ! OBB tree tool instance
  MBOrientedBoxTreeTool* m_hObbTree;

  // ! intersected surface geometry list
  std::vector<MBEntityHandle> m_vhInterSurf;

  // ! intersected facet list
  std::vector<MBEntityHandle> m_vhInterFacet;

  // ! overlapped surface list
  std::map<MBEntityHandle, int> m_mhOverlappedSurf;

  std::map<MBEntityHandle, MBEntityHandle>  m_mRootSets;

  double m_minCoord[3], m_maxCoord[3];
#endif
};

inline void EBMesher::increase_box(double box_increase)
{
  m_boxIncrease = box_increase;
}

inline void EBMesher::use_whole_geom(int use)
{
  m_bUseWholeGeom = use;
}

inline void EBMesher::use_mesh_geometry(bool use)
{
  m_bUseMeshGeom = use;
}
}

#endif

  
