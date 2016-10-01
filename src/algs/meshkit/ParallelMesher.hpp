// It does parallel meshing
// Read and distribute geometry in parallel by iGeom
// Inserts shared vertex/edge/surface entity mesh operations and
// exchange operations to graph
// shared entities are meshed and exchanged, then volume mesh is
// performed at last
// 

#ifndef PARALLELMESHER_HPP
#define PARALLELMESHER_HPP

#include "meshkit/MeshScheme.hpp"
#include "moab/ParallelComm.hpp"
#include "TDParallel.hpp"

namespace MeshKit
{
class ModelEnt;

using namespace moab;

enum PARALLEL_OP_TYPE
{
  MESH_VERTEX = 0,
  MESH_EDGE,
  MESH_INTER_SURF,
  MESH_NINTER_SURF,
  EXCHANGE_VERTEX,
  EXCHANGE_EDGE,
  SEND_POST_SURF_MESH,
  RECV_SURF_MESH,
  MESH_VOLUME
};

class ParallelMesher : public MeshScheme
{
public:

  ParallelMesher(MKCore *mkcore, const MEntVector &me_vec);

  virtual ~ParallelMesher();

  /**\brief Get class name */
  static const char* name() 
    { return "ParallelMesher"; }

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
  static bool can_mesh(ModelEnt *me)
    { return canmesh_region(me); }
  
  /**\brief Get list of mesh entity types that can be generated.
   *\return array terminated with \c moab::MBMAXTYPE
   */
  static const moab::EntityType* output_types();

  /** \brief Return the mesh entity types operated on by this scheme
   * \return array terminated with \c moab::MBMAXTYPE
   */
  virtual const moab::EntityType* mesh_types_arr() const
    { return output_types(); }

  /**\brief Setup is a no-op, but must be provided since it's pure virtual
   */
  virtual void setup_this();

  /**\ The only setup/execute function we need, since meshing vertices is trivial
   */
  virtual void execute_this();

  /** \brief set mesh size
   * \param size mesh size
   * \param interval # of interval
   */
  void set_mesh_size(double size, int interval);

  // debug function
  void print_mesh();

  // debug function
  void print_geom_info(ModelEnt* me, const int dim,
                       const bool local);
  
private:

  /** \brief Construct a MeshOp that can generate input type elements
   * \param type type of parallel operation
   * \return A MeshOp that can generate input type elements
   */
  MeshOp* get_mesher(PARALLEL_OP_TYPE type);

  /** \brief insert mesh operation for parallel meshing
   * \param type type of parallel operation
   * \param after if insert after this mesh operation
   * \return A MeshOp that can generate input type elements
   */
  void add_parallel_mesh_op(PARALLEL_OP_TYPE type, bool after = false);

  void check_partition(TDParallel* td_par, ModelEnt* me, int dim);

  ParallelComm* m_mpcomm; // mesh parallel communication

  unsigned int m_rank;

  iBase_TagHandle m_mPuniqueIDTag;
  
  std::vector< MEntSet > m_sEntity;
};
}

#endif
