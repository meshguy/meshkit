// It exchanges entities and sets
// calls moab::ParallelComm::exchange_owned_meshs
// 

#ifndef PAR_EXCHANGE_MESH_HPP
#define PAR_EXCHANGE_MESH_HPP

#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/ParallelComm.hpp"

namespace MeshKit
{
using namespace moab;

class ParExchangeMesh : public MeshScheme
{
public:

  ParExchangeMesh(MKCore *mkcore, const MEntVector &me_vec);

  virtual ~ParExchangeMesh();

  /**\brief Get class name */
  static const char* name() 
    { return "ParExchangeMesh"; }

  /**\brief Function returning whether this scheme can mesh entities of t
   *        the specified dimension.
   *\param dim entity dimension
   */
  static bool can_mesh(iBase_EntityType dim)
    { return iBase_ALL_TYPES == dim; }
   
  /** \brief Function returning whether this scheme can mesh the specified entity
   * 
   * Used by MeshOpFactory to find scheme for an entity.
   * \param me ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
  static bool can_mesh(ModelEnt *me)
    { return (me->dimension() == 4); }
  
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

  // debug function
  void print_mesh();

private:

  ParallelComm* m_mpcomm; // mesh parallel communication

  unsigned int m_rank, m_proc_size;

  iBase_TagHandle m_mPuniqueIDTag;

  std::vector< Range* > m_shared_entities;

  std::vector< unsigned int > m_shared_procs;

  std::vector< MPI_Request > m_recv_reqs, m_recv_remoteh_reqs;

  int get_shared_list(const int to_proc);
};
}

#endif
