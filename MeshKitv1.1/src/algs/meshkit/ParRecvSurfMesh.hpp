#ifndef PAR_RECV_SURF_MESH_HPP
#define PAR_RECV_SURF_MESH_HPP

#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "moab/ParallelComm.hpp"

namespace MeshKit
{
using namespace moab;

class ParRecvSurfMesh : public MeshScheme
{
public:

  ParRecvSurfMesh(MKCore *mkcore, const MEntVector &me_vec);

  virtual ~ParRecvSurfMesh();

  /**\brief Get class name */
  static const char* name() 
    { return "ParRecvSurfMesh"; }

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

private:

  ParallelComm* m_mpcomm; // mesh parallel communication

  unsigned int m_rank;

  std::vector< Range* > m_shared_entities;

  std::vector< unsigned int > m_shared_procs;

  int get_shared_list(const int to_proc);
};
}

#endif
