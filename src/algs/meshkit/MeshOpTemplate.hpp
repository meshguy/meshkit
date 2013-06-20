#ifndef MESHKIT_MESH_OP_TEMPLATE_HPP
#define MESHKIT_MESH_OP_TEMPLATE_HPP

#include <cassert>
#include <string>
#include <vector>
#include <set>

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"

#include "CESets.hpp"
#include "meshkit/LocalTag.hpp"
#include "meshkit/Matrix.hpp"

#include "meshkit/iMesh.hpp"
#include "meshkit/iGeom.hpp"

namespace MeshKit {

class MKCore;

class MeshOpTemplate : public MeshScheme
{
public:
  /* \brief Constructor
   *
   * Create a new MeshOpTemplate instance
   * \param impl the iGeom instance handle for the Geom
   */
  MeshOpTemplate(MKCore *mkcore, const MEntVector &me_vec);

  /* \brief Destructor
   */
  virtual ~MeshOpTemplate();

  /**\brief Get class name */
  static const char* name();

  /**\brief Function returning whether this scheme can mesh entities of t
   *        the specified dimension.
   *\param dim entity dimension
   */
  static bool can_mesh(iBase_EntityType dim);

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
  virtual const moab::EntityType* mesh_types_arr() const;

  /** \brief Re-implemented here so we can check topological dimension of model_ent
   * \param model_ent ModelEnt being added
   */
  virtual bool add_modelent(ModelEnt *model_ent);

  //! Setup is a no-op, but must be provided since it's pure virtual
  virtual void setup_this();

  //! The only setup/execute function we need, since meshing vertices is trivial
  virtual void execute_this();

  /* \brief set size
   *
   * \param vector with x y z size
   */
  void set_size(const Vector<3> &dx);

  void tag_copied_sets(const char **tag_names, const char **tag_vals,
                       const int num_tags);
private:
  iGeom *igeomImpl;
  Vector<3> s_x;
};

inline const char* MeshOpTemplate::name()
{
  return "MeshOpTemplate";
}

inline bool MeshOpTemplate::can_mesh(iBase_EntityType)
{
  // Given just a dimension, MeshOpTemplate can't do anything since it doesn't know
  // what to copy.
  return false;
}

inline bool MeshOpTemplate::can_mesh(ModelEnt *)
{
  return true;
}

inline const moab::EntityType* MeshOpTemplate::mesh_types_arr() const
{
  return output_types();
}

} // namespace MeshKit
#endif
