#ifndef MESHKIT_NG_TRI_MESHER_HPP
#define MESHKIT_NG_TRI_MESHER_HPP

#include "meshkit/iGeom.hpp"
#include <set>
#include <vector>
#include "meshkit/MeshScheme.hpp"
#include "moab/Interface.hpp"

/** \file NGTriMesher.hpp
 */

#include "meshkit/MeshScheme.hpp"

namespace MeshKit
{

class MKCore;

/** \class NGTriMesher NGTriMesher.hpp "meshkit/NGTriMesher.hpp"
 * \brief The interface to the triangle meshing algorithm that is part of the netgen tet mesher
 *
 * This class implements the interface to tri mesher algorithm that is included in netgen.  The
 * algorithm is closesly tied to the .sat file format, so it cannot be used for arbitrary
 * geometry.
 */
class NGTriMesher : public MeshScheme
{
public:
    /** \brief Constructor
     * \param mk_core MKCore instance
     * \param me_vec ModelEnts this mesher will be applied to
     */
  NGTriMesher(MKCore *mk_core, const MEntVector &me_vec);

    /** \brief Destructor
     */
  ~NGTriMesher();

    /** \brief Setup function for this mesher, simply calls setup_boundary
     */
  virtual void setup_this();

    /** \brief Execute the mesher
     */
  virtual void execute_this();

    /** \brief Static list of mesh types treated by this scheme
     */
  static moab::EntityType meshTps[];


  /**\brief Get class name */
  static const char* name()
    { return "NGTriMesher"; }

  /**\brief Function returning whether this scheme can mesh entities of
   *        the specified dimension.
   *\param dim entity dimension
   */
  static bool can_mesh(iBase_EntityType dim)
    { return iBase_FACE == dim; }

  /** \brief Function returning whether this scheme can mesh the specified entity
   *
   * Used by MeshOpFactory to find scheme for an entity.
   * \param me ModelEnt being queried
   * \return If true, this scheme can mesh the specified ModelEnt
   */
  static bool can_mesh(ModelEnt *me)
    { return canmesh_face(me); }

  /**\brief Get list of mesh entity types that can be generated.
   *\return array terminated with \c moab::MBMAXTYPE
   */
  static const moab::EntityType* output_types()
    { return meshTps; }

  /** \brief Return the mesh entity types operated on by this scheme
   * \return array terminated with \c moab::MBMAXTYPE
   */
  virtual const moab::EntityType* mesh_types_arr() const
    { return output_types(); }

};

} // namespace MeshKit

#endif // MESHKIT_NG_TRI_MESHER_HPP
