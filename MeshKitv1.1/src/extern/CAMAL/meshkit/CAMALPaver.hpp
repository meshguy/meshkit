#ifndef __CAMALPAVER_HPP
#define __CAMALPAVER_HPP

#include "meshkit/iGeom.hpp"
#include <set>
#include <vector>
#include "meshkit/MeshScheme.hpp"
#include "moab/Interface.hpp"

/** \file CAMALPaver.hpp
 */

#include "meshkit/MeshScheme.hpp"

namespace MeshKit
{

class MKCore;
    
/** \class CAMALPaver CAMALPaver.hpp "meshkit/CAMALPaver.hpp"
 * \brief The MeshKit-based surface evaluator for CAMAL meshing algorithms
 *
 * This class implements the interface to the CAMAL paver.
 */
class CAMALPaver : public MeshScheme
{
public:
    /** \brief Constructor
     * \param mk_core MKCore instance
     * \param me_vec ModelEnts this mesher will be applied to
     */
  CAMALPaver(MKCore *mk_core, const MEntVector &me_vec);

    /** \brief Destructor
     */
  ~CAMALPaver();

    /** \brief Setup function for this mesher, simply calls setup_boundary
     */
  virtual void setup_this();

    /** \brief Execute the mesher (calling CAMAL mesher on this surface)
     */
  virtual void execute_this();
	
    /** \brief Static list of mesh types treated by this scheme
     */
  static moab::EntityType meshTps[];
  
  
  /**\brief Get class name */
  static const char* name() 
    { return "CAMALPaver"; }

  /**\brief Function returning whether this scheme can mesh entities of t
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
  
private:

};

} // namespace MeshKit

#endif
