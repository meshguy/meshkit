#ifndef MESHKIT_TETGEN_TET_MESHER_HPP
#define MESHKIT_TETGEN_TET_MESHER_HPP

#include "meshkit/iGeom.hpp"
#include <set>
#include <vector>
#include "meshkit/MeshScheme.hpp"
#include "moab/Interface.hpp"

/** \file TGTetMesher.hpp
 */

#include "meshkit/MeshScheme.hpp"

class CMLTetMesher;

namespace MeshKit
{

class MKCore;
    
/** \class TGTetMesher TGTetMesher.hpp "meshkit/TGTetMesher.hpp"
 * \brief The wrapper for the Tetgen tet mesher
 *
 * This class calls the Tetgen tet mesher and converts mesh to/from the required format for that mesher.
 */
class TGTetMesher : public MeshScheme
{
public:
    /** \brief Constructor
     * \param mk_core MKCore instance
     * \param me_vec ModelEnts this mesher will be applied to
     */
  TGTetMesher(MKCore *mk_core, const MEntVector &me_vec);

    /** \brief Destructor
     */
  ~TGTetMesher();

    /** \brief Setup function for this mesher, simply calls setup_boundary
     */
  virtual void setup_this();

    /** \brief Execute the mesher (calling Tetgen mesher on this volume)
     */
  virtual void execute_this();
	
    /** \brief Static variable for registering this meshop
     */
  static bool meshopRegistered;
  
    /** \brief Static list of geometry types treated by this scheme
     */
  static iBase_EntityType geomTps[];
  
    /** \brief Static list of mesh types treated by this scheme
     */
  static moab::EntityType meshTps[];
  
  
  /**\brief Get class name */
  static const char* name() 
    { return "TGTetMesher"; }

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
  static const moab::EntityType* output_types()
    { return meshTps; }

  /** \brief Return the mesh entity types operated on by this scheme
   * \return array terminated with \c moab::MBMAXTYPE
   */
  virtual const moab::EntityType* mesh_types_arr() const
    { return output_types(); }

  /** \brief This variable controls tegen command line switches
   * \
   */
  char* tet_switch;
  
private:

    /** \brief Construct a MeshOp that can generate triangle elements
     * \return A MeshOp that can generate tri elements
     */
  MeshOp *get_tri_mesher();
};

} // namespace MeshKit

#endif
