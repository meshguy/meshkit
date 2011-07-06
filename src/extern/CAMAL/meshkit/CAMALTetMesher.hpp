#ifndef MESHKIT_CAMAL_TET_MESHER_HPP
#define MESHKIT_CAMAL_TET_MESHER_HPP

#include "meshkit/iGeom.hpp"
#include <set>
#include <vector>
#include "meshkit/MeshScheme.hpp"
#include "moab/Interface.hpp"

/** \file CAMALTetMesher.hpp
 */

#include "meshkit/MeshScheme.hpp"

class CMLTetMesher;

namespace MeshKit
{

class MKCore;
    
/** \class CAMALTetMesher CAMALTetMesher.hpp "meshkit/CAMALTetMesher.hpp"
 * \brief The wrapper for the CAMAL tet mesher
 *
 * This class calls the CAMAL tet mesher and converts mesh to/from the required format for that mesher.
 */
class CAMALTetMesher : public MeshScheme
{
public:
    /** \brief Constructor
     * \param mk_core MKCore instance
     * \param me_vec ModelEnts this mesher will be applied to
     */
  CAMALTetMesher(MKCore *mk_core, const MEntVector &me_vec);

    /** \brief Destructor
     */
  ~CAMALTetMesher();

    /** \brief Setup function for this mesher, simply calls setup_boundary
     */
  virtual void setup_this();

    /** \brief Execute the mesher (calling CAMAL mesher on this volume)
     */
  virtual void execute_this();
	
    /** \brief Static list of mesh types treated by this scheme
     */
  static moab::EntityType meshTps[];
  
  
  /**\brief Get class name */
  static const char* name() 
    { return "CAMALTetMesher"; }

  /**\brief Function returning whether this scheme can mesh entities of t
   *        the specified dimension.
   *\param dim entity dimension
   */
  static bool can_mesh(iBase_EntityType dim)
    { return iBase_REGION == dim; }

  /** \brief Function returnign whether this scheme can mesh the specified entity
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
  
private:

    /** \brief Construct a MeshOp that can generate triangle elements
     * \return A MeshOp that can generate tri elements
     */
  MeshOp *get_tri_mesher();

  /** \brief print debug information and save input boundary mesh for CAMAL library
   */
  void print_debug(ModelEnt *me, std::vector<double> &coords,
                   moab::Range &bdy_vrange,
                   std::vector<moab::EntityHandle> &bdy,
                   std::vector<int> &group_sizes,
                   std::vector<int> &bdy_ids);
  
    /** \brief CAMAL mesher object called by this mesher
     */
  CMLTetMesher *cmlTetMesher;
};

} // namespace MeshKit

#endif
