#ifndef __CAMALTRIADV_HPP
#define __CAMALTRIADV_HPP

#include "meshkit/iGeom.hpp"
#include <set>
#include <vector>
#include "meshkit/MeshScheme.hpp"
#include "moab/Interface.hpp"

/** \file CAMALTriAdvance.hpp
 */

#include "meshkit/MeshScheme.hpp"

class CMLTriAdv;

namespace MeshKit
{

class MKCore;
    
/** \class CAMALSurfEval CAMALSurfEval.hpp "meshkit/CAMALSurfEval.hpp"
 * \brief The MeshKit-based surface evaluator for CAMAL meshing algorithms
 *
 * This class implements surface evaluation functions required by CAMAL algorithms in terms of functions
 * on ModelEnt objects (and, in some cases, functions in iGeom).
 */
class CAMALTriAdvance : public MeshScheme
{
public:
    /** \brief Constructor
     * \param mk_core MKCore instance
     * \param me_vec ModelEnts this mesher will be applied to
     */
  CAMALTriAdvance(MKCore *mk_core, const MEntVector &me_vec);

    /** \brief Destructor
     */
  ~CAMALTriAdvance();

    /** \brief Return the moab entity types produced by this mesher
     * \param tps Moab types returned
     */
  void mesh_types(std::vector<moab::EntityType> &tps);

    /** \brief Setup function for this mesher, simply calls setup_boundary
     */
  virtual void setup_this();

    /** \brief Execute the mesher (calling CAMAL mesher on this surface)
     */
  virtual void execute_this();

    /** \brief Construct a mesher of this type
     * \param mkcore MKCore associated with the mesher
     * \param me_vec ModelEnts to which this mesher will be applied
     */
  static MeshOp *factory(MKCore *mkcore, const MEntVector &me_vec);
	
    /** \brief Static variable for registering this meshop
     */
  static bool meshopRegistered;
  
    /** \brief Static list of geometry types treated by this scheme
     */
  static iBase_EntityType geomTps[];
  
    /** \brief Static list of mesh types treated by this scheme
     */
  static moab::EntityType meshTps[];
  
private:

    /** \brief CAMAL mesher object called by this mesher
     */
  CMLTriAdv *cmlTriAdv;
};

} // namespace MeshKit

#endif
