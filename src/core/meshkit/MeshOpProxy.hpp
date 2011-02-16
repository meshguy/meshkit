#ifndef MESHKIT_MESH_OP_PROXY
#define MESHKIT_MESH_OP_PROXY

#include "meshkit/Types.hpp"
#include "iBase.h"

namespace MeshKit 
{

class MeshOp;
class MKCore;

/**\brief Polymorphic behavior for MeshOp classes (as opposed to instances)
 *
 * Each registered MeshOp sub-class must register via an implementation
 * of this interface.  An implementation of this class provides the necessary
 * information about the class as a whole.
 */
class MeshOpProxy
{
  public:
    /**\brief Factory method */
    virtual MeshOp* create( MKCore* core, const MEntVector &vec ) = 0;
    /**\brief String name of MeshOp class */
    virtual const char* name() const = 0;
    /**\brief Get array of mesh entity types produced by MeshOp instances
     *
     * \return an array of \c MOAB::EntityType terminated by the
     * \c MOAB::MBMAXTYPE value.
     */
    virtual const moab::EntityType* output_types() const = 0;
    /**\brief Check if class can be used to mesh entites of the specified dimension */
    virtual bool can_mesh( iBase_EntityType dimension ) const = 0;
    /**\brief Check if class can be used to mesh passed entity */
    virtual bool can_mesh( ModelEnt* entity ) const = 0;
};

} // namespace MeshKit

#endif // #ifdef MESHKIT_MESH_OP_PROXY
