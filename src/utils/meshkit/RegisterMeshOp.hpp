
/**\file RegisterMeshOp.hpp
 * Declare \c RegisterMeshOp and requisite utility templates
 */

#ifndef MESHKIT_REGISTER_MESH_OP_HPP
#define MESHKIT_REGISTER_MESH_OP_HPP

#include "meshkit/MeshOpProxy.hpp"
#include "meshkit/MKCore.hpp"

namespace MeshKit {

/**\brief Utility class for registering MeshOps with MKCore
 *
 * This template class provides two related functionalities:
 * - It provides an implementation of MKCore::MeshOpProxy for the class
 * - It handles registration with MKCore
 *
 * To use this class, simply declare a global variable templated
 * with the type of the MeshOp to be registered.
 *
 *\param TYPE The MeshOp sub-class to register
 */
template <class TYPE> 
class RegisterMeshOp : public MeshOpProxy
{
  public:

    /** \brief Register a new MeshOp factory
     */
    RegisterMeshOp()
      { MKCore::register_meshop(this); }

    /**\brief Implementation of factory method */
    MeshOp* create( MKCore* core, const MEntVector& vec )
      { TYPE *t = new TYPE(core, vec); t->set_name(name()); return t;}
    
    /**\brief Call TYPE::name() to get name for class */
    const char* name() const
      { return TYPE::name(); }
      
    /**\brief Call TYPE::output_types 
     *
     * Returns an array of \c MOAB::EntityType values, terminated
     * by a value of \c MOAB::MBMAXTYPE .
     */
    const moab::EntityType* output_types() const
      { return TYPE::output_types(); }
      
    /**\brief Call TYPE::can_mesh(iBase_EntityType)
     *
     * Ask if the class is capable of meshing an entity of the
     * specified dimension.
     */
    bool can_mesh( iBase_EntityType dimension ) const
      { return TYPE::can_mesh(dimension); }
    
    /**\brief Implementation of can_mesh() method for proxy of \c TYPE
     *
     *\return \c TYPE::can_mesh() if \c HAVE_CAN_MESH is true.
     *        \c default MeshOpProxy::can_mesh() value if \c HAVE_CAN_MESH is false.
     */
    bool can_mesh( ModelEnt* entity ) const
      { return TYPE::can_mesh(entity); }
};

} // namespace MeshKit

#endif 

