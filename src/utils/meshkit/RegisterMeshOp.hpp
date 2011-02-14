
/**\file RegisterMeshOp.hpp
 * Declare \c RegisterMeshOp and requisite utility templates
 */

#ifndef MESHKIT_REGISTER_MESH_OP_HPP
#define MESHKIT_REGISTER_MESH_OP_HPP

#include "meshkit/MKCore.hpp"

namespace MeshKit {

/**\brief Utility template for RegisterMeshOp */
template <class TYPE, bool HAVE_CAN_MESH>
class MeshOpCanMesh 
{};

/**\brief Utility template for RegisterMeshOp */
template <class TYPE>
class MeshOpCanMesh<TYPE,true> 
{
  public:
  static bool can_mesh( const MKCore::MeshOpProxy& p, ModelEnt* entity )
    { return TYPE::can_mesh(entity); }
};

/**\brief Utility template for RegisterMeshOp */
template <class TYPE>
class MeshOpCanMesh<TYPE,false> 
{
  public:
  static bool can_mesh( const MKCore::MeshOpProxy& p, ModelEnt* entity )
    { return p.MKCore::MeshOpProxy::can_mesh(entity); }
};

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
 *\param HAVE_CAN_MESH Specify as true if class contains a static can_mesh()
 *          method that should be used to implement MKCore::MeshOpProxy::can_mesh()
 */
template <class TYPE, bool HAVE_CAN_MESH> 
class RegisterMeshOp : public MKCore::MeshOpProxy
{
  public:
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tp The iBase_EntityType (or dimension) operated on by this MeshOp type
     * \param tp The MOAB entity type produced by this MeshOp
     */
    RegisterMeshOp(const char* op_name,
                   iBase_EntityType model_tp, 
                   moab::EntityType tp)
      { MKCore::register_meshop( op_name, &model_tp, 1, &tp, 1, this ); }
    
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tp The iBase_EntityType (or dimension) operated on by this MeshOp type
     * \param tp The iMesh entity topology produced by this MeshOp
     */
    RegisterMeshOp(const char *op_name, 
                   iBase_EntityType model_tp, 
                   iMesh::EntityTopology tp)
      { MKCore::register_meshop( op_name, &model_tp, 1, &tp, 1, this ); }
  
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tps The iBase_EntityType's (or dimensions) operated on by this MeshOp type
     * \param num_mtps Number of model entity types in model_tps
     * \param tps The MOAB entity types operated on by this MeshOp
     * \param num_tps Number of entity types in tps
     */
    RegisterMeshOp(const char *op_name, 
                   iBase_EntityType *model_tps, int num_mtps,
                   moab::EntityType *tps, int num_tps)
      { MKCore::register_meshop( op_name, model_tps, num_mtps, tps, num_tps, this ); }
  
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tps The iBase_EntityType's (or dimensions) operated on by this MeshOp type
     * \param num_mtps Number of model entity types in model_tps
     * \param tps The iMesh entity types operated on by this MeshOp
     * \param num_tps Number of entity types in tps
     */
    RegisterMeshOp(const char *op_name, 
                   iBase_EntityType *model_tps, int num_mtps,
                   iMesh::EntityTopology *tps, int num_tps)
      { MKCore::register_meshop( op_name, model_tps, num_mtps, tps, num_tps, this ); }

    /**\brief Implementation of factory method */
    MeshOp* create( MKCore* core, const MEntVector& vec )
      { return new TYPE(core, vec); }
    
    /**\brief Implementation of can_mesh() method for proxy of \c TYPE
     *
     *\return \c TYPE::can_mesh() if \c HAVE_CAN_MESH is true.
     *        \c default MeshOpProxy::can_mesh() value if \c HAVE_CAN_MESH is false.
     */
    bool can_mesh( ModelEnt* entity )
      { return MeshOpCanMesh<TYPE,HAVE_CAN_MESH>::can_mesh(*this,entity); }
};

} // namespace MeshKit

#endif 

