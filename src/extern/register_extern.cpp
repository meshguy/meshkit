#include "meshkit/RegisterMeshOp.hpp"
#ifdef HAVE_CAMAL
#  include "meshkit/CAMALTriAdvance.hpp"
#  include "meshkit/CAMALTetMesher.hpp"
#endif

namespace MeshKit {

/**\brief Dummy function to force load from static library */
extern int register_extern_mesh_ops() { return 1; }

/**\brief Register a MeshOp type during initialization
 *\param NAME Class name for MeshOp subclass
 */
#define REGISTER_MESH_OP(NAME) \
  RegisterMeshOp<NAME> NAME ## _GLOBAL_PROXY

#ifdef HAVE_CAMAL
  REGISTER_MESH_OP(CAMALTriAdvance);
  REGISTER_MESH_OP(CAMALTetMesher);
#endif


} // namespace MeshKit

