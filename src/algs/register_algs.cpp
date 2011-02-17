#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/EBMesher.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/SCDMesh.hpp"
#include "meshkit/VertexMesher.hpp"

namespace MeshKit {

/**\brief Dummy function to force load from static library */
extern int register_algs_mesh_ops() { return 1; }

/**\brief Register a MeshOp type during initialization
 *\param NAME Class name for MeshOp subclass
 */
#define REGISTER_MESH_OP(NAME) \
  RegisterMeshOp<NAME> NAME ## _GLOBAL_PROXY

REGISTER_MESH_OP(CopyMesh);
REGISTER_MESH_OP(EBMesher);
REGISTER_MESH_OP(EdgeMesher);
REGISTER_MESH_OP(ExtrudeMesh);
REGISTER_MESH_OP(OneToOneSwept);
REGISTER_MESH_OP(SCDMesh);
REGISTER_MESH_OP(VertexMesher);


} // namespace MeshKit
