#include "meshkit/RegisterMeshOp.hpp"
#ifdef HAVE_CAMAL_TRIADV
#  include "meshkit/CAMALTriAdvance.hpp"
#endif
#ifdef HAVE_CAMAL_TET
#  include "meshkit/CAMALTetMesher.hpp"
#endif
#ifdef HAVE_CAMAL_PAVER
#  include "meshkit/CAMALPaver.hpp"
#endif
#ifdef HAVE_NETGEN
#  include "meshkit/NGTetMesher.hpp"
#endif
#ifdef HAVE_MESQUITE
#  include "meshkit/MesquiteOpt.hpp"
#endif
#ifdef HAVE_TRIANGLE
#  include "meshkit/TriangleMesher.hpp"
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

#ifdef HAVE_CAMAL_TRIADV
  REGISTER_MESH_OP(CAMALTriAdvance);
#endif
#ifdef HAVE_CAMAL_PAVER
  REGISTER_MESH_OP(CAMALPaver);
#endif
#ifdef HAVE_CAMAL_TET
  REGISTER_MESH_OP(CAMALTetMesher);
#endif

#endif

#ifdef HAVE_TRIANGLE
  REGISTER_MESH_OP(TriangleMesher);
#endif

#ifdef HAVE_NETGEN
  REGISTER_MESH_OP(NGTetMesher);
#endif    

#ifdef HAVE_MESQUITE
  REGISTER_MESH_OP(MesquiteOpt);
#endif 

} // namespace MeshKit

