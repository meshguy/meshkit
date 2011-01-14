#ifndef MKTYPES
#define MKTYPES

#include <vector>
#include <set>
namespace MeshKit {

    class MeshOp;
    
    enum Firmness {SOFT, LIMP, HARD};

    enum MeshedState {
        NO_MESH = 0,
        BOUNDARY_MESH,
        SOME_MESH,
        COMPLETE_MESH,
        REFINED_MESH,
        POST_MESH
    };
    
} // namespace MeshKit

#endif
