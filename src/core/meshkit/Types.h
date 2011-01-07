#ifndef MKTYPES
#define MKTYPES

#include <vector>
#include <set>
namespace MeshKit {

    class MeshOp;
    class ModelEnt;
    
    enum Firmness {SOFT, LIMP, HARD};

    enum {MK_SUCCESS, MK_FAILURE};
    
    typedef std::vector<MeshOp*> MOVector;

    typedef std::vector<ModelEnt*> MEVector;
    typedef std::set<ModelEnt*> MESet;
    
} // namespace MeshKit

#endif
