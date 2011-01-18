#ifndef MKTYPES
#define MKTYPES

#include <vector>
#include <set>
#include <map>

#include "moab/Range.hpp"

namespace MeshKit {

    class MeshOp;
    class ModelEnt;
    
      /** \brief Type used to store a vector of ModelEnt*'s
       */
    typedef std::vector<ModelEnt*> MEVector;

      /** \brief Type used to store a set of ModelEnt*'s
       */
    typedef std::set<ModelEnt*> MESet;

      /** \brief Type used to store pairs of ModelEnt* and moab::Range, used to associate partial meshes
       * with the ModelEnt's they resolve
       */
    typedef std::map<ModelEnt*, moab::Range> MESelection;

      /** \brief Type used to store a set of MeshOp*'s
       */
    typedef std::vector<MeshOp*> MOVector;

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
