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
    typedef std::vector<ModelEnt*> MEntVector;

      /** \brief Type used to store a set of ModelEnt*'s
       */
    typedef std::set<ModelEnt*> MEntSet;

      /** \brief Type used to store pairs of ModelEnt* and moab::Range, used to associate partial meshes
       * with the ModelEnt's they resolve
       */
    typedef std::map<ModelEnt*, moab::Range> MEntSelection;

      /** \brief Type used to store a set of MeshOp*'s
       */
    typedef std::vector<MeshOp*> MOpVector;

    enum Firmness {DEFAULT, SOFT, HARD};

    enum MeshedState {
        NO_MESH = 0,
        BOUNDARY_MESH,
        PARTIAL_MESH,
        COMPLETE_MESH,
        REFINED_MESH,
        POST_MESH
    };

    enum BooleanType {
        INTERSECT, UNION
    };
    
} // namespace MeshKit

#endif
