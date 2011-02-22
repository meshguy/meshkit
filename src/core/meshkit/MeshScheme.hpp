#ifndef MESHKIT_MESHSCHEME_HPP
#define MESHKIT_MESHSCHEME_HPP

/** \file MeshScheme.hpp
 */
#include "meshkit/Types.hpp"
#include "meshkit/MeshOp.hpp"

namespace MeshKit {

class MKCore;
    
/** \class MeshScheme MeshScheme.hpp "meshkit/MeshScheme.hpp"
 * \brief A class that generates mesh.
 *
 * A MeshScheme generates mesh on one or more ModelEnt objects.
 */
class MeshScheme : public MeshOp
{
public:

    /** \brief Full constructor
     * \param mkcore MKCore instance to which this scheme instance is associated
     * \param me_vec MEntVector of model entities
     */
  MeshScheme(MKCore *mkcore,
             const MEntVector &me_vec = MEntVector());
  
    /** \brief Copy constructor
     * \param mesh_scheme Object being copied
     */
  MeshScheme(const MeshScheme &mesh_scheme);
  
    //! Destructor
  virtual ~MeshScheme();

    //! Constrain all edges bounding ModelEnts to be even
  void constrain_even();
  
private:
  
};

inline MeshScheme::MeshScheme(MKCore *mkcore,
                              const MEntVector &me_vec)
        : MeshOp(mkcore, me_vec)
{}

  //! Copy constructor
inline MeshScheme::MeshScheme(const MeshScheme &mesh_scheme)
        : MeshOp(mesh_scheme)
{}

inline MeshScheme::~MeshScheme()
{}

} // namespace MeshKit

#endif
