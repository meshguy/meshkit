#ifndef MESHSCHEME
#define MESHSCHEME

/** \file MeshScheme.hpp
 */
#include "meshkit/Types.h"
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

    //! Setup function, called in reverse order before execute
  virtual void setup();

    //! Execute function, called in forward order after setup
  virtual void execute();

    //! This function is pure virtual, all derived classes must define
  virtual void setup_this()=0;

    //! This function is pure virtual, all derived classes must define
  virtual void execute_this()=0;

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
