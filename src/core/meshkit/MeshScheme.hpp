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
     * \param me_vec MEVector of model entities
     */
  MeshScheme(MKCore *mkcore,
             const MEVector &me_vec = MEVector());
  
    /** \brief Copy constructor
     * \param mesh_scheme Object being copied
     */
  MeshScheme(const MeshScheme &mesh_scheme);
  
    /** \brief operator=
     * \param mesh_scheme Object being copied
     */
  virtual MeshScheme &operator=(const MeshScheme &mesh_scheme);
  
    //! Destructor
  virtual ~MeshScheme();

    /** \brief Return what types of mesh entities this algorithm generates; pure virtual so every scheme must define them
     * \param mesh_types Types handled by this meshop
     */
  virtual void mesh_types(std::vector<moab::EntityType> &mesh_types)=0;
  
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
                              const MEVector &me_vec)
        : MeshOp(mkcore, me_vec)
{}

  //! Copy constructor
inline MeshScheme::MeshScheme(const MeshScheme &mesh_scheme)
        : MeshOp(mesh_scheme)
{}

  //! operator=
inline MeshScheme &MeshScheme::operator=(const MeshScheme &mesh_scheme)
{
  mkCore = mesh_scheme.mk_core();
  opName = mesh_scheme.get_name();
  meSelection = mesh_scheme.me_selection();
  return *this;
}

inline MeshScheme::~MeshScheme()
{}

} // namespace MeshKit

#endif
