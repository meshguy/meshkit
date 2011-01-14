#ifndef MESHSCHEME
#define MESHSCHEME

/** \file MeshScheme.hpp
 */
#include "meshkit/Types.h"
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"

namespace MeshKit {
    
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
     * \param mesh_size SizingFunction instance to be used by this meshing algorithm
     */
  MeshScheme(MKCore *mkcore,
             const MEVector &me_vec = MEVector(),
             const SizingFunction &mesh_size = SizingFunction(NULL, 1.0));
  
    //! Copy constructor
  MeshScheme(const MeshScheme &mesh_scheme);
  
    //! Copy constructor
  virtual MeshScheme &operator=(const MeshScheme &mesh_scheme);
  
    //! Destructor
  virtual ~MeshScheme();

    //! Get a reference to this scheme's sizing function
  virtual const SizingFunction &size() const;
  
    //! Set a reference to this scheme's sizing function
  virtual void size(const SizingFunction &mesh_size);
  
    //! Return what types of mesh entities this algorithm generates; pure virtual so every scheme must define them
  virtual void mesh_types(std::vector<moab::EntityType> &mesh_types)=0;
  
private:

    //! Sizing function
  SizingFunction meshSize;
  
};

inline MeshScheme::MeshScheme(MKCore *mkcore,
                              const MEVector &me_vec,
                              const SizingFunction &mesh_size)
        : MeshOp(mkcore, me_vec), meshSize(mesh_size)
{}

  //! Copy constructor
inline MeshScheme::MeshScheme(const MeshScheme &mesh_scheme)
        : MeshOp(mesh_scheme), meshSize(mesh_scheme.size())
{}

  //! Copy constructor
inline MeshScheme &MeshScheme::operator=(const MeshScheme &mesh_scheme)
{
  mkCore = mesh_scheme.mk_core();
  opName = mesh_scheme.name();
  meSelection = mesh_scheme.me_selection();
  meshSize = mesh_scheme.size();
  return *this;
}

inline MeshScheme::~MeshScheme()
{}

inline const SizingFunction &MeshScheme::size() const
{
  return meshSize;
}

inline void MeshScheme::size(const SizingFunction &mesh_size)
{
  meshSize = mesh_size;
}

} // namespace MeshKit

#endif
