#ifndef MESHSCHEME
#define MESHSCHEME

/** \file MeshScheme.hpp
 */
#include "meshkit/Types.h"
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/SizingFunction.hpp"

namespace MeshKit {
    
class ModelEnt;

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
  MeshScheme(MKCore &mkcore,
             const MEVector &me_vec = MEVector(),
             const SizingFunction &mesh_size = SizingFunction(NULL, 1.0)) throw(Error);
  
    //! Copy constructor
  MeshScheme(const MeshScheme &mesh_scheme) throw(Error);
  
    //! Copy constructor
  virtual MeshScheme &operator=(const MeshScheme &mesh_scheme) throw(Error);
  
    //! Destructor
  virtual ~MeshScheme() throw(Error);

    //! Add a ModelEnt from this scheme's MEVector
  void add_modelent(ModelEnt &model_ent) throw(Error);

    //! Remove a ModelEnt from this scheme's MEVector
  void remove_modelent(ModelEnt &model_ent) throw(Error);

    //! Get the mkCore instance
  MKCore *mk_core() const;

    //! Get a reference to this scheme's sizing function
  virtual const SizingFunction &size() const;
  
    //! Set a reference to this scheme's sizing function
  virtual void size(const SizingFunction &mesh_size) throw(Error);
  
    //! Get scheme name
  virtual std::string name() const;

    //! Get the entity list for this scheme
  virtual const MEVector &get_entities() const;

private:

    //! MeshKit instance
  MKCore *mkCore;
  
    //! Local name for this scheme
  std::string schemeName;
  
    //! List of ModelEnt objects to which this scheme is applied
  MEVector meVector;

    //! Sizing function
  SizingFunction meshSize;
  
};

inline MeshScheme::MeshScheme(MKCore &mkcore,
                              const MEVector &me_vec,
                              const SizingFunction &mesh_size) throw(Error)
        : MeshOp(mkcore), meVector(me_vec), meshSize(mesh_size)
{}

    //! Copy constructor
inline MeshScheme::MeshScheme(const MeshScheme &mesh_scheme) throw(Error)
  : MeshOp(mesh_scheme.mk_core()), schemeName(mesh_scheme.name()), 
    meVector(mesh_scheme.get_entities()), meshSize(mesh_scheme.size())
{}

    //! Copy constructor
inline MeshScheme &MeshScheme::operator=(const MeshScheme &mesh_scheme) throw(Error)
{
  mkCore = mesh_scheme.mk_core();
  schemeName = mesh_scheme.name();
  meVector = mesh_scheme.get_entities();
  meshSize = mesh_scheme.size();
  return *this;
}

MeshScheme::~MeshScheme() throw(Error)
{}

inline void MeshScheme::add_modelent(ModelEnt &model_ent) throw(Error)
{
  meVector.push_back(model_ent);
}
  
inline void MeshScheme::remove_modelent(ModelEnt &model_ent) throw(Error)
{
  MEVector::iterator vit = std::find(meVector.begin(), meVector.end(), model_ent);
  if (vit != meVector.end()) meVector.erase(vit);
  else throw Error(MK_FAILURE, "Entity not found in ModelEnt."); 
}

inline std::string MeshScheme::name() const
{
  return schemeName;
}

inline const SizingFunction &MeshScheme::size() const
{
  return meshSize;
}

inline void MeshScheme::size(const SizingFunction &mesh_size) throw(Error)
{
  meshSize = mesh_size;
}

inline const MEVector &MeshScheme::get_entities() const
{
  return meVector;
}


} // namespace MeshScheme

#endif
