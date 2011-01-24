#ifndef MESHOP_HPP
#define MESHOP_HPP

#include "meshkit/Types.h"
#include "meshkit/MKCore.hpp"
#include "meshkit/GraphNode.hpp"
#include "moab/Types.hpp"
#include <vector>

namespace MeshKit {

class ModelEnt;
    
/** \class MeshOp MeshOp.hpp "meshkit/MeshOp.hpp"
 * \brief An operation that operates on mesh data
 *
 * The class encapsulating setup/execute mesh operations on collections of entities.
 * MeshOp objects derive from GraphNode, which encapsulates all graph-related
 * operations on a MeshOp.
 */
class MeshOp : public GraphNode
{
public:

    //! Copy constructor
  MeshOp(const MeshOp &mesh_op);
  
    //! Bare constructor
  MeshOp(MKCore *mkcore, const MEVector &me_vec = MEVector());

    //! Destructor
  virtual ~MeshOp();
  
    /** \brief Add a ModelEnt to this operation's MESelection
     * The MESelection is a map, and can only have a given entity once.
     * \param model_ent ModelEnt being added
     * \return Returns true if model_ent was actually added, false otherwise
     */
  virtual bool add_modelent(ModelEnt *model_ent);

    /** \brief Removes a ModelEnt from this operation's MESelection
     * \param model_ent ModelEnt being removed
     * \return Returns true if model_ent was in the map and was actually removed, false if model_ent wasn't in the map
     */
  virtual bool remove_modelent(ModelEnt *model_ent);

    //! Return a reference to the MESelection list
  virtual MESelection &me_selection();

    //! Return a const reference to the MESelection list
  virtual const MESelection &me_selection() const;

    //! Get the associated MKCore object; this applies a dynamic_cast to the parent's MKGraph member
  MKCore *mk_core() const;

    /** \brief Return what types of mesh entities this algorithm generates; pure virtual so every scheme must define them
     * \param mesh_types Types handled by this meshop
     */
  virtual void mesh_types(std::vector<moab::EntityType> &mesh_types);
  
protected:
    //! MESelection that stores what this operation generated or otherwise worked on
  MESelection meSelection;
  
private:

};

inline const MESelection &MeshOp::me_selection() const
{
  return meSelection;
}

inline MESelection &MeshOp::me_selection() 
{
  return meSelection;
}

inline MKCore *MeshOp::mk_core() const
{
  return dynamic_cast<MKCore*>(mkGraph);
}

inline void MeshOp::mesh_types(std::vector<moab::EntityType> &mesh_types) 
{}

} // namespace MeshKit

#endif

  
