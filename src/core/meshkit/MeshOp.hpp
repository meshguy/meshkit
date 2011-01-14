#ifndef MESHOP
#define MESHOP

#include "meshkit/Types.h"
#include "meshkit/Error.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/GraphNode.hpp"
#include "iGeom.h"
#include "moab/Interface.hpp"
#include <vector>

namespace MeshKit {

class MKCore;
    
typedef std::vector<MeshOp*> MOVector;

/** \class MeshOp MeshOp.hpp "meshkit/MeshOp.hpp"
 * \brief An operation that operates on mesh data
 *
 * The purpose of this class is to enable mesh-related operations to be added to
 * an operation graph.  Once the graph has been created, it can be executed by calling
 * the execute() member function.
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
  
    //! Get operation name
  virtual std::string name() const;

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

    //! Get the associated MKCore object
  MKCore *mk_core() const;

protected:
    //! Local name for this operation
  std::string opName;
  
    //! MESelection that stores what this operation generated or otherwise worked on
  MESelection meSelection;
  
    //! MKCore associated with this MeshOp
  MKCore *mkCore;
  
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
  return mkCore;
}

}

#endif

  
