#ifndef MESHOP_HPP
#define MESHOP_HPP

#include "meshkit/Types.h"
#include "meshkit/Error.hpp"
#include "iGeom.h"
#include "moab/Interface.hpp"
#include "lemon/list_graph.h"
#include <vector>

namespace MeshKit {

class MKCore;
class ModelEnt;
    
/** \class MeshOp MeshOp.hpp "meshkit/MeshOp.hpp"
 * \brief An operation that operates on mesh data
 *
 * The purpose of this class is to enable mesh-related operations to be added to
 * an operation graph.  Once the graph has been created, it can be executed by calling
 * the execute() member function.
 */
class MeshOp
{
public:

    //! Copy constructor
  MeshOp(const MeshOp &mesh_op);
  
    //! Bare constructor
  MeshOp(MKCore *mkcore, const MEVector &me_vec = MEVector());

    //! Destructor
  virtual ~MeshOp();
  
    //! Get operation name
  virtual std::string get_name() const;

    //! Get operation name
  virtual void set_name(std::string new_name);

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

    //! Get the graph node corresponding to this MeshOp
  lemon::ListDigraph::Node op_node() const;

    //! Return an iterator over incoming graph edges
  lemon::ListDigraph::InArcIt in_arcs();
  
    //! Return an iterator over outgoing graph edges
  lemon::ListDigraph::OutArcIt out_arcs();

    /** \brief Return the MeshOp at the start of a graph edge
     * \param arc Edge being queried
     * \return MeshOp corresponding to graph node at start of edge
     */
  MeshOp *source(lemon::ListDigraph::Arc arc);

    /** \brief Return the MeshOp at the end of a graph edge
     * \param arc Edge being queried
     * \return MeshOp corresponding to graph node at end of edge
     */
  MeshOp *target(lemon::ListDigraph::Arc arc);
  
    /** \brief Return what types of mesh entities this algorithm generates; pure virtual so every scheme must define them
     * \param mesh_types Types handled by this meshop
     */
  virtual void mesh_types(std::vector<moab::EntityType> &mesh_types)=0;
  
    //! Setup function, called in reverse order before execute
  virtual void setup();

    //! Execute function, called in forward order after setup
  virtual void execute();
  
    //! Pure virtual, derived class must define
  virtual void setup_this()=0;

    //! Pure virtual, derived class must define
  virtual void execute_this()=0;

protected:
    //! Local name for this operation
  std::string opName;
  
    //! MESelection that stores what this operation generated or otherwise worked on
  MESelection meSelection;
  
    //! MKCore associated with this MeshOp
  MKCore *mkCore;

    //! The graph node associated with this MeshOP
  lemon::ListDigraph::Node opNode;

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

    //! Get the graph node corresponding to this MeshOp
inline lemon::ListDigraph::Node MeshOp::op_node() const
{
  return opNode;
}

inline std::string MeshOp::get_name() const
{
  return opName;
}

inline void MeshOp::set_name(std::string new_name)
{
  opName = new_name;
}

}

#endif

  
