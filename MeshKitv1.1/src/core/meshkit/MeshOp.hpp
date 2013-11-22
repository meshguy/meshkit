#ifndef MESHKIT_MESHOP_HPP
#define MESHKIT_MESHOP_HPP

#include "meshkit/Types.hpp"
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
  MeshOp(MKCore *mkcore, const MEntVector &me_vec = MEntVector());

    //! Destructor
  virtual ~MeshOp();
  
    /** \brief Add a ModelEnt to this operation's MEntSelection
     * The MEntSelection is a map, and can only have a given entity once.
     * \param model_ent ModelEnt being added
     * \return Returns true if model_ent was actually added, false otherwise
     */
  virtual bool add_modelent(ModelEnt *model_ent);

    /** \brief Removes a ModelEnt from this operation's MEntSelection
     * \param model_ent ModelEnt being removed
     * \return Returns true if model_ent was in the map and was actually removed, false if model_ent wasn't in the map
     */
  virtual bool remove_modelent(ModelEnt *model_ent);

    //! Return a reference to the MEntSelection list
  virtual MEntSelection &me_selection();

    //! Return a const reference to the MEntSelection list
  virtual const MEntSelection &me_selection() const;

    //! Get the associated MKCore object; this applies a dynamic_cast to the parent's MKGraph member
  MKCore *mk_core() const;
  
    /** \brief Return the mesh entity types operated on by this scheme
     * \return array terminated with \c moab::MBMAXTYPE
     */
  virtual const moab::EntityType* mesh_types_arr() const = 0;

    /** \brief Return what types of mesh entities this algorithm generates; pure virtual so every scheme must define them
     * \param mesh_types Types handled by this meshop
     */
  void mesh_types(std::vector<moab::EntityType> &mesh_types);
  
    /** \brief Check that bounding entities have an assigned MeshOp, and create them for ones that don't
     *
     * Uses default MeshOp for a given dimension from MeshOpFactory.  If there isn't a registered MeshOp for
     * the dimension requested, this function throws an exception with mode MK_MESHOP_NOT_FOUND.
     */
  void setup_boundary();

    /** \brief Helper function for meshop registration, returns true if specified ModelEnt is a vertex
     * \param model_ent Model entity being evaluated
     * \return True if model_ent has dimension() == 0
     */
  static bool canmesh_vertex(ModelEnt *model_ent);
  
    /** \brief Helper function for meshop registration, returns true if specified ModelEnt is an edge
     * \param model_ent Model entity being evaluated
     * \return True if model_ent has dimension() == 1
     */
  static bool canmesh_edge(ModelEnt *model_ent);
  
    /** \brief Helper function for meshop registration, returns true if specified ModelEnt is a face
     * \param model_ent Model entity being evaluated
     * \return True if model_ent has dimension() == 2
     */
  static bool canmesh_face(ModelEnt *model_ent);
  
    /** \brief Helper function for meshop registration, returns true if specified ModelEnt is a region
     * \param model_ent Model entity being evaluated
     * \return True if model_ent has dimension() == 3
     */
  static bool canmesh_region(ModelEnt *model_ent);

  /**
   * brief  mesh-based geometry model ents based on previous ops
   * it will be used for Camal advancing front meshers
   */
  void create_model_ents_from_previous_ops();

protected:
    //! MEntSelection that stores what this operation generated or otherwise worked on
  MEntSelection mentSelection;
  
private:

};

inline const MEntSelection &MeshOp::me_selection() const
{
  return mentSelection;
}

inline MEntSelection &MeshOp::me_selection() 
{
  return mentSelection;
}

inline MKCore *MeshOp::mk_core() const
{
  return dynamic_cast<MKCore*>(mkGraph);
}

} // namespace MeshKit

#endif

  
