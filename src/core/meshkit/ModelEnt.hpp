#ifndef MESHKIT_MODELENT_HPP
#define MESHKIT_MODELENT_HPP

#include "meshkit/iGeom.hpp"
#include "meshkit/Types.hpp"
#include "moab/Interface.hpp"
#include <cassert>
#include <vector>
#include <set>
#include <map>

namespace MeshKit {

class ModelEnt;
class MKCore;
class IAVariable;

/** \class ModelEnt ModelEnt.hpp "meshkit/ModelEnt.hpp"
 * \brief The class used in MeshKit for referring to model entities and mesh associated with them.
 *
 * Instances of this class are references to the entities stored in the geometry and mesh
 * databases.  Destruction/construction of ModelEnt instances does not cause destruction/construction
 * of the corresponding database entities (one exception to this statement is that an argument to the
 * constructor can request that the mesh entity set be created for the corresponding geometry entity).
 *
 * A ModelEnt instance keeps a reference to the MKCore object to which it is associated; through that
 * instance it can get access to the geometry/mesh/relation database instances.
 * \nosubgrouping
 */
class ModelEnt
{
public:

    /** \name Constructor/destructor
     */
    /**@{*/

    /** \brief Constructor; mesh entity can be missing, in which case it's retrieved or created
     *
     * \param mk MeshKit instance
     * \param geom_ent Geometry entity corresponding to this ModelEntity
     * \param igeom iGeom instance
     * \param mesh_ent Mesh entity set corresponding to this ModelEntity
     * \param sizing_index Sizing function index (from MBCore) to be used to mesh this entity
     */
  ModelEnt(MKCore *mk,
           iGeom::EntityHandle geom_ent,
           int geom_index,
           moab::EntityHandle mesh_ent = NULL,
           int mesh_index = -1,
           int irel_index = -1,
           int sizing_index = -1);

    /** \brief Constructor; mesh entity can be missing, in which case it's retrieved or created
     *
     * \param mk MeshKit instance
     * \param geom_set Geometry entity set corresponding to this ModelEntity
     * \param igeom iGeom instance
     * \param mesh_ent Mesh entity set corresponding to this ModelEntity
     * \param sizing_index Sizing function index (from MBCore) to be used to mesh this entity
     */
  ModelEnt(MKCore *mk,
           iGeom::EntitySetHandle geom_set,
           int geom_index,
           moab::EntityHandle mesh_ent = NULL,
           int mesh_index = -1,
           int irel_index = -1,
           int sizing_index = -1);

    /** \brief Destructor
     */
  virtual ~ModelEnt();

    /**@}*/

    /** \name Topology
     */
    /**@{*/

    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
  void children(std::vector<MEntVector> &child_ents) const;

    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
  void children(std::vector<MEntVector> &child_ents);

    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
  void children(std::vector<std::vector<moab::EntityHandle> > &child_ents) const;

    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
  void children(std::vector<std::vector<iGeom::EntityHandle> > &child_ents) const;

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function appends adjacencies to the input vector.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim, MEntVector &adjs) const;

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function appends adjacencies to the input vector.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim, std::vector<iGeom::EntityHandle> &adjs) const;

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function appends adjacencies to the input vector.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim, std::vector<moab::EntityHandle> &adjs) const;

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function adds adjacencies to the input range.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim, moab::Range &adjs) const;

    /** \brief Return a shared entity of specified dimension
     *
     * If no shared entities are found, NULL is returned.  If more than one are found,
     * an Error is thrown with MK_MULTIPLE_FOUND as the code.
     * \param ent2 Other entity
     * \param to_dim Dimension of shared entity
     * \return Shared entity
     */
  ModelEnt *shared_entity(ModelEnt *ent2, int to_dim);
  
    /** \brief Get adjacent entities, with specified boolean on results
     * \param from_ents Entities whose adjacencies are being queried
     * \param to_dim Dimension of adjacencies requested
     * \param to_ents Adjacent entities
     * \param op_type Boolean type, intersect or union
     * \param only_to_ents if is true, get adjacencies from to_ents
     */
  void get_adjs_bool(MEntVector &from_ents,
                     int to_dim,
                     MEntVector &to_ents,
                     BooleanType op_type,
                     bool only_to_ents = true);
    /**@}*/

    /** \name Geometric evaluation
     */
    /**@{*/

    //! Return the topological dimension of this model entity
  int dimension() const;
  
    //! Measure of this entity, -DBL_MAX for vertices, or length, area, or volume
  double measure() const;
  
    //! Similar as measure, but based on mesh 
  double measure_discrete() const;
  
    /** \brief Evaluate entity at a point
     *
     * If any of the evaluation results are not wanted, pass NULL for them.  At least one
     * must be provided.  This function throws an error if called on a vertex or volume and
     * non-NULL is passed for direction or curvature.  Direction is the tangent for edges,
     * normal for surfaces.
     * \param x Start point x
     * \param y Start point y
     * \param z Start point z
     * \param close Closest point
     * \param direction Tangent or normal returned for edges or surfaces
     * \param curvature1 1st curvature for edges or surfaces
     * \param curvature2 2nd curvature for edges or surfaces
     */
  void evaluate(double x, double y, double z, 
                double *close = NULL,
                double *direction = NULL,
                double *curvature1 = NULL,
                double *curvature2 = NULL) const;
  
    /** \brief Evaluate entity at a point, based on mesh data
     *
     * If any of the evaluation results are not wanted, pass NULL for them.  At least one
     * must be provided.  This function throws an error if called on a vertex or volume and
     * non-NULL is passed for direction or curvature.  Direction is the tangent for edges,
     * normal for surfaces.
     * \param x Start point x
     * \param y Start point y
     * \param z Start point z
     * \param close Closest point
     * \param direction Tangent or normal returned for edges or surfaces
     * \param curvature1 1st curvature for edges or surfaces
     * \param curvature2 2nd curvature for edges or surfaces
     */
  void evaluate_discrete(double x, double y, double z, 
                         double *close = NULL,
                         double *direction = NULL,
                         double *curvature1 = NULL,
                         double *curvature2 = NULL) const;
  
    //! Get the id of the geometry/mesh entity, if any 
  int id() const;
  
    /**@}*/

    /** \name Mesh
     */
    /**@{*/

    /** \brief Create mesh set for this ModelEnt
     *
     * This function creates the entity set for a geometry entity, and tags that set
     * according to convention (GEOM_DIMENSION and GLOBAL_ID).  It also sets the relation
     * with the corresponding geometry entity.  The value of flag controls whether a
     * vector- or set-based entity set is created:
     * flag = 1 : create a vector-based set
     * flag = 0 : create a set-based set
     * flag = -1: check dimension of iGeomEnt; if 1, vector-based set, otherwise set-based set
     * \param ordered_flag Ordered flag
     */
  void create_mesh_set(int ordered_flag = -1);

    /** \brief Commit mesh to a model entity
     *
     * Takes the input mesh entities, adds them to the entity set for this model entity,
     * and (if both-type relation on the mesh side) sets the relations to the corresponding
     * geometry entity.
     * \param mesh_ents Mesh entities being assigned to this model entity
     * \param mstate The meshed state after this mesh is added
     */
  void commit_mesh(moab::Range &mesh_ents,
                   MeshedState mstate);
  
    /** \brief Commit mesh to a model entity
     *
     * Takes the input mesh entities, adds them to the entity set for this model entity,
     * and (if both-type relation on the mesh side) sets the relations to the corresponding
     * geometry entity.
     * \param mesh_ents Pointer to mesh entities being assigned to this model entity
     * \param num_ents Number of mesh entities in list
     * \param mstate The meshed state after this mesh is added
     */
  void commit_mesh(moab::EntityHandle *mesh_ents,
                   int num_ents,
                   MeshedState mstate);

    /** \brief Return the mesh on this entity, of the specified dimension or (if dim=-1) all dimensions
     *
     * This function will return the mesh entities on bounding entities too, if requested.  This is useful
     * for example for assembling the nodes bounding a surface, as some of those nodes are on model
     * vertices.  If requesting mesh on a periodic edge, end vertex will appear on both ends of the returned
     * vector.
     * \param dim Dimension of entities requested
     * \param ments Entities
     * \param bdy_too If true, returns (dim)-dimensional entities from bounding entities too
     */
  void get_mesh(int dim,
                std::vector<moab::EntityHandle> &ments,
                bool bdy_too = false);

    /** \brief Return the mesh bounding this entity, their senses, and optionally the loops/shells
     *
     * In the case where vertices are requested (dim=0), vertices on end of loops are *not* repeated
     * \param dim Dimension of boundary entities requested
     * \param bdy Boundary entities
     * \param senses Senses of boundary entities
     * \param group_sizes If non-NULL, pointer to vector where group sizes will be returned
     */
  void boundary(int dim,
                std::vector<moab::EntityHandle> &bdy,
                std::vector<int> *senses = NULL,
                std::vector<int> *group_sizes = NULL);

    /** \brief Return the model entities bounding this entity, their senses, and optionally the loops/shells
     * \param dim Dimension of boundary entities requested
     * \param elements Boundary entities
     * \param senses Senses of boundary entities
     * \param group_sizes If non-NULL, pointer to vector where group sizes will be returned
     */
  void boundary(int dim,
                MEntVector &elements,
                std::vector<int> *senses,
                std::vector<int> *group_sizes = NULL);

    /** \brief Return mesh bounding this model entity in a Range (unordered wrt model entity)
     *
     * This function returns boundary entities in a Range, which means they're not guaranteed
     * to be ordered consistently with respect to the model entity orientation.  The input range
     * is appended to (i.e. it is not emptied first).
     * \param dim Dimension of requested boundary entities
     * \param ents Boundary entities requested
     */
  void boundary(int dim, 
                moab::Range &ents) const;

    /** \brief From a vector of entity handles, return connectivity as integer ids, and coords of vertices if requested
     *
     * The entity handles are ordered and made unique, assigned ids, then the id vector is assembled
     * in the same order as the input handle vector.  Optionally this function returns the moab tag 
     * used to assign ids, and the moab range used to get unique ids
     * \param ents Vector of entities whose ids to assign
     * \param senses Senses of entities in ents; if NULL, assume forward sense, or no sense if ents are vertices
     * \param tagh Tag handle to use for ids; if zero, local tag is used then deleted
     * \param ents_as_ids Vector of ids that get assigned, ordered same as connectivity of ents (or ents itself if they're vertices)
     * \param coords Coordinates of vertices corresponding to ents_as_ids
     * \param verts_range Vertices in ents, as a range; if non-NULL, populated range is returned; this range is
     *       cleared in this function
     * \param start_index Starting index for vertex numbering
     */
  void get_indexed_connect_coords(std::vector<moab::EntityHandle> &ents,
                                  std::vector<int> *senses,
                                  moab::Tag tagh,
                                  std::vector<int> &ents_as_ids,
                                  std::vector<double> &coords,
                                  moab::Range *verts_range = NULL,
                                  int start_index = 0);
  
    /** \brief Get meshed state
     * \return Meshed state
     */
  MeshedState get_meshed_state();

    /** \brief Set the meshed state
     * \param mstate
     */
  void set_meshed_state(MeshedState mstate);

    /**@}*/

    /** \name Mesh intervals, sizing functions
     */

    /**@{*/

    /** \brief Get sizing function index
     * \return Returns -1 if it has not been set for this entity
     */
  int sizing_function_index() const;

    /** \brief Set sizing function index
     * \param index Sizing function index being set
     * \param children_too If true, sizing function index is propagated to child (bounding) entities
     */
  void sizing_function_index(int index, bool children_too = true);

    /** \brief Get mesh interval size, if any
     * Returns -1 if no size set on this entity.  If intervals are set and this is a model edge, returns computed size.
     * \return Interval size for this ModelEnt.
     */
  double mesh_interval_size() const;

    //! Get intervals
  int mesh_intervals() const;

    /** \brief Set intervals
     * \param ints Set value
     */
  void mesh_intervals(int ints);

    //! Get firmness
  Firmness interval_firmness() const;

    /** \brief Set firmness
     * \param firm Set value
     */
  void interval_firmness(Firmness firm);

    //! Get constrainEven
  bool constrain_even() const;
  
    /** \brief Set constrainEven
     * \param flag Set value
     */
  void constrain_even(bool flag);

    /**@}*/

    /** \name MeshOps
     */

    /**@{*/

    /** \brief Add a MeshOp that points to this ModelEnt
     *
     * \param meshop MeshOp to add
     */
  void add_meshop(MeshOp *meshop);
  
    /** \brief Remove a MeshOp that pointed to this ModelEnt
     *
     * \param meshop MeshOp to remove
     */
  void remove_meshop(MeshOp *meshop);
  
    /** \brief Get MeshOps pointing to this ModelEnt
     *
     * \param meshops MeshOps returned
     */
  void get_meshops(std::vector<MeshOp*> &meshops);

  /** \brief decide if there is a meshop defined for this model entity
       *
       * return true if the list of meshops is not empty
       */
  bool is_meshops_list_empty()
  {
    return meshOps.empty();
  }

    /**@}*/

    /** \name Interfaces, entities
     */

    /**@{*/

    //! Get the MKCore object
  MKCore *mk_core() const;
  
    //! Get geometry entity handle
  iGeom::EntityHandle geom_handle() const;

    //! Get geometry entity handle for a given mesh entity set
  iGeom::EntityHandle geom_handle(moab::EntityHandle ment) const;

    /** \brief Get iGeom instance associated with this model entity
     * \return Returns 0 by default
     */
  iGeom * igeom_instance() const ;

    /** \brief Get moab instance associated with this model entity
     * \return Returns 0 by default
     */
  moab::Interface * moab_instance() const ;

    //! Get mesh entity set handle
  moab::EntityHandle mesh_handle() const;

    //! Get mesh entity set handle for a given geometry entity
  moab::EntityHandle mesh_handle(iGeom::EntityHandle gent) const;

    //! Get mesh entity set handle for a given geometry entity
  moab::EntityHandle mesh_handle(iGeom::EntitySetHandle gent) const;

  int iGeomIndex()
   {return igeomIndex;}

  int iRelPairIndex()
   {return irelIndex;}

  IAVariable *ia_var() {return iaVariable;}
  const IAVariable *ia_var() const {return iaVariable;}
  void ia_var(IAVariable *var) {iaVariable = var;}

private:

    /** \brief Return the next entity in the loop, using winding number
     * \param this_edge Edge next to one being requested
     * \param this_sense Sense of this_edge
     * \param tmp_adjs Optional vector of candidates
     * \return Next edge in loop
     */
  ModelEnt *next_winding(ModelEnt *this_edge, 
                            int this_sense, 
                            MEntVector &tmp_adjs);
  
    /** \brief Initialize mesh set parents/children to parallel geometric topology
     */
  void init_parents_children();

    /** \brief Initialize mesh set contents to parallel geometry set contents
     */
  void init_group_contents();
  
    //! Set senses tag on moabEntSet wrt adjacent entities
  void set_downward_senses();
  
    //! MeshKit instance to which this model entity is associated
  MKCore *mkCore;

    //! iGeom instance index associated with this model entity
  int igeomIndex;

    //! MOAB instance index associated with this model entity
  int meshIndex;

    //! iRel pair index associated with this model entity
  int irelIndex;

    //! Geometry entity for this model entity
  iGeom::EntityHandle iGeomEnt;
  
    //! Geometry set for this model entity
  iGeom::EntitySetHandle iGeomSet;
  
    //! Mesh entity set for this model entity
  moab::EntityHandle moabEntSet;
  
    //! Sizing function associated with this model entity
  int sizingFunctionIndex;
  
    //! Mesh intervals for this model entity
  int meshIntervals;
  
    //! Mesh interval firmness for this model entity
  Firmness intervalFirmness;

    //! Constrain mesh intervals to be even
  bool constrainEven;
  
    //! Meshed state of this entity
  MeshedState meshedState;

    //! MeshOps pointing to this entity
  std::vector<MeshOp*> meshOps;

    //! interval assignment variable
  IAVariable *iaVariable;
};

inline MKCore *ModelEnt::mk_core() const
{
  return mkCore;
}

inline iGeom::EntityHandle ModelEnt::geom_handle() const
{
  return iGeomEnt;
}

inline moab::EntityHandle ModelEnt::mesh_handle() const
{
  return moabEntSet;
}

inline void ModelEnt::children(std::vector<std::vector<moab::EntityHandle> > &child_ents) const 
{
  std::vector<MEntVector> tmp_vec;
  children(tmp_vec);
  std::vector<moab::EntityHandle> tmp_vec2;
  for (std::vector<MEntVector>::iterator vit1 = tmp_vec.begin(); vit1 != tmp_vec.end(); vit1++) {
    tmp_vec2.clear();
    for (MEntVector::iterator vit2 = (*vit1).begin(); vit2 != (*vit1).end(); vit2++)
      tmp_vec2.push_back((*vit2)->mesh_handle());
    child_ents.push_back(tmp_vec2);
  }
}

    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
inline void ModelEnt::children(std::vector<std::vector<iGeom::EntityHandle> > &child_ents) const
{
  std::vector<MEntVector> tmp_vec;
  children(tmp_vec);
  std::vector<iGeom::EntityHandle> tmp_vec2;
  for (std::vector<MEntVector>::iterator vit1 = tmp_vec.begin(); vit1 != tmp_vec.end(); vit1++) {
    tmp_vec2.clear();
    for (MEntVector::iterator vit2 = (*vit1).begin(); vit2 != (*vit1).end(); vit2++)
      tmp_vec2.push_back((*vit2)->geom_handle());
    child_ents.push_back(tmp_vec2);
  }
}

inline int ModelEnt::sizing_function_index() const 
{
  return sizingFunctionIndex;
}

    //! Get intervals
inline int ModelEnt::mesh_intervals() const 
{
  return meshIntervals;
}

//! Set intervals
inline void ModelEnt::mesh_intervals(int ints)
{
  meshIntervals = ints;
}

//! Get firmness
inline Firmness ModelEnt::interval_firmness() const 
{
  return intervalFirmness;
}

//! Set firmness
inline void ModelEnt::interval_firmness(Firmness firm) 
{
  intervalFirmness = firm;
}

inline bool ModelEnt::constrain_even() const 
{
  return constrainEven;
}

inline void ModelEnt::constrain_even(bool flag) 
{
  constrainEven = flag;
}

inline MeshedState ModelEnt::get_meshed_state() 
{
  return meshedState;
}

inline void ModelEnt::set_meshed_state(MeshedState mstate) 
{
  meshedState = mstate;
}

inline void ModelEnt::evaluate_discrete(double x, double y, double z, 
                                        double *close,
                                        double *direction,
                                        double *curvature1,
                                        double *curvature2) const
{
  evaluate(x, y, z, close, direction, curvature1, curvature2);
}

      /** \brief Add a MeshOp that points to this ModelEnt
       * \param meshop MeshOp to add
       */
inline void ModelEnt::add_meshop(MeshOp *meshop) 
{
  assert(std::find(meshOps.begin(), meshOps.end(), meshop) == meshOps.end());
  meshOps.push_back(meshop);
}
  
    /* \brief Remove a MeshOp that pointed to this ModelEnt
     * \param meshop MeshOp to remove
     */
inline void ModelEnt::remove_meshop(MeshOp *meshop)
{
  assert(std::find(meshOps.begin(), meshOps.end(), meshop) != meshOps.end());
  meshOps.erase(std::remove(meshOps.begin(), meshOps.end(), meshop), meshOps.end());
}
  
    /* \brief Get MeshOps pointing to this ModelEnt
     * \param meshop MeshOps returned
     */
inline void ModelEnt::get_meshops(std::vector<MeshOp*> &meshops) 
{
  meshops.resize(meshOps.size());
  std::copy(meshOps.begin(), meshOps.end(), meshops.begin());
}

// TODO: unimplemented functions (these are just here because SWIG expects some
// implementation!
inline void ModelEnt::get_adjacencies(int dim, std::vector<iGeom::EntityHandle> &adjs) const {}
inline void ModelEnt::get_adjacencies(int dim, std::vector<moab::EntityHandle> &adjs) const {}
inline void ModelEnt::get_adjacencies(int dim, moab::Range &adjs) const {}
inline void ModelEnt::children(std::vector<MEntVector> &child_ents) const {}
inline void ModelEnt::children(std::vector<MEntVector> &child_ents) {}

} // namespace meshkit

#endif
