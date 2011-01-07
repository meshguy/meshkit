#ifndef MODELENT
#define MODELENT

#include "meshkit/Types.h"
#include "meshkit/SizingFunction.hpp"
#include "iGeom.h"
#include "moab/Interface.hpp"
#include <vector>

namespace MeshKit {

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
     * \param mesh_ent Mesh entity set corresponding to this ModelEntity
     * \param mesh_size Sizing function to be used to mesh this entity
     */
  virtual ModelEnt(MKCore *mk,
                   iGeom::EntityHandle &geom_ent,
                   moab::EntityHandle &mesh_ent = 0,
                   const SizingFunction &mesh_size = SizingFunction(1.0)) throw(Error);
    /** \brief Destructor
     */
  virtual ~ModelEnt() throw(Error);

    /**@}*/

    /** \name Topology
     */
    /**@{*/

    /** \brief Return parents as ModelEnts
     *
     * \param parent_ents Parent entities returned
     */
  void parents(MEVector &parent_ents) const throw(Error);

    /** \brief Return parents as mesh entity sets
     *
     * \param parent_ents Parent entities returned
     */
  void parents(std::vector<moab::EntityHandle> &parent_ents) const throw(Error);

    /** \brief Return parents as GeomEnts
     *
     * \param parent_ents Parent entities returned
     */
  void parents(std::vector<iGeom::EntityHandle> &parent_ents) const throw(Error);
  
    /** \brief Return children as ModelEnts
     *
     * \param child_ents Child entities returned
     */
  void children(MEVector &child_ents) const throw(Error);

    /** \brief Return children as mesh entity sets
     *
     * \param child_ents Child entities returned
     */
  void children(std::vector<moab::EntityHandle> &child_ents) const throw(Error);

    /** \brief Return children as GeomEnts
     *
     * \param child_ents Child entities returned
     */
  void children(std::vector<iGeom::EntityHandle> &child_ents) const throw(Error);
  
    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
  void children(std::vector<MEVector> &child_ents) const throw(Error);

    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
  void children(std::vector<std::vector<moab::EntityHandle> > &child_ents) const throw(Error);

    /** \brief Return children as vectors of vectors, e.g. for loops or shells
     *
     * No ordered flag, since that's implied by definition
     * \param child_ents Child entities returned
     */
  void children(std::vector<std::vector<iGeom::EntityHandle> > &child_ents) const throw(Error);

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function appends adjacencies to the input vector.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim = -1, MEVector &adjs) const throw(Error);

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function appends adjacencies to the input vector.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim = -1, std::vector<iGeom::EntityHandle> &adjs) const throw(Error);

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function appends adjacencies to the input vector.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim = -1, std::vector<moab::EntityHandle> &adjs) const throw(Error);

    /** \brief Return adjacencies to this model entity of the input dimension 
     *
     * This function adds adjacencies to the input range.
     * \param dim Dimension of entities requested; if -1, all adjacent entities are returned
     * \param adjs Adjacent entities 
     */
  void get_adjacencies(int dim = -1, moab::Range &adjs) const throw(Error);

    /**@}*/

    /** \name Geometric evaluation
     */
    /**@{*/

    //! Return the topological dimension of this model entity
  int dimension() const throw(Error);
  
    //! Measure of this entity, -DBL_MAX for vertices, or length, area, or volume
  double measure() const throw(Error);
  
    //! Similar as measure, but based on mesh 
  double measure_discrete() const throw(Error);
  
    //! Closest point to input
  unit_vector closest(double x, double y, double z) const throw(Error);

    //! Similar to closest_point, but based on mesh
  unit_vector closest_discrete(double x, double y, double z) const throw(Error);

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
     * \param flag Ordered flag
     */
  void create_mesh_set(int flag = -1) throw(Error);

    /** \brief Commit mesh to a model entity
     *
     * Takes the input mesh entities, adds them to the entity set for this model entity,
     * and (if both-type relation on the mesh side) sets the relations to the corresponding
     * geometry entity.
     * \param mesh_ents Mesh entities being assigned to this model entity
     */
  void commit_mesh(moab::Range &mesh_ents) throw(Error);
  
    /** \brief Commit mesh to a model entity
     *
     * Takes the input mesh entities, adds them to the entity set for this model entity,
     * and (if both-type relation on the mesh side) sets the relations to the corresponding
     * geometry entity.
     * \param mesh_ents Mesh entities being assigned to this model entity
     */
  void commit_mesh(std::vector<moab::EntityHandle> &mesh_ents) throw(Error);
  
    /** \brief Return mesh bounding this model entity
     *
     * This function returns boundary entities in a list, ordered consistently with
     * the underlying geometry entity's natual direction.
     * \param dim Dimension of requested boundary entities
     * \param mesh_ents Ordered list of boundary entities returned
     */
  void boundary(int dim, std::vector<moab::EntityHandle> &mesh_ents) const throw(Error);

    /** \brief Return mesh bounding this model entity, distinguished by orientation
     *
     * This function returns boundary entities in two lists, the first with entities oriented
     * consistently with the underlying model entity, and the other with entities oriented opposite
     * the model entity.
     * \param dim Dimension of requested boundary entities
     * \param forward_ents Boundary entities with order consistent with model entity
     * \param reverse_ents Boundary entities with order opposite that of model entity
     */
  void boundary(int dim, 
                std::vector<moab::EntityHandle> &forward_ents,
                std::vector<moab::EntityHandle> &reverse_ents) const throw(Error);

    /** \brief Return mesh bounding this model entity in a Range (unordered wrt model entity)
     *
     * This function returns boundary entities in a Range, which means they're not guaranteed
     * to be ordered consistently with respect to the model entity orientation.  The input range
     * is appended to (i.e. it is not emptied first).
     * \param dim Dimension of requested boundary entities
     * \param ents Boundary entities requested
     */
  void boundary(int dim, 
                moab::Range &ents) const throw(Error);

    /**@}*/

    /** \name Member get/set
     */

    //! Get the MKCore object
  MKCore *mk_core() const;
  
    /**@{*/
    //! Get geometry entity handle
  iGeom::EntityHandle geom_handle() const;

    //! Get (const) geometry entity handle
  const iGeom::EntityHandle geom_handle() const;

    //! Get mesh entity set handle
  moab::EntityHandle mesh_handle() const;

    //! Get (const) mesh entity set handle
  const moab::EntityHandle mesh_handle() const;

    //! Get sizing function
  const SizingFunction &sizing_function() const;

    //! Set sizing function
  void sizing_function(const SizingFunction &sizing);

    //! Get intervals
  int mesh_intervals() const;

    //! Set intervals
  void mesh_intervals(int ints);

    //! Get firmness
  Firmness interval_firmness() const;

    //! Set firmness
  void interval_firmness(Firmness firm);
  
    /**@}*/

private:
    //! MeshKit instance to which this model entity is associated
  MKCore *mkCore;

    //! Geometry entity for this model entity
  iGeom::EntityHandle iGeomEnt;
  
    //! Mesh entity set for this model entity
  moab::EntityHandle moabEntSet;
  
    //! Sizing function associated with this model entity
  SizingFunction sizingFunction;
  
    //! Mesh intervals for this model entity
  int meshIntervals;
  
    //! Mesh interval firmness for this model entity
  Firmness intervalFirmness;
};

inline MKCore *ModelEnt::mk_core() const
{
  return mkCore;
}

inline iGeom::EntityHandle ModelEnt::geom_handle() const
{
  return iGeomEnt;
}

inline const iGeom::EntityHandle ModelEnt::geom_handle() const
{
  return iGeomEnt;
}

inline moab::EntityHandle ModelEnt::mesh_handle() const
{
  return moabEntSet;
}

inline const moab::EntityHandle ModelEnt::mesh_handle() const
{
  return moabEntSet;
}

inline const SizingFunction &ModelEnt::sizing_function() const 
{
  return sizingFunction;
}

inline void ModelEnt::sizing_function(const SizingFunction &sizing) 
{
  sizingFunction = sizing;
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

} // namespace meshkit
