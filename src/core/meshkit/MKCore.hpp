#ifndef MKCORE
#define MKCORE

/** \file MKCore.hpp
 */
#include "meshkit/Error.hpp"
#include "meshkit/MeshOp.hpp"
#include "iGeom.hh"
#include "moab/Interface.hpp"
#include "iRel.hh"
#include <vector>

namespace MeshKit {

/** \class MKCore MKCore.hpp "meshkit/MKCore.hpp"
 * \brief The core MeshKit instance
 *
 * The MKCore object stores MeshKit-wide data like the geometry/mesh/relations instances,
 * and is the object through which other MeshKit class objects are accessed.  Since it is
 * a child class of MeshOp, the MKCore instance can also be a node in the MeshOp graph.  By
 * convention, the MKCore instance serves as the root of the this (directed) graph.
 *
 * If the MKCore constructor is called with no arguments, then instances of the geometry, mesh,
 * and relations databases are created inside the constructor; these instances are then deleted
 * in the MKCore destructor.  If this is not the desired behavior, either pass in non-NULL
 * instances to the MKCore constructor, or re-set the iCreatedIgeom, iCreatedMoab,
 * and/or iCreatedIrel flags in this class.
 */

class MKCore : public MeshOp
{
public:

    /** \brief Constructor
     * \param igeom iGeom instance
     * \param moab MOAB instance
     * \param irel iRel instance
     * \param construct_missing_ifaces If true, constructs the interfaces whose handles are passed in NULL
     */
  MKCore(iGeom *igeom = NULL, 
         moab::Interface *moab = NULL, 
         iRel *irel = NULL,
         bool construct_missing_ifaces = true) throw(Error);
  
    //! destructor
  ~MKCore() throw(Error);

    //! initialize, creating missing geom/mesh/rel interfaces if requested
  void init(bool construct_missing_ifaces) throw(Error);
  
    /** \brief Load a geometry model from a file
     * \param filename The file to load
     * \param options File options to be passed to the load function
     */
  void load_geometry(const char *filename, const char *options = NULL) throw(Error);

    /** \brief Load a mesh model from a file
     * \param filename The file to load
     * \param options File options to be passed to the load function
     */
  void load_mesh(const char *filename, const char *options = NULL) throw(Error);

    /** \brief Save a geometry model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     */
  void save_geometry(const char *filename, const char *options = NULL) throw(Error);

    /** \brief Save a mesh model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     */
  void save_mesh(const char *filename, const char *options = NULL) throw(Error);

    /** \brief populate mesh/relations data for geometric entities
     */
  void populate_mesh() throw(Error);

    /** \brief Get model entities of a given dimension
     * \param dim Dimension of entities to get
     * \param model_ents The list these entities get appended to
     */
  void get_entities_by_dimension(int dim, MEVector &model_ents) throw(Error);

    /** \brief Get all model entities
     * \param model_ents The list these entities get appended to
     */
  void get_entities_by_handle(MEVector &model_ents) throw(Error);

    //! Geometry api instance
  iGeom *iGeomInstance;
  
    //! Mesh api instance
  moab::Interface *moabInstance;
  
    //! IREL api instance
  iRel *iRelInstance;
  
    //! Tag used to associate geometry entities with model entities
  iGeom::TagHandle igeomModelTag;
  
    //! Tag used to associate mesh entities with model entities
  moab::Tag imeshModelTag;

private:
  bool iCreatedIgeom, iCreatedMoab, iCreatedIrel;

    //! Model entity sets, in array by topological dimension
  moab::Range modelEnts[4];
  
};

}

#endif

  
