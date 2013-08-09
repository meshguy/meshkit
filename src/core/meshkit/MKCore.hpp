#ifndef MESHKIT_MKCORE_HPP
#define MESHKIT_MKCORE_HPP

/** \file MKCore.hpp
 */
#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MKGraph.hpp"
#include "meshkit/iGeom.hpp"
#include "moab/Interface.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/iRel.hpp"
#include "lemon/list_graph.h"
#include <vector>

namespace MeshKit {

      //! Forward declare since we store a vector of these
class SizingFunction;

      //! Single instance of VertexMesher
class VertexMesher;

      //! Single instance of EBMesher
class EBMesher;

      //! Returned by some MKCore functions
class MeshOpProxy;

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

class MKCore : public MKGraph
{
public:

    /** \name Constructor/destructor */

    /**@{*/

    /** \brief Constructor
     * \param igeom iGeom instance
     * \param moab MOAB instance
     * \param imesh iMesh instance; if non-NULL, should use/point to moab parameter
     * \param irel iRel instance
     * \param construct_missing_ifaces If true, constructs the interfaces whose handles are passed in NULL
     */
  MKCore(iGeom *igeom = NULL, 
         moab::Interface *moab = NULL, 
         iMesh *imesh = NULL, 
         iRel *irel = NULL,
         bool construct_missing_ifaces = true);
  
    //! destructor
  ~MKCore();

    /**@}*/

    /** \brief Register a new MeshOp factory
     * \param proxy class-specific (as opposed to instance-specific) polymorphic methods
     */
  static void register_meshop(MeshOpProxy* proxy);
  
    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return MeshOpProxy for the corresponding MeshOp type
     */
  static MeshOpProxy* meshop_proxy(const char *op_name);
  
    /** \brief Return the MeshOp type at the sepcified index
     * \param index Operation index requested
     * \return MeshOpProxy for the corresponding MeshOp type
     */
  static MeshOpProxy* meshop_proxy(unsigned index);
  
    /**\brief Number of registered MeshOps
     *
     *\return One more than largest possible MeshOp index
     */
  static unsigned num_meshops();
  
    /** \brief Return the MeshOp index given the name
     * \param op_name Operation name requested
     * \return OpInfo index for the corresponding MeshOp type
     */
  static unsigned int meshop_index(const char *op_name);

    /** \name Non-static MeshOp construction and query */

    /**@{*/

    /** \brief Make the specified MeshOp name the default for the given dimension(s)
     * 
     * If the specified MeshOp cannot produce entities of the specified dimension, an error is
     * thrown with type MK_BAD_INPUT.
     * \param op_name MeshOp name being set
     * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
     */
  void set_default_meshop(const char *op_name, unsigned short dims);
  
    /** \brief Make the specified MeshOp the default for the given dimension(s)
     * 
     * If the specified MeshOp cannot produce entities of the specified dimension, an error is
     * thrown with type MK_BAD_INPUT.
     * \param op_index MeshOp index being set
     * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
     */
  void set_default_meshop(unsigned op_index, unsigned short dims);
  
    /** \brief Make the specified MeshOp the default for the given dimension(s)
     * 
     * If the specified MeshOp cannot produce entities of the specified dimension, an error is
     * thrown with type MK_BAD_INPUT.
     * \param mesh_op MeshOp being set
     * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
     */
  void set_default_meshop(MeshOpProxy* mesh_op, unsigned short dims);
  
    /**\brief Get default MeshOp for dimension 
     *
     * Get the default meshop for the specified dimension of entity.
     * Throws excpetion if no there are no registered MeshOps that
     * can mesh entities of the passed dimension.
     *
     *\param dimension Dimension of entity to be meshed.
     *\return pointer to proxy for MeshOp 
     */
  MeshOpProxy* get_default_meshop( unsigned dimension );
  
    /** \brief Return MeshOp types that can generate or operate on the specified mesh entity type
     * \param tp Entity type requested
     * \param ops MeshOp types returned
     */
  void meshop_by_mesh_type(moab::EntityType tp, std::vector<MeshOpProxy*> &ops);
    
    /** \brief Return MeshOp types that can operate on geometric entities of specified dimension
     * \param dim Entity dimension requested
     * \param ops MeshOp types returned
     */
  void meshop_by_dimension(int dim, std::vector<MeshOpProxy*> &ops);
    
    /** \brief Return MeshOp types that can mesh the specified ModelEnt
     * \param ent ModelEnt* requested
     * \param ops MeshOp types returned
     */
  void meshop_by_modelent(ModelEnt * const ent, std::vector<MeshOpProxy*> &ops);
    
    /** \brief Construct a new MeshOp of the specified OpInfo type
     * \param op_info OpInfo for the MeshOp being requested
     * \param me_vec MEntVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(MeshOpProxy* op_info, const MEntVector &me_vec = MEntVector());
    
    /** \brief Construct a new MeshOp of the specified name
     * \param op_name MeshOp name being requested
     * \param me_vec MEntVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(std::string op_name, const MEntVector &me_vec = MEntVector());
    
    /** \name MeshOp construction and query */

    /**@{*/

    /** \brief Construct the default type of MeshOp for the specified dimension
     * \param dim Dimension requested
     * \param me_vec MEntVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(unsigned int dim, const MEntVector &me_vec = MEntVector());
    
    /**@}*/

    /** \name Loading/saving models */

    /**@{*/

    /** \brief Load a geometry model from a file, and populate mesh entity sets
     * \param geom_filename The file from which to load geometry
     * \param mesh_filename The file from which to load mesh
     * \param geom_options File options to be passed to the load_geometry function
     * \param mesh_options File options to be passed to the load_mesh function
     * \param geom_index Index of geometry instance to use
     * \param mesh_index Index of mesh instance to use
     * \param irel_index Index of irel instance to use
     * \param relate_too If true, infers all geometry-mesh relations after mesh load
     * \param populate_too If true, calls populate_modelents after load and relate
     */
  void load_geometry_mesh(const char *geom_filename, 
                          const char *mesh_filename, 
                          const char *geom_options = NULL, 
                          const char *mesh_options = NULL, 
                          int geom_index = 0,
                          int mesh_index = 0,
                          int irel_index = 0,
                          bool relate_too = true,
                          bool populate_too = true);

    /** \brief Load a geometry model from a file, and populate mesh entity sets
     * \param filename The file to load
     * \param options File options to be passed to the load function
     * \param geom_index Index of geometry instance to use
     * \param mesh_index Index of mesh instance to use
     * \param irel_index Index of irel instance to use
     * \param relate_too If true, infers all geometry-mesh relations after load
     * \param populate_too If true, calls populate_modelents after load and relate
     */
  void load_geometry(const char *filename, const char *options = NULL, 
                     int geom_index = 0,
                     int mesh_index = 0,
                     int irel_index = 0,
                     bool relate_too = false,
                     bool populate_too = true);

    /** \brief Load a mesh model from a file
     * \param filename The file to load
     * \param options File options to be passed to the load function
     * \param geom_index Index of geometry instance to use
     * \param mesh_index Index of mesh instance to use
     * \param irel_index Index of irel instance to use
     * \param relate_too If true, infers all geometry-mesh relations after load
     * \param populate_too If true, calls populate_modelents after load and relate
     */
  void load_mesh(const char *filename, const char *options = NULL,
                 int geom_index = 0,
                 int mesh_index = 0,
                 int irel_index = 0,
                 bool relate_too = false,
                 bool populate_too = true);

    /** \brief Save a geometry model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     * \param index Index of geometry instance to use
     */
  void save_geometry(const char *filename, const char *options = NULL,
                     int index = 0);

    /** \brief Save a mesh model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     * \param index Index of mesh instance to use
     */
  void save_mesh(const char *filename, const char *options = NULL,
                 int index = 0);

  /** \brief Save part of a model mesh to a file
      * \param filename The file to save
      * \param ments vector of model entities to be exported
      * \param options File options to be passed to the save function
      * \param index Index of mesh instance to use
      */
  void save_mesh_from_model_ents(const char *filename, MEntVector & ments, const char *options = NULL,
                  int index = 0);

    /** \brief Populate model entities for geometry and mesh
     * Pass -1 for geom, mesh, or irel index to disable checking for entities in those
     * interfaces (geom, mesh) or skip the check for relations (irel)
     * \param geom_index Index of geometry instance to use
     * \param mesh_index Index of mesh instance to use
     * \param irel_index Index of irel instance to use
     * \param use_irel flag to create relations or not
     */
  void populate_model_ents(int geom_index = 0, 
                           int mesh_index = 0, 
                           int irel_index = 0, bool only_geom = false);

  /**
   *  \ brief  create model entities based on a set root for mesh based geometry
   *  \param modelRootSet root model for a complete topology model (for GeomTopoTool)
   *  \param geometry flag true for FBiGeom geometry, false for simple mesh ModelEnts
   */
  void create_mbg_model_entities(moab::EntityHandle modelRootSet, bool geometry);

    /**@}*/

    /** \name Entity access */

    /**@{*/

    /** \brief Get model entities of a given dimension
     * \param dim Dimension of entities to get
     * \param model_ents The list these entities get appended to
     * \param igindx the geometry instance index of the retrieved model ents
     *   (default -1 means no check on igeom instance index)
     */
  void get_entities_by_dimension(int dim, MEntVector &model_ents,
      int igindx = -1);

    /** \brief Get all model entities
     * \param model_ents The list these entities get appended to
     */
  void get_entities_by_handle(MEntVector &model_ents);

    /**@}*/

    /** \name Interface instances */

    /**@{*/

    /** \brief Return the iGeom instance pointer
     * \param index Index of desired iGeom, default to first
     */
  iGeom *igeom_instance(unsigned index = 0);
  
    /** \brief Return the index of iGeom instance added
     * \param igeom iGeom interface to be added
     */
  unsigned int add_igeom_instance(iGeom * igeom);

    /** \brief Return the MOAB instance pointer
     * \param index Index of desired moab instance, default to first
     */
  moab::Interface *moab_instance(unsigned index = 0);
  
    /** \brief Return the iMesh instance pointer
     * \param index Index of desired iMesh, default to first
     */
  iMesh *imesh_instance(unsigned index = 0);
  
    /** \brief Return the iRel instance pointer
     */
  iRel *irel_instance(unsigned index = 0);

    /** \brief Return the index of irel_pair added
       * \param pair iRel::PairHandle pointer to be added
      */
  inline unsigned int add_irel_pair(iRel::PairHandle * pair);

    /** \brief Return the iRel pair handle used to relate geometry/mesh entities
     */
  iRel::PairHandle *irel_pair(unsigned index = 0);

  /** \brief Return the iRel pair handle used to relate geometry sets to mesh entity sets
     */
  iRel::PairHandle *group_set_pair(unsigned index = 0);

    /** \brief Return the (iGeom) tag used to relate geometry entities to ModelEnts
     */
  iGeom::TagHandle igeom_model_tag(unsigned index = 0);

    /** \brief Return the (MOAB) tag used to relate mesh entities to ModelEnts
     */
  moab::Tag moab_model_tag(unsigned index = 0);

    /** \brief Return the (MOAB) geometry dimension tag
     */
  moab::Tag moab_geom_dim_tag(unsigned index = 0);

    /** \brief Return the (MOAB) global id tag
     */
  moab::Tag moab_global_id_tag(unsigned index = 0);

    /**@}*/

    /** \name Singleton MeshOps */

    /**@{*/

    /** \brief Get the (single) VertexMesher instance
     * \return VertexMesher for this MKCore
     */
  VertexMesher *vertex_mesher() const;
  
    /** \brief Set the (single) VertexMesher instance
     * \param vm VertexMesher for this MKCore
     */
  void vertex_mesher(VertexMesher *vm);

    /** \brief Get the (single) EBMesher instance
     * \return EBMesher for this MKCore
     */
  EBMesher *eb_mesher() const;
  
    /** \brief Set the (single) EBMesher instance
     * \param ebm EBMesher for this MKCore
     */
  void eb_mesher(EBMesher *ebm);

    /**@}*/

    /** \name Sizing functions */

    /**@{*/

    /** \brief Get sizing function by index
     * If the requested index is outside the range of SizingFunction's currently registered,
     * throws an Error.
     * \param index Index of sizing function requested
     * \return SizingFunction* to requested sizing function, NULL of no SizingFunction with that index
     */
  SizingFunction *sizing_function(int index);
  
    /** \brief Get sizing function by size
     * If there is no sizing function with that size and create_if_missing is true, one is constructed and registered with MKCore.
     * \param size Requested size
     * \param create_if_missing If true and no sizing function exists with the specified size, one is created.
     * \return SizingFunction* to requested sizing function, NULL if no SizingFunction with that size
     */
  SizingFunction *sizing_function(double size, bool create_if_missing);
  
    /** \brief Add sizing function to those managed by MeshKit
     *
     * The argument to this function is a SizingFunction*; once added, it is MKCore's
     * responsibility to delete this SizingFunction.  Applications can tell MKCore to delete
     * a given SizingFunction (e.g. if it requires lots of memory) by calling delete_sizing_function.
     * \param sf SizingFunction* to be added
     * \return Index of sizing function in MKCore's list of SizingFunction's
     */
  int add_sizing_function(SizingFunction *sf);

    /** \brief Remove and, optionally, delete sizing function
     *
     * This function removes the referenced sizing function from MKCore's list (setting the
     * corresponding SizingFunction* to NULL, to keep others at the same index position).  
     * Throws an Error if requested sizing function is NULL.
     * \param index Index of SizingFunction to be removed
     * \param delete_too If true, deletes the SizingFunction object too
     */
  void remove_sizing_function(int index, bool delete_too = true);
  
#if 1
  /** \brief interpret mesh-based geometry
       *
       * This function interprets the mesh from a moab db as mesh based geometry
       * it will create a new FBiGeom instance and populate the model ents
       * accordingly, with the geometry implied from moab sets and tags
       * it returns the index of FBiGeom in iGeomInstances vector
       */

  int initialize_mesh_based_geometry(iBase_EntitySetHandle modelSet);

  /** \brief remove mesh-based geometry instance
     *
     * This function removes the FBiGeom instance and clears some memory
     * (will delete the smooth tags, among other things), but it should preserve
     * the geom topo model intact
     */
  void remove_mesh_based_geometry(int index);
#endif
  /** \brief retrieve number of iGeom instances
     *
     * This function returns number of igeom instances, to be used
     * by some mesh operations.
     */
  int number_of_igeom_instances()
  {
    return (int)iGeomInstances.size();
  }
  /** \brief delete model entities
       *
       * This function clears memory occupied by model entities, and
       * empties the vectors.
       * more careful with the lifecycles.
       */
  void delete_model_entities();
  /** \brief delete all
   *
   * This function clears memory occupied by model entities,
   * deletes graph,
   * deletes mesh and geometry
   */
  void delete_all();
    /**@}*/

private:
    //! initialize, creating missing geom/mesh/rel interfaces if requested
  void init(bool construct_missing_ifaces);

    //! Initialize the opsByDim array
  void init_opsbydim();

    //! Geometry api instance
  std::vector<iGeom *> iGeomInstances;
  
    //! MOAB instance
  std::vector<moab::Interface *> moabInstances;
  
    //! iMesh api instance, for use in iRel
  std::vector<iMesh *> iMeshInstances;
  
    //! IREL api instance
  std::vector<iRel *> iRelInstances;

    //! iRel pair handle used to relate geometry/mesh entities
  std::vector<iRel::PairHandle *> iRelPairs;
  
    //! iRel pair handle used to relate geometry groups to mesh entity sets
  std::vector<iRel::PairHandle *> groupSetPairs;
  
    //! Tag used to associate geometry entities with model entities
  std::vector<iGeom::TagHandle> iGeomModelTags;
  
    //! Tag used to associate mesh entities with model entities
  std::vector<moab::Tag> moabModelTags;

    //! Tag used to associate existing mesh entities with model entities
  std::vector<moab::Tag> moabGeomDimTags;

    //! Tag used to associate existing mesh entities with model entities
  std::vector<moab::Tag> moabIDTags;

    //! If true, the corresponding interfaces will get deleted from the destructor
  std::vector<bool> iCreatedIgeoms, iCreatedMoabs, iCreatedImeshs, iCreatedIrels;

    //! Model entities, in array by topological dimension
  MEntVector modelEnts[5];
  
    //! (Single) VertexMesher scheme for this MKCore
  VertexMesher *vertexMesher;

      //! (Single) EBMesher scheme for this MKCore
  EBMesher *ebMesher;
  
    //! Default algorithms for each dimension
  MeshOpProxy* defaultMeshOps[4];
  
    //! SizingFunction vector
  std::vector<SizingFunction*> sizingFunctions;
};

inline iGeom *MKCore::igeom_instance(unsigned index) 
{
  if (iGeomInstances.size() <= index)
    throw Error(MK_BAD_INPUT, "No instance of that index.");
  
  return iGeomInstances[index];
}

inline moab::Interface *MKCore::moab_instance(unsigned index)
{
  if (moabInstances.size() <= index)
    throw Error(MK_BAD_INPUT, "No instance of that index.");
  
  return moabInstances[index];
}

inline iMesh *MKCore::imesh_instance(unsigned index) 
{
  if (iMeshInstances.size() <= index)
    throw Error(MK_BAD_INPUT, "No instance of that index.");
  
  return iMeshInstances[index];
}

inline iRel *MKCore::irel_instance(unsigned index)
{
  if (iRelInstances.size() <= index)
    throw Error(MK_BAD_INPUT, "No instance of that index.");

  return iRelInstances[index];
}

inline iRel::PairHandle *MKCore::irel_pair(unsigned index)
{
  if (iRelPairs.size() <= index)
    throw Error(MK_BAD_INPUT, "No pair of that index.");

  return iRelPairs[index];
}

inline unsigned int MKCore::add_irel_pair(iRel::PairHandle * pair)
{
  iRelPairs.push_back(pair);
  return iRelPairs.size()-1;
}

inline iRel::PairHandle *MKCore::group_set_pair(unsigned index)
{
  if (groupSetPairs.size() <= index)
    throw Error(MK_BAD_INPUT, "No pair of that index.");

  return groupSetPairs[index];
}

inline iGeom::TagHandle MKCore::igeom_model_tag(unsigned index)
{
  if (iGeomModelTags.size() <= index)
    throw Error(MK_BAD_INPUT, "No tag of that index.");

  return iGeomModelTags[index];
}

inline moab::Tag MKCore::moab_model_tag(unsigned index)
{
  if (moabModelTags.size() <= index)
    throw Error(MK_BAD_INPUT, "No tag of that index.");

  return moabModelTags[index];
}

inline moab::Tag MKCore::moab_geom_dim_tag(unsigned index)
{
  if (moabGeomDimTags.size() <= index)
    throw Error(MK_BAD_INPUT, "No tag of that index.");

  return moabGeomDimTags[index];
}

inline moab::Tag MKCore::moab_global_id_tag(unsigned index)
{
  if (moabIDTags.size() <= index)
    throw Error(MK_BAD_INPUT, "No tag of that index.");

  return moabIDTags[index];
}

inline void MKCore::vertex_mesher(VertexMesher *vm) 
{
  vertexMesher = vm;
}

inline VertexMesher *MKCore::vertex_mesher() const
{
  return vertexMesher;
}

inline EBMesher *MKCore::eb_mesher() const 
{
  return ebMesher;
}

inline void MKCore::eb_mesher(EBMesher *ebm) 
{
  ebMesher = ebm;
}

inline SizingFunction *MKCore::sizing_function(int index) 
{
    // don't check for NULL here 'cuz sometimes we just want to know there isn't one
    // with that index
  if (index >= (int)sizingFunctions.size())
    throw Error(MK_BAD_INPUT, "Sizing function index outside range of valid indices.");
  else if (index == -1)
    return NULL;

  return sizingFunctions[index];
}
  
} // namespace MeshKit

#endif

  
