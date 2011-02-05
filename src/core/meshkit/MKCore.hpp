#ifndef MKCORE_HPP
#define MKCORE_HPP

/** \file MKCore.hpp
 */
#include "meshkit/Types.h"
#include "meshkit/Error.hpp"
#include "meshkit/MKGraph.hpp"
#include "iGeom.hh"
#include "moab/Interface.hpp"
#include "iMesh.hh"
#include "iRel.hh"
#include "lemon/list_graph.h"
#include <vector>

// outside the namespace 'cuz it just is
class MBiMesh;

namespace MeshKit {

      //! Forward declare since we store a vector of these
class SizingFunction;

      //! MKCore keeps a single instance of this
class VertexMesher;

//! SCDMesh keeps a single instance of this
class SCDMesh;

class EBMesher;

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

    /** \brief Constructor
     * \param igeom iGeom instance
     * \param moab MOAB instance
     * \param mbi MBiMesh instance; if non-NULL, should use/point to moab parameter
     * \param irel iRel instance
     * \param construct_missing_ifaces If true, constructs the interfaces whose handles are passed in NULL
     */
  MKCore(iGeom *igeom = NULL, 
         moab::Interface *moab = NULL, 
         MBiMesh *mbi = NULL, 
         iRel *irel = NULL,
         bool construct_missing_ifaces = true);
  
    //! destructor
  ~MKCore();

    //! initialize, creating missing geom/mesh/rel interfaces if requested
  void init(bool construct_missing_ifaces);

    //! Function for constructing MeshOp instances, provided by registrants
  typedef MeshOp* (*meshop_factory_t)(MKCore *, const MEntVector &vec);

    //! Function returning whether a given MeshOp can mesh the specified ModelEnt
  typedef bool (*meshop_canmesh_t)(ModelEnt *);

    /** \struct OpInfo MeshOpFactory.hpp "meshkit/MeshOpFactory.hpp"
     * \brief Struct used to store information about each MeshOp type registered with MeshOpFactory
     */
  struct OpInfo
  {
    std::string opName;
    unsigned short opIndex;
    std::vector<iBase_EntityType> modelEntTypes;
    std::vector<moab::EntityType> meshEntTypes;
    meshop_factory_t opFactory;
    meshop_canmesh_t opCanMesh;
  };

    /** \brief A map from MeshOp name to index in the OpInfos vector
     */
  typedef std::map<std::string, unsigned short> OpNameMap;
  
  class MeshOpFactory 
  {
   public:
    std::vector<OpInfo> registeredOps;
    OpNameMap opNameMap;
  };

  static MeshOpFactory *op_factory();

    //! Initialize the opsByDim array
  void init_opsbydim();

    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tp The iBase_EntityType (or dimension) operated on by this MeshOp type
     * \param tp The MOAB entity type produced by this MeshOp
     * \param meshop The (static) factory function producing instances of this MeshOp type
     * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
     * \return Returns true if successful, false otherwise
     */
  static bool register_meshop(const char *op_name, 
                              iBase_EntityType model_tp, moab::EntityType tp, 
                              meshop_factory_t meshop, meshop_canmesh_t canmesh = NULL);
  
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tp The iBase_EntityType (or dimension) operated on by this MeshOp type
     * \param tp The iMesh entity topology produced by this MeshOp
     * \param meshop The (static) factory function producing instances of this MeshOp type
     * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
     * \return Returns true if successful, false otherwise
     */
  static bool register_meshop(const char *op_name, 
                              iBase_EntityType model_tp, iMesh::EntityTopology tp, 
                              meshop_factory_t meshop, meshop_canmesh_t canmesh = NULL);
  
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tps The iBase_EntityType's (or dimensions) operated on by this MeshOp type
     * \param num_mtps Number of model entity types in model_tps
     * \param tps The MOAB entity types operated on by this MeshOp
     * \param num_tps Number of entity types in tps
     * \param meshop The (static) factory function producing instances of this MeshOp type
     * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
     * \return Returns true if successful, false otherwise
     */
  static bool register_meshop(const char *op_name, 
                              iBase_EntityType *model_tps, int num_mtps,
                              moab::EntityType *tps, int num_tps,
                              meshop_factory_t meshop, meshop_canmesh_t canmesh = NULL);
  
    /** \brief Register a new MeshOp factory
     * \param op_name The name by which this type of MeshOp can be requested
     * \param model_tps The iBase_EntityType's (or dimensions) operated on by this MeshOp type
     * \param num_mtps Number of model entity types in model_tps
     * \param tps The iMesh entity types operated on by this MeshOp
     * \param num_tps Number of entity types in tps
     * \param meshop The (static) factory function producing instances of this MeshOp type
     * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
     * \return Returns true if successful, false otherwise
     */
  static bool register_meshop(const char *op_name, 
                              iBase_EntityType *model_tps, int num_mtps,
                              iMesh::EntityTopology *tps, int num_tps,
                              meshop_factory_t meshop, meshop_canmesh_t canmesh);
  
    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return OpInfo for the corresponding MeshOp type
     */
  static OpInfo meshop_info(const char *op_name);
  
    /** \brief Return the MeshOp index given the name
     * \param op_name Operation name requested
     * \return OpInfo index for the corresponding MeshOp type
     */
  static unsigned int meshop_index(const char *op_name);

    /** \brief Make the specified MeshOp name the default for the given dimension(s)
     * 
     * If the specified MeshOp cannot produce entities of the specified dimension, an error is
     * thrown with type MK_BAD_INPUT.
     * \param op_name MeshOp name being set
     * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
     */
  void set_default_meshop(const char *op_name, unsigned short dims);
  
    /** \brief Make the specified MeshOp name the default for the given dimension(s)
     * 
     * If the specified MeshOp cannot produce entities of the specified dimension, an error is
     * thrown with type MK_BAD_INPUT.
     * \param op_index MeshOp index being set
     * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
     */
  void set_default_meshop(unsigned short op_index, unsigned short dims);
  
    /** \brief Return MeshOp types that can operate on the specified entity type
     * \param tp Entity type requested
     * \param ops MeshOp types returned
     */
  void meshop_by_mesh_type(moab::EntityType tp, std::vector<OpInfo> &ops);
    
    /** \brief Return MeshOp types that can operate on mesh of specified dimension
     * \param dim Entity dimension requested
     * \param ops MeshOp types returned
     */
  void meshop_by_dimension(int dim, std::vector<OpInfo> &ops);
    
    /** \brief Return MeshOp types that can mesh the specified ModelEnt
     * \param ent ModelEnt* requested
     * \param ops MeshOp types returned
     */
  void meshop_by_modelent(ModelEnt * const ent, std::vector<OpInfo> &ops);
    
    /** \brief Construct a new MeshOp of the specified OpInfo type
     * \param op_info OpInfo for the MeshOp being requested
     * \param me_vec MEntVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(OpInfo &op_info, const MEntVector &me_vec = MEntVector());
    
    /** \brief Construct a new MeshOp of the specified name
     * \param op_name MeshOp name being requested
     * \param me_vec MEntVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(std::string op_name, const MEntVector &me_vec = MEntVector());
    
    /** \brief Construct the default type of MeshOp for the specified dimension
     * \param dim Dimension requested
     * \param me_vec MEntVector of entities the operation applies to
     * \return Pointer to new MeshOp constructed
     */
  MeshOp *construct_meshop(unsigned int dim, const MEntVector &me_vec = MEntVector());
    
    /** \brief Load a geometry model from a file, and populate mesh entity sets
     * \param filename The file to load
     * \param options File options to be passed to the load function
     * \param populate_too If true, calls populate_mesh after load
     */
  void load_geometry(const char *filename, const char *options = NULL, 
                     bool populate_too = true);

    /** \brief Load a mesh model from a file
     * \param filename The file to load
     * \param options File options to be passed to the load function
     */
  void load_mesh(const char *filename, const char *options = NULL);

    /** \brief Save a geometry model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     */
  void save_geometry(const char *filename, const char *options = NULL);

    /** \brief Save a mesh model to a file
     * \param filename The file to save
     * \param options File options to be passed to the save function
     */
  void save_mesh(const char *filename, const char *options = NULL);

    /** \brief Populate mesh/relations data for geometric entities
     */
  void populate_mesh();

    /** \brief Get model entities of a given dimension
     * \param dim Dimension of entities to get
     * \param model_ents The list these entities get appended to
     */
  void get_entities_by_dimension(int dim, MEntVector &model_ents);

    /** \brief Get all model entities
     * \param model_ents The list these entities get appended to
     */
  void get_entities_by_handle(MEntVector &model_ents);

    /** \brief Return the iGeom instance pointer
     */
  iGeom *igeom_instance();
  
    /** \brief Return the MOAB instance pointer
     */
  moab::Interface *moab_instance();
  
    /** \brief Return the iMesh instance pointer
     */
  MBiMesh *mb_imesh();
  
    /** \brief Return the iRel instance pointer
     */
  iRel *irel_instance();

    /** \brief Return the iRel pair handle used to relate geometry/mesh entities
     */
  iRel::PairHandle *irel_pair();

    /** \brief Return the iRel pair handle used to relate geometry sets to mesh entity sets
     */
  iRel::PairHandle *group_set_pair();

    /** \brief Return the (iGeom) tag used to relate geometry entities to ModelEnts
     */
  iGeom::TagHandle igeom_model_tag();

    /** \brief Return the (MOAB) tag used to relate mesh entities to ModelEnts
     */
  moab::Tag moab_model_tag();

    /** \brief Return the (MOAB) geometry dimension tag
     */
  moab::Tag moab_geom_dim_tag();

    /** \brief Return the (MOAB) global id tag
     */
  moab::Tag moab_global_id_tag();

    /** \brief Get the (single) VertexMesher instance
     * \return VertexMesher for this MKCore
     */
  VertexMesher *vertex_mesher() const;
  
    /** \brief Set the (single) VertexMesher instance
     * \param vm VertexMesher for this MKCore
     */
  void vertex_mesher(VertexMesher *vm);

   /** \brief Get the (single) SCDMesh instance                           
    * \return SCDMesh for this MKCore                                     
    */
  SCDMesh *scd_mesh() const;

   /** \brief Set the (single) SCDMesh instance                           
    * \param vm SCDMesh for this MKCore                                   
    */
  void scd_mesh(SCDMesh *scd);

    /** \brief Get the (single) EBMesher instance
     * \return EBMesher for this MKCore
     */
  EBMesher *eb_mesher() const;
  
    /** \brief Set the (single) EBMesher instance
     * \param vm EBMesher for this MKCore
     */
  void eb_mesher(EBMesher *ebm);

    /** \brief Get sizing function by index
     * If the requested index is outside the range of SizingFunction's currently registered,
     * throws an Error.
     * \param index Index of sizing function requested
     * \return SizingFunction* to requested sizing function, NULL of no SizingFunction with that index
     */
  SizingFunction *sizing_function(int index);
  
    /** \brief Add sizing function to those managed by MeshKit
     *
     * The argument to this function is a SizingFunction*; once added, it is MKCore's
     * responsibility to delete this SizingFunction.  Applications can tell MKCore to delete
     * a given SizingFunction (e.g. if it requires lots of memory) by calling delete_sizing_function.
     * \param sf SizingFunction* to be added
     * \return Index of sizing function in MKCore's list of SizingFunction's
     */
  int add_sizing_function(SizingFunction *sf);

    /** \brief Delete sizing function
     *
     * This function removes the referenced sizing function from MKCore's list (setting the
     * corresponding SizingFunction* to NULL, to keep others at the same index position).  
     * Throws an Error if requested sizing function is NULL.
     * \param index Index of SizingFunction to be removed
     */
  void remove_sizing_function(int index);
  
private:
    //! Geometry api instance
  iGeom *iGeomInstance;
  
    //! MOAB instance
  moab::Interface *moabInstance;
  
    //! iMesh api instance, for use in iRel
  MBiMesh *mbImesh;
  
    //! IREL api instance
  iRel *iRelInstance;

    //! iRel pair handle used to relate geometry/mesh entities
  iRel::PairHandle *iRelPair;
  
    //! iRel pair handle used to relate geometry groups to mesh entity sets
  iRel::PairHandle *groupSetPair;
  
    //! Tag used to associate geometry entities with model entities
  iGeom::TagHandle iGeomModelTag;
  
    //! Tag used to associate mesh entities with model entities
  moab::Tag moabModelTag;

    //! Tag used to associate existing mesh entities with model entities
  moab::Tag moabGeomDimTag;

    //! Tag used to associate existing mesh entities with model entities
  moab::Tag moabIDTag;

    //! If true, the corresponding interfaces will get deleted from the destructor
  bool iCreatedIgeom, iCreatedMoab, iCreatedMbimesh, iCreatedIrel;

    //! Model entities, in array by topological dimension
  MEntVector modelEnts[5];

  static MeshOpFactory *opFactory;

    //! OpInfo's by dimension they can produce; first for a given dimension is the default
  std::vector<unsigned short> *opsByDim;

    //! Number of ops checked, cached to avoid redundant checking of registered ops
  unsigned int numOpsByDim;
  
    //! (Single) VertexMesher scheme for this MKCore
  VertexMesher *vertexMesher;

    //! (Single) SCDMesh scheme for this MKCore
  SCDMesh *scdMesh;

      //! (Single) EBMesher scheme for this MKCore
  EBMesher *ebMesher;
  
    //! SizingFunction vector
  std::vector<SizingFunction*> sizingFunctions;
};

inline iGeom *MKCore::igeom_instance() 
{
  return iGeomInstance;
}

inline moab::Interface *MKCore::moab_instance()
{
  return moabInstance;
}

inline MBiMesh *MKCore::mb_imesh()
{
  return mbImesh;
}

inline iRel *MKCore::irel_instance()
{
  return iRelInstance;
}

inline iRel::PairHandle *MKCore::irel_pair()
{
  return iRelPair;
}

inline iRel::PairHandle *MKCore::group_set_pair()
{
  return groupSetPair;
}

inline iGeom::TagHandle MKCore::igeom_model_tag()
{
  return iGeomModelTag;
}

inline moab::Tag MKCore::moab_model_tag()
{
  return moabModelTag;
}

inline moab::Tag MKCore::moab_geom_dim_tag()
{
  return moabGeomDimTag;
}

inline moab::Tag MKCore::moab_global_id_tag()
{
  return moabIDTag;
}


/** \brief Register a new MeshOp factory
 * \param op_name The name by which this type of MeshOp can be requested
 * \param tp The MOAB entity type operated on by this MeshOp
 * \param meshop The (static) factory function producing instances of this MeshOp type
 * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
 */
inline bool MKCore::register_meshop(const char *op_name, 
                                    iBase_EntityType model_tp, iMesh::EntityTopology tp, 
                                    meshop_factory_t meshop, meshop_canmesh_t canmesh) 
{
  return register_meshop(op_name, model_tp, iMesh::mb_topology_table[tp], meshop, canmesh);
}
  
/** \brief Register a new MeshOp factory
 * \param op_name The name by which this type of MeshOp can be requested
 * \param tps The MOAB entity types operated on by this MeshOp
 * \param num_tps Number of entity types in tps
 * \param meshop The (static) factory function producing instances of this MeshOp type
 * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
 */
inline bool MKCore::register_meshop(const char *op_name, 
                                    iBase_EntityType *model_tps, int num_mtps,
                                    iMesh::EntityTopology *tps, int num_tps,
                                    meshop_factory_t meshop, meshop_canmesh_t canmesh)
{
  std::vector<moab::EntityType> tmp_tps(num_tps);
  for (int i = 0; i < num_tps; i++) tmp_tps[i] = iMesh::mb_topology_table[tps[i]];
  return register_meshop(op_name, model_tps, num_mtps, &tmp_tps[0], num_tps, meshop, canmesh);
}

inline VertexMesher *MKCore::vertex_mesher() const 
{
  return vertexMesher;
}

inline void MKCore::vertex_mesher(VertexMesher *vm) 
{
  vertexMesher = vm;
}

inline SCDMesh *MKCore::scd_mesh() const
{
  return scdMesh;
}
  
inline void MKCore::scd_mesh(SCDMesh *scd)
{
  scdMesh = scd;
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

  
