#include "meshkit/MeshOpFactory.hpp"
#include "meshkit/MKCore.hpp"
#include "moab/CN.hpp"
#include "moab/Interface.hpp"

namespace MeshKit 
{
    
MeshOpFactory *MeshOpFactory::instance_ = NULL;

/** \brief Register a new MeshOp factory
 * \param op_name The name by which this type of MeshOp can be requested
 * \param tp The MOAB entity type operated on by this MeshOp
 * \param meshop The (static) factory function producing instances of this MeshOp type
 * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
 */
bool MeshOpFactory::register_meshop(const char *op_name, moab::EntityType tp, 
                                    meshop_factory_t meshop, meshop_canmesh_t canmesh) 
{
  OpInfoMap::iterator oit = registeredOps.find(std::string(op_name));
  if (oit != registeredOps.end()) 
    throw Error(MK_ALREADY_DEFINED, "A MeshOp with this name has already been registered.");
  
  OpInfo &oi = registeredOps[std::string(op_name)];
  oi.opName = std::string(op_name);
  oi.opEntTypes.push_back(tp);
  oi.opFactory = meshop;
  oi.opCanMesh = canmesh;

  return true;
}
  
/** \brief Register a new MeshOp factory
 * \param op_name The name by which this type of MeshOp can be requested
 * \param tps The MOAB entity types operated on by this MeshOp
 * \param num_tps Number of entity types in tps
 * \param meshop The (static) factory function producing instances of this MeshOp type
 * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
 */
bool MeshOpFactory::register_meshop(const char *op_name, moab::EntityType *tps, int num_tps,
                                    meshop_factory_t meshop, meshop_canmesh_t canmesh)
{
  OpInfoMap::iterator oit = registeredOps.find(std::string(op_name));
  if (oit != registeredOps.end()) 
    throw Error(MK_ALREADY_DEFINED, "A MeshOp with this name has already been registered.");
  
  OpInfo &oi = registeredOps[std::string(op_name)];
  oi.opName = std::string(op_name);
  std::copy(tps, tps+num_tps, oi.opEntTypes.begin());
  oi.opFactory = meshop;
  oi.opCanMesh = canmesh;

  return true;
}
  
MeshOpFactory::MeshOpFactory(MKCore *mk_core, bool create_if_missing) 
        : mkCore(mk_core), iCreatedCore(false)
{
  if (!mk_core && create_if_missing) {
    mkCore = new MKCore();
    iCreatedCore = true;
  }
}

    //! Destructor; virtual because applications may want to substitute their own factory
MeshOpFactory::~MeshOpFactory() 
{}

    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     */
MeshOpFactory::OpInfo MeshOpFactory::meshop_info(const char *op_name) 
{
  OpInfoMap::iterator oit = registeredOps.find(op_name);
  if (oit != registeredOps.end()) return oit->second;
  else throw Error(MK_NOT_FOUND, "A MeshOp with that name was not found.");
}
  
/** \brief Return MeshOp types that can operate on the specified entity type
 * \param tp Entity type requested
 * \param ops MeshOp types returned
 */
void MeshOpFactory::meshop_by_type(moab::EntityType tp, std::vector<OpInfo> &ops) 
{
  for (OpInfoMap::iterator oi = registeredOps.begin(); oi != registeredOps.end(); oi++) {
    if (std::find(oi->second.opEntTypes.begin(), oi->second.opEntTypes.end(), tp) != oi->second.opEntTypes.end())
      ops.push_back(oi->second);
  }
}
    
  /** \brief Return MeshOp types that can operate on mesh of specified dimension
   * \param dim Entity dimension requested
   * \param ops MeshOp types returned
   */
void MeshOpFactory::meshop_by_dimension(int dim, std::vector<OpInfo> &ops) 
{
  for (OpInfoMap::iterator oit = registeredOps.begin(); oit != registeredOps.end(); oit++) {
    for (std::vector<moab::EntityType>::iterator vit = oit->second.opEntTypes.begin();
         vit != oit->second.opEntTypes.end(); vit++) {
      if (moab::CN::Dimension(*vit) == dim)
        ops.push_back(oit->second);
    }
  }
}
    
    /** \brief Return MeshOp types that can mesh the specified ModelEnt
     * \param ent ModelEnt* requested
     * \param ops MeshOp types returned
     */
void MeshOpFactory::meshop_by_modelent(ModelEnt * const ent, std::vector<OpInfo> &ops) 
{
  for (OpInfoMap::iterator oit = registeredOps.begin(); oit != registeredOps.end(); oit++) {
    if (oit->second.opCanMesh && oit->second.opCanMesh(ent))
      ops.push_back(oit->second);
  }
}
    
/** \brief Construct a new MeshOp of the specified name
 * \param op_name MeshOp name being requested
 * \param me_vec Model entity vector to which this operation applies
 * \return Pointer to new MeshOp constructed
 */
MeshOp *MeshOpFactory::construct_meshop(std::string op_name, const MEVector &me_vec) 
{
  OpInfoMap::iterator oit = registeredOps.find(op_name);
  if (oit != registeredOps.end()) return oit->second.opFactory(mkCore, me_vec);
}

/** \brief Find an existing MeshOp in the graph, starting from the root
 * \param op_name MeshOp name being requested
 * \return Pointer to MeshOp found, NULL if not found
 */
MeshOp *MeshOpFactory::find_meshop(std::string op_name) 
{
  return NULL;
}

void MeshOpFactory::destroy_instance(bool dont_destroy_core) 
{
  if (!dont_destroy_core && iCreatedCore && mkCore) 
    delete mkCore;
  
  delete this;
}
    
  
    
} // namespace MeshKit


  
