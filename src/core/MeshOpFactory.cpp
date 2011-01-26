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
  OpNameMap::iterator oit = opNameMap.find(std::string(op_name));
  if (oit != opNameMap.end()) 
    throw Error(MK_ALREADY_DEFINED, "A MeshOp with this name has already been registered.");
  
  OpInfo oi = {std::string(op_name), registeredOps.size(), std::vector<moab::EntityType>(1, tp),
               meshop, canmesh};
  opsByDim[moab::CN::Dimension(tp)].push_back(registeredOps.size());
  opNameMap[oi.opName] = registeredOps.size();
  
  registeredOps.push_back(oi);

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
  OpNameMap::iterator oit = opNameMap.find(std::string(op_name));
  if (oit != opNameMap.end()) 
    throw Error(MK_ALREADY_DEFINED, "A MeshOp with this name has already been registered.");
  
  OpInfo oi = {std::string(op_name), registeredOps.size(), std::vector<moab::EntityType>(tps, tps+num_tps),
               meshop, canmesh};
  for (int i = 0; i < num_tps; i++)
    opsByDim[moab::CN::Dimension(tps[i])].push_back(registeredOps.size());
  opNameMap[oi.opName] = registeredOps.size();
  
  registeredOps.push_back(oi);

  return true;
}
  
MeshOpFactory::MeshOpFactory(MKCore *mk_core, bool create_if_missing) 
        : mkCore(mk_core), iCreatedCore(false)
{
  if (!mk_core) {
    mkCore = new MKCore();
    iCreatedCore = true;
  }
}

    //! Destructor; virtual because applications may want to substitute their own factory
MeshOpFactory::~MeshOpFactory() 
{}

    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return OpInfo for the corresponding MeshOp type
     */
MeshOpFactory::OpInfo MeshOpFactory::meshop_info(const char *op_name) 
{
  OpNameMap::iterator oit = opNameMap.find(op_name);
  if (oit != opNameMap.end()) {
    assert(registeredOps.size() > oit->second);
    return registeredOps[oit->second];
  }
  else throw Error(MK_NOT_FOUND, "A MeshOp with that name was not found.");
}
  
    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return OpInfo index for the corresponding MeshOp type
     */
unsigned int MeshOpFactory::meshop_index(const char *op_name) 
{
  OpNameMap::iterator oit = opNameMap.find(op_name);
  if (oit != opNameMap.end()) {
    assert(registeredOps.size() > oit->second);
    return oit->second;
  }
  else throw Error(MK_NOT_FOUND, "A MeshOp with that name was not found.");
}
  
  /** \brief Make the specified MeshOp name the default for the given dimension(s)
   * 
   * If the specified MeshOp cannot produce entities of the specified dimension, an error is
   * thrown with type MK_BAD_INPUT.
   * \param op_name MeshOp name being set
   * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
   */
void MeshOpFactory::set_default_meshop(const char *op_name, unsigned short dims) 
{
    // get the meshop index
  OpNameMap::iterator oit = opNameMap.find(op_name);
  if (oit == opNameMap.end()) throw Error(MK_NOT_FOUND, "A MeshOp with that name was not found.");

  set_default_meshop(oit->second, dims);
}
    
  /** \brief Make the specified MeshOp name the default for the given dimension(s)
   * 
   * If the specified MeshOp cannot produce entities of the specified dimension, an error is
   * thrown with type MK_BAD_INPUT.
   * \param op_index MeshOp index being set
   * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
   */
void MeshOpFactory::set_default_meshop(unsigned short op_index, unsigned short dims) 
{
    // check the specified dimension(s) against the types the meshop can mesh
  std::vector<unsigned short> good_dims;
  for (unsigned short i = 0; i < 4; i++) {
    if (dims & (1u << i)) {
        // verify this op can generate elements with dimension i; just check for the index in the vector
      if (std::find(opsByDim[i].begin(), opsByDim[i].end(), op_index) == opsByDim[i].end())
        throw Error(MK_BAD_INPUT, "Specified MeshOp type cannot generate elements of specified dimension.");
      else good_dims.push_back(i);
    }
  }
  
    // if we're here, need to move specified op to front of list for all good_dims
  for (std::vector<unsigned short>::iterator vit = good_dims.begin(); vit != good_dims.end(); vit++) {
    std::vector<unsigned short>::iterator vit2 = std::find(opsByDim[*vit].begin(), opsByDim[*vit].end(), op_index);
    assert(vit2 != opsByDim[*vit].end());
      // don't need to do anything if it's already the default for this dim
    if (vit2 == opsByDim[*vit].begin()) continue;
    
      // else switch with 1st element
    *vit2 = *opsByDim[*vit].begin();
    *opsByDim[*vit].begin() = op_index;
  }
}

/* \brief Return MeshOp types that can operate on the specified entity type
 * \param tp Entity type requested
 * \param ops MeshOp types returned
 */
void MeshOpFactory::meshop_by_type(moab::EntityType tp, std::vector<OpInfo> &ops) 
{
  for (std::vector<OpInfo>::iterator oi = registeredOps.begin(); oi != registeredOps.end(); oi++) {
    if (std::find(oi->opEntTypes.begin(), oi->opEntTypes.end(), tp) != oi->opEntTypes.end())
      ops.push_back(*oi);
  }
}
    
  /** \brief Return MeshOp types that can operate on mesh of specified dimension
   * \param dim Entity dimension requested
   * \param ops MeshOp types returned
   */
void MeshOpFactory::meshop_by_dimension(int dim, std::vector<OpInfo> &ops) 
{
  for (std::vector<unsigned short>::iterator vit = opsByDim[dim].begin(); vit != opsByDim[dim].end();
       vit++)
    ops.push_back(registeredOps[*vit]);
}
    
    /** \brief Return MeshOp types that can mesh the specified ModelEnt
     * \param ent ModelEnt* requested
     * \param ops MeshOp types returned
     */
void MeshOpFactory::meshop_by_modelent(ModelEnt * const ent, std::vector<OpInfo> &ops) 
{
  for (std::vector<OpInfo>::iterator oit = registeredOps.begin(); oit != registeredOps.end(); oit++) {
    if (oit->opCanMesh && oit->opCanMesh(ent))
      ops.push_back(*oit);
  }
}
    
/** \brief Construct a new MeshOp of the specified name
 * \param op_name MeshOp name being requested
 * \param me_vec Model entity vector to which this operation applies
 * \return Pointer to new MeshOp constructed
 */
MeshOp *MeshOpFactory::construct_meshop(std::string op_name, const MEVector &me_vec) 
{
  OpNameMap::iterator oit = opNameMap.find(op_name);
  if (oit != opNameMap.end()) return registeredOps[oit->second].opFactory(mkCore, me_vec);
}

/** \brief Construct a new MeshOp of the specified name
 * \param op_name MeshOp name being requested
 * \param me_vec Model entity vector to which this operation applies
 * \return Pointer to new MeshOp constructed
 */
MeshOp *MeshOpFactory::construct_meshop(unsigned int dim, const MEVector &me_vec) 
{
  if (opsByDim[dim].empty()) throw Error(MK_MESHOP_NOT_FOUND, "No default MeshOp for that dimension.");
  else return construct_meshop(registeredOps[opsByDim[dim][0]], me_vec);
}

MeshOp *MeshOpFactory::construct_meshop(OpInfo &info, const MEVector &me_vec) 
{
  return registeredOps[info.opIndex].opFactory(mkCore, me_vec);
}

void MeshOpFactory::destroy_instance(bool dont_destroy_core) 
{
  if (!dont_destroy_core && iCreatedCore && mkCore) 
    delete mkCore;
  
  delete this;
}
    
  
    
} // namespace MeshKit


  
