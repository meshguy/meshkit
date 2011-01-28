#include "iGeom.hh"
#include "iRel.hh"
#include "moab/Core.hpp"
#include "moab/CN.hpp"
#include "MBiMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/NoOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "lemon/bfs.h"
#include "lemon/adaptors.h"

namespace MeshKit 
{
    
std::vector<MKCore::OpInfo> MKCore::registeredOps;
MKCore::OpNameMap MKCore::opNameMap;

MKCore::MKCore(iGeom *igeom, moab::Interface *moab, MBiMesh *mbi, iRel *irel,
               bool construct_missing_ifaces) 
        : iGeomInstance(igeom), moabInstance(moab), mbImesh(mbi), iRelInstance(irel),
          iRelPair(NULL), groupSetPair(0), iGeomModelTag(0), moabModelTag(0),
          iCreatedIgeom(false), iCreatedMoab(false), iCreatedMbimesh(false), iCreatedIrel(false),
          opsByDim(NULL), numOpsByDim(0)
{
    // leave initialization of root/leaf nodes to hear (and not in MKGraph), so that we have an MKCore
    // to pass to MeshOp's constructor
    // make the leaf/root nodes, link them with an edge; don't need to initialize map to MeshOp, since
    // by default that's NULL
  rootNode = new NoOp(this);
  leafNode = new NoOp(this);
  mkGraph.addArc(rootNode->get_node(), leafNode->get_node());

  init(construct_missing_ifaces);
}

MKCore::~MKCore() 
{
  int err;
  if (iCreatedIrel)
    delete iRelInstance;

  if (iCreatedIgeom)
    delete iGeomInstance;
  
  if (iCreatedMbimesh)
    delete mbImesh;

  if (iCreatedMoab)
    delete moabInstance;
  
  for (std::vector<SizingFunction*>::iterator vit = sizingFunctions.begin(); vit != sizingFunctions.end(); vit++)
    if (*vit) delete *vit;
  sizingFunctions.clear();
}

void MKCore::init(bool construct_missing_ifaces) 
{
  iBase_ErrorType err;

  if (!iGeomInstance && construct_missing_ifaces) {
    iGeomInstance = new iGeom();
    if (!iGeomInstance) 
      throw Error(MK_FAILURE, "Failure creating iGeom instance.");
    else iCreatedIgeom = true;
  }
  
  if (!moabInstance && construct_missing_ifaces) {
    moabInstance = new moab::Core();
    if (!moabInstance) throw Error(MK_FAILURE, "Failure creating MOAB instance.");
    else iCreatedMoab = true;
  }
  
  if (!mbImesh && construct_missing_ifaces) {
    mbImesh = new MBiMesh(moabInstance);
    if (!mbImesh) throw Error(MK_FAILURE, "Failure creating MBiMesh instance.");
    else iCreatedMbimesh = true;
  }
  
  if (!iRelInstance && construct_missing_ifaces) {
    iRelInstance = new iRel();
    if (!iRelInstance) 
      throw Error(MK_FAILURE, "Failure creating iRel instance.");
    else iCreatedIrel = true;
  }

  if (!iRelPair) {
    err = iRelInstance->createPair(iGeomInstance->instance(), iRel::ENTITY, iRel::IGEOM_IFACE,
                                   mbImesh, iRel::SET, iRel::IMESH_IFACE, iRelPair);
    IBERRCHK(err, "Failure to create relation pair.");
      // don't need to keep track of whether I created the pair, since it'll be deleted anyway when
      // the iRel instance is deleted.

      // FIXME: need a better scheme for finding any existing relation pairs or inferring them from 
      // imported model(s)
  }
  
  if (!groupSetPair) {
    err = iRelInstance->createPair(iGeomInstance->instance(), iRel::SET, iRel::IGEOM_IFACE,
                                   mbImesh, iRel::SET, iRel::IMESH_IFACE, groupSetPair);
    IBERRCHK(err, "Failure to create relation pair.");
  }
  
  err = iGeomInstance->createTag("__MKModelEntity", sizeof(MeshKit::ModelEnt*), iBase_BYTES,
                                 iGeomModelTag);
  IBERRCHK(err, "Failure to create MKModelEnt tag in iGeom.");

  ModelEnt *null_me = NULL;
  moab::ErrorCode rval = moabInstance->tag_create("__MKModelEntity", sizeof(MeshKit::ModelEnt*), moab::MB_TAG_SPARSE,
                                                  moab::MB_TYPE_OPAQUE, moabModelTag, &null_me);
  if (moab::MB_SUCCESS != rval) 
    IBERRCHK(iBase_FAILURE, "Failure to create MKModelEnt tag in iMesh.");

}

void MKCore::populate_mesh() 
{
    // populate mesh entity sets for geometric entities, relate them through iRel, and construct 
    // ModelEnts for them; also handle geometry groups

  std::vector<iGeom::EntityHandle> ents;
  std::vector<ModelEnt*> new_ents;
  iBase_ErrorType err;
  ModelEnt *this_me;
  
  for (int dim = iBase_VERTEX; dim <= iBase_REGION; dim++) {
    // get geometry entities
    ents.clear();
    err = iGeomInstance->getEntities(iGeomInstance->getRootSet(), (iBase_EntityType)dim, ents);
    IBERRCHK(err, "Failed to get entities from iGeom.");

    for (std::vector<iGeom::EntityHandle>::iterator eit = ents.begin(); eit != ents.end(); eit++) {
        // get the modelent
      this_me = NULL;
      err = iGeomInstance->getData(*eit, iGeomModelTag, &this_me);
      if (NULL == this_me || iBase_TAG_NOT_FOUND == err) {
          // construct a new ModelEnt and set the geom ent to point to it
        this_me = new ModelEnt(this, *eit);
        err = iGeomInstance->setData(*eit, iGeomModelTag, &this_me);
        IBERRCHK(err, "Failed to set iGeom ModelEnt tag.");
      }
      
        // check for a mesh ent, and populate one if there is none
      if (!this_me->mesh_handle()) {
        this_me->create_mesh_set();
        new_ents.push_back(this_me);
      }

        // save the mesh set handle here
      modelEnts[dim].push_back(this_me);
    }
  }
  
  for (std::vector<ModelEnt*>::iterator vit = new_ents.begin(); vit != new_ents.end(); vit++) 
    (*vit)->set_senses();

  std::vector<iGeom::EntitySetHandle> gsets;
  err = iGeomInstance->getEntSets(iGeomInstance->getRootSet(), -1, gsets);
  IBERRCHK(err, "Failed to get entity sets from iGeom.");
  for (std::vector<iGeom::EntitySetHandle>::iterator vit = gsets.begin();
       vit != gsets.end(); vit++) {
    this_me = NULL;
    err = iGeomInstance->getEntSetData(*vit, iGeomModelTag, &this_me);
    if (NULL == this_me || iBase_TAG_NOT_FOUND == err) {
        // construct a new ModelEnt and set the geom ent to point to it
      this_me = new ModelEnt(this, *vit);
      err = iGeomInstance->setEntSetData(*vit, iGeomModelTag, &this_me);
      IBERRCHK(err, "Failed to set iGeom ModelEnt tag.");
    }
      
      // check for a mesh ent, and populate one if there is none
    if (!this_me->mesh_handle()) {
      this_me->create_mesh_set();
    }
  }
}

void MKCore::load_geometry(const char *filename, const char *options, bool populate_too) 
{
  iBase_ErrorType err = iGeomInstance->load(filename, options);
  IBERRCHK(err, "Failed to load geometry model.");

  if (populate_too) populate_mesh();
}

void MKCore::load_mesh(const char *filename, const char *options) 
{
  moab::ErrorCode rval = moabInstance->load_file(filename, NULL, options);
  MBERRCHK(rval, "Failed to load mesh.");
}

void MKCore::save_geometry(const char *filename, const char *options) 
{
  iBase_ErrorType err = iGeomInstance->save(filename, options);
  IBERRCHK(err, "Failed to save geometry model.");
}

void MKCore::save_mesh(const char *filename, const char *options)
{
  moab::ErrorCode rval = moabInstance->write_file(filename, NULL, options);
  MBERRCHK(rval, "Failed to save mesh.");
}

void MKCore::get_entities_by_dimension(int dim, MEntVector &model_ents) 
{
  int start = dim, end = dim;
  if (iBase_ALL_TYPES == dim) {
    start = 0;
    end = iBase_REGION;
  }
  
  std::vector<ModelEnt*> tmp_ents;
  for (dim = start; dim <= end; dim++) {
    std::copy(modelEnts[dim].begin(), modelEnts[dim].end(), model_ents.end());
  }
}

void MKCore::get_entities_by_handle(MEntVector &model_ents) 
{
  get_entities_by_dimension(iBase_ALL_TYPES, model_ents);
}

int MKCore::add_sizing_function(SizingFunction *sf) 
{
  sizingFunctions.push_back(sf);
  return sizingFunctions.size()-1;
}

void MKCore::remove_sizing_function(int index) 
{
  if (index >= (int)sizingFunctions.size() || !sizingFunctions[index]) 
    throw Error(MK_BAD_INPUT, "No sizing function with that index.");

  sizingFunctions[index] = NULL;
}

/** \brief Register a new MeshOp factory
 * \param op_name The name by which this type of MeshOp can be requested
 * \param tp The MOAB entity type operated on by this MeshOp
 * \param meshop The (static) factory function producing instances of this MeshOp type
 * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
 */
bool MKCore::register_meshop(const char *op_name, moab::EntityType tp, 
                                    meshop_factory_t meshop, meshop_canmesh_t canmesh) 
{
  OpNameMap::iterator oit = opNameMap.find(std::string(op_name));
  if (oit != opNameMap.end()) 
    throw Error(MK_ALREADY_DEFINED, "A MeshOp with this name has already been registered.");
  
  OpInfo oi = {std::string(op_name), registeredOps.size(), std::vector<moab::EntityType>(1, tp),
               meshop, canmesh};
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
bool MKCore::register_meshop(const char *op_name, moab::EntityType *tps, int num_tps,
                                    meshop_factory_t meshop, meshop_canmesh_t canmesh)
{
  OpNameMap::iterator oit = opNameMap.find(std::string(op_name));
  if (oit != opNameMap.end()) 
    throw Error(MK_ALREADY_DEFINED, "A MeshOp with this name has already been registered.");
  
  OpInfo oi = {std::string(op_name), registeredOps.size(), std::vector<moab::EntityType>(tps, tps+num_tps),
               meshop, canmesh};
  opNameMap[oi.opName] = registeredOps.size();
  
  registeredOps.push_back(oi);

  return true;
}
  
    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return OpInfo for the corresponding MeshOp type
     */
MKCore::OpInfo MKCore::meshop_info(const char *op_name) 
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
unsigned int MKCore::meshop_index(const char *op_name) 
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
void MKCore::set_default_meshop(const char *op_name, unsigned short dims) 
{
    // get the meshop index
  OpNameMap::iterator oit = opNameMap.find(op_name);
  if (oit == opNameMap.end()) throw Error(MK_NOT_FOUND, "A MeshOp with that name was not found.");

  set_default_meshop(oit->second, dims);
}

  /** \brief Initialize the opsByDim vectors from (static) registeredOps
   */
void MKCore::init_opsbydim() 
{
    // compare cached size to registered size, return if it hasn't changed since we last checked
  if (opsByDim && numOpsByDim == registeredOps.size()) return;
  
  if (!opsByDim) opsByDim = new std::vector<unsigned short>[4];
  
    // reset numOpsByDim, which we'll increment next
  opsByDim = 0;

    // check each registered op
  for (std::vector<OpInfo>::iterator vit = registeredOps.begin(); vit != registeredOps.end(); vit++) {
      // cache dimension so we don't check the same dim twice
    int last_dim = -1;
    opsByDim++;
    for (std::vector<moab::EntityType>::iterator vit2 = (*vit).opEntTypes.begin(); 
         vit2 != (*vit).opEntTypes.end(); vit2++) {
      if (moab::CN::Dimension(*vit2) == last_dim) continue;
      last_dim = moab::CN::Dimension(*vit2);
      if (std::find(opsByDim[last_dim].begin(), opsByDim[last_dim].end(), (*vit).opIndex) == 
          opsByDim[last_dim].end()) opsByDim[last_dim].push_back((*vit).opIndex);
    }
  }
}
    
  /** \brief Make the specified MeshOp name the default for the given dimension(s)
   * 
   * If the specified MeshOp cannot produce entities of the specified dimension, an error is
   * thrown with type MK_BAD_INPUT.
   * \param op_index MeshOp index being set
   * \param dims Bitmask, where 2^x indicates that this MeshOp should be the default for dimension x 
   */
void MKCore::set_default_meshop(unsigned short op_index, unsigned short dims) 
{
  if (!opsByDim) init_opsbydim();

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
void MKCore::meshop_by_type(moab::EntityType tp, std::vector<OpInfo> &ops) 
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
void MKCore::meshop_by_dimension(int dim, std::vector<OpInfo> &ops) 
{
  if (!opsByDim) init_opsbydim();

  for (std::vector<unsigned short>::iterator vit = opsByDim[dim].begin(); vit != opsByDim[dim].end();
       vit++)
    ops.push_back(registeredOps[*vit]);
}
    
    /** \brief Return MeshOp types that can mesh the specified ModelEnt
     * \param ent ModelEnt* requested
     * \param ops MeshOp types returned
     */
void MKCore::meshop_by_modelent(ModelEnt * const ent, std::vector<OpInfo> &ops) 
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
MeshOp *MKCore::construct_meshop(std::string op_name, const MEntVector &me_vec) 
{
  OpNameMap::iterator oit = opNameMap.find(op_name);
  if (oit != opNameMap.end()) return registeredOps[oit->second].opFactory(this, me_vec);
}

/** \brief Construct a new MeshOp of the specified name
 * \param op_name MeshOp name being requested
 * \param me_vec Model entity vector to which this operation applies
 * \return Pointer to new MeshOp constructed
 */
MeshOp *MKCore::construct_meshop(unsigned int dim, const MEntVector &me_vec) 
{
  if (!opsByDim) init_opsbydim();
  
  if (opsByDim[dim].empty()) throw Error(MK_MESHOP_NOT_FOUND, "No default MeshOp for that dimension.");
  else return construct_meshop(registeredOps[opsByDim[dim][0]], me_vec);
}

MeshOp *MKCore::construct_meshop(OpInfo &info, const MEntVector &me_vec) 
{
  return registeredOps[info.opIndex].opFactory(this, me_vec);
}

} // namespace meshkit
