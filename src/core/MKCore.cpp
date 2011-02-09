#include "meshkit/iGeom.hh"
#include "meshkit/iMesh.hh"
#include "meshkit/iRel.hh"
#include "moab/Core.hpp"
#include "MBiMesh.hpp"
#include "moab/CN.hpp"
#include "MBTagConventions.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/NoOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "lemon/bfs.h"
#include "lemon/adaptors.h"

namespace MeshKit 
{
    
MKCore::MeshOpFactory *MKCore::opFactory = NULL;

MKCore::MeshOpFactory *MKCore::op_factory()
{ 
  if (!opFactory) opFactory = new MeshOpFactory;
  return opFactory;
}

MKCore::MKCore(iGeom *igeom, moab::Interface *moab, iMesh *imesh, iRel *irel,
               bool construct_missing_ifaces) 
        : opsByDim(NULL), numOpsByDim(0), vertexMesher(NULL)
{

  if (igeom) {
    iGeomInstances.push_back(igeom);
    iCreatedIgeoms.push_back(false);
  }
  
  if (moab) {
    moabInstances.push_back(moab);
    iCreatedMoabs.push_back(false);
  }
  
  if (imesh) {
    iMeshInstances.push_back(imesh);
    iCreatedImeshs.push_back(false);
  }
  
  if (irel) {
    iRelInstances.push_back(irel);
    iCreatedIrels.push_back(false);
  }

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
  for (unsigned int i = 0; i < iRelInstances.size(); i++) 
    if (iCreatedIrels[i]) delete iRelInstances[i];

  for (unsigned int i = 0; i < iGeomInstances.size(); i++) 
    if (iCreatedIgeoms[i]) delete iGeomInstances[i];

  for (unsigned int i = 0; i < moabInstances.size(); i++) 
    if (iCreatedMoabs[i]) delete moabInstances[i];

  for (unsigned int i = 0; i < iMeshInstances.size(); i++) 
    if (iCreatedImeshs[i]) delete iMeshInstances[i];

  for (std::vector<SizingFunction*>::iterator vit = sizingFunctions.begin(); vit != sizingFunctions.end(); vit++)
    if (*vit) delete *vit;
  sizingFunctions.clear();
}

void MKCore::init(bool construct_missing_ifaces) 
{
  iBase_ErrorType err;

  if (iGeomInstances.empty() && construct_missing_ifaces) {
    iGeomInstances.push_back(new iGeom());
    iCreatedIgeoms.push_back(true);
  }
  
  if (moabInstances.empty() && construct_missing_ifaces) {
    moabInstances.push_back(new moab::Core());
    iCreatedMoabs.push_back(true);
  }
  
  if (iMeshInstances.empty() && construct_missing_ifaces) {
    iMeshInstances.push_back(new iMesh((iMesh_Instance)new MBiMesh(moabInstances[0])));
    iCreatedImeshs.push_back(true);
  }
  
  if (iRelInstances.empty() && construct_missing_ifaces) {
    iRelInstances.push_back(new iRel());
    iCreatedIrels.push_back(true);
  }

  if (iRelPairs.empty() && !iRelInstances.empty() && !iGeomInstances.empty() && !iMeshInstances.empty()) {
    iRelPairs.resize(1);
    err = iRelInstances[0]->createPair(iGeomInstances[0]->instance(), iRel::ENTITY, iRel::IGEOM_IFACE,
                                       iMeshInstances[0]->instance(), iRel::SET, iRel::IMESH_IFACE, iRelPairs[0]);
    IBERRCHK(err, "Failure to create relation pair.");
      // don't need to keep track of whether I created the pair, since it'll be deleted anyway when
      // the iRel instance is deleted.

      // FIXME: need a better scheme for finding any existing relation pairs or inferring them from 
      // imported model(s)
  }
  
  if (groupSetPairs.empty() && !iRelInstances.empty() && !iGeomInstances.empty() && !iMeshInstances.empty()) {
    groupSetPairs.resize(1);
    err = iRelInstances[0]->createPair(iGeomInstances[0]->instance(), iRel::SET, iRel::IGEOM_IFACE,
                                       iMeshInstances[0]->instance(), iRel::SET, iRel::IMESH_IFACE, groupSetPairs[0]);
    IBERRCHK(err, "Failure to create relation pair.");
  }

  if (iGeomModelTags.empty()) {
    iGeomModelTags.resize(1);
    err = iGeomInstances[0]->createTag("__MKModelEntity", sizeof(MeshKit::ModelEnt*), iBase_BYTES,
                                       iGeomModelTags[0]);
    IBERRCHK(err, "Failure to create MKModelEnt tag in iGeom.");
  }

  moab::ErrorCode rval;
  if (moabModelTags.empty()) {
    ModelEnt *null_me = NULL;
    moabModelTags.resize(1);
    rval = moabInstances[0]->tag_create("__MKModelEntity", sizeof(MeshKit::ModelEnt*), moab::MB_TAG_SPARSE,
                                        moab::MB_TYPE_OPAQUE, moabModelTags[0], &null_me);
    if (moab::MB_SUCCESS != rval && moab::MB_ALREADY_ALLOCATED != rval) 
      MBERRCHK(rval, "Failure to create MKModelEnt tag in iMesh.");
  }

  if (moabGeomDimTags.empty()) {
      // moab geometry dimension tag
    moabGeomDimTags.resize(1);
    rval = moabInstances[0]->tag_create(GEOM_DIMENSION_TAG_NAME, sizeof(int), moab::MB_TAG_SPARSE,
                                        moab::MB_TYPE_INTEGER, moabGeomDimTags[0], 0, true);
    if (moab::MB_SUCCESS != rval && moab::MB_ALREADY_ALLOCATED != rval) 
      MBERRCHK(rval, "Failure to create geometry dimension tag in iMesh.");
  }
  
  // moab global id tag
  if (moabIDTags.empty()) {
    moabIDTags.resize(1);
    rval = moabInstances[0]->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), moab::MB_TAG_DENSE,
                                        moab::MB_TYPE_INTEGER, moabIDTags[0], 0, true);
    if (moab::MB_SUCCESS != rval && moab::MB_ALREADY_ALLOCATED != rval) 
      MBERRCHK(rval, "Failure to create global id tag in iMesh.");
  }
}

void MKCore::populate_mesh(int index) 
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
    err = igeom_instance(index)->getEntities(igeom_instance(index)->getRootSet(), (iBase_EntityType)dim, ents);
    IBERRCHK(err, "Failed to get entities from iGeom.");

    for (std::vector<iGeom::EntityHandle>::iterator eit = ents.begin(); eit != ents.end(); eit++) {
        // get the modelent
      this_me = NULL;
      err = igeom_instance(index)->getData(*eit, iGeomModelTags[index], &this_me);
      if (NULL == this_me || iBase_TAG_NOT_FOUND == err) {
          // construct a new ModelEnt and set the geom ent to point to it
        this_me = new ModelEnt(this, *eit);
        err = igeom_instance(index)->setData(*eit, iGeomModelTags[index], &this_me);
        IBERRCHK(err, "Failed to set iGeom ModelEnt tag.");
      }
      
        // check for a mesh ent, and populate one if there is none
      if (!this_me->mesh_handle()) {
	if (!this_me->exist_mesh_set()) {
	  this_me->create_mesh_set();
	}
	new_ents.push_back(this_me);
      }

        // save the model entity here
      modelEnts[dim].push_back(this_me);
    }
  }
  
  for (std::vector<ModelEnt*>::iterator vit = new_ents.begin(); vit != new_ents.end(); vit++) 
    (*vit)->set_senses();

  std::vector<iGeom::EntitySetHandle> gsets;
  err = igeom_instance(index)->getEntSets(igeom_instance(index)->getRootSet(), -1, gsets);
  IBERRCHK(err, "Failed to get entity sets from iGeom.");
  for (std::vector<iGeom::EntitySetHandle>::iterator vit = gsets.begin();
       vit != gsets.end(); vit++) {
    this_me = NULL;
    err = igeom_instance(index)->getEntSetData(*vit, iGeomModelTags[index], &this_me);
    if (NULL == this_me || iBase_TAG_NOT_FOUND == err) {
        // construct a new ModelEnt and set the geom ent to point to it
      this_me = new ModelEnt(this, *vit);
      err = igeom_instance(index)->setEntSetData(*vit, iGeomModelTags[index], &this_me);
      IBERRCHK(err, "Failed to set iGeom ModelEnt tag.");
    }
      
      // check for a mesh ent, and populate one if there is none
    if (!this_me->mesh_handle()) {
      this_me->create_mesh_set(index);
    }
  }
}

    void MKCore::load_geometry(const char *filename, const char *options, int index, bool populate_too) 
{
  iBase_ErrorType err = igeom_instance(index)->load(filename, options);
  IBERRCHK(err, "Failed to load geometry model.");

  if (populate_too) populate_mesh(index);
}

void MKCore::load_mesh(const char *filename, const char *options, int index) 
{
  moab::ErrorCode rval = moabInstances[index]->load_file(filename, NULL, options);
  MBERRCHK(rval, "Failed to load mesh.");
}

void MKCore::save_geometry(const char *filename, const char *options, int index) 
{
  iBase_ErrorType err = igeom_instance(index)->save(filename, options);
  IBERRCHK(err, "Failed to save geometry model.");
}

void MKCore::save_mesh(const char *filename, const char *options, int index)
{
  moab::ErrorCode rval = moab_instance(index)->write_file(filename, NULL, options);
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
    std::copy(modelEnts[dim].begin(), modelEnts[dim].end(), std::back_inserter(model_ents));
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

SizingFunction *MKCore::sizing_function(double size, bool create_if_missing) 
{
  for (unsigned int i = 0; i < sizingFunctions.size(); i++)
    if (sizingFunctions[i]->size() == size) return sizingFunctions[i];
    
    // if we got here, either create one or return NULL
  if (!create_if_missing) return NULL;
  
  return new SizingFunction(this, -1, size);
}
  
/** \brief Register a new MeshOp factory
 * \param op_name The name by which this type of MeshOp can be requested
 * \param tp The MOAB entity type operated on by this MeshOp
 * \param meshop The (static) factory function producing instances of this MeshOp type
 * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
 */
bool MKCore::register_meshop(const char *op_name, 
                             iBase_EntityType model_tp, moab::EntityType tp,
                             meshop_factory_t meshop, meshop_canmesh_t canmesh) 
{
  OpNameMap::iterator oit = op_factory()->opNameMap.find(std::string(op_name));
  if (oit != op_factory()->opNameMap.end()) 
    throw Error(MK_ALREADY_DEFINED, "%s, line %d: A MeshOp with the name \"%s\" has already been registered.", 
                __FILE__, __LINE__, op_name);
  
  OpInfo oi = {std::string(op_name), op_factory()->registeredOps.size(), 
               std::vector<iBase_EntityType>(1, model_tp), std::vector<moab::EntityType>(1, tp),
               meshop, canmesh};
  op_factory()->opNameMap[oi.opName] = op_factory()->registeredOps.size();
  
  op_factory()->registeredOps.push_back(oi);

  return true;
}
  
/** \brief Register a new MeshOp factory
 * \param op_name The name by which this type of MeshOp can be requested
 * \param tps The MOAB entity types operated on by this MeshOp
 * \param num_tps Number of entity types in tps
 * \param meshop The (static) factory function producing instances of this MeshOp type
 * \param canmesh If provided, a static function that returns whether the MeshOp can mesh a given ModelEnt
 */
bool MKCore::register_meshop(const char *op_name, 
                             iBase_EntityType *model_tps, int num_mtps,
                             moab::EntityType *tps, int num_tps,
                             meshop_factory_t meshop, meshop_canmesh_t canmesh)
{
  OpNameMap::iterator oit = op_factory()->opNameMap.find(std::string(op_name));
  if (oit != op_factory()->opNameMap.end()) 
    throw Error(MK_ALREADY_DEFINED, "%s, line %d: A MeshOp with the name \"%s\" has already been registered.", 
                __FILE__, __LINE__, op_name);
  
  OpInfo oi = {std::string(op_name), op_factory()->registeredOps.size(), 
               std::vector<iBase_EntityType>(model_tps, model_tps+num_mtps),
               std::vector<moab::EntityType>(tps, tps+num_tps),
               meshop, canmesh};
  op_factory()->opNameMap[oi.opName] = op_factory()->registeredOps.size();
  
  op_factory()->registeredOps.push_back(oi);

  return true;
}
  
    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return OpInfo for the corresponding MeshOp type
     */
MKCore::OpInfo MKCore::meshop_info(const char *op_name) 
{
  OpNameMap::iterator oit = op_factory()->opNameMap.find(op_name);
  if (oit != op_factory()->opNameMap.end()) {
    assert(op_factory()->registeredOps.size() > oit->second);
    return op_factory()->registeredOps[oit->second];
  }
  else throw Error(MK_NOT_FOUND, "%s, line %d: A MeshOp with name %s was not found.", __FILE__, __LINE__, op_name);
}
  
    /** \brief Return the MeshOp type with the given name
     * \param op_name Operation name requested
     * \return OpInfo index for the corresponding MeshOp type
     */
unsigned int MKCore::meshop_index(const char *op_name) 
{
  OpNameMap::iterator oit = op_factory()->opNameMap.find(op_name);
  if (oit != op_factory()->opNameMap.end()) {
    assert(op_factory()->registeredOps.size() > oit->second);
    return oit->second;
  }
  else throw Error(MK_NOT_FOUND, "%s, line %d: A MeshOp with name %s was not found.", __FILE__, __LINE__, op_name);
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
  OpNameMap::iterator oit = op_factory()->opNameMap.find(op_name);
  if (oit == op_factory()->opNameMap.end()) 
    throw Error(MK_NOT_FOUND, "%s, line %d: A MeshOp with name %s was not found.", __FILE__, __LINE__, op_name);

  set_default_meshop(oit->second, dims);
}

  /** \brief Initialize the opsByDim vectors from (static) op_factory()->registeredOps
   */
void MKCore::init_opsbydim() 
{
    // compare cached size to registered size, return if it hasn't changed since we last checked
  if (opsByDim && numOpsByDim == op_factory()->registeredOps.size()) return;
  
  if (!opsByDim) opsByDim = new std::vector<unsigned short>[4];
  
    // reset numOpsByDim, which we'll increment next
  numOpsByDim = 0;

    // check each registered op
  for (std::vector<OpInfo>::iterator vit = op_factory()->registeredOps.begin(); vit != op_factory()->registeredOps.end(); vit++) {
      // cache dimension so we don't check the same dim twice
    numOpsByDim++;
    for (std::vector<iBase_EntityType>::iterator vit2 = (*vit).modelEntTypes.begin(); 
         vit2 != (*vit).modelEntTypes.end(); vit2++) {
      int this_dim = (int)(*vit2);
      if (std::find(opsByDim[this_dim].begin(), opsByDim[this_dim].end(), (*vit).opIndex) == 
          opsByDim[this_dim].end()) opsByDim[this_dim].push_back((*vit).opIndex);
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
void MKCore::meshop_by_mesh_type(moab::EntityType tp, std::vector<OpInfo> &ops) 
{
  for (std::vector<OpInfo>::iterator oi = op_factory()->registeredOps.begin(); oi != op_factory()->registeredOps.end(); oi++) {
    if (std::find(oi->meshEntTypes.begin(), oi->meshEntTypes.end(), tp) != oi->meshEntTypes.end())
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
    ops.push_back(op_factory()->registeredOps[*vit]);
}
    
    /** \brief Return MeshOp types that can mesh the specified ModelEnt
     * \param ent ModelEnt* requested
     * \param ops MeshOp types returned
     */
void MKCore::meshop_by_modelent(ModelEnt * const ent, std::vector<OpInfo> &ops) 
{
  for (std::vector<OpInfo>::iterator oit = op_factory()->registeredOps.begin(); oit != op_factory()->registeredOps.end(); oit++) {
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
  OpNameMap::iterator oit = op_factory()->opNameMap.find(op_name);
  if (oit != op_factory()->opNameMap.end()) return op_factory()->registeredOps[oit->second].opFactory(this, me_vec);
  else 
    throw Error(MK_MESHOP_NOT_FOUND, "%s, line %d: Didn't find a MeshOp with name \"%s\".", 
                __FILE__, __LINE__, op_name.c_str());
}

/** \brief Construct a new MeshOp of the specified name
 * \param op_name MeshOp name being requested
 * \param me_vec Model entity vector to which this operation applies
 * \return Pointer to new MeshOp constructed
 */
MeshOp *MKCore::construct_meshop(unsigned int dim, const MEntVector &me_vec) 
{
  if (!opsByDim) init_opsbydim();
  
  if (opsByDim[dim].empty()) 
    throw Error(MK_MESHOP_NOT_FOUND, "%s, line %d: No default MeshOp for dimension %d.", __FILE__, __LINE__, dim);
  else return construct_meshop(op_factory()->registeredOps[opsByDim[dim][0]], me_vec);
}

MeshOp *MKCore::construct_meshop(OpInfo &info, const MEntVector &me_vec) 
{
  return op_factory()->registeredOps[info.opIndex].opFactory(this, me_vec);
}

} // namespace meshkit
