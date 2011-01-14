#include "iGeom.hh"
#include "iRel.hh"
#include "moab/Core.hpp"
#include "MBiMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"

namespace MeshKit 
{
    
MKCore::MKCore(iGeom *igeom, moab::Interface *moab, MBiMesh *mbi, iRel *irel,
               bool construct_missing_ifaces) 
        : iGeomInstance(igeom), moabInstance(moab), mbImesh(mbi), iRelInstance(irel),
          iRelPair(NULL), iGeomModelTag(0), moabModelTag(0),
          iCreatedIgeom(false), iCreatedMoab(false), iCreatedMbimesh(false), iCreatedIrel(false)

{
  init(construct_missing_ifaces);

    // add myself as the first parent and first child
  add_child(this);
}

MKCore::~MKCore() 
{
  int err;
  if (iCreatedIgeom)
    delete iGeomInstance;
  
  if (iCreatedMoab)
    delete moabInstance;
  
  if (iCreatedIrel)
    delete iRelInstance;

  if (iCreatedMbimesh)
    delete mbImesh;
  
    // be careful about how parent/child ops are handled, to avoid infinite loop
  if (!nodeParents.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<GraphNode*>::reverse_iterator vit = nodeParents.rbegin(); vit != nodeParents.rend(); vit++) {
      if (*vit == this) continue;
      else (*vit)->remove_child(this);
    }
  }
  
  if (!nodeChildren.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<GraphNode*>::reverse_iterator vit = nodeChildren.rbegin(); vit != nodeChildren.rend(); vit++) {
      if (*vit == this) continue;
      else remove_child(*vit);
    }
  }

  assert(nodeParents.size() == 1 && nodeParents[0] == this &&
         nodeChildren.size() == 1 && nodeChildren[0] == this);

    // now clear them so we don't try to remove them in the parent class dtor
  nodeParents.clear();
  nodeChildren.clear();
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
                                   mbImesh, iRel::ENTITY, iRel::IMESH_IFACE, iRelPair);
    IBERRCHK(err, "Failure to create relation pair.");
      // don't need to keep track of whether I created the pair, since it'll be deleted anyway when
      // the iRel instance is deleted.

      // FIXME: need a better scheme for finding any existing relation pairs or inferring them from 
      // imported model(s)
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
    // ModelEnts for them

  std::vector<iGeom::EntityHandle> ents;
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
      if (!this_me->mesh_handle()) this_me->create_mesh_set();

        // save the mesh set handle here
      modelEnts[dim].insert(this_me->mesh_handle());
    }
  }
}

void MKCore::load_geometry(const char *filename, const char *options) 
{
  iBase_ErrorType err = iGeomInstance->load(filename, options);
  IBERRCHK(err, "Failed to load geometry model.");
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

void MKCore::get_entities_by_dimension(int dim, MEVector &model_ents) 
{
  int start = dim, end = dim;
  if (iBase_ALL_TYPES == dim) {
    start = 0;
    end = iBase_REGION;
  }
  
  std::vector<ModelEnt*> tmp_ents;
  for (dim = start; dim <= end; dim++) {
    tmp_ents.resize(modelEnts[dim].size());
    moab::ErrorCode rval = moabInstance->tag_get_data(moabModelTag, modelEnts[dim], &tmp_ents[0]);
    MBERRCHK(rval, "Failed to get ModelEnt tag.");
    std::copy(tmp_ents.begin(), tmp_ents.end(), model_ents.end());
  }
}

void MKCore::get_entities_by_handle(MEVector &model_ents) 
{
  get_entities_by_dimension(iBase_ALL_TYPES, model_ents);
}

bool MKCore::is_root() 
{
  return true;
}

} // namespace meshkit
