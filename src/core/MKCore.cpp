#include "meshkit/MKCore.hpp"

namespace meshkit 
{
    
MKCore::MKCore(iGeom *igeom, moab::Interface *moab, iRel *irel,
               bool construct_missing_ifaces) throw(Error) 
: MeshOp(this), iGeomInstance(igeom), moabInstance(moab), iRelInstance(irel)
{
  init(construct_missing_ifaces);

    // add myself as the first parent and first child
  add_child(this);
}

void MKCore::~MKCore() throw(Error) 
{
  int err;
  if (iCreatedIgeom)
    delete iGeomInstance;
  
  if (iCreatedMoab)
    delete moabInstance;
  
  if (iCreatedIrel)
    delete iRelInstance;

    // be careful about how parent/child ops are handled, to avoid infinite loop
  if (!opParents.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<MeshOp*>::reverse_iterator vit = opParents.rbegin(); vit != opParents.rend(); vit++) {
      if (*vit == this) continue;
      else vit->remove_child(this);
    }
  }
  
  if (!opChildren.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<MeshOp*>::reverse_iterator vit = opChildren.rbegin(); vit != opChildren.rend(); vit++) {
      if (*vit == this) continue;
      else remove_child(*vit);
    }
  }

  assert(opParents.size() == 1 && opParents[0] == this &&
         opChildren.size() == 1 && opChildren[0] == this);

    // now clear them so we don't try to remove them in the parent class dtor
  opParents.clear();
  opChildren.clear();
}

void MKCore::init(bool construct_missing_ifaces) throw(Error) 
{
  iBase::Error err;

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
  
  if (!iRelInstance && construct_missing_ifaces) {
    iRelInstance = new iRel();
    if (iRelInstance) 
      throw Error(MK_FAILURE, "Failure creating iRel instance.");
    else iCreatedIrel = true;
  }
  
  err = iGeomInstance.createTag("__MKModelEntity", sizeof(MeshKit::ModelEnt*), iBase_BYTES,
                                igeomModelTag);
  MKERRCHK(err, "Failure to create MKModelEnt tag in iGeom.");

  moab::ErrorCode rval = Moab->create_tag("__MKModelEntity", sizeof(MeshKit::ModelEnt*), MB_TAG_SPARSE,
                                          MB_TYPE_OPAQUE, imeshModelTag);
  if (MB_SUCCESS != rval) 
    MKERRCHK(MK_FAILURE, "Failure to create MKModelEnt tag in iMesh.");
}

void MKCore::populate_mesh() throw(Error) 
{
    // populate mesh entity sets for geometric entities, relate them through iRel, and construct 
    // ModelEnts for them

  std::vector<iGeom::EntityHandle> ents;
  Error err;
  ModelEnt *this_me;
  
  for (int dim = iBase_VERTEX; dim <= iBase_REGION; dim++) {
    // get geometry entities
    ents.clear();
    err = iGeomInstance.getEntities(iGeom.rootSet(), dim, ents);
    MKERRCHK(err, "Failed to get entities from iGeom.");

    for (std::vector<iGeom::EntityHandle>::iterator eit = ents.begin(); eit != ents.end(); eit++) {
        // get the modelent
      this_me = NULL;
      err = iGeom.getData(*eit, igeomModelTag, &this_me);
      if (NULL == this_me || iBase_TAG_NOT_FOUND == err) {
          // construct a new ModelEnt and set the geom ent to point to it
        this_me = new ModelEnt(*this, *eit);
        err = iGeom.setData(*eit, igeomModelTag, &this_me);
        MKERRCHK(err, "Failed to set iGeom ModelEnt tag.");
      }
      
        // check for a mesh ent, and populate one if there is none
      if (!this_me->mesh_handle()) this_me->create_mesh_set();

        // save the mesh set handle here
      modelEnts[dim].insert(this_me->mesh_handle());
    }
  }
}

void MKCore::load_geometry(const char *filename, const char *options) throw(Error) 
{
  iBase::Error err = iGeomInstance.load(filename, options);
  MKERRCHK(err, "Failed to load geometry model.");
}

void MKCore::load_mesh(const char *filename, const char *options) throw(Error) 
{
  moab::ErrorCode rval = moabInstance->load_file(filename, NULL, options);
  MBERRCHK(rval, "Failed to load mesh.");
}

void MKCore::save_geometry(const char *filename, const char *options) throw(Error) 
{
  iBase::Error err = iGeomInstance.save(filename, options);
  MKERRCHK(err, "Failed to save geometry model.");
}

void MKCore::save_mesh(const char *filename, const char *options) throw(Error)
{
  moab::ErrorCode rval = moabInstance->write_file(filename, NULL, options);
  MBERRCHK(rval, "Failed to save mesh.");
}

void MKCore::get_entities_by_dimension(int dim, MEVector &model_ents) throw(Error) 
{
  int start = dim, end = dim;
  if (iBase_ALL_TYPES == dim) {
    start = 0;
    end = iBase_REGION;
  }
  
  std::vector<ModelEnt*> tmp_ents;
  for (dim = start; dim <= end; dim++) {
    tmp_ents.resize(modelEnts[dim].size());
    moab::ErrorCode rval = moabInstance->tag_get_data(imeshModelTag, modelEnts[dim], &tmp_ents[0]);
    MBERRCHK(rval, "Failed to get ModelEnt tag.");
    std::copy(tmp_ents.begin(), tmp_ents.end(), model_ents.end());
  }
}

void MKCore::get_entities_by_handle(MEVector &model_ents) throw(Error) 
{
  get_entities_by_dimension(iBase_ALL_TYPES, model_ents);
}

} // namespace meshkit
