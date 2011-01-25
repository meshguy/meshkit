#include "iGeom.hh"
#include "iRel.hh"
#include "moab/Core.hpp"
#include "MBiMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/NoOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "lemon/bfs.h"
#include "lemon/adaptors.h"

namespace MeshKit 
{
    
MKCore::MKCore(iGeom *igeom, moab::Interface *moab, MBiMesh *mbi, iRel *irel,
               bool construct_missing_ifaces) 
        : iGeomInstance(igeom), moabInstance(moab), mbImesh(mbi), iRelInstance(irel),
          iRelPair(NULL), iGeomModelTag(0), moabModelTag(0),
          iCreatedIgeom(false), iCreatedMoab(false), iCreatedMbimesh(false), iCreatedIrel(false)
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
    std::copy(modelEnts[dim].begin(), modelEnts[dim].end(), model_ents.end());
  }
}

void MKCore::get_entities_by_handle(MEVector &model_ents) 
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

} // namespace meshkit
