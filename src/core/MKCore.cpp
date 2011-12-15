#include "meshkit/iGeom.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/iRel.hpp"
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
#include "MeshOpSet.hpp"
#include "meshkit/MeshOpProxy.hpp"
#include "meshkit/VertexMesher.hpp"
#if HAVE_FBIGEOM
#include "meshkit/FBiGeom.hpp"
// this is defined in moab library...
extern void FBiGeom_newGeomFromMesh( iMesh_Instance mesh, iBase_EntitySetHandle set,
                          const char *options, FBiGeom_Instance *geom,
                          int *err, int options_len);
#endif
namespace MeshKit 
{

// Declare these here rather than in a common header because 
// we don't want to introduce any API dependences on anything
// outside of core/ beyond what is absolutely necessary (these
// two functions).  And it is unlikely that the signature of 
// these functions will ever change.
int register_algs_mesh_ops();
int register_extern_mesh_ops();

// Call functions to register MeshOp classes.  These two variables
// are checked in the MKCore constructor soley to ensure that they
// are not optimized away.
const int have_algs_mesh_ops = register_algs_mesh_ops();
const int have_extern_mesh_ops = register_extern_mesh_ops();

bool MeshOpProxy::can_mesh( ModelEnt* entity ) const
  { return true; }

MKCore::MKCore(iGeom *igeom, moab::Interface *moab, iMesh *imesh, iRel *irel,
               bool construct_missing_ifaces) 
        : vertexMesher(NULL)
{
    // This error should never happen.  The primary purpose of this
    // check is to use the variables such that they are not optimized
    // away, thus ensuring that registration actually happens.
  if (!have_algs_mesh_ops || !have_extern_mesh_ops)
    throw Error(MK_MESHOP_NOT_FOUND, "Registration of mesh opts did not happen.");

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
  //mkGraph.addArc(rootNode->get_node(), leafNode->get_node());

  init(construct_missing_ifaces);
  
  for (int i = 0; i < 4; ++i)
    defaultMeshOps[i] = 0;
}

MKCore::~MKCore() 
{
    // delete the graphnodes here, since they point back to me and depend on me being still here
  clear_graph();
  
  // delete the model entities, otherwise they would leak memory
  delete_model_entities();

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
    err = iRelInstances[0]->createPair(iGeomInstances[0]->instance(),
                                       iRel::ENTITY, iRel::IGEOM_IFACE,
                                       iRel::ACTIVE,
                                       iMeshInstances[0]->instance(),
                                       iRel::SET, iRel::IMESH_IFACE,
                                       iRel::ACTIVE,
                                       iRelPairs[0]);
    IBERRCHK(err, "Failure to create relation pair.");
      // don't need to keep track of whether I created the pair, since it'll be deleted anyway when
      // the iRel instance is deleted.

      // FIXME: need a better scheme for finding any existing relation pairs or inferring them from 
      // imported model(s)
  }
  
  if (groupSetPairs.empty() && !iRelInstances.empty() && !iGeomInstances.empty() && !iMeshInstances.empty()) {
    groupSetPairs.resize(1);
    err = iRelInstances[0]->createPair(iGeomInstances[0]->instance(),
                                       iRel::SET, iRel::IGEOM_IFACE,
                                       iRel::ACTIVE,
                                       iMeshInstances[0]->instance(),
                                       iRel::SET, iRel::IMESH_IFACE,
                                       iRel::ACTIVE,
                                       groupSetPairs[0]);
    IBERRCHK(err, "Failure to create relation pair.");
  }

  if (iGeomModelTags.empty()) {
    iGeomModelTags.resize(1);
    // change the name of ModelEntity pointer tag, for iGeom, to
    // differentiate from the tag with the same name in Moab
    // for mesh-based geometry, the model ent tag from moab will conflict
    //  with igeom tag !!! "__MKModelEntityGeo" != "__MKModelEntity"
    err = iGeomInstances[0]->createTag("__MKModelEntityGeo", sizeof(MeshKit::ModelEnt*), iBase_BYTES,
                                       iGeomModelTags[0]);
    IBERRCHK(err, "Failure to create MKModelEnt tag in iGeom.");
  }

  moab::ErrorCode rval;
  if (moabModelTags.empty()) {
    ModelEnt *null_me = NULL;
    moabModelTags.resize(1);
    //rval = moabInstances[0]->tag_create("__MKModelEntity", sizeof(MeshKit::ModelEnt*), moab::MB_TAG_SPARSE,
    //                                    moab::MB_TYPE_OPAQUE, moabModelTags[0], &null_me);
    rval = moabInstances[0]->tag_get_handle("__MKModelEntity", sizeof(MeshKit::ModelEnt*), moab::MB_TYPE_OPAQUE, 
                                            moabModelTags[0], moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT, &null_me);
    if (moab::MB_SUCCESS != rval) 
      MBERRCHK(rval, moab_instance());
  }

  if (moabGeomDimTags.empty()) {
      // moab geometry dimension tag
    moabGeomDimTags.resize(1);
    //rval = moabInstances[0]->tag_create(GEOM_DIMENSION_TAG_NAME, sizeof(int), moab::MB_TAG_SPARSE,
    //                                    moab::MB_TYPE_INTEGER, moabGeomDimTags[0], 0, true);
    rval = moabInstances[0]->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, moabGeomDimTags[0], 
                                            moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
    if (moab::MB_SUCCESS != rval) 
      MBERRCHK(rval, moab_instance());
  }
  
  // moab global id tag
  if (moabIDTags.empty()) {
    moabIDTags.resize(1);
    rval = moabInstances[0]->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, moabIDTags[0], 
                                            moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
    if (moab::MB_SUCCESS != rval) 
      MBERRCHK(rval, moab_instance());
  }
}

unsigned int MKCore::add_igeom_instance(iGeom * igeom)
{
  iBase_ErrorType err;
  iGeomInstances.push_back(igeom);
  unsigned int index = iGeomInstances.size()-1;
  iGeomModelTags.resize(index+1);
  // add the iGeomModelTag
  // it could exist already
  iGeom::TagHandle mkmodeltag;
  err = igeom->getTagHandle("__MKModelEntityGeo", mkmodeltag);
  if (iBase_SUCCESS == err)
  {
    iGeomModelTags[index]= mkmodeltag;
  }
  else
  {
    err = igeom->createTag("__MKModelEntityGeo", sizeof(MeshKit::ModelEnt*), iBase_BYTES,
                                             iGeomModelTags[index]);
    IBERRCHK(err, "Failure to create MKModelEnt tag in iGeom.");
  }


  return index;
}

void MKCore::populate_model_ents(int geom_index, 
                                 int mesh_index,
                                 int irel_index, bool only_geom)
{
    // go through the geometry and mesh models (including geometry groups), and make sure they 
    // all have corresponding model entities; this 

  std::vector<iGeom::EntityHandle> ents;
  std::vector<ModelEnt*> new_ents;
  iBase_ErrorType err;
  ModelEnt *this_me;

  if (-1 != geom_index) {
    for (int dim = iBase_VERTEX; dim <= iBase_REGION; dim++) {
        // get geometry entities
      ents.clear();
      err = igeom_instance(geom_index)->getEntities(igeom_instance(geom_index)->getRootSet(), 
                                                    (iBase_EntityType)dim, ents);
      IBERRCHK(err, "Failed to get entities from iGeom.");

      for (std::vector<iGeom::EntityHandle>::iterator eit = ents.begin(); eit != ents.end(); eit++) {
          // get the modelent
        this_me = NULL;
        err = igeom_instance(geom_index)->getData(*eit, iGeomModelTags[geom_index], &this_me);
        if (NULL == this_me || iBase_TAG_NOT_FOUND == err) {
            // construct a new ModelEnt and set the geom ent to point to it
          this_me = new ModelEnt(this, *eit, geom_index, NULL, mesh_index, irel_index);
          modelEnts[dim].push_back(this_me);

          //check whether there is the mesh related to this model entity
          if (irel_index != -1) {
            iMesh::EntitySetHandle msets;	
            err = irel_pair(irel_index)->getEntSetRelation(this_me->geom_handle(), 0, msets);
            IBERRCHK(err, "Failed to get imesh entity set from model entity.");
            int num = -1;
            err = imesh_instance(mesh_index)->getNumOfType(msets, (iBase_EntityType)dim, num);
            IBERRCHK(err, "Failed to get the mesh entities.");
            if (num > 0) {//there is the mesh related to this model entities
              this_me->set_meshed_state(COMPLETE_MESH);
            }
          }
        }
      }
    }
    
    std::vector<iGeom::EntitySetHandle> gsets;
    err = igeom_instance(geom_index)->getEntSets(igeom_instance(geom_index)->getRootSet(), -1, gsets);
    IBERRCHK(err, "Failed to get entity sets from iGeom.");
    for (std::vector<iGeom::EntitySetHandle>::iterator vit = gsets.begin();
         vit != gsets.end(); vit++) {
      this_me = NULL;
      err = igeom_instance(geom_index)->getEntSetData(*vit, iGeomModelTags[geom_index], &this_me);
      if (NULL == this_me || iBase_TAG_NOT_FOUND == err) {
          // construct a new ModelEnt and set the geom ent to point to it
        this_me = new ModelEnt(this, *vit, geom_index, NULL, mesh_index, irel_index);
        modelEnts[4].push_back(this_me);
      }
    }
  }

  if (-1 != mesh_index && !only_geom) {
    moab::Tag dum_tag = moab_geom_dim_tag(mesh_index);
    moab::Range geom_sets;
    moab::ErrorCode rval = moab_instance(mesh_index)->get_entities_by_type_and_tag(0, MBENTITYSET, &dum_tag, NULL, 1,
                                                                                   geom_sets);
    MBERRCHK(rval, moab_instance(mesh_index));
    for (moab::Range::iterator rit = geom_sets.begin(); rit != geom_sets.end(); rit++) {
        // get the modelent
      this_me = NULL;
      rval = moab_instance(mesh_index)->tag_get_data(moabModelTags[mesh_index], &(*rit), 1, &this_me);
      if (NULL == this_me || MB_TAG_NOT_FOUND == rval) {
            // construct a new ModelEnt and set the geom ent to point to it
        this_me = new ModelEnt(this, (iGeom::EntityHandle)NULL, geom_index, *rit, mesh_index, irel_index);
        int dim = this_me->dimension();
        modelEnts[dim].push_back(this_me);
      }
    }
  }
  // if there are no model ents after parsing through geometry entities or mesh sets, create at least one
  // model entity, with the existing mesh inside, and with a dimension decided upon the highest degree
  // this is useful for loading a mesh for qslim, if there are no sets available to start with
  if (-1 != mesh_index) {
  // first test if there are no ments so far
    int nments = 0;
    for (int k=0; k<5; k++)
    {
      nments += modelEnts[k].size();
    }
    if (nments > 0)
      return; // do not force creation of a model ent, we have at least one and we should be fine
    // put all ents we have so far in one set
    // we know that none of the sets we have so far is a geo set (we would have gotten at least
    // one model ent)
    // so, get all sets first
    moab::ErrorCode rval;
    for (int dimension=3; dimension >=1; dimension--)
    {
      moab::Range ents;
      rval = moab_instance(mesh_index)->get_entities_by_dimension(0, dimension, ents);
      MBERRCHK(rval, moab_instance(mesh_index));
      if (ents.empty())
        continue;
      // put all ents in one moab set, and create one ModelEnt with it, of given dimension
      moab::EntityHandle newSet;
      rval = moab_instance(mesh_index)->create_meshset(moab::MESHSET_SET, newSet);
      MBERRCHK(rval, moab_instance(mesh_index));
      rval = moab_instance(mesh_index)->add_entities(newSet, ents);
      MBERRCHK(rval, moab_instance(mesh_index));
      // also, set the geom dimension
      rval = moab_instance(mesh_index)->tag_set_data(moab_geom_dim_tag(mesh_index), &newSet, 1, &dimension);
      MBERRCHK(rval, moab_instance(mesh_index));
      this_me = new ModelEnt(this, (iGeom::EntityHandle)NULL, geom_index, newSet, mesh_index, irel_index);
      modelEnts[dimension].push_back(this_me);
      return;// get out if we created one model ent
    }
  }
}

void MKCore::load_geometry_mesh(const char *geom_filename, 
                                const char *mesh_filename, 
                                const char *geom_options,
                                const char *mesh_options,
                                int geom_index, 
                                int mesh_index, 
                                int irel_index, 
                                bool relate_too,
                                bool populate_too) 
{
  if (geom_filename) {
    iBase_ErrorType err = igeom_instance(geom_index)->load(geom_filename, geom_options);
    IBERRCHK(err, "Failed to load geometry model.");
  }
  
  if (mesh_filename) {
    moab::ErrorCode rval = moabInstances[mesh_index]->load_file(mesh_filename, NULL, mesh_options);
    MBERRCHK(rval, moab_instance());
  }
  
  if (relate_too) {
    assert(irel_pair(irel_index));
    iRel::Error err = irel_pair(irel_index)->inferAllRelations();
    IBERRCHK(err, "Failed to infer relations.");
  }
  
  if (populate_too) {
    populate_model_ents(geom_index, mesh_index, irel_index);
  }
}

void MKCore::load_geometry(const char *filename, const char *options, int geom_index, 
                           int mesh_index, int irel_index, bool relate_too, bool populate_too) 
{
  iBase_ErrorType err = igeom_instance(geom_index)->load(filename, options);
  IBERRCHK(err, "Failed to load geometry model.");

  if (relate_too) {
    assert(irel_pair(irel_index));
    iRel::Error err = irel_pair(irel_index)->inferAllRelations();
    IBERRCHK(err, "Failed to infer relations.");
  }
  
  if (populate_too) populate_model_ents(geom_index, mesh_index, irel_index);
}

void MKCore::load_mesh(const char *filename, const char *options, 
                       int geom_index, int mesh_index, int irel_index, 
                       bool relate_too, bool populate_too) 
{
  moab::ErrorCode rval = moabInstances[mesh_index]->load_file(filename, NULL, options);
  MBERRCHK(rval, moab_instance(mesh_index));

  if (relate_too) {
    assert(irel_pair(irel_index));
    iRel::Error err = irel_pair(irel_index)->inferAllRelations();
    IBERRCHK(err, "Failed to infer relations.");
  }
  
  if (populate_too) populate_model_ents(geom_index, mesh_index, irel_index);
}

void MKCore::save_geometry(const char *filename, const char *options, int geom_index) 
{
  iBase_ErrorType err = igeom_instance(geom_index)->save(filename, options);
  IBERRCHK(err, "Failed to save geometry model.");
}

void MKCore::save_mesh(const char *filename, const char *options, int index)
{
  moab::ErrorCode rval = moab_instance(index)->write_file(filename, NULL, options);
  MBERRCHK(rval, moab_instance(index));
}

void MKCore::save_mesh_from_model_ents(const char *filename,
    MEntVector & ments, const char *options, int index)
{
  // first, extract moab ent sets from the ments vector
  moab::Range output_sets;
  for (unsigned int i=0; i<ments.size(); i++)
  {
    moab::EntityHandle mset;
    if ((mset = ments[i]->mesh_handle()))
      output_sets.insert(mset);
  }
  moab::ErrorCode rval = moab_instance(index)->write_file(
      filename, NULL, options, output_sets);
  MBERRCHK(rval, moab_instance(index));
}

void MKCore::get_entities_by_dimension(int dim, MEntVector &model_ents,
    int igindx)
{
  int start = dim, end = dim;
  if (iBase_ALL_TYPES == dim) {
    start = 0;
    end = iBase_REGION;
  }
  for (dim = start; dim <= end; dim++)
  {
    if (igindx < 0)
    {

        std::copy(modelEnts[dim].begin(), modelEnts[dim].end(), std::back_inserter(model_ents));
    }
    else
    {
      // retrieve only the model ents with the desired igeom index
      for (unsigned int i=0; i<modelEnts[dim].size(); i++)
      {
        ModelEnt * ment = modelEnts[dim][i];
        if (ment->iGeomIndex() == igindx)
        {
          model_ents.push_back(ment);
        }
      }
    }
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

void MKCore::remove_sizing_function(int index, bool delete_too) 
{
  if (index >= (int)sizingFunctions.size() || !sizingFunctions[index]) 
    throw Error(MK_BAD_INPUT, "No sizing function with that index.");

  if (delete_too) delete sizingFunctions[index];
  
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
  
void MKCore::register_meshop(MeshOpProxy* proxy) 
{
  MeshOpSet::instance().register_mesh_op(proxy);
}
  
MeshOpProxy* MKCore::meshop_proxy(const char *op_name) 
{
  return MeshOpSet::instance().mesh_op(op_name);
}
  
MeshOpProxy* MKCore::meshop_proxy(unsigned index) 
{
  return MeshOpSet::instance().mesh_op(index);
}

unsigned MKCore::num_meshops()
{
  return MeshOpSet::instance().mesh_ops().size();
}
  
unsigned int MKCore::meshop_index(const char *op_name) 
{
  return MeshOpSet::instance().index(op_name);
}
  
void MKCore::set_default_meshop(const char *op_name, unsigned short dims) 
{
  set_default_meshop( meshop_proxy( op_name ), dims );
}
  
void MKCore::set_default_meshop(unsigned index, unsigned short dims) 
{
  set_default_meshop( meshop_proxy( index ), dims );
}

MeshOpProxy* MKCore::get_default_meshop( unsigned dimension )
{
    // If default hasn't been explicity set, set it to the first
    // available algorithm
  if (!defaultMeshOps[dimension]) {
    const MeshOpSet::OpList& list = MeshOpSet::instance().mesh_ops( dimension );
    if (list.empty()) 
      throw Error(MK_NOT_FOUND, "No MeshOp available for dimension %u", dimension);
    defaultMeshOps[dimension] = list.front();
  }
  
  return defaultMeshOps[dimension];
}
    
void MKCore::set_default_meshop(MeshOpProxy* mesh_op, unsigned short dims) 
{
    // check the specified dimension(s) against the types the meshop can mesh
  for (unsigned i = 0; i < 4; ++i)
    if ((dims & (1u << i)) && !mesh_op->can_mesh(iBase_EntityType(i)))
      throw Error(MK_BAD_INPUT, "Specified MeshOp type cannot generate elements of specified dimension.");
  
    // set as default for specified dimensions
  for (unsigned i = 0; i < 4; ++i)
    if (dims & (1u << i))
      defaultMeshOps[i] = mesh_op;
}

void MKCore::meshop_by_mesh_type(moab::EntityType tp, std::vector<MeshOpProxy*> &ops) 
{
  unsigned dim = moab::CN::Dimension(tp);
  const MeshOpSet::OpList& list = MeshOpSet::instance().mesh_ops( dim );
  for (MeshOpSet::iterator i = list.begin(); i != list.end(); ++i) {
    const moab::EntityType* list = (*i)->output_types();
    for (int j = 0; list[j] != moab::MBMAXTYPE; ++j) {
      if (list[j] == tp) {
        ops.push_back(*i);
        break;
      }
    }
  }
}
    
void MKCore::meshop_by_dimension(int dim, std::vector<MeshOpProxy*> &ops) 
{
  const MeshOpSet::OpList& list = MeshOpSet::instance().mesh_ops( dim );
  std::copy( list.begin(), list.end(), std::back_inserter(ops) );
}
    
void MKCore::meshop_by_modelent(ModelEnt * const ent, std::vector<MeshOpProxy*> &ops) 
{
  const MeshOpSet::OpList& list = MeshOpSet::instance().mesh_ops( ent->dimension() );
  for (MeshOpSet::iterator i = list.begin(); i != list.end(); ++i) 
    if ((*i)->can_mesh(ent))
      ops.push_back(*i);
}
    
MeshOp *MKCore::construct_meshop(std::string op_name, const MEntVector &me_vec) 
{
  return construct_meshop( meshop_proxy(op_name.c_str()), me_vec );
}

MeshOp *MKCore::construct_meshop(unsigned int dim, const MEntVector &me_vec) 
{
  if (dim == 0) {
    if (!vertexMesher) {
      vertexMesher = new VertexMesher(this, me_vec);
      vertexMesher->set_name("VertexMesher");
    }
    else 
      for (MEntVector::const_iterator i = me_vec.begin(); i != me_vec.end(); ++i)
        vertex_mesher()->add_modelent( *i );
    return vertex_mesher();
  }
  else 
    return construct_meshop( get_default_meshop(dim), me_vec );
}

MeshOp *MKCore::construct_meshop(MeshOpProxy* proxy, const MEntVector &me_vec) 
{
  return proxy->create( this, me_vec );
}
#if HAVE_FBIGEOM
// this method takes the model from moab modelRootSet (0 by default means all)
// also, assumes the first moab instance from MKCore, and the first iMesh instance
int MKCore::initialize_mesh_based_geometry(moab::EntityHandle modelRootSet)
{
  // initialize the FBiGeom, in a different way
  iMesh * imeshMK = imesh_instance(); // this is MeshKit::iMesh class
  iMesh_Instance mesh = imeshMK->instance();
  FBiGeom_Instance geom;

  iBase_EntitySetHandle root_set = reinterpret_cast<iBase_EntitySetHandle>(modelRootSet);
  std::string opts("SMOOTH;");
     // new constructor
  int err;
  // this does the initialization of smoothing (heavy duty computation)
  FBiGeom_newGeomFromMesh(mesh, root_set, opts.c_str(), &geom, &err, opts.length());
  if (err!=0)
    return -1; // error
  FBiGeom * fbigeom = new FBiGeom(geom);

  int index = add_igeom_instance((iGeom*)fbigeom);
  // use this index to populate
  populate_model_ents(index, 0, -1, true);
  return index;

}

void MKCore::remove_mesh_based_geometry(int index)
{
  // first, check if the index does make sense
  if (index <0 || index >= (int)iGeomInstances.size())
    throw Error(MK_BAD_INPUT, "Specified index too big for FBiGeom.");
  FBiGeom * fbIGeom = reinterpret_cast<FBiGeom*>(iGeomInstances[index]);
  unsigned int sz = iGeomInstances.size();
  delete fbIGeom;
  fbIGeom = NULL;
  for (unsigned int i=index; i<sz-1; i++)
  {
    iGeomInstances[i] = iGeomInstances[i+1];
  }
  iGeomInstances.resize(sz-1);
}
#endif

void MKCore::delete_model_entities()
{

  std::vector<ModelEnt*> tmp_ents;
  for (int dim = 0; dim <= iBase_REGION; dim++) {
    unsigned int sz = modelEnts[dim].size();
    for (unsigned int i=0; i<sz; i++)
       delete modelEnts[dim][i];
    modelEnts[dim].clear();
  }
}

void MKCore::create_mbg_model_entities(moab::EntityHandle modelRootSet, bool geometry)
{
  if (geometry)
  {
    //
  }
  else
  {
    moab::Tag dum_tag = moab_geom_dim_tag();// default 0, no issues here
    moab::Range geom_sets;
    // the enclosing set is the only difference from populate method above
    moab::ErrorCode rval = moab_instance()->get_entities_by_type_and_tag(modelRootSet, MBENTITYSET, &dum_tag, NULL, 1,
                                                                                   geom_sets);
    MBERRCHK(rval, moab_instance());
    for (moab::Range::iterator rit = geom_sets.begin(); rit != geom_sets.end(); rit++) {
      // get the modelent

          // construct a new ModelEnt and set the geom ent to point to it
      // assume mesh index 0, as always, and no relations yet
      ModelEnt * this_me = new ModelEnt(this, (iGeom::EntityHandle)NULL, 0, *rit, 0, -1);
      int dim = this_me->dimension();// will get it from moab set, actually
      modelEnts[dim].push_back(this_me);
    }
  }
}
} // namespace meshkit
