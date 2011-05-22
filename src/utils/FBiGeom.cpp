
/** \file FBiGeom.cpp;
 *
 */

#include "meshkit/FBiGeom.hpp"
#include "meshkit/MKCore.hpp"

namespace MeshKit {

FBiGeom::FBiGeom( MKCore * mk, bool smooth):
    iGeom(true), _fbEngine(NULL), _smooth(smooth) , _mk(mk), _index(0)
{
  // use the first moab instance, as a basis for
  //_mbImpl = _mk->moabInstances[0];
  //_index = _mk->add_igeom_instance((iGeom*)this);
}

FBiGeom::~FBiGeom()
{
  if (_fbEngine) {
    if (_smooth)
      _fbEngine->delete_smooth_tags();
    delete _fbEngine;
  }
}
iGeom::Error FBiGeom::Init()
{
  // what will happen here:
  // create facet-based engine, smooth option
  _fbEngine= new moab::FBEngine(_mk->moab_instance(), NULL, true); //
  if (!_fbEngine)
    return (Error)1;
  _index = _mk->add_igeom_instance((iGeom*)this);

  _fbEngine->Init();// this will trigger the initialization
  // of smoothing
  // now, we need to create the ModelEnts, use populate_mesh !!!
  _mk->populate_model_ents(meshkit_index(), 0, -1);// no irel !!
  return (Error)0; // success
}

inline iGeom::Error
FBiGeom::load( const char* file_name,
             const char* options )
{
  //moab::ErrorCode rval = _mbImpl->load_file(file_name, 0, options);
  // will load the file as a moab database, and populate it as needed
  moab::ErrorCode rval = _mk->moab_instance()->load_file(file_name, 0, options);
  if (rval!=moab::MB_SUCCESS)
    return (Error)1;

  // the file that is loaded in the moab database
  // right now, assume only one file is read into
  // for multiple files, we will have to add them to different sets
  // in the moab database
  Init();
  return (Error)0;
}
// include here the iBaseVirtual methods that will be overloaded
// start copy iBaseVirtual

iGeom::EntitySetHandle FBiGeom::getRootSet()
{
  // the root set in MOAB is always null
  // we will return the set from FBEngine
  moab::EntityHandle modelSet;
  _fbEngine->getRootSet(&modelSet);
  return (iGeom::EntitySetHandle)modelSet;
}

iGeom::Error FBiGeom::getEntType( EntityHandle handle, EntityType& type_out )
{
  int type;
  moab::ErrorCode rval = _fbEngine-> getEntType( (moab::EntityHandle)handle, &type);
  type_out = (EntityType)type;
  if (rval == moab::MB_SUCCESS)
    return (Error)0;
  return (Error)1;
}


#if 0

PFX(Base)::Error
PFX(Base)::createEntSet( bool is_list, EntitySetHandle& handle_out )
{
  int err;
  PFX(_createEntSet)( mInstance, is_list, &handle_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::destroyEntSet( EntitySetHandle handle )
{
  int err;
  PFX(_destroyEntSet)( mInstance, handle, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::isList( EntitySetHandle handle, bool& is_list )
{
  int err, result;
  PFX(_isList)( mInstance, handle, &result, &err );
  is_list = (result != 0);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getNumEntSets( EntitySetHandle set, int num_hops, int& num_sets_out )
{
  int err;
  PFX(_getNumEntSets)( mInstance, set, num_hops, &num_sets_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getEntSets( EntitySetHandle set, int num_hops,
                   std::vector<EntitySetHandle>& contained_sets_out )
{
  int err, count;
  PFX(_getNumEntSets)( mInstance, set, num_hops, &count, &err );
  if (iBase_SUCCESS != err)
    return (Error)err;
  contained_sets_out.resize(count);
  int alloc = contained_sets_out.size(), size;
  EntitySetHandle* ptr = &contained_sets_out[0];
  PFX(_getEntSets)( mInstance, set, num_hops, &ptr, &alloc, &size, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::addEntToSet( EntityHandle entity, EntitySetHandle set )
{
  int err;
  PFX(_addEntToSet)( mInstance, entity,set, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::rmvEntFromSet( EntityHandle entity, EntitySetHandle set )
{
  int err;
  PFX(_rmvEntFromSet)( mInstance, entity,set, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::addEntArrToSet( const EntityHandle* entity_handles,
                             int entity_handles_size,
                             EntitySetHandle entity_set )
{
  int err;
  PFX(_addEntArrToSet)( mInstance, entity_handles, entity_handles_size, entity_set, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::rmvEntArrFromSet( const EntityHandle* entity_handles,
                               int entity_handles_size,
                               EntitySetHandle entity_set )
{
  int err;
  PFX(_rmvEntArrFromSet)( mInstance, entity_handles, entity_handles_size, entity_set, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::addEntSet( EntitySetHandle to_add, EntitySetHandle add_to )
{
  int err;
  PFX(_addEntSet)( mInstance, to_add, add_to, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::rmvEntSet( EntitySetHandle to_rmv, EntitySetHandle rmv_from )
{
  int err;
  PFX(_rmvEntSet)( mInstance, to_rmv, rmv_from, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::isEntContained( EntitySetHandle set, EntityHandle ent, bool& contained_out )
{
  int err, result;
  PFX(_isEntContained)( mInstance, set, ent, &result, &err );
  contained_out = (result != 0);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::isEntArrContained( EntitySetHandle containing_set,
                                const EntityHandle* entity_handles,
                                int num_entity_handles,
                                bool* is_contained_out )
{
  int err, *ptr = 0, alloc = 0, size = 0;
  PFX(_isEntArrContained)( mInstance, containing_set,
                           entity_handles, num_entity_handles,
                           &ptr, &alloc, &size, &err );
  if (iBase_SUCCESS != err)
    return (Error)err;
  for (int i = 0; i < num_entity_handles; ++i)
    is_contained_out[i] = (ptr[i] != 0);
  free(ptr);
  return iBase_SUCCESS;
}

PFX(Base)::Error
PFX(Base)::isEntSetContained( EntitySetHandle containing_set,
                          EntitySetHandle contained_set,
                          bool& contained_out )
{
  int err, result;
  PFX(_isEntSetContained)( mInstance, containing_set, contained_set, &result, &err );
  contained_out = (result != 0);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::addPrntChld( EntitySetHandle parent, EntitySetHandle child )
{
  int err;
  PFX(_addPrntChld)( mInstance, parent, child, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::rmvPrntChld( EntitySetHandle parent, EntitySetHandle child )
{
  int err;
  PFX(_rmvPrntChld)( mInstance, parent, child, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::isChildOf( EntitySetHandle parent, EntitySetHandle child, bool& is_child_out )
{
  int err, result;
  PFX(_isChildOf)( mInstance, parent, child, &result, &err );
  is_child_out = (result != 0);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getNumChld( EntitySetHandle parent, int num_hops, int& num_child_out )
{
  int err;
  PFX(_getNumChld)( mInstance, parent, num_hops, &num_child_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getNumPrnt( EntitySetHandle child, int num_hops, int& num_parent_out )
{
  int err;
  PFX(_getNumPrnt)( mInstance, child, num_hops, &num_parent_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getChldn( EntitySetHandle parent, int num_hops,
                 std::vector<EntitySetHandle>& children_out )
{
  int err, count;
  PFX(_getNumChld)( mInstance, parent, num_hops, &count, &err );
  if (iBase_SUCCESS != err)
    return (Error)err;
  children_out.resize(count);
  int alloc = children_out.size(), size;
  EntitySetHandle* ptr = &children_out[0];
  PFX(_getEntSets)( mInstance, parent, num_hops, &ptr, &alloc, &size, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getPrnts( EntitySetHandle child, int num_hops,
                 std::vector<EntitySetHandle>& parents_out )
{
  int err, count;
  PFX(_getNumPrnt)( mInstance, child, num_hops, &count, &err );
  if (iBase_SUCCESS != err)
    return (Error)err;
  parents_out.resize(count);
  int alloc = parents_out.size(), size;
  EntitySetHandle* ptr = &parents_out[0];
  PFX(_getEntSets)( mInstance, child, num_hops, &ptr, &alloc, &size, &err );
  return (Error)err;
}


PFX(Base)::Error
PFX(Base)::subtract( EntitySetHandle set1, EntitySetHandle set2,
                 EntitySetHandle& result_set_out )
{
  int err;
  PFX(_subtract)( mInstance, set1, set1, &result_set_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::intersect( EntitySetHandle set1, EntitySetHandle set2,
                  EntitySetHandle& result_set_out )
{
  int err;
  PFX(_intersect)( mInstance, set1, set1, &result_set_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::unite( EntitySetHandle set1, EntitySetHandle set2,
              EntitySetHandle& result_set_out )
{
  int err;
  PFX(_unite)( mInstance, set1, set1, &result_set_out, &err );
  return (Error)err;
}


PFX(Base)::Error
PFX(Base)::createTag( const char* tag_name,
                  int tag_num_type_values,
                  TagValueType tag_type,
                  TagHandle& tag_handle_out )
{
  int err;
  PFX(_createTag)( mInstance, tag_name, tag_num_type_values, tag_type,
                   &tag_handle_out, &err, strlen(tag_name) );
  return (Error)err;
}


PFX(Base)::Error
PFX(Base)::destroyTag( TagHandle tag_handle, bool forced )
{
  int err;
  PFX(_destroyTag)( mInstance, tag_handle, forced, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getTagName( TagHandle tag_handle, std::string& name_out )
{
  int err;
  char buffer[1024];
  memset( buffer, 0, sizeof(buffer) );
  PFX(_getTagName)( mInstance, tag_handle, buffer, &err, sizeof(buffer) );
  name_out = buffer;
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getTagSizeValues( TagHandle tag_handle, int& size_out )
{
  int err;
  PFX(_getTagSizeValues)( mInstance, tag_handle, &size_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getTagSizeBytes( TagHandle tag_handle, int& size_out )
{
  int err;
  PFX(_getTagSizeBytes)( mInstance, tag_handle, &size_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getTagHandle( const char* name, TagHandle& handle_out )
{
  int err;
  PFX(_getTagHandle)( mInstance, name, &handle_out, &err, strlen(name) );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getTagType( TagHandle tag_handle, TagValueType& type_out )
{
  int err, result;
  PFX(_getTagType)( mInstance, tag_handle, &result, &err );
  type_out = (TagValueType)result;
  return (Error)err;
}


PFX(Base)::Error
PFX(Base)::setEntSetData( EntitySetHandle set_handle,
                      TagHandle tag_handle,
                      const void* tag_value )
{
  int err, size = 1;
  PFX(_getTagSizeBytes)( mInstance, tag_handle, &size, &err );
  PFX(_setEntSetData)( mInstance, set_handle, tag_handle,
                       (const char*)tag_value, size, &err);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setEntSetIntData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         int value )
{
  int err;
  PFX(_setEntSetIntData)( mInstance, set_handle, tag_handle, value, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setEntSetDblData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         double value )
{
  int err;
  PFX(_setEntSetDblData)( mInstance, set_handle, tag_handle, value, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setEntSetEHData( EntitySetHandle set_handle,
                        TagHandle tag_handle,
                        EntityHandle value )

{
  int err;
  PFX(_setEntSetEHData)( mInstance, set_handle, tag_handle, value, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getEntSetData( EntitySetHandle set_handle,
                      TagHandle tag_handle,
                      void* tag_value_out )
{
  int err, alloc = std::numeric_limits<int>::max(), size;
  char* ptr = (char*)tag_value_out;
  PFX(_getEntSetData)( mInstance, set_handle, tag_handle,
                      &ptr, &alloc, &size, &err);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getEntSetIntData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         int& value_out )
{
  int err;
  PFX(_getEntSetIntData)( mInstance, set_handle, tag_handle, &value_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getEntSetDblData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         double& value_out )
{
  int err;
  PFX(_getEntSetDblData)( mInstance, set_handle, tag_handle, &value_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getEntSetEHData( EntitySetHandle set_handle,
                        TagHandle tag_handle,
                        EntityHandle& value_out )

{
  int err;
  PFX(_getEntSetEHData)( mInstance, set_handle, tag_handle, &value_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getAllEntSetTags( EntitySetHandle set,
                         std::vector<TagHandle>& tags_out )
{
  if (tags_out.capacity() == 0)
    tags_out.resize( 32 );
  else
    tags_out.resize( tags_out.capacity() );

  int err, alloc = tags_out.size(), size = 0;
  TagHandle* ptr = &tags_out[0];
  PFX(_getAllEntSetTags)( mInstance, set, &ptr, &alloc, &size, &err );
  tags_out.resize(size);

  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = tags_out.size();
    ptr = &tags_out[0];
    PFX(_getAllEntSetTags)( mInstance, set, &ptr, &alloc, &size, &err );
  }

  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getAllTags( EntityHandle entity,
                   std::vector<TagHandle>& tags_out )

{
  if (tags_out.capacity() == 0)
    tags_out.resize( 32 );
  else
    tags_out.resize( tags_out.capacity() );

  int err, alloc = tags_out.size(), size = 0;
  TagHandle* ptr = &tags_out[0];
  PFX(_getAllTags)( mInstance, entity, &ptr, &alloc, &size, &err );
  tags_out.resize(size);

  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = tags_out.size();
    ptr = &tags_out[0];
    PFX(_getAllTags)( mInstance, entity, &ptr, &alloc, &size, &err );
  }

  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::rmvEntSetTag( EntitySetHandle set, TagHandle tag )
{
  int err;
  PFX(_rmvEntSetTag)( mInstance, set, tag, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::rmvTag( EntityHandle entity, TagHandle tag )
{
  int err;
  PFX(_rmvTag)( mInstance, entity, tag, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::rmvArrTag( const EntityHandle* handles, int size, TagHandle tag )
{
  int err;
  PFX(_rmvArrTag)( mInstance, handles, size, tag, &err );
  return (Error)err;
}


PFX(Base)::Error
PFX(Base)::getArrData( const EntityHandle* entity_handles,
                   int entity_handles_size,
                   TagHandle tag_handle,
                   void* tag_values_out )
{
  int err, alloc = std::numeric_limits<int>::max(), size;
  char* ptr = (char*)tag_values_out;
  PFX(_getArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                    &ptr, &alloc, &size, &err);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getIntArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      int* tag_values_out )
{
  int err, alloc = std::numeric_limits<int>::max(), size;
  PFX(_getIntArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                       &tag_values_out, &alloc, &size, &err);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getDblArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      double* tag_values_out )
{
  int err, alloc = std::numeric_limits<int>::max(), size;
  PFX(_getDblArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                       &tag_values_out, &alloc, &size, &err);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getEHArrData( const EntityHandle* entity_handles,
                     int entity_handles_size,
                     TagHandle tag_handle,
                     EntityHandle* tag_values_out )

{
  int err, alloc = std::numeric_limits<int>::max(), size;
  PFX(_getEHArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                      &tag_values_out, &alloc, &size, &err);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setArrData( const EntityHandle* entity_handles,
                   int entity_handles_size,
                   TagHandle tag_handle,
                   const void* tag_values )
{
  int err, size = 1;
  PFX(_getTagSizeBytes)( mInstance, tag_handle, &size, &err );
  PFX(_setArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                    (const char*)tag_values, size*entity_handles_size,
                    &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setIntArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      const int* tag_values )
{
  int err, size = 1;
  PFX(_getTagSizeValues)( mInstance, tag_handle, &size, &err );
  PFX(_setIntArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                       tag_values, size*entity_handles_size, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setDblArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      const double* tag_values )
{
  int err, size = 1;
  PFX(_getTagSizeValues)( mInstance, tag_handle, &size, &err );
  PFX(_setDblArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                       tag_values, size*entity_handles_size, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setEHArrData( const EntityHandle* entity_handles,
                     int entity_handles_size,
                     TagHandle tag_handle,
                     const EntityHandle* tag_values )
{
  int err, size = 1;
  PFX(_getTagSizeValues)( mInstance, tag_handle, &size, &err );
  PFX(_setEHArrData)( mInstance, entity_handles, entity_handles_size, tag_handle,
                      tag_values, size*entity_handles_size, &err );
  return (Error)err;
}



PFX(Base)::Error
PFX(Base)::setData( EntityHandle entity_handle,
                TagHandle tag_handle,
                const void* tag_value )
{
  int err, size = 1;
  PFX(_getTagSizeBytes)( mInstance, tag_handle, &size, &err );
  PFX(_setData)( mInstance, entity_handle, tag_handle,
                 (const char*)tag_value, size, &err);
  return (Error)err;
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setIntData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   int value )
{
  int err;
  PFX(_setIntData)( mInstance, entity_handle, tag_handle, value, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setDblData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   double value )
{
  int err;
  PFX(_setDblData)( mInstance, entity_handle, tag_handle, value, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::setEHData( EntityHandle entity_handle,
                  TagHandle tag_handle,
                  EntityHandle value )

{
  int err;
  PFX(_setEHData)( mInstance, entity_handle, tag_handle, value, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getData( EntityHandle entity_handle,
                TagHandle tag_handle,
                void* tag_value_out )
{
  int err, alloc = std::numeric_limits<int>::max(), size;
  char* ptr = (char*)tag_value_out;
  PFX(_getData)( mInstance, entity_handle, tag_handle,
                &ptr, &alloc, &size, &err);
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getIntData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   int& value_out )
{
  int err;
  PFX(_getIntData)( mInstance, entity_handle, tag_handle, &value_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getDblData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   double& value_out )
{
  int err;
  PFX(_getDblData)( mInstance, entity_handle, tag_handle, &value_out, &err );
  return (Error)err;
}

PFX(Base)::Error
PFX(Base)::getEHData( EntityHandle entity_handle,
                  TagHandle tag_handle,
                  EntityHandle& value_out )
{
  int err;
  PFX(_getEHData)( mInstance, entity_handle, tag_handle, &value_out, &err );
  return (Error)err;
}
#endif

// end copy iBaseVirtual

// no inlining anymore
iGeom::Error FBiGeom::getEntities( EntitySetHandle set,
                    EntityType type,
                    std::vector<EntityHandle>& entities_out )
{
  // ErrorCode getEntities(EntityHandle root_set, int ent_type, Range & gentities);
  moab::Range gentities;
  moab::ErrorCode rval = _fbEngine->getEntities((moab::EntityHandle)set,  (int) type, gentities);
  if (rval!=moab::MB_SUCCESS)
    return (Error)1;
  entities_out.clear();
  entities_out.resize(gentities.size());
  moab::Range::iterator it = gentities.begin();
  for (unsigned int i=0; i<gentities.size(); i++, it++)
    entities_out[i] = (iGeom::EntityHandle)*it;
  return (Error)0;
}

iGeom::Error FBiGeom::getEntAdj( EntityHandle handle,
                  EntityType type_requested,
                  std::vector<EntityHandle>& adj_entities_out )
{
  moab::Range adjEnts;
  moab::ErrorCode rval = _fbEngine->getEntAdj((moab::EntityHandle)handle,
      (int)type_requested, adjEnts  );
  if (rval==moab::MB_SUCCESS)
  {
    adj_entities_out.clear();
    adj_entities_out.resize(adjEnts.size());
    moab::Range::iterator it = adjEnts.begin();
    for (unsigned int i=0; i<adjEnts.size(); i++, it++)
      adj_entities_out[i] = (iGeom::EntityHandle)*it;
    return (Error)0;
  }
  else
    return (Error)1;
}
iGeom::Error FBiGeom::getEgFcSense( EntityHandle edge, EntityHandle face, int& sense )
{
  moab::EntityHandle mbedge = (moab::EntityHandle) edge;
  moab::EntityHandle mbface = (moab::EntityHandle) face;
  moab::ErrorCode rval = _fbEngine->getEgFcSense( mbedge, mbface, sense );
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}

iGeom::Error FBiGeom::measure( const EntityHandle* entities,
                              int entities_size,
                              double* measures )
{
  const moab::EntityHandle * moab_entities = ( moab::EntityHandle *) entities;

  moab::ErrorCode rval = _fbEngine->measure(moab_entities, entities_size, measures);
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}

iGeom::Error FBiGeom::getEntNrmlSense( EntityHandle face, EntityHandle region,
        int& sense )
{
  moab::EntityHandle mbregion = (moab::EntityHandle) region;
  moab::EntityHandle mbface = (moab::EntityHandle) face;
  moab::ErrorCode rval = _fbEngine->getEgFcSense(  mbface, mbregion, sense );
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}

iGeom::Error FBiGeom::getEntNrmlXYZ( EntityHandle entity,
                                    double x, double y, double z,
                                    double& i, double& j, double& k )
{
  moab::EntityHandle mbentity = (moab::EntityHandle) entity;
  moab::ErrorCode rval = _fbEngine->getEntNrmlXYZ( mbentity,
      x,   y,  z,
       &i,  &j,  &k ); // some inconsistency here, we use pointers, not refs
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}

iGeom::Error FBiGeom::getVtxCoord( EntityHandle vertex,
                                  double& x, double& y, double& z )
{
  moab::EntityHandle mbentity = (moab::EntityHandle) vertex;
  moab::ErrorCode rval = _fbEngine->getVtxCoord(mbentity, &x, &y, &z);
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}

iGeom::Error FBiGeom::getEntClosestPt( EntityHandle entity,
                                      double near_x, double near_y, double near_z,
                                      double& on_x, double& on_y, double& on_z )
{
  moab::EntityHandle mbentity = (moab::EntityHandle) entity;
  moab::ErrorCode rval = _fbEngine->getEntClosestPt(mbentity,
      near_x,  near_y,   near_z,
        &on_x,  &on_y, &on_z );
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}
iGeom::Error FBiGeom::getEgEvalXYZ( EntityHandle edge,
                                   double x, double y, double z,
                                   double& on_x, double& on_y, double& on_z,
                                   double& tngt_i, double& tngt_j, double& tngt_k,
                                   double& cvtr_i, double& cvtr_j, double& cvtr_k )
{
  moab::ErrorCode rval = _fbEngine->getEgEvalXYZ( (moab::EntityHandle) edge,
                                   x, y, z,
                                   on_x, on_y, on_z,
                                   tngt_i, tngt_j, tngt_k,
                                   cvtr_i,cvtr_j, cvtr_k );
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}
iGeom::Error FBiGeom::getFcEvalXYZ( EntityHandle face,
                                   double x, double y, double z,
                                   double& on_x, double& on_y, double& on_z,
                                   double& nrml_i, double& nrml_j, double& nrml_k,
                                   double& cvtr1_i, double& cvtr1_j, double& cvtr1_k,
                                   double& cvtr2_i, double& cvtr2_j, double& cvtr2_k )
{
  /*
  moab::ErrorCode rval = _fbEngine->getFcEvalXYZ( (moab::EntityHandle) face,
                            x,  y, z,
                            on_x, on_y, on_z,
                            nrml_i, nrml_j, nrml_k,
                            cvtr1_i, cvtr1_j, cvtr1_k,
                            cvtr2_i, cvtr2_j,  cvtr2_k );*/
  // camal really does not use curvatures
  // the most it is calling for normals and for closest point
  // return all curvatures = 0 for the time being, because we
  // know camal is not requesting them

  cvtr1_i=cvtr1_j=cvtr1_k= cvtr2_i= cvtr2_j=  cvtr2_k=0.;
  // first get closest point, then normal, separately
  moab::ErrorCode rval = _fbEngine->getEntClosestPt((moab::EntityHandle) face,
       x, y, z,
         &on_x,  &on_y, &on_z );
  if (rval != moab::MB_SUCCESS)
    return (Error)1;
  rval = _fbEngine->getEntNrmlXYZ( (moab::EntityHandle) face,
        x,   y,  z,
        &nrml_i, &nrml_j, &nrml_k ); // some inconsistency here, we use pointers, not refs
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;

}

iGeom::Error FBiGeom::getEgVtxSense( EntityHandle edge, EntityHandle vtx1, EntityHandle vtx2, int& sense )
{
  moab::ErrorCode rval = _fbEngine->getEgVtxSense((moab::EntityHandle) edge,
      (moab::EntityHandle) vtx1, (moab::EntityHandle) vtx2,  sense);
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}

iGeom::Error FBiGeom::getEntURange( EntityHandle edge,
                                   double& u_min, double& u_max )
{
  moab::ErrorCode rval = _fbEngine->getEntURange((moab::EntityHandle) edge,
               u_min,  u_max );
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}
iGeom::Error FBiGeom::getEntUtoXYZ( EntityHandle edge, double u,
                                   double& x, double& y, double& z )
{
  moab::ErrorCode rval = _fbEngine->getEntUtoXYZ( (moab::EntityHandle) edge,   u,
                                    x,  y,   z );
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}
iGeom::Error FBiGeom::isEntAdj( EntityHandle entity1, EntityHandle entity2,
        bool& adjacent_out )
{
  moab::ErrorCode rval = _fbEngine->isEntAdj( (moab::EntityHandle) entity1,
      (moab::EntityHandle) entity2,
       adjacent_out );
  if (rval == moab::MB_SUCCESS)
    return (Error) 0;
  else
    return (Error)1;
}

iGeom::Error FBiGeom::getTagHandle( const char* name, TagHandle& handle_out )
{
  moab::Tag handle;
  moab::ErrorCode rval = _fbEngine->moab_instance()->tag_get_handle( name, handle );

  if (rval==moab::MB_SUCCESS)
  {
    handle_out = (TagHandle) handle;
    return (Error)0;
  }
  else
    return (Error)1;
}

iGeom::Error FBiGeom::createTag( const char* tag_name,
                  int tag_num_type_values,
                  TagValueType tag_type,
                  TagHandle& tag_handle_out )
{
  moab::Tag handle;
  // this is not implemented yet
  moab::ErrorCode rval = _fbEngine->createTag(tag_name,
      tag_num_type_values, tag_type, handle );
  if (rval == moab::MB_SUCCESS)
  {
    tag_handle_out = (TagHandle) handle;
    return (Error)0;
  }
  return (Error)1;
}

iGeom::Error FBiGeom::getData( EntityHandle entity_handle,
                TagHandle tag_handle,
                void* tag_value_out )
{
  // we know this is really a moab thingy
  return getArrData(  &entity_handle, 1, tag_handle, tag_value_out );

}

iGeom::Error FBiGeom::getIntData( EntityHandle entity_handle,
                           TagHandle tag_handle,
                           int& value_out )
{
  // we know this is really a moab thingy
   return getArrData(  &entity_handle, 1, tag_handle, (void*)&value_out );
}
iGeom::Error FBiGeom::getArrData( const EntityHandle* entity_handles,
                   int entity_handles_size,
                   TagHandle tag_handle,
                   void* tag_values_out )
{
  moab::Tag tagh=(moab::Tag)tag_handle;
  const moab::EntityHandle *ents = reinterpret_cast<const moab::EntityHandle *>(entity_handles);
  moab::ErrorCode rval =  _fbEngine->getArrData( ents,
      entity_handles_size, tagh, tag_values_out);
  if (rval == moab::MB_SUCCESS)
    return (Error)0;

  return (Error)1;
}

iGeom::Error FBiGeom::setData( EntityHandle entity_handle,
                TagHandle tag_handle,
                const void* tag_value )
{
  // we know this is really a moab thingy
  return setArrData(  &entity_handle, 1, tag_handle, tag_value );
}

iGeom::Error FBiGeom::setArrData( const EntityHandle* entity_handles,
                   int entity_handles_size,
                   TagHandle tag_handle,
                   const void* tag_values )
{
  moab::Tag tagh=(moab::Tag)tag_handle;
  const moab::EntityHandle *ents = reinterpret_cast<const moab::EntityHandle *>(entity_handles);
  moab::ErrorCode rval =  _fbEngine->setArrData( ents,
      entity_handles_size, tagh, tag_values);
  if (rval == moab::MB_SUCCESS)
    return (Error)0;

  return (Error)1;

}

iGeom::Error FBiGeom::getFaceType( EntityHandle face, std::string& type )
{
  type = "nonplanar";
  return (Error)0;
}

iGeom::Error FBiGeom::isEntParametric( EntityHandle entity, bool& parametric )
{
  parametric = false;
  return (Error)0;
}

iGeom::Error FBiGeom::getEntBoundBox( EntityHandle entity,
                       double& min_x, double& min_y, double& min_z,
                       double& max_x, double& max_y, double& max_z )
{
  moab::ErrorCode rval =  _fbEngine->getEntBoundBox((moab::EntityHandle) entity,
      &min_x, &min_y, &min_z,
      &max_x, &max_y, &max_z);
  if (rval == moab::MB_SUCCESS)
    return (Error)0;

  return (Error)1;
}

iGeom::Error FBiGeom::isEntPeriodic( EntityHandle entity, bool& in_u, bool& in_v )
{

  in_u = false; // false
  in_v = false; // false
  EntityType tp_ent;
  Error err = getEntType(entity, tp_ent);

  if (iBase_SUCCESS != err) return err;
  if (tp_ent == iBase_FACE)
    return (Error)0; // we do not have parametric surfaces , so it cannot be periodic
                     // in u or v
  // if it is an edge/curve, it could be periodic in u, if there is
  // only one vertex as neighbor
  if (tp_ent == iBase_EDGE)
  {
    std::vector<EntityHandle> adj_entities_out;
    getEntAdj( entity, iBase_VERTEX, adj_entities_out);
    if (adj_entities_out.size() == 1)
    {
      in_u = true; //
    }
  }
  return (Error)0; // no error
}
#if 0

inline iGeom::Error
iGeom::save( const char* file_name,
             const char* options )
{
  int err, len = options ? strlen(options) : 0;
  iGeom_save( mInstance, file_name, options, &err, strlen(file_name), len );
  return (Error)err;
}

//inline iGeom::StorageOrder 
//iGeom::getDfltStorage()
//{
//  int err, order;
//  iGeom_getDfltStorage( mInstance, &order, &err );
//  return (iBase_SUCCESS == err) ? (StorageOrder)order : iBase_UNDETERMINED;
//}
  
inline iGeom::Error
iGeom::getNumOfType( EntitySetHandle set, EntityType type, int& count_out  )
{
  int err;
  iGeom_getNumOfType( mInstance, set, type, &count_out, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::deleteEnt( EntityHandle handle )
{
  int err;
  iGeom_deleteEnt( mInstance, handle, &err );
  return (Error)err;
}

//inline iGeom::Error
//iGeom::deleteEntArr( const EntityHandle* entity_handles, int num_handles )
//{
//  int err;
//  iGeom_deleteEntArr( mInstance, entity_handles, num_handles, &err );
//  return (Error)err;
//}

//inline iGeom::Error
//iGeom::getAdjEntities( EntitySetHandle set,
//                       EntityType type_requestor,
//                       EntityType type_requested,
//                       std::vector<EntityHandle>& adj_entity_handles,
//                       std::vector<int>& offset )
//{
//  std::vector<EntityHandle> entities;
//  Error err = getEntities( set, type_requestor, entities );
//  if (iBase_SUCCESS != err)
//    return err;
//  
//  offset.resize( entities.size() + 1 );
//  return getArrAdj( &entities[0], entities.size(), type_requested,
//                    adj_entity_handles, &offset[0] );
//}
    
inline iGeom::Error 
iGeom::initEntIter( EntitySetHandle set,
                    EntityType requested_type,
                    iGeom::EntIter& iter )
{
  int err;
  iter.mInstance = mInstance;
  iGeom_initEntIter( mInstance, set, requested_type,
                     &iter.mHandle, &err );
  return (Error)err;
}

inline iGeom::Error 
iGeom::initEntArrIter( EntitySetHandle set,
                       EntityType requested_type,
                       int requested_array_size,
                       iGeom::EntArrIter& iter )
{
  int err;
  iter.mInstance = mInstance;
  iter.mSize = requested_array_size;
  iGeom_initEntArrIter( mInstance, set, requested_type,
                        requested_array_size, &iter.mHandle, &err );
  return (Error)err;
}

inline 
iGeom::EntArrIter::~EntArrIter()
{
  int err;
  if (mHandle != 0) {
    iGeom_endEntArrIter( mInstance, mHandle, &err );
    mHandle = 0;
  }
}

inline 
iGeom::EntIter::~EntIter()
{
  int err;
  if (mHandle != 0) {
    iGeom_endEntIter( mInstance, mHandle, &err );
    mHandle = 0;
  }
}

inline iGeom::Error
iGeom::EntArrIter::getNext( EntityHandle* entity_handles,
                            int& size_out,
                            bool& has_more_data_out )
{
  int err, alloc = mSize, has_data;
  iGeom_getNextEntArrIter( mInstance, mHandle, &entity_handles, &alloc,
                           &size_out, &has_data, &err );
  has_more_data_out = (has_data != 0);
  return (Error)err;
}

inline iGeom::Error
iGeom::EntIter::getNext( EntityHandle& handle_out, bool& has_more_data_out )
{
  int err, has_data;
  iGeom_getNextEntIter( mInstance, mHandle, &handle_out, &has_data, &err );
  has_more_data_out = (has_data != 0);
  return (Error)err;
}

inline iGeom::Error
iGeom::EntArrIter::reset()
{
  int err;
  iGeom_resetEntArrIter( mInstance, mHandle, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::EntIter::reset()
{
  int err;
  iGeom_resetEntIter( mInstance, mHandle, &err );
  return (Error)err;
}
 
inline iGeom::Error
iGeom::getEntType( EntityHandle handle, EntityType& type_out )
{
  int err, result;
  iGeom_getEntType( mInstance, handle, &result, &err );
  type_out = (EntityType)result;
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrType( const EntityHandle* entity_handles,
                   int entity_handles_size,
                   EntityType* types_out )
{
  int err, alloc = entity_handles_size, junk, *ptr;
  std::vector<int> storage;
  if (sizeof(EntityType) == sizeof(int))
    ptr = reinterpret_cast<int*>(types_out);
  else {
    storage.resize( entity_handles_size );
    ptr = &storage[0];
  }
  
  iGeom_getArrType( mInstance, entity_handles, entity_handles_size,
                    &ptr, &alloc, &junk, &err );
  
  if (sizeof(EntityType) != sizeof(int))
    for (int i = 0; i < entity_handles_size; ++i)
      types_out[i] = (EntityType)storage[i];
      
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntAdj( EntityHandle handle, 
                  EntityType type_requested,
                  std::vector<EntityHandle>& adj_entities_out )                  
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  EntityHandle* ptr = &adj_entities_out[0];
  iGeom_getEntAdj( mInstance, handle, type_requested, 
                   &ptr, &alloc, &size, &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iGeom_getEntAdj( mInstance, handle, type_requested, 
                     &ptr, &alloc, &size, &err );
  }
  
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrAdj( const EntityHandle* entity_handles,
                  int entity_handles_size,
                  EntityType type_requested,
                  std::vector<EntityHandle>& adj_entities_out,
                  int* offsets_out )
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12*entity_handles_size);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  int off_alloc = entity_handles_size+1, junk;
  EntityHandle* ptr = &adj_entities_out[0];
  iGeom_getArrAdj( mInstance, entity_handles, entity_handles_size, 
                   type_requested,
                   &ptr, &alloc, &size,
                   &offsets_out, &off_alloc, &junk,
                   &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iGeom_getArrAdj( mInstance, entity_handles, entity_handles_size, 
                     type_requested,
                     &ptr, &alloc, &size,
                     &offsets_out, &off_alloc, &junk,
                     &err );
  }
  
  return (Error)err;
}

                   
inline iGeom::Error
iGeom::getEnt2ndAdj( EntityHandle handle, 
                     EntityType bridge_entity_type,
                     EntityType type_requested,
                     std::vector<EntityHandle>& adj_entities_out )                  
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  EntityHandle* ptr = &adj_entities_out[0];
  iGeom_getEnt2ndAdj( mInstance, handle, bridge_entity_type, type_requested, 
                      &ptr, &alloc, &size, &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iGeom_getEnt2ndAdj( mInstance, handle, bridge_entity_type, type_requested, 
                       &ptr, &alloc, &size, &err );
  }
  
  return (Error)err;
}

inline iGeom::Error
iGeom::getArr2ndAdj( const EntityHandle* entity_handles,
                     int entity_handles_size,
                     EntityType order_key,
                     EntityType type_requested,
                     std::vector<EntityHandle>& adj_entities_out,
                     int* offsets_out )
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12*entity_handles_size);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  int off_alloc = entity_handles_size+1, junk;
  EntityHandle* ptr = &adj_entities_out[0];
  iGeom_getArr2ndAdj( mInstance, entity_handles, entity_handles_size, 
                      order_key, type_requested,
                      &ptr, &alloc, &size,
                      &offsets_out, &off_alloc, &junk,
                      &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iGeom_getArr2ndAdj( mInstance, entity_handles, entity_handles_size, 
                        order_key, type_requested,
                        &ptr, &alloc, &size,
                        &offsets_out, &off_alloc, &junk,
                        &err );
  }
  
  return (Error)err;
}


inline iGeom::Error
iGeom::getBoundBox( double& min_x, double& min_y, double& min_z,
                    double& max_x, double& max_y, double& max_z )
{
  int err;
  iGeom_getBoundBox( mInstance, &min_x, &min_y, &min_z, &max_x, &max_y, &max_z, &err );
  return (Error)err;
}

    
inline int
iGeom::getParametric()
{
  int err, result;
  iGeom_getParametric( mInstance, &result, &err );
  return result;
}

    
inline iGeom::Error
iGeom::isEntAdj( EntityHandle entity1, EntityHandle entity2, bool& adjacent_out )
{
  int err, result;
  iGeom_isEntAdj( mInstance, entity1, entity2, &result, &err );
  adjacent_out = (result != 0);
  return (Error)err;
}

    
inline iGeom::Error
iGeom::isArrAdj( const EntityHandle* entities1, 
                 const EntityHandle* entities2,
                 int num_entity_pairs,
                 int* is_adj_out )
{
  int err, alloc = num_entity_pairs, size = 0;
  iGeom_isArrAdj( mInstance, entities1, num_entity_pairs, entities2, num_entity_pairs, 
            &is_adj_out, &alloc, &size, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::getEntClosestPt( EntityHandle entity,
                        double near_x, double near_y, double near_z,
                        double& on_x, double& on_y, double& on_z )
{
  int err;
  iGeom_getEntClosestPt( mInstance, entity, near_x, near_y, near_z,
                         &on_x, &on_y, &on_z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrClosestPt( const EntityHandle* handles, 
                        int handles_size,
                        StorageOrder order,
                        const double* near_coordinates,
                        int near_coordinates_size,
                        double* on_coordinates )
{
  int err, alloc = std::max(near_coordinates_size,3*handles_size), size = 0;
  iGeom_getArrClosestPt( mInstance, 
                         handles, handles_size, order,
                         near_coordinates, near_coordinates_size,
                         &on_coordinates, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getEntNrmlXYZ( EntityHandle entity,
                      double x, double y, double z,
                      double& i, double& j, double& k )
{
  int err;
  iGeom_getEntNrmlXYZ( mInstance, entity, x, y, z, &i, &j, &k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrNrmlXYZ( const EntityHandle* entities,
                      int entities_size,
                      StorageOrder order,
                      const double* xyz,
                      int xyz_size,
                      double* ijk )
{
  int err, alloc = std::max(xyz_size,3*entities_size), size = 0;
  iGeom_getArrNrmlXYZ( mInstance, 
                       entities, entities_size,
                       order,
                       xyz, xyz_size,
                       &ijk, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getEntNrmlPlXYZ( EntityHandle entity,
                        double x, double y, double z,
                        double& on_x, double& on_y, double& on_z,
                        double& i, double& j, double& k )
{
  int err;
  iGeom_getEntNrmlPlXYZ( mInstance, entity, x, y, z, &on_x, &on_y, &on_z, &i, &j, &k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrNrmlPlXYZ( const EntityHandle* entities,
                        int entities_size,
                        StorageOrder order,
                        const double* near_xyz,
                        int near_xyz_size,
                        double* on_xyz,
                        double* nrml_ijk )
{
  int err, alloc = std::max(near_xyz_size,3*entities_size), size = 0;
  iGeom_getArrNrmlPlXYZ( mInstance, 
                   entities, entities_size,
                   order,
                   near_xyz, near_xyz_size,
                   &on_xyz, &alloc, &size,
                   &nrml_ijk, &alloc, &size,
                   &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getEntTgntXYZ( EntityHandle entity,
                      double x, double y, double z,
                      double& i, double& j, double& k )
{
  int err;
  iGeom_getEntTgntXYZ( mInstance, entity, x, y, z, &i, &j, &k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrTgntXYZ( const EntityHandle* entities,
                      int entities_size,
                      StorageOrder order,
                      const double* xyz,
                      int xyz_size,
                      double* ijk )
{
  int err, alloc = std::max(xyz_size,3*entities_size), size = 0;
  iGeom_getArrTgntXYZ( mInstance, 
                 entities, entities_size,
                 order,
                 xyz, xyz_size,
                 &ijk, &alloc, &size, 
                 &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getFcCvtrXYZ( EntityHandle face,
                     double x, double y, double z,
                     double& i1, double& j1, double& k1,
                     double& i2, double& j2, double& k2 )
{
  int err;
  iGeom_getFcCvtrXYZ( mInstance, face, x, y, z, &i1, &j1, &k1, &i2, &j2, &k2, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEgCvtrXYZ( EntityHandle edge,
                     double x, double y, double z,
                     double& i, double& j, double& k )
{
  int err;
  iGeom_getEgCvtrXYZ( mInstance, edge, x, y, z, &i, &j, &k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntArrCvtrXYZ( const EntityHandle* entities,
                         int entities_size,
                         StorageOrder order,
                         const double* xyz,
                         int xyz_size,
                         double* cvtr_1,
                         double* cvtr_2 )
{
  int err, alloc = std::max(xyz_size,3*entities_size), size = 0;
  iGeom_getEntArrCvtrXYZ( mInstance, 
                          entities, entities_size,
                          order, xyz, xyz_size,
                          &cvtr_1, &alloc, &size,
                          &cvtr_2, &alloc, &size, 
                          &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getEgEvalXYZ( EntityHandle edge,
                     double x, double y, double z,
                     double& on_x, double& on_y, double& on_z,
                     double& tngt_i, double& tngt_j, double& tngt_k,
                     double& cvtr_i, double& cvtr_j, double& cvtr_k )
{
  int err;
  iGeom_getEgEvalXYZ( mInstance, edge, x, y, z, &on_x, &on_y, &on_z,
                      &tngt_i, &tngt_j, &tngt_k, 
                      &cvtr_i, &cvtr_j, &cvtr_k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getFcEvalXYZ( EntityHandle face,
                     double x, double y, double z,
                     double& on_x, double& on_y, double& on_z,
                     double& tngt_i, double& tngt_j, double& tngt_k,
                     double& cvtr1_i, double& cvtr1_j, double& cvtr1_k,
                     double& cvtr2_i, double& cvtr2_j, double& cvtr2_k )
{
  int err;
  iGeom_getFcEvalXYZ( mInstance, face, x, y, z, &on_x, &on_y, &on_z,
                      &tngt_i, &tngt_j, &tngt_k, 
                      &cvtr1_i, &cvtr1_j, &cvtr1_k, 
                      &cvtr2_i, &cvtr2_j, &cvtr2_k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrEgEvalXYZ( const EntityHandle* edges,
                        int edges_size,
                        StorageOrder order,
                        const double* near_coords,
                        int near_coords_size,
                        double* on_coords,
                        double* tangent,
                        double* curvature )
{
  int err, alloc = std::max(near_coords_size,3*edges_size), size = 0;
  iGeom_getArrEgEvalXYZ( mInstance, edges, edges_size, order,
                         near_coords, near_coords_size,
                         &on_coords, &alloc, &size,
                         &tangent, &alloc, &size,
                         &curvature, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrFcEvalXYZ( const EntityHandle* faces,
                        int faces_size,
                        StorageOrder order,
                        const double* near_coords,
                        int near_coords_size,
                        double* on_coords,
                        double* tangent,
                        double* curvature1,
                        double* curvature2 )
{
  int err, alloc = std::max(near_coords_size,3*faces_size), size = 0;
  iGeom_getArrFcEvalXYZ( mInstance, faces, faces_size, order,
                         near_coords, near_coords_size,
                         &on_coords, &alloc, &size,
                         &tangent, &alloc, &size,
                         &curvature1, &alloc, &size,
                         &curvature2, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntBoundBox( EntityHandle entity,
                       double& min_x, double& min_y, double& min_z,
                       double& max_x, double& max_y, double& max_z )
{
  int err;
  iGeom_getEntBoundBox( mInstance, entity, &min_x, &min_y, &min_z, &max_x, &max_y, &max_z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrBoundBox( const EntityHandle* entities,
                       int entities_size,
                       StorageOrder order,
                       double* min_corners,
                       double* max_corners )
{
  int err, alloc = 3*entities_size, size = 0, order_int = order;
  iGeom_getArrBoundBox( mInstance, entities, entities_size,
                        order_int,
                        &min_corners, &alloc, &size,
                        &max_corners, &alloc, &size, 
                        &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getVtxCoord( EntityHandle vertex, double& x, double& y, double& z )
{
  int err;
  iGeom_getVtxCoord( mInstance, vertex, &x, &y, &z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getVtxArrCoords( const EntityHandle* vertices,
                        int vertices_size,
                        StorageOrder order,
                        double* coords )
{
  int err, alloc = vertices_size, size = 0;
  iGeom_getVtxArrCoords( mInstance, vertices, vertices_size, order,
                         &coords, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getPntRayIntsct( double x, double y, double z,
                        double i, double j, double k,
                        StorageOrder order,
                        std::vector<EntityHandle>& entities_out,
                        std::vector<double>& points_out,
                        std::vector<double>& params_out )
{
  int err, count;
  Error err2 = getNumOfType( getRootSet(), iBase_ALL_TYPES, count );
  if (err2 != iBase_SUCCESS)
    return err2;
  
  entities_out.resize( count );
  points_out.resize( 3*count );
  params_out.resize( 2*count );
  int entities_alloc = entities_out.size(), entities_size = 0;
  int   points_alloc =   points_out.size(),   points_size = 0;
  int   params_alloc =   params_out.size(),   params_size = 0;
  EntityHandle* entities_ptr = &entities_out[0];
  double      *   points_ptr = &  points_out[0];
  double      *   params_ptr = &  params_out[0];
  
  iGeom_getPntRayIntsct( mInstance, x, y, z, i, j, k, 
                   &entities_ptr, &entities_alloc, &entities_size,
                   order,
                   &points_ptr, &points_alloc, &points_size,
                   &params_ptr, &params_alloc, &params_size,
                   &err );
  entities_out.resize( entities_size );
  points_out.resize( points_size );
  params_out.resize( params_size );
  if (err == iBase_BAD_ARRAY_SIZE || err == iBase_BAD_ARRAY_DIMENSION) {
    entities_alloc = entities_out.size();
      points_alloc =   points_out.size();
      params_alloc =   params_out.size();
    entities_ptr = &entities_out[0];
      points_ptr = &  points_out[0];
      params_ptr = &  params_out[0];
    iGeom_getPntRayIntsct( mInstance, x, y, z, i, j, k, 
                     &entities_ptr, &entities_alloc, &entities_size,
                     order,
                     &points_ptr, &points_alloc, &points_size,
                     &params_ptr, &params_alloc, &params_size,
                     &err );
  }
  
  return (Error)err;
}

    
inline iGeom::Error 
iGeom::getPntClsf( double x, double y, double z, EntityHandle& handle_out )
{
  int err;
  iGeom_getPntClsf( mInstance, x, y, z, &handle_out, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getPntArrClsf( StorageOrder order,
                      const double* coords,
                      int coords_size,
                      EntityHandle* entities_out )
{
  int err, alloc = coords_size/3, size = 0;
  iGeom_getPntArrClsf( mInstance, order, coords, coords_size, 
                 &entities_out, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error iGeom::getSense(EntityHandle ent, EntityHandle wrt_ent, int &sense) 
{
  EntityType tp_wrt, tp_ent;
  Error err = getEntType(wrt_ent, tp_wrt);
  if (iBase_SUCCESS != err) return err;
  err = getEntType(ent, tp_ent);
  if (iBase_SUCCESS != err) return err;
  if (tp_wrt-tp_ent != 1) return iBase_FAILURE;
  switch (tp_wrt) {
    case iBase_REGION:
        return getEntNrmlSense(ent, wrt_ent, sense);
        break;
    case iBase_FACE:
        return getEgFcSense(ent, wrt_ent, sense);
    case iBase_EDGE:
        return getEgVtxSense(wrt_ent, ent, ent, sense);
    case iBase_VERTEX:
        return iBase_FAILURE;
  }
  return iBase_FAILURE;
}
  
inline iGeom::Error iGeom::getArrSense(const EntityHandle *ent, int num_ents, EntityHandle wrt_ent, int *sense) 
{
  EntityType tp_wrt, tp_ent;
  Error err = getEntType(wrt_ent, tp_wrt);
  if (iBase_SUCCESS != err) return err;
  err = getEntType(ent[0], tp_ent);
  if (iBase_SUCCESS != err) return err;
  if (tp_wrt-tp_ent != 1) return iBase_FAILURE;
  std::vector<EntityHandle> dum_wrts(num_ents, wrt_ent);
  switch (tp_wrt) {
    case iBase_REGION:
        return getArrNrmlSense(ent, num_ents, &dum_wrts[0], num_ents, sense);
        break;
    case iBase_FACE:
        return getEgFcArrSense(ent, num_ents, &dum_wrts[0], num_ents, sense);
        break;
    case iBase_EDGE:
        return getEgVtxArrSense(&dum_wrts[0], num_ents, ent, num_ents, ent, num_ents, sense);
        break;
    case iBase_VERTEX:
        return iBase_FAILURE;
        break;
  }
  return iBase_FAILURE;
}
    
inline iGeom::Error
iGeom::getEntNrmlSense( EntityHandle face, EntityHandle region, int& sense )
{
  int err;
  iGeom_getEntNrmlSense( mInstance, face, region, &sense, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEgFcSense( EntityHandle edge, EntityHandle face, int& sense )
{
  int err;
  iGeom_getEgFcSense( mInstance, edge, face, &sense, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEgVtxSense( EntityHandle edge, EntityHandle vtx1, EntityHandle vtx2, int& sense )
{
  int err;
  iGeom_getEgVtxSense( mInstance, edge, vtx1, vtx2, &sense, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getArrNrmlSense( const EntityHandle* faces, int faces_size,
                        const EntityHandle* vols,  int vols_size,
                        int* senses_out )
{
  int err, alloc = std::max(vols_size, faces_size), size = 0;
  iGeom_getArrNrmlSense( mInstance, 
                         faces, faces_size, 
                         vols, vols_size,
                         &senses_out, &alloc, &size, 
                         &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEgFcArrSense( const EntityHandle* edges, int edges_size,
                        const EntityHandle* faces, int faces_size,
                        int* senses_out )
{
  int err, alloc = std::max(edges_size, faces_size), size = 0;
  iGeom_getEgFcArrSense( mInstance, 
                         edges, edges_size,
                         faces, faces_size,
                         &senses_out, &alloc, &size,
                         &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEgVtxArrSense( const EntityHandle* edges, int edges_size,
                         const EntityHandle* vertices1, int vertices1_size,
                         const EntityHandle* vertices2, int vertices2_size,
                         int* senses_out )
{
  int err, alloc = std::max(vertices1_size,std::max(vertices2_size,edges_size)), size = 0;
  iGeom_getEgVtxArrSense( mInstance, 
                          edges, edges_size,
                          vertices1, vertices1_size,
                          vertices2, vertices2_size,
                          &senses_out, &alloc, &size, 
                          &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::measure( const EntityHandle* entities,
                int entities_size,
                double* measures )
{
  int err, alloc = entities_size, size = 0;
  iGeom_measure( mInstance, entities, entities_size, &measures, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getFaceType( EntityHandle face, std::string& type )
{
  char buffer[1024];
  int err, len = sizeof(buffer);
  iGeom_getFaceType( mInstance, face, buffer, &err, &len );
  type = std::string( buffer, len );
  return (Error)err;
}

inline iGeom::Error
iGeom::isEntParametric( EntityHandle entity, bool& parametric )
{
  int err, result;
  iGeom_isEntParametric( mInstance, entity, &result, &err );
  parametric = (result != 0);
  return (Error)err;
}

inline iGeom::Error
iGeom::isArrParametric( const EntityHandle* entities,
                        int entities_size,
                        int* is_parametric )
{
  int err, alloc = entities_size, size = 1;
  iGeom_isArrParametric( mInstance, entities, entities_size, 
                   &is_parametric, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getEntUVtoXYZ( EntityHandle face,
                      double u, double v,
                      double& x, double& y, double& z )
{
  int err;
  iGeom_getEntUVtoXYZ( mInstance, face, u, v, &x, &y, &z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntUtoXYZ( EntityHandle edge, double u,
                     double& x, double& y, double& z )
{
  int err;
  iGeom_getEntUtoXYZ( mInstance, edge, u, &x, &y, &z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrUVtoXYZ( const EntityHandle* faces,
                      int faces_size,
                      StorageOrder order,
                      const double* uv,
                      int uv_size,
                      double* xyz )
{
  int err, alloc = std::max(3*uv_size/2,3*faces_size), size = 0;
  iGeom_getArrUVtoXYZ( mInstance, faces, faces_size, order, uv, uv_size, 
                 &xyz, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrUtoXYZ( const EntityHandle* edges,
                     int edges_size,
                     const double* u,
                     int u_size,
                     StorageOrder order,
                     double* xyz )
{
  int err, alloc = std::max(3*u_size,3*edges_size), size = 0;
  iGeom_getArrUtoXYZ( mInstance, edges, edges_size, u, u_size, order, 
                      &xyz, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getEntXYZtoUV( EntityHandle face,
                      double x, double y, double z,
                      double& u, double& v )
{
  int err;
  iGeom_getEntXYZtoUV( mInstance, face, x, y, z, &u, &v, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntXYZtoU( EntityHandle edge,
                     double x, double y, double z,
                     double& u )
{
  int err;
  iGeom_getEntXYZtoU( mInstance, edge, x, y, z, &u, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrXYZtoUV( const EntityHandle* faces,
                      int faces_size,
                      StorageOrder order,
                      const double* coords,
                      int coords_size,
                      double* uv )
{
  int err, alloc = std::max(2*coords_size/3,2*faces_size), size = 0;
  iGeom_getArrXYZtoUV( mInstance, faces, faces_size, order, coords, coords_size,
                 &uv, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrXYZtoU( const EntityHandle* edges,
                     int edges_size,
                     StorageOrder order,
                     const double* coords,
                     int coords_size,
                     double* u )
{
  int err, alloc = std::max(coords_size/3,edges_size), size = 0;
  iGeom_getArrXYZtoU( mInstance, edges, edges_size, order, coords, coords_size,
                &u, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntXYZtoUVHint( EntityHandle face,
                          double x, double y, double z,
                          double& u, double& v )
{
  int err;
  iGeom_getEntXYZtoUVHint( mInstance, face, x, y, z, &u, &v, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrXYZtoUVHint( const EntityHandle* faces,
                          int faces_size,
                          StorageOrder order,
                          const double* coords,
                          int coords_size,
                          double* uv )
{
  int err, alloc = std::max(2*coords_size/3,2*faces_size), size = 0;
  iGeom_getArrXYZtoUVHint( mInstance, faces, faces_size, order, coords, coords_size,
                           &uv, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getEntNrmlUV( EntityHandle face,
                     double u, double v,
                     double& i, double& j, double& k )
{
  int err;
  iGeom_getEntNrmlUV( mInstance, face, u, v, &i, &j, &k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrNrmlUV( const EntityHandle* faces,
                     int faces_size,
                     StorageOrder order,
                     const double* uv,
                     int uv_size,
                     double* normals )
{
  int err, alloc = std::max(3*uv_size/2,3*faces_size), size = 0;
  iGeom_getArrNrmlUV( mInstance, faces, faces_size,order, uv, uv_size,
                &normals, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntTgntU( EntityHandle edge,
                    double u, 
                    double& i, double& j, double& k )
{
  int err;
  iGeom_getEntTgntU( mInstance, edge, u, &i, &j, &k, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrTgntU( const EntityHandle* edges,
                    int edges_size,
                    StorageOrder order,
                    const double* u,
                    int u_size,
                    double* normals )
{
  int err, alloc = std::max(3*u_size,3*edges_size), size = 0;
  iGeom_getArrTgntU( mInstance, edges, edges_size, order, u, u_size,
               &normals, &alloc, &size, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::getEnt1stDrvt( EntityHandle handle,
                      double u, double v,
                      double& du_i, double& du_j, double& du_k,
                      double& dv_i, double& dv_j, double& dv_k )
{
  int err, du_alloc = 3, dv_alloc = 3, du_size = 0, dv_size = 0;
  double du[3], dv[3], *du_ptr = du, *dv_ptr = dv;
  iGeom_getEnt1stDrvt( mInstance, handle, u, v, 
                       &du_ptr, &du_alloc, &du_size,
                       &dv_ptr, &dv_alloc, &dv_size, &err );
  du_i = du[0];
  du_j = du[1];
  du_k = du[2];
  dv_i = dv[0];
  dv_j = dv[1];
  dv_k = dv[2];
  return (Error)err;
}

inline iGeom::Error
iGeom::getEnt2ndDrvt( EntityHandle handle,
                      double u, double v,
                      double& duu_i, double& duu_j, double& duu_k,
                      double& dvv_i, double& dvv_j, double& dvv_k,
                      double& duv_i, double& duv_j, double& duv_k )
{
  int err, uu_alloc = 3, uv_alloc = 3, vv_alloc = 3, uu_size = 0, uv_size = 0, vv_size = 0;
  double uu[3], uv[3], vv[3], *uu_ptr = uu, *vv_ptr = vv, *uv_ptr = uv;
  iGeom_getEnt2ndDrvt( mInstance, handle, u, v, 
                       &uu_ptr, &uu_alloc, &uu_size,
                       &vv_ptr, &vv_alloc, &vv_size,
                       &uv_ptr, &uv_alloc, &uv_size, 
                       &err ); 
  duu_i = uu[0];
  duu_j = uu[1];
  duu_k = uu[2];
  dvv_i = vv[0];
  dvv_j = vv[1];
  dvv_k = vv[2];
  duv_i = uv[0];
  duv_j = uv[1];
  duv_k = uv[2];
  return (Error)err;
}

inline iGeom::Error
iGeom::getArr1stDrvt( const EntityHandle* entities,
                      int entities_size,
                      StorageOrder order,
                      const double* uv,
                      int uv_size,
                      double* dvtr_u,
                      double* dvtr_v )
{
  int err, allocu = std::max(3*uv_size/2,3*entities_size), sizeu = 0,
           allocv = allocu, sizev = 0;
  std::vector<int> offset1(std::max(uv_size/2,entities_size)+1);
  std::vector<int> offset2(std::max(uv_size/2,entities_size)+1);
  int alloc1 = offset1.size(), size1 = 0, *ptr1 = &offset1[0];
  int alloc2 = offset2.size(), size2 = 0, *ptr2 = &offset2[0];
  iGeom_getArr1stDrvt( mInstance, entities, entities_size, 
                       order, uv, uv_size, 
                       &dvtr_u, &allocu, &sizeu,
                       &ptr1, &alloc1, &size1,
                       &dvtr_v, &allocv, &sizev,
                       &ptr2, &alloc2, &size2,
                       &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArr2ndDrvt( const EntityHandle* entities,
                      int entities_size,
                      StorageOrder order,
                      const double* uv,
                      int uv_size,
                      double* dvtr_uu,
                      double* dvtr_vv,
                      double* dvtr_uv )
{
  int err, allocuu = std::max(3*uv_size/2,3*entities_size), sizeuu = 0,
           allocvv = allocuu, sizevv = 0, allocuv = allocuu, sizeuv = 0;
  std::vector<int> offset1(std::max(uv_size/2,entities_size)+1);
  std::vector<int> offset2(std::max(uv_size/2,entities_size)+1);
  std::vector<int> offset3(std::max(uv_size/2,entities_size)+1);
  int alloc1 = offset1.size(), size1 = 0, *ptr1 = &offset1[0];
  int alloc2 = offset2.size(), size2 = 0, *ptr2 = &offset2[0];
  int alloc3 = offset3.size(), size3 = 0, *ptr3 = &offset3[0];
  iGeom_getArr2ndDrvt( mInstance, entities, entities_size, 
                       order, uv, uv_size,
                       &dvtr_uu, &allocuu, &sizeuu, 
                       &ptr1, &alloc1, &size1,
                       &dvtr_vv, &allocvv, &sizevv, 
                       &ptr2, &alloc2, &size2,
                       &dvtr_uv, &allocuv, &sizeuv, 
                       &ptr3, &alloc3, &size3,
                       &err );
  return (Error)err;
}

    

inline iGeom::Error
iGeom::getFcCvtrUV( EntityHandle face, double u, double v,
                    double& i1, double& j1, double& k1,
                    double& i2, double& j2, double& k2 )
{
  int err;
  iGeom_getFcCvtrUV( mInstance, face, u, v, &i1, &j1, &k1, &i2, &j2, &k2, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getFcArrCvtrUV( const EntityHandle* faces, int faces_size,
                       StorageOrder order,
                       const double* uv, int uv_size,
                       double* cvtr1, double* cvtr2 )
{
  int err, alloc = std::max(3*uv_size/2, 3*faces_size), size = 0;
  iGeom_getFcArrCvtrUV( mInstance, faces, faces_size, order, uv, uv_size, 
                  &cvtr1, &alloc, &size,
                  &cvtr2, &alloc, &size, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::isEntPeriodic( EntityHandle entity, bool& in_u, bool& in_v )
{
  int err, u, v;
  iGeom_isEntPeriodic( mInstance, entity, &u, &v, &err );
  in_u = (u != 0);
  in_v = (v != 0);
  return (Error)err;
}

inline iGeom::Error
iGeom::isArrPeriodic( const EntityHandle* entities,
                      int entities_size,
                      int* in_uv )
{
  int err, alloc = 2*entities_size, size = 0;
  iGeom_isArrPeriodic( mInstance, entities, entities_size, 
                       &in_uv, &alloc, &size, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::isFcDegenerate( EntityHandle face, bool& is_degenerate )
{
  int err, result;
  iGeom_isFcDegenerate( mInstance, face, &result, &err );
  is_degenerate = (result != 0);
  return (Error)err;
}

inline iGeom::Error
iGeom::isFcArrDegenerate( const EntityHandle* faces,
                          int faces_size,
                          int* degenerate )
{
  int err, alloc = faces_size, size = 0;
  iGeom_isFcArrDegenerate( mInstance, faces, faces_size, &degenerate, &alloc, &size, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::getTolerance( int& type_out, double& tolerance_out )
{
  int err;
  iGeom_getTolerance( mInstance, &type_out, &tolerance_out, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntTolerance( EntityHandle entity, double& tolerance )
{
  int err;
  iGeom_getEntTolerance( mInstance, entity, &tolerance, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrTolerance( const EntityHandle* entities,
                        int entities_size,
                        double* tolerances )
{
  int err, alloc = entities_size, size = 0;
  iGeom_getArrTolerance( mInstance, entities, entities_size, 
                         &tolerances, &alloc, &size, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::getEntUVRange( EntityHandle face,
                      double& u_min, double& v_min,
                      double& u_max, double& v_max )
{
  int err;
  iGeom_getEntUVRange( mInstance, face, &u_min, &v_min, &u_max, &v_max, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getEntURange( EntityHandle edge, double& u_min, double& u_max )
{
  int err;
  iGeom_getEntURange( mInstance, edge, &u_min, &u_max, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrUVRange( const EntityHandle* faces,
                      int faces_size,
                      StorageOrder order,
                      double* uv_min,
                      double* uv_max )
{
  int err, alloc = faces_size, size = 0;
  iGeom_getArrUVRange( mInstance, faces, faces_size, order,
                 &uv_min, &alloc, &size,
                 &uv_max, &alloc, &size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrURange( const EntityHandle* edges,
                     int edges_size,
                     double* u_min,
                     double* u_max )
{
  int err, alloc = edges_size, size = 0;
  iGeom_getArrURange( mInstance, edges, edges_size, 
                &u_min, &alloc, &size,
                &u_max, &alloc, &size, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::getEntUtoUV( EntityHandle edge, 
                    EntityHandle face,
                    double edge_u,
                    double& face_u, double& face_v )
{
  int err;
  iGeom_getEntUtoUV( mInstance, edge, face, edge_u, &face_u, &face_v, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getVtxToUV( EntityHandle vertex,
                   EntityHandle face,
                   double& u, double& v )
{
  int err;
  iGeom_getVtxToUV( mInstance, vertex, face, &u, &v, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getVtxToU( EntityHandle vertex, EntityHandle edge, double& u )
{
  int err;
  iGeom_getVtxToU( mInstance, vertex, edge, &u, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getArrUtoUV( const EntityHandle* edges, int edges_size,
                    const EntityHandle* faces, int faces_size,
                    const double* edge_u, int edge_u_size,
                    StorageOrder order,
                    double* face_uv )
{
  int err, alloc = std::max(edge_u_size, std::max(edges_size, faces_size));
  int size = 0;
  iGeom_getArrUtoUV( mInstance, 
               edges, edges_size,
               faces, faces_size,
               edge_u, edge_u_size,
               order, 
               &face_uv, &alloc, &size,
               &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getVtxArrToUV( const EntityHandle* vertices, int vertices_size,
                      const EntityHandle* faces, int faces_size,
                      StorageOrder order,
                      double* face_uv )
{
  int err, alloc = std::max(vertices_size, faces_size), size = 0;
  iGeom_getVtxArrToUV( mInstance, 
                       vertices, vertices_size,
                       faces, faces_size, 
                       order,
                       &face_uv, &alloc, &size,
                       &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::getVtxArrToU( const EntityHandle* vertices, int vertices_size,
                     const EntityHandle* edges, int edges_size,
                     double* edge_u )
{
  int err, alloc = std::max(vertices_size, edges_size), size = 0;
  iGeom_getVtxArrToU( mInstance, 
                      vertices, vertices_size,
                      edges, edges_size,
                      &edge_u, &alloc, &size,
                      &err );
  return (Error)err;
}

    
    
inline iGeom::Error
iGeom::deleteAll()
{
  int err;
  iGeom_deleteAll( mInstance, &err );
  return (Error)err;
}


inline iGeom::Error
iGeom::copyEnt( EntityHandle source, EntityHandle& copy )
{
  int err;
  iGeom_copyEnt( mInstance, source, &copy, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::createSphere( double radius, EntityHandle& sphere )
{
  int err;
  iGeom_createSphere( mInstance, radius, &sphere, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::createPrism( double height, int num_sides,
                    double maj_radius, double min_radius,
                    EntityHandle& prism )
{
  int err;
  iGeom_createPrism( mInstance, height, num_sides, maj_radius, min_radius,
                     &prism, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::createBrick( double x, double y, double z, EntityHandle& brick )
{
  int err;
  iGeom_createBrick( mInstance, x, y, z, &brick, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::createCylinder( double height, double maj_rad, double min_rad,
                       EntityHandle& cylinder )
{
  int err;
  iGeom_createCylinder( mInstance, height, maj_rad, min_rad, &cylinder, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::createTorus( double maj_rad, double min_rad, EntityHandle& torus )
{
  int err;
  iGeom_createTorus( mInstance, maj_rad, min_rad, &torus, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::moveEnt( EntityHandle entity, double x, double y, double z )
{
  int err;
  iGeom_moveEnt( mInstance, entity, x, y, z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::rotateEnt( EntityHandle entity, double angle,
                  double axis_x, double axis_y, double axis_z )
{
  int err;
  iGeom_rotateEnt( mInstance, entity, angle, axis_x, axis_y, axis_z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::reflectEnt( EntityHandle entity, 
                   double norm_x, double norm_y, double norm_z )
{
  int err;
  iGeom_reflectEnt( mInstance, entity, norm_x, norm_y, norm_z, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::scaleEnt( EntityHandle entity, 
                 double x_factor, double y_factor, double z_factor )
{
  int err;
  iGeom_scaleEnt( mInstance, entity, x_factor, y_factor, z_factor, &err );
  return (Error)err;
}

  
inline iGeom::Error
iGeom::uniteEnts( const EntityHandle* entities,
                  int entities_size,
                  EntityHandle& result_entity )
{
  int err;
  iGeom_uniteEnts( mInstance, entities, entities_size, &result_entity, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::subtractEnts( EntityHandle blank, EntityHandle tool, 
                     EntityHandle& result )
{
  int err;
  iGeom_subtractEnts( mInstance, blank, tool, &result, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::intersectEnts( EntityHandle entity1, EntityHandle entity2,
                      EntityHandle& result )
{
  int err;
  iGeom_intersectEnts( mInstance, entity1, entity2, &result, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::sectionEnt( EntityHandle entity,
                   double plane_x, double plane_y, double plane_z,
                   double offset,
                   bool reverse,
                   EntityHandle& result )
{
  int err;
  iGeom_sectionEnt( mInstance, entity, plane_x, plane_y, plane_z, offset, reverse, &result, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::sweepEntAboutAxis( EntityHandle entity,
                          double angle,
                          double axis_x, 
                          double axis_y,
                          double axis_z,
                          EntityHandle& swept_entity )
{
  int err;
  iGeom_sweepEntAboutAxis( mInstance, entity, angle, axis_x, axis_y, axis_z,
                     &swept_entity, &err );
  return (Error)err;
}

    
inline iGeom::Error
iGeom::imprintEnts( const EntityHandle* entities, int entities_size )
{
  int err;
  iGeom_imprintEnts( mInstance, entities, entities_size, &err );
  return (Error)err;
}

inline iGeom::Error
iGeom::mergeEnts( const EntityHandle* entities, int entities_size, double tolerance )
{
  int err;
  iGeom_mergeEnts( mInstance, entities, entities_size, tolerance, &err );
  return (Error)err;
}
#endif

} // end namespace MeshKit
