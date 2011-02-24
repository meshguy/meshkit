

#define ITAPS_PREFIX iGeom
#include "meshkit/iBaseVirtual.hpp"
#define PFX__(A,B) A ## B
#define PFX_(A,B) PFX__(A,B)
#define PFX(B) PFX_( ITAPS_PREFIX, B )

 /*PFX(_Instance) PFX(Base)::instance()
{
  return mInstance;
}*/
    
PFX(Base)::Error
PFX(Base)::getErrorType() 
{
  /*int err1, err2;
  PFX(_getErrorType)( mInstance, &err1, &err2 );
  if (iBase_SUCCESS != err2)
    return (Error)err2;
  else
    return (Error)err1;*/
  return (Error)0;
}

std::string
PFX(Base)::getDescription()
{
  std::vector<char> buffer(1024);
  int err;
  PFX(_getDescription)( mInstance, &buffer[0], &err, buffer.size() );
  if (iBase_SUCCESS != err)
    return std::string();
  else
    return std::string(&buffer[0]);
}




PFX(Base)::EntitySetHandle
PFX(Base)::getRootSet()
{
  int err;
  EntitySetHandle result;
  PFX(_getRootSet)( mInstance, &result, &err );
  return iBase_SUCCESS == err ? result : 0;
}



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

#undef PFX__
#undef PFX_
#undef PFX

#undef ITAPS_PREFIX
