#ifndef IBASEVIRTUALH
#define IBASEVIRTUALH


#include "iGeom.h"

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <string>

#ifndef IBSH
#define IBSH(h) reinterpret_cast<iBase_EntitySetHandle>(h)
#endif    

#ifndef IBSHR
#define IBSHR(h) reinterpret_cast<iBase_EntitySetHandle&>(h)
#endif    

#ifndef IBEH
#define IBEH(h) reinterpret_cast<iBase_EntityHandle>(h)
#endif    

#ifndef IBEHR
#define IBEHR(h) reinterpret_cast<iBase_EntityHandle&>(h)
#endif    

#define PFX__(A,B) A ## B
#define PFX_(A,B) PFX__(A,B)
#define PFX(B) PFX_( ITAPS_PREFIX, B )

class PFX(Base) {
  protected:
    PFX(_Instance) mInstance;

  public:  
    typedef iBase_EntitySetHandle EntitySetHandle;
    typedef iBase_EntityHandle EntityHandle;
    typedef iBase_TagHandle TagHandle;
    typedef iBase_ErrorType Error;
    typedef iBase_EntityType EntityType;
    typedef iBase_StorageOrder StorageOrder;
    typedef iBase_TagValueType TagValueType;

    PFX(_Instance) instance() { return mInstance ;}

    virtual std::string getDescription();

    Error getErrorType();
    
    virtual EntitySetHandle getRootSet();

    virtual Error createEntSet( bool is_list, EntitySetHandle& handle_out );
    virtual Error destroyEntSet( EntitySetHandle handle );
    virtual Error isList( EntitySetHandle handle, bool& is_list );
    
    virtual Error getNumEntSets( EntitySetHandle set, int num_hops, int& num_sets_out );
    virtual Error getEntSets( EntitySetHandle set, int num_hops,
                             std::vector<EntitySetHandle>& contained_sets_out );
    
    virtual Error addEntToSet( EntityHandle entity, EntitySetHandle set );
    virtual Error rmvEntFromSet( EntityHandle entity, EntitySetHandle set );
    
    virtual Error addEntArrToSet( const EntityHandle* entity_handles,
                                 int entity_handles_size,
                                 EntitySetHandle entity_set );
    virtual Error rmvEntArrFromSet( const EntityHandle* entity_handles,
                                   int entity_handles_size,
                                   EntitySetHandle entity_set );
    
    virtual Error addEntSet( EntitySetHandle to_add, EntitySetHandle add_to );
    virtual Error rmvEntSet( EntitySetHandle to_rmv, EntitySetHandle rmv_from );
    
    virtual Error isEntContained( EntitySetHandle set, EntityHandle ent, bool& contained_out );
    virtual Error isEntArrContained( EntitySetHandle containing_set,
                                    const EntityHandle* entity_handles,
                                    int num_entity_handles,
                                    bool* is_contained_out );
    virtual Error isEntSetContained( EntitySetHandle containing_set,
                                    EntitySetHandle contained_set, 
                                    bool& contained_out );
    
    virtual Error addPrntChld( EntitySetHandle parent, EntitySetHandle child );
    virtual Error rmvPrntChld( EntitySetHandle parent, EntitySetHandle child );
    virtual Error isChildOf( EntitySetHandle parent, EntitySetHandle child, bool& is_child_out );
    virtual Error getNumChld( EntitySetHandle parent, int num_hops, int& num_child_out );
    virtual Error getNumPrnt( EntitySetHandle child, int num_hops, int& num_parent_out );
    virtual Error getChldn( EntitySetHandle parent, int num_hops,
                           std::vector<EntitySetHandle>& children_out );
    virtual Error getPrnts( EntitySetHandle child, int num_hops,
                           std::vector<EntitySetHandle>& parents_out );
    
    virtual Error subtract( EntitySetHandle set1, EntitySetHandle set2,
                           EntitySetHandle& result_set_out );
    virtual Error intersect( EntitySetHandle set1, EntitySetHandle set2,
                            EntitySetHandle& result_set_out );
    virtual Error unite( EntitySetHandle set1, EntitySetHandle set2,
                        EntitySetHandle& result_set_out );

    virtual Error createTag( const char* tag_name,
                            int tag_num_type_values,
                            TagValueType tag_type,
                            TagHandle& tag_handle_out );
    
    virtual Error destroyTag( TagHandle tag_handle, bool forced );
    virtual Error getTagName( TagHandle tag_handle, std::string& name_out );
    virtual Error getTagSizeValues( TagHandle tag_handle, int& size_out );
    virtual Error getTagSizeBytes( TagHandle tag_handle, int& size_out );
    virtual Error getTagHandle( const char* name, TagHandle& handle_out );
    virtual Error getTagType( TagHandle tag_handle, TagValueType& type_out );
    
    virtual Error setEntSetData( EntitySetHandle set_handle,
                                TagHandle tag_handle,
                                const void* tag_value );
    virtual Error setEntSetIntData( EntitySetHandle set_handle,
                                   TagHandle tag_handle,
                                   int value );
    virtual Error setEntSetDblData( EntitySetHandle set_handle,
                                   TagHandle tag_handle,
                                   double value );
    virtual Error setEntSetEHData( EntitySetHandle set_handle,
                                  TagHandle tag_handle,
                                  EntityHandle value );
    
    virtual Error getEntSetData( EntitySetHandle set_handle,
                                TagHandle tag_handle,
                                void* tag_value_out );
    virtual Error getEntSetIntData( EntitySetHandle set_handle,
                                   TagHandle tag_handle,
                                   int& value_out );
    virtual Error getEntSetDblData( EntitySetHandle set_handle,
                                   TagHandle tag_handle,
                                   double& value_out );
    virtual Error getEntSetEHData( EntitySetHandle set_handle,
                                  TagHandle tag_handle,
                                  EntityHandle& value_out );
    
    virtual Error getAllEntSetTags( EntitySetHandle set,
                                   std::vector<TagHandle>& tags_out );
    virtual Error getAllTags( EntityHandle entity,
                             std::vector<TagHandle>& tags_out );
                                   
    virtual Error rmvEntSetTag( EntitySetHandle set, TagHandle tag );
    virtual Error rmvTag( EntityHandle entity, TagHandle tag );
    virtual Error rmvArrTag( const EntityHandle* handles, int size, TagHandle tag );
    
    virtual Error getArrData( const EntityHandle* entity_handles,
                             int entity_handles_size,
                             TagHandle tag_handle,
                             void* tag_values_out );
    virtual Error getIntArrData( const EntityHandle* entity_handles,
                                int entity_handles_size,
                                TagHandle tag_handle,
                                int* tag_values_out );
    virtual Error getDblArrData( const EntityHandle* entity_handles,
                                int entity_handles_size,
                                TagHandle tag_handle,
                                double* tag_values_out );
    virtual Error getEHArrData( const EntityHandle* entity_handles,
                               int entity_handles_size,
                               TagHandle tag_handle,
                               EntityHandle* tag_values_out );
    
    virtual Error setArrData( const EntityHandle* entity_handles,
                             int entity_handles_size,
                             TagHandle tag_handle,
                             const void* tag_values );
    virtual Error setIntArrData( const EntityHandle* entity_handles,
                                int entity_handles_size,
                                TagHandle tag_handle,
                                const int* tag_values );
    virtual Error setDblArrData( const EntityHandle* entity_handles,
                                int entity_handles_size,
                                TagHandle tag_handle,
                                const double* tag_values );
    virtual Error setEHArrData( const EntityHandle* entity_handles,
                               int entity_handles_size,
                               TagHandle tag_handle,
                               const EntityHandle* tag_values );
    
    
    virtual Error setData( EntityHandle entity_handle,
                          TagHandle tag_handle,
                          const void* tag_value );
    virtual Error setIntData( EntityHandle entity_handle,
                             TagHandle tag_handle,
                             int value );
    virtual Error setDblData( EntityHandle entity_handle,
                             TagHandle tag_handle,
                             double value );
    virtual Error setEHData( EntityHandle entity_handle,
                            TagHandle tag_handle,
                            EntityHandle value );
    
    virtual Error getData( EntityHandle entity_handle,
                          TagHandle tag_handle,
                          void* tag_value_out );
    virtual Error getIntData( EntityHandle entity_handle,
                             TagHandle tag_handle,
                             int& value_out );
    virtual Error getDblData( EntityHandle entity_handle,
                             TagHandle tag_handle,
                             double& value_out );
    virtual Error getEHData( EntityHandle entity_handle,
                            TagHandle tag_handle,
                            EntityHandle& value_out );
};

#undef PFX__
#undef PFX_
#undef PFX

// IBASEVIRTUALH
#endif
