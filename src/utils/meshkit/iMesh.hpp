#ifndef MESHKIT_ITAPS_MESH_HPP
#define MESHKIT_ITAPS_MESH_HPP

/** \file iMesh.hpp
 */

#include "moab/EntityType.hpp"
#include "iMesh.h"
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <string>

#define ITAPS_PREFIX iMesh
#include "meshkit/iBase.hpp"
#undef ITAPS_PREFIX

/** \class iMesh
 * \brief C++ interface to ITAPS iMesh interface
 *
 * This class is a simple wrapper for the ITAPS iGeom interface.  The primary benefit to using this class
 * instead of iMesh directly is that lists of handles are passed as std::vectors instead of pointers to
 * handle arrays.  This file includes both declaration and definition of all iGeom class functions, i.e.
 * all functions are inlined.  The class can be constructed and destructed in the standard C++ way; the
 * implementation of those functions call into the standard iMesh C functions newMesh and dtor.
 *
 * For complete documentation of these functions, see http://www.itaps.org/software/iMesh_html/index.html.
 */
class iMesh : public iMeshBase {
  public:

  /* map to MB's entity type from TSTT's entity topology */
  static const moab::EntityType mb_topology_table[iMesh_ALL_TOPOLOGIES+1];

    typedef iMesh_EntityTopology EntityTopology;
    typedef iBase_AdjacencyCost AdjacencyCost;
    
    explicit inline iMesh( const char* options = 0 );
    inline iMesh( iMesh_Instance imesh );
    
    inline ~iMesh();
    
    inline Error load( EntitySetHandle set,
                       const char* file_name,
                       const char* options = 0 );

    inline Error save( EntitySetHandle set,
                       const char* file_name,
                       const char* options = 0 );
    
    inline int getGeometricDimension() const;
    
    inline Error setGeometricDimension( int dim );
    
    inline StorageOrder getDfltStorage() const;
    
    typedef AdjacencyCost (*AdjTableType)[4];
    inline const AdjTableType getAdjTable();
    
    inline Error getNumOfType( EntitySetHandle set, EntityType type, int& count_out  ) const;
    
    inline Error getNumOfTopo( EntitySetHandle set, EntityTopology topo, int& count_out ) const;
    
    inline bool optimize();
    
    inline Error getEntities( EntitySetHandle set,
                              EntityType type,
                              EntityTopology topo,
                              std::vector<EntityHandle>& entities_out ) const;
    
    inline Error getVtxCoord( EntityHandle vertex, double& x, double& y, double& z ) const;
    inline Error setVtxCoord( EntityHandle vertex, double x, double y, double z );
    
    
    inline Error getVtxArrCoords( const EntityHandle* vertex_handles,
                                  int vertex_handles_size,
                                  StorageOrder storage_order,
                                  double* coords_out ) const;
    
    inline Error setVtxArrCoords( const EntityHandle* vertex_handles,
                                  int vertex_handles_size,
                                  StorageOrder storage_order,
                                  const double* new_coords );
    
    inline Error createVtx( double x, double y, double z, EntityHandle& vertex_out );
    
    inline Error createVtxArr( int num_verts, 
                               StorageOrder storage_order,
                               const double* new_coords,
                               EntityHandle* new_vertex_handles_out );
    
    inline Error createEnt( EntityTopology topology,
                            const EntityHandle* lower_order_entity_handles,
                            int lower_order_entity_handles_size,
                            EntityHandle& new_entity_handle_out );
    
    inline Error createEntArr( EntityTopology new_entity_topology,
                               const EntityHandle* lower_order_entity_handles,
                               int lower_order_entity_handles_size,
                               EntityHandle* new_entity_handles_out );
    
    inline Error deleteEnt( EntityHandle handle );
    
    inline Error deleteEntArr( const EntityHandle* entity_handles, int num_handles );
    
    inline Error getAdjEntities( EntitySetHandle set,
                                 EntityType type_requestor,
                                 EntityTopology topo_requestor,
                                 EntityType type_requested,
                                 std::vector<EntityHandle>& adj_entity_handles,
                                 std::vector<int>& offset ) const;
    
/** \class EntArrIter iMesh.hpp "iMesh.hpp"
 * \brief Class for iterating over %iMesh entity arrays.
 */
    class EntArrIter {
      private:
        friend class iMesh;
        iBase_EntityArrIterator mHandle;
        iMesh_Instance mInstance;
        int mSize;
      public:
        EntArrIter() : mHandle(0), mInstance(0), mSize(0) {}
        inline ~EntArrIter();
        inline Error getNext( EntityHandle* entity_handles_out,
                              int& size_out,
                              bool& has_more_data_out );
        inline Error reset();
    };
    
/** \class EntIter iMesh.hpp "iMesh.hpp"
 * \brief Class for iterating over %iMesh entities.
 */
    class EntIter {
      private:
        friend class iMesh;
        iBase_EntityIterator mHandle;
        iMesh_Instance mInstance;
      public:
        EntIter() : mHandle(0), mInstance(0) {}
        inline ~EntIter();
        inline Error getNext( EntityHandle& entity_handle_out,
                              bool& has_more_data_out );
        inline Error reset();
    };
    
    inline Error initEntIter( EntitySetHandle set,
                              EntityType requested_type,
                              EntityTopology requested_topo,
                              EntIter& iter );
    inline Error initEntArrIter( EntitySetHandle set,
                                 EntityType requested_type,
                                 EntityTopology requested_topo,
                                 int requested_array_size,
                                 EntArrIter& iter );
    
    inline Error getEntTopo( EntityHandle handle, EntityTopology& topo_out ) const;
    inline Error getEntArrTopo( const EntityHandle* entity_handles,
                                int entity_handles_Size,
                                EntityTopology* topos_out ) const;
     
    inline Error getEntType( EntityHandle handle, EntityType& type_out ) const;
    inline Error getEntArrType( const EntityHandle* entity_handles,
                                int entity_handles_Size,
                                EntityType* types_out ) const;
    
    inline Error getEntAdj( EntityHandle handle, 
                            EntityType type_requested,
                            std::vector<EntityHandle>& adj_entities_out ) const;                   
    inline Error getEntArrAdj( const EntityHandle* entity_handles,
                               int entity_handles_size,
                               EntityType type_requested,
                               std::vector<EntityHandle>& adjacent_entity_handles_out,
                               int* offsets_out ) const ;
                       
   inline Error getEnt2ndAdj( EntityHandle handle, 
                              EntityType bridge_entity_type,
                              EntityType type_requested,
                              std::vector<EntityHandle> adj_entities_out ) const;                   
   inline Error getEntArr2ndAdj( const EntityHandle* entity_handles,
                                  int entity_handles_size,
                                  EntityType order_key,
                                  EntityType type_requested,
                                  std::vector<EntityHandle>& adjacent_entity_handles_out,
                                  int* offsets_out ) const;

    inline Error getAdjEntIndices( EntitySetHandle set,
                                   EntityType requestor_type,
                                   EntityTopology requestor_topo,
                                   EntityType type_requested,
                                   std::vector<EntityHandle>& entity_handles_out,
                                   std::vector<EntityHandle>& adj_entity_handles_out,
                                   std::vector<int>& adj_entity_indices_out,
                                   std::vector<int>& offsets_out ) const;
    
  private:
    bool iMeshInstanceOwner;
    AdjacencyCost adjTable[4][4];
    Error adjTableErr;
    
    void cacheAdjTable();
    
      // prohibit copying
    iMesh( const iMesh& ) {}
    void operator=(const iMesh&) {}
};

inline void
iMesh::cacheAdjTable()
{
  int err;
  if (sizeof(int) == sizeof(AdjacencyCost)) {
    int* ptr = (int*)&adjTable[0][0];
    int size = 16, alloc = 16;
    iMesh_getAdjTable( mInstance, &ptr, &alloc, &size, &err );
    adjTableErr = (Error)err;
  }
  else {
    int data[16];
    int* ptr = data;
    int size = 16, alloc = 16;
    iMesh_getAdjTable( mInstance, &ptr, &alloc, &size, &err );
    adjTableErr = (Error)err;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        adjTable[i][j] = (AdjacencyCost)data[4*i+1];
  }
}

inline 
iMesh::iMesh( const char* options )
{
  int err, len = options ? strlen(options) : 0;
  iMesh_newMesh( options, &mInstance, &err, len );
  if (iBase_SUCCESS != err) {
    mInstance = 0;
    iMeshInstanceOwner = false;
  }
  else {
    iMeshInstanceOwner = true;
    cacheAdjTable();
  }
}

inline 
iMesh::iMesh( iMesh_Instance instance )
  : iMeshInstanceOwner(false)
{
  mInstance = instance;
  cacheAdjTable();
}

inline iMesh::~iMesh()
{
  if (iMeshInstanceOwner) {
    int err;
    iMesh_dtor( mInstance, &err );
  }
}

inline iMesh::Error
iMesh::load( EntitySetHandle set,
             const char* file_name,
             const char* options )
{
  int err, len = options ? strlen(options) : 0;
  iMesh_load( mInstance, set, file_name, options, &err, strlen(file_name), len );
  return (Error)err;
}


inline iMesh::Error
iMesh::save( EntitySetHandle set,
             const char* file_name,
             const char* options )
{
  int err, len = options ? strlen(options) : 0;
  iMesh_save( mInstance, set, file_name, options, &err, strlen(file_name), len );
  return (Error)err;
}


inline int 
iMesh::getGeometricDimension() const
{
  int err, result;
  iMesh_getGeometricDimension( mInstance, &result, &err );
  return iBase_SUCCESS == err ? result : -err;
}


inline iMesh::Error
iMesh::setGeometricDimension( int dim )
{
  int err;
  iMesh_setGeometricDimension( mInstance, dim, &err );
  return (Error)err;
}

inline iMesh::StorageOrder 
iMesh::getDfltStorage() const
{
  int err, order;
  iMesh_getDfltStorage( mInstance, &order, &err );
  return (iBase_SUCCESS == err) ? (StorageOrder)order : iBase_BLOCKED;
}

inline const iMesh::AdjTableType
iMesh::getAdjTable() 
{
  return (iBase_SUCCESS == adjTableErr) ? adjTable : 0;
}
  
inline iMesh::Error
iMesh::getNumOfType( EntitySetHandle set, EntityType type, int& count_out  ) const
{
  int err;
  iMesh_getNumOfType( mInstance, set, type, &count_out, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::getNumOfTopo( EntitySetHandle set, EntityTopology topo, int& count_out ) const
{
  int err;
  iMesh_getNumOfTopo( mInstance, set, topo, &count_out, &err );
  return (Error)err;
}

inline bool 
iMesh::optimize()
{
  int err, result;
  iMesh_optimize( mInstance, &result, &err );
  return (iBase_SUCCESS == err) ? !!result : true;
}

inline iMesh::Error
iMesh::getEntities( EntitySetHandle set,
                    EntityType type,
                    EntityTopology topo,
                    std::vector<EntityHandle>& entities_out ) const
{
    // if input vect has no allocated space, allocate some so
    // we don't accidentally ask the impl to allocate an array
  if (entities_out.capacity() == 0) {
    Error err2;
    int count;
    if (topo == iMesh_ALL_TOPOLOGIES)
      err2 = getNumOfType( set, type, count );
    else
      err2 = getNumOfTopo( set, topo, count );
    if (err2 != iBase_SUCCESS)
      return err2;
    entities_out.resize( count );
  }
  
    // try getting results using whatever space input vector has allocated
  int err, size = 0, alloc = entities_out.capacity();
  entities_out.resize( entities_out.capacity() );
  EntityHandle* ptr = &entities_out[0];
  iMesh_getEntities( mInstance, set, type, topo, &ptr, &alloc, &size, &err );
  entities_out.resize(size);
  
    // if input vector was too small, try again with increased size
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = entities_out.size();
    ptr = &entities_out[0];
    iMesh_getEntities( mInstance, set, type, topo, &ptr, &alloc, &size, &err );
  }
  
  return (Error)err;
}

inline iMesh::Error
iMesh::getVtxCoord( EntityHandle vertex, double& x, double& y, double& z ) const
{
  int err;
  iMesh_getVtxCoord( mInstance, vertex, &x, &y, &z, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::setVtxCoord( EntityHandle vertex, double x, double y, double z )
{
  int err;
  iMesh_setVtxCoord( mInstance, vertex, x, y, z, &err );
  return (Error)err;
}


inline iMesh::Error
iMesh::getVtxArrCoords( const EntityHandle* vertex_handles,
                        int vertex_handles_size,
                        StorageOrder storage_order,
                        double* coords_out ) const
{
  int err, alloc = 3*vertex_handles_size, junk, order = storage_order;
  iMesh_getVtxArrCoords( mInstance, vertex_handles, vertex_handles_size,
                         order, &coords_out, &alloc, &junk, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::setVtxArrCoords( const EntityHandle* vertex_handles,
                        int vertex_handles_size,
                        StorageOrder storage_order,
                        const double* new_coords )
{
  int err, dim = getGeometricDimension();
  iMesh_setVtxArrCoords( mInstance, vertex_handles, vertex_handles_size, 
                         storage_order, new_coords, dim*vertex_handles_size,
                         &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::createVtx( double x, double y, double z, EntityHandle& vertex_out )
{
  int err;
  iMesh_createVtx( mInstance, x, y, z, &vertex_out, &err );
  return (Error)err;
}


inline iMesh::Error
iMesh::createVtxArr( int num_verts, 
                     StorageOrder storage_order,
                     const double* new_coords,
                     EntityHandle* new_vertex_handles_out )
{
  int err, alloc = num_verts, junk, dim = getGeometricDimension();
  iMesh_createVtxArr( mInstance, num_verts, storage_order, new_coords, 
                      dim * num_verts, &new_vertex_handles_out,
                      &alloc, &junk, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::createEnt( EntityTopology topology,
                  const EntityHandle* entities,
                  int entities_size,
                  EntityHandle& new_entity_out )
{
  int err, status;
  iMesh_createEnt( mInstance, topology, entities, entities_size, &new_entity_out, &status, &err );
  if (err == iBase_SUCCESS) 
    if (status == iBase_ALREADY_EXISTED || status == iBase_CREATION_FAILED)
      err = iBase_ENTITY_CREATION_ERROR;
  return (Error)err;
}

inline iMesh::Error
iMesh::createEntArr( EntityTopology new_entity_topology,
                     const EntityHandle* lower_order_handles,
                     int lower_order_handles_size,
                     EntityHandle* new_entity_handles_out )
{
  std::vector<int> status(lower_order_handles_size, iBase_CREATION_FAILED);
  int err, alloc = lower_order_handles_size, size = 0;
  int status_alloc = status.size(), junk, *stat_ptr = &status[0];
  
  iMesh_createEntArr( mInstance, new_entity_topology, 
                      lower_order_handles, lower_order_handles_size, 
                      &new_entity_handles_out, &alloc, &size,
                      &stat_ptr, &status_alloc, &junk, &err );
  
  if (iBase_SUCCESS == err) {
    for (int i = 0; i < size; ++i)
      if (status[i] != iBase_NEW && status[i] != iBase_CREATED_DUPLICATE)
        err = iBase_ENTITY_CREATION_ERROR;
  }
  
  if (iBase_SUCCESS != err) {
    int err2, w = 0;
    for (int r = 0; r < size; ++r) 
      if (status[r] == iBase_NEW || status[r] == iBase_CREATED_DUPLICATE) 
        new_entity_handles_out[w++] = new_entity_handles_out[r];
    iMesh_deleteEntArr( mInstance, new_entity_handles_out, w, &err2 );
  }
  
  return (Error)err;
}

inline iMesh::Error
iMesh::deleteEnt( EntityHandle handle )
{
  int err;
  iMesh_deleteEnt( mInstance, handle, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::deleteEntArr( const EntityHandle* entity_handles, int num_handles )
{
  int err;
  iMesh_deleteEntArr( mInstance, entity_handles, num_handles, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::getAdjEntities( EntitySetHandle set,
                       EntityType type_requestor,
                       EntityTopology topo_requestor,
                       EntityType type_requested,
                       std::vector<EntityHandle>& adj_entity_handles,
                       std::vector<int>& offset ) const
{
  std::vector<EntityHandle> entities;
  Error err = getEntities( set, type_requestor, topo_requestor, entities );
  if (iBase_SUCCESS != err)
    return err;
  
  offset.resize( entities.size() + 1 );
  return getEntArrAdj( &entities[0], entities.size(), type_requested,
                       adj_entity_handles, &offset[0] );
}
    
inline iMesh::Error 
iMesh::initEntIter( EntitySetHandle set,
                    EntityType requested_type,
                    EntityTopology requested_topo,
                    iMesh::EntIter& iter )
{
  int err;
  iter.mInstance = mInstance;
  iMesh_initEntIter( mInstance, set, requested_type, requested_topo,
                     &iter.mHandle, &err );
  return (Error)err;
}

inline iMesh::Error 
iMesh::initEntArrIter( EntitySetHandle set,
                       EntityType requested_type,
                       EntityTopology requested_topo,
                       int requested_array_size,
                       iMesh::EntArrIter& iter )
{
  int err;
  iter.mInstance = mInstance;
  iter.mSize = requested_array_size;
  iMesh_initEntArrIter( mInstance, set, requested_type, requested_topo,
                        requested_array_size, &iter.mHandle, &err );
  return (Error)err;
}

inline 
iMesh::EntArrIter::~EntArrIter()
{
  int err;
  if (mHandle != 0) {
    iMesh_endEntArrIter( mInstance, mHandle, &err );
    mHandle = 0;
  }
}

inline 
iMesh::EntIter::~EntIter()
{
  int err;
  if (mHandle != 0) {
    iMesh_endEntIter( mInstance, mHandle, &err );
    mHandle = 0;
  }
}

inline iMesh::Error
iMesh::EntArrIter::getNext( EntityHandle* entity_handles,
                            int& size_out,
                            bool& has_more_data_out )
{
  int err, alloc = mSize, has_data;
  iMesh_getNextEntArrIter( mInstance, mHandle, &entity_handles, &alloc,
                           &size_out, &has_data, &err );
  has_more_data_out = (has_data != 0);
  return (Error)err;
}

inline iMesh::Error
iMesh::EntIter::getNext( EntityHandle& handle_out, bool& has_more_data_out )
{
  int err, has_data;
  iMesh_getNextEntIter( mInstance, mHandle, &handle_out, &has_data, &err );
  has_more_data_out = (has_data != 0);
  return (Error)err;
}

inline iMesh::Error
iMesh::EntArrIter::reset()
{
  int err;
  iMesh_resetEntArrIter( mInstance, mHandle, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::EntIter::reset()
{
  int err;
  iMesh_resetEntIter( mInstance, mHandle, &err );
  return (Error)err;
}

inline iMesh::Error
iMesh::getEntTopo( EntityHandle handle, EntityTopology& topo_out ) const
{
  int err, result;
  iMesh_getEntTopo( mInstance, handle, &result, &err );
  topo_out = (EntityTopology)result;
  return (Error)err;
}

inline iMesh::Error
iMesh::getEntArrTopo( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      EntityTopology* topos_out ) const 
{
  int err, alloc = entity_handles_size, junk, *ptr;
  std::vector<int> storage;
  if (sizeof(EntityTopology) == sizeof(int))
    ptr = reinterpret_cast<int*>(topos_out);
  else {
    storage.resize( entity_handles_size );
    ptr = &storage[0];
  }
  
  iMesh_getEntArrTopo( mInstance, entity_handles, entity_handles_size,
                       &ptr, &alloc, &junk, &err );
  
  if (sizeof(EntityTopology) != sizeof(int))
    for (int i = 0; i < entity_handles_size; ++i)
      topos_out[i] = (EntityTopology)storage[i];
      
  return (Error)err;
}
 
inline iMesh::Error
iMesh::getEntType( EntityHandle handle, EntityType& type_out ) const
{
  int err, result;
  iMesh_getEntType( mInstance, handle, &result, &err );
  type_out = (EntityType)result;
  return (Error)err;
}

inline iMesh::Error
iMesh::getEntArrType( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      EntityType* types_out ) const
{
  int err, alloc = entity_handles_size, junk, *ptr;
  std::vector<int> storage;
  if (sizeof(EntityType) == sizeof(int))
    ptr = reinterpret_cast<int*>(types_out);
  else {
    storage.resize( entity_handles_size );
    ptr = &storage[0];
  }
  
  iMesh_getEntArrType( mInstance, entity_handles, entity_handles_size,
                       &ptr, &alloc, &junk, &err );
  
  if (sizeof(EntityType) != sizeof(int))
    for (int i = 0; i < entity_handles_size; ++i)
      types_out[i] = (EntityType)storage[i];
      
  return (Error)err;
}

inline iMesh::Error
iMesh::getEntAdj( EntityHandle handle, 
                  EntityType type_requested,
                  std::vector<EntityHandle>& adj_entities_out ) const
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  EntityHandle* ptr = &adj_entities_out[0];
  iMesh_getEntAdj( mInstance, handle, type_requested, 
                   &ptr, &alloc, &size, &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iMesh_getEntAdj( mInstance, handle, type_requested, 
                     &ptr, &alloc, &size, &err );
  }
  
  return (Error)err;
}

inline iMesh::Error
iMesh::getEntArrAdj( const EntityHandle* entity_handles,
                           int entity_handles_size,
                           EntityType type_requested,
                           std::vector<EntityHandle>& adj_entities_out,
                           int* offsets_out ) const
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12*entity_handles_size);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  int off_alloc = entity_handles_size+1, junk;
  EntityHandle* ptr = &adj_entities_out[0];
  iMesh_getEntArrAdj( mInstance, entity_handles, entity_handles_size, 
                      type_requested,
                      &ptr, &alloc, &size,
                      &offsets_out, &off_alloc, &junk,
                      &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iMesh_getEntArrAdj( mInstance, entity_handles, entity_handles_size, 
                        type_requested,
                        &ptr, &alloc, &size,
                        &offsets_out, &off_alloc, &junk,
                        &err );
  }
  
  return (Error)err;
}

                   
inline iMesh::Error
iMesh::getEnt2ndAdj( EntityHandle handle, 
                     EntityType bridge_entity_type,
                     EntityType type_requested,
                     std::vector<EntityHandle> adj_entities_out ) const
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  EntityHandle* ptr = &adj_entities_out[0];
  iMesh_getEnt2ndAdj( mInstance, handle, bridge_entity_type, type_requested, 
                      &ptr, &alloc, &size, &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iMesh_getEnt2ndAdj( mInstance, handle, bridge_entity_type, type_requested, 
                       &ptr, &alloc, &size, &err );
  }
  
  return (Error)err;
}

inline iMesh::Error
iMesh::getEntArr2ndAdj( const EntityHandle* entity_handles,
                        int entity_handles_size,
                        EntityType order_key,
                        EntityType type_requested,
                        std::vector<EntityHandle>& adj_entities_out,
                        int* offsets_out ) const
{
  if (adj_entities_out.capacity() == 0) 
    adj_entities_out.resize(12*entity_handles_size);
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
  
  int err, alloc = adj_entities_out.size(), size = 0;
  int off_alloc = entity_handles_size+1, junk;
  EntityHandle* ptr = &adj_entities_out[0];
  iMesh_getEntArr2ndAdj( mInstance, entity_handles, entity_handles_size, 
                         order_key, type_requested,
                         &ptr, &alloc, &size,
                         &offsets_out, &off_alloc, &junk,
                         &err );
  adj_entities_out.resize(size);
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
    alloc = adj_entities_out.size();
    ptr = &adj_entities_out[0];
    iMesh_getEntArr2ndAdj( mInstance, entity_handles, entity_handles_size, 
                           order_key, type_requested,
                           &ptr, &alloc, &size,
                           &offsets_out, &off_alloc, &junk,
                           &err );
  }
  
  return (Error)err;
}


inline iMesh::Error
iMesh::getAdjEntIndices( EntitySetHandle set,
                         EntityType requestor_type,
                         EntityTopology requestor_topo,
                         EntityType type_requested,
                         std::vector<EntityHandle>& entities_out,
                         std::vector<EntityHandle>& adj_entities_out,
                         std::vector<int>& adj_indices_out,
                         std::vector<int>& offsets_out ) const
{
    // if input vects have no allocated space, allocate some so
    // we don't accidentally ask the impl to allocate an array
  if (entities_out.capacity() == 0) {
    Error err2;
    int count;
    if (requestor_topo == iMesh_ALL_TOPOLOGIES)
      err2 = getNumOfType( set, requestor_type, count );
    else
      err2 = getNumOfTopo( set, requestor_topo, count );
    if (err2 != iBase_SUCCESS)
      return err2;
    entities_out.resize( count );
  }
  else
    entities_out.resize( entities_out.capacity() );
    
  if (adj_entities_out.capacity() == 0)
    adj_entities_out.resize( 12*entities_out.size() );
  else
    adj_entities_out.resize( adj_entities_out.capacity() );
    
  if (adj_indices_out.capacity() == 0)
    adj_indices_out.resize( 12*entities_out.size() );
  else
    adj_indices_out.resize( adj_indices_out.capacity() );
    
  if (offsets_out.capacity() == 0)
    offsets_out.resize( entities_out.size() + 1 );
  else
    offsets_out.resize( offsets_out.capacity() );

  int err;
  int     entities_size = 0,     entities_alloc =     entities_out.size();
  int adj_entities_size = 0, adj_entities_alloc = adj_entities_out.size();
  int  adj_indices_size = 0,  adj_indices_alloc =  adj_indices_out.size();
  int      offsets_size = 0,      offsets_alloc =      offsets_out.size();
  EntityHandle*     entities_ptr = &    entities_out[0];
  EntityHandle* adj_entities_ptr = &adj_entities_out[0];
  int         *  adj_indices_ptr = & adj_indices_out[0];
  int         *      offsets_ptr = &     offsets_out[0];
  
  iMesh_getAdjEntIndices( mInstance, 
                          set, requestor_type, requestor_topo,
                          type_requested,
                          &    entities_ptr, &    entities_alloc, &    entities_size,
                          &adj_entities_ptr, &adj_entities_alloc, &adj_entities_size,
                          & adj_indices_ptr, & adj_indices_alloc, & adj_indices_size,
                          &     offsets_ptr, &     offsets_alloc, &     offsets_size,
                          &err );
  
      entities_out.resize(     entities_size );
  adj_entities_out.resize( adj_entities_size );
   adj_indices_out.resize(  adj_indices_size );
       offsets_out.resize(      offsets_size );
  
  if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
        entities_alloc =     entities_out.size();
    adj_entities_alloc = adj_entities_out.size();
     adj_indices_alloc =  adj_indices_out.size();
         offsets_alloc =      offsets_out.size();
        entities_ptr = &    entities_out[0];
    adj_entities_ptr = &adj_entities_out[0];
     adj_indices_ptr = & adj_indices_out[0];
         offsets_ptr = &     offsets_out[0];

    iMesh_getAdjEntIndices( mInstance, 
                            set, requestor_type, requestor_topo,
                            type_requested,
                            &    entities_ptr, &    entities_alloc, &    entities_size,
                            &adj_entities_ptr, &adj_entities_alloc, &adj_entities_size,
                            & adj_indices_ptr, & adj_indices_alloc, & adj_indices_size,
                            &     offsets_ptr, &     offsets_alloc, &     offsets_size,
                            &err );
  }
  
  return (Error)err;
}

#endif

