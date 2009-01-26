#include "MKUtils.hpp"
#include "iMesh.h"
#include <vector>
#include <iostream>
#include <cstring>

#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

//! assign a global id space, for largest-dimension or all entities (and
//! in either case for vertices too)
int MKUtils::assign_global_ids(iBase_EntitySetHandle this_set,
                               const int dimension, 
                               const int start_id,
                               const bool largest_dim_only,
                               const bool parallel,
                               const char *tag_name) 
{
  int err;
  std::vector<unsigned char> pstatus;
  iBase_EntityHandle *entities[4] = {NULL, NULL, NULL, NULL};
  int nprocs = 1, nrank = 0;
    // use nrank in init so compiler doesn't warn about unused var
  int entities_alloc[4] = {nrank, nrank, nrank, nrank}, entities_size[4];
  std::vector<int> num_elements(nprocs*4);

  for (int dim = 0; dim <= dimension; dim++) {
    if (dim == 0 || !largest_dim_only || dim == dimension) {
      iMesh_getEntities(imeshImpl, this_set, dim, iMesh_ALL_TOPOLOGIES,
                        &entities[dim], &entities_alloc[dim], &entities_size[dim],
                        &err); 
      ERRORR("Failed to get entities in assign_global_ids.", err);
    }
  }
  

#ifdef MPI
  if (parallel) {
      // need to filter out non-locally-owned entities!!!
    err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    ERRORR("Couldn't get number of procs.", err);
    err = MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
    ERRORR("Couldn't get proc rank.", err);

    pstatus.resize(entities[dim].size());
    err = mbImpl->tag_get_data(pstatus_tag(), entities[dim], &pstatus[0]);
    ERRORR("Failed to get pstatus in assign_global_ids.", err);
    
    int local_num_elements[4];
    MBRange dum_range;
    MBRange::iterator rit;
    unsigned int i;
    for (rit = entities[dim].begin(); rit != entities[dim].end(); rit++, i++)
      if (pstatus[i] & PSTATUS_NOT_OWNED)
        dum_range.insert(*rit);
    entities[dim] = entities[dim].subtract(dum_range);
    
    local_num_elements[dim] = entities[dim].size();
  
      // communicate numbers
    if (nprocs > 1) {
      int retval = MPI_Alltoall(local_num_elements, 4, MPI_INT,
                                &num_elements[0], nprocs*4, 
                                MPI_INT, MPI_COMM_WORLD);
      if (0 != retval) {
        err = iBase_FAILURE;
        return;
      }
    }
  }
  else
#endif
    for (int dim = 0; dim < 4; dim++) num_elements[dim] = entities_size[dim];
  
    // my entities start at one greater than total_elems[d]
  int total_elems[4] = {start_id, start_id, start_id, start_id};

#ifdef MPI
  if (parallel) {
    for (int dim = 0; dim < 4; dim++) {
      if (largest_dim_only && dim != dimension) continue;
        
      for (unsigned int proc = 0; proc < nrank; proc++) {
        total_elems[dim] += num_elements[4*proc + dim];
      }
    }
  }
    
#endif

    //assign global ids now
  iBase_TagHandle gid_tag;
  iMesh_getTagHandle(imeshImpl, tag_name, &gid_tag, &err, 9);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(imeshImpl, tag_name, 1, iBase_INTEGER, 
                    &gid_tag, &err, strlen(tag_name));
    ERRORR("Couldn't create global id tag", err);
  }
  
  for (int dim = 0; dim < 4; dim++) {
    if (!entities[dim] || 0 == entities_size[dim]) continue;
    num_elements.resize(entities_size[dim]);
    
    for (int i = 0; i < entities_size[dim]; i++)
      num_elements[i] = total_elems[dim]++;
    
    iMesh_setIntArrData(imeshImpl, entities[dim], entities_size[dim], gid_tag, 
                        &num_elements[0], entities_size[dim], &err); 
    ERRORR("Failed to set global id tag in assign_global_ids.", err);
  }
  
  return iBase_SUCCESS;
}

