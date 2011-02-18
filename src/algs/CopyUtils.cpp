#include "CopyUtils.hpp"

#include "ArrayManager.hpp"
#include "MKException.hpp"
#include "SimpleArray.hpp"
#include "MBCN.h"
#include <iMesh_extensions.h>
#include <vector>

void iMesh_getStructure(iMesh_Instance instance, iBase_EntitySetHandle set,
                        iBase_EntityHandle **ents,
                        int *ents_allocated,
                        int *ents_size,
                        iBase_EntityHandle **unique_adj,
                        int *unique_adj_allocated,
                        int *unique_adj_size,
                        int **indices,
                        int *indices_allocated,
                        int *indices_size,
                        int **offsets,
                        int *offsets_allocated,
                        int *offsets_size,
                        int *err)
{
  // 1) Get source entities, making sure verts are first
  int num;
  iMesh_getNumOfTypeRec(instance, set, iBase_ALL_TYPES, true, &num, err);
  if (*err != iBase_SUCCESS) return;

  ALLOC_CHECK_ARRAY(ents, num);
  ALLOC_CHECK_ARRAY(offsets, num+1);

  iBase_EntityHandle *block = *ents;
  int block_alloc = *ents_allocated, block_size, num_verts = 0;
  for (int t = iMesh_POINT; t < iMesh_ALL_TOPOLOGIES && block_alloc; ++t) {
    iMesh_getEntitiesRec(instance, set, iBase_ALL_TYPES, t, true,
                         &block, &block_alloc, &block_size, err);
    if (*err != iBase_SUCCESS) return;

    block_alloc -= block_size;
    block += block_size;
    if (t == iMesh_POINT)
      num_verts = block_size;
  }

  // 2) Get verts adjacent to all source entitites (verts are adj to themselves)
  std::vector<iBase_EntityHandle> all_adj(*ents, *ents+num_verts);

  // first, fill the vertex-vertex adjacencies
  for (int i = 0; i < num_verts; ++i)
    (*offsets)[i] = i;

  iBase_EntityHandle *tmp_adj = NULL;
  int tmp_adj_alloc = 0, tmp_adj_size;
  int *tmp_off = *offsets + num_verts;
  int tmp_off_alloc = *offsets_allocated - num_verts, tmp_off_size;
  iMesh_getEntArrAdj(instance, *ents+num_verts, *ents_size-num_verts,
                     iBase_VERTEX, &tmp_adj, &tmp_adj_alloc, &tmp_adj_size,
                     &tmp_off, &tmp_off_alloc, &tmp_off_size, err);
  if (*err != iBase_SUCCESS) return;

  // shift all the offsets to account for vertices
  for(int i = num_verts; i < *offsets_size; ++i)
    (*offsets)[i] += num_verts;

  all_adj.reserve(all_adj.size() + tmp_adj_size);
  all_adj.insert(all_adj.end(), tmp_adj, tmp_adj+tmp_adj_size);
  free(tmp_adj);

  // 3) Get unique adjacent vertices and offsets
  // TODO: this might put unncessary restrictions on the size of the input array
  ALLOC_CHECK_ARRAY(unique_adj, all_adj.size());
  ALLOC_CHECK_ARRAY(indices, all_adj.size());

  std::copy(all_adj.begin(), all_adj.end(), *unique_adj);
  std::sort(*unique_adj, *unique_adj+*unique_adj_size);
  *unique_adj_size = std::unique(*unique_adj, *unique_adj+*unique_adj_size) -
    *unique_adj;

  for (size_t i = 0; i < all_adj.size(); ++i) {
    (*indices)[i] = std::lower_bound(*unique_adj, *unique_adj+*unique_adj_size,
                                     all_adj[i]) - *unique_adj;
  }

  KEEP_ARRAY(ents);
  KEEP_ARRAY(unique_adj);
  KEEP_ARRAY(indices);
  KEEP_ARRAY(offsets);
}

void connect_the_dots(iMesh_Instance imeshImpl, iBase_EntityHandle *ents,
                      int ents_size, iBase_TagHandle local_tag, int *indices,
                      int *offsets, iBase_EntityHandle *verts)
{
  int err;

  SimpleArray<int> topos;
  iMesh_getEntArrTopo(imeshImpl, ents, ents_size, ARRAY_INOUT(topos), &err);
  check_error(imeshImpl, err);
  
  // scan forward to first non-vertex
  int pos = 0;
  while (pos < topos.size() && iMesh_POINT == topos[pos])
    pos++;
  if (pos == topos.size()) return;

  // for each run of same size & type
  std::vector<iBase_EntityHandle> connect, new_ents;
  std::vector<int> status;
  int begin, end = pos;
  while (end < ents_size) {
    // get next run; end points to start of *next* element,
    // or ents_size if no elems left
    begin = end++;

    int topo = topos[begin];
    int vtx_per_ent = offsets[end] - offsets[begin];
    while (end < ents_size &&
           topos[end] == topo &&
           offsets[end+1] - offsets[end] == vtx_per_ent)
      end++;
    int num_ents = end - begin;

    int mbcn_type;
    int num_corner_verts;
    iMesh_MBCNType(topo, &mbcn_type);
    MBCN_VerticesPerEntity(mbcn_type, &num_corner_verts);

    // build vector of vtx handles
    connect.resize(vtx_per_ent * num_ents);
    for (size_t i = 0; i < connect.size(); i++)
      connect[i] = verts[indices[offsets[begin] + i]];

    // create entities
    new_ents.resize(num_ents);

    if (num_corner_verts == vtx_per_ent) {
      status.resize(num_ents);

      iBase_EntityHandle *new_ents_ptr = &new_ents[0];
      int new_ents_alloc = num_ents, new_ents_size;
      int *status_ptr = &status[0];
      int status_alloc = num_ents, status_size;
      iMesh_createEntArr(imeshImpl, topo, &connect[0], connect.size(),
                         &new_ents_ptr, &new_ents_alloc, &new_ents_size,
                         &status_ptr, &status_alloc, &status_size, &err);
      check_error(imeshImpl, err);
    }
    else {
      // use single-entity function in this case, entity might have higher-order
      // nodes (imesh fcn doesn't have argument for # entities)
      for (int i = 0; i < num_ents; i++) {
        int status;
        iMesh_createEnt(imeshImpl, topo, &connect[i*vtx_per_ent],
                        vtx_per_ent, &new_ents[i], &status, &err);
        check_error(imeshImpl, err);
      }
    }

    // set the local copy tags
    iMesh_setEHArrData(imeshImpl, &ents[begin], num_ents, local_tag, 
                       &new_ents[0], num_ents, &err);
    check_error(imeshImpl, err);
  }
}
