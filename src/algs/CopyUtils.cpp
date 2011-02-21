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
