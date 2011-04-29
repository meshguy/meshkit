#include "CopyUtils.hpp"

#include <iMesh_extensions.h>
#include <vector>
#include <algorithm>

#define CHKERR(err) do {                        \
  if ((err) != iBase_SUCCESS)                   \
    return iBase_ErrorType(err);                \
  } while(false)


// TODO: eliminate this function
iMesh::Error iMesh_getStructure(iMesh_Instance instance,
                                iBase_EntitySetHandle set,
                                std::vector<iBase_EntityHandle> &ents,
                                std::vector<iBase_EntityHandle> &unique_adj,
                                std::vector<int> &indices,
                                std::vector<int> &offsets)
{
  // 1) Get source entities, making sure verts are first
  int num;
  int err;
  iMesh_getNumOfTypeRec(instance, set, iBase_ALL_TYPES, true, &num, &err);
  CHKERR(err);

  ents.resize(num);
  offsets.resize(num+1);

  iBase_EntityHandle *block = &ents[0];
  int block_alloc = ents.size(), block_size, num_verts = 0;
  for (int t = iMesh_POINT; t < iMesh_ALL_TOPOLOGIES && block_alloc; ++t) {
    iMesh_getEntitiesRec(instance, set, iBase_ALL_TYPES, t, true,
                         &block, &block_alloc, &block_size, &err);
    CHKERR(err);

    block_alloc -= block_size;
    block += block_size;
    if (t == iMesh_POINT)
      num_verts = block_size;
  }

  // 2) Get verts adjacent to all source entitites (verts are adj to themselves)
  std::vector<iBase_EntityHandle> all_adj(ents.begin(), ents.begin()+num_verts);

  // first, fill the vertex-vertex adjacencies
  for (int i = 0; i < num_verts; ++i)
    offsets[i] = i;

  iBase_EntityHandle *tmp_adj = NULL;
  int tmp_adj_alloc = 0, tmp_adj_size;
  int *tmp_off = &offsets[num_verts];
  int tmp_off_alloc = offsets.size() - num_verts, tmp_off_size;
  iMesh_getEntArrAdj(instance, &ents[num_verts], ents.size()-num_verts,
                     iBase_VERTEX, &tmp_adj, &tmp_adj_alloc, &tmp_adj_size,
                     &tmp_off, &tmp_off_alloc, &tmp_off_size, &err);
  CHKERR(err);

  // shift all the offsets to account for vertices
  for(unsigned int i = num_verts; i < offsets.size(); ++i)
    offsets[i] += num_verts;

  all_adj.reserve(all_adj.size() + tmp_adj_size);
  all_adj.insert(all_adj.end(), tmp_adj, tmp_adj+tmp_adj_size);
  free(tmp_adj);

  // 3) Get unique adjacent vertices and offsets
  unique_adj.resize(all_adj.size());
  indices.resize(all_adj.size());
  std::copy(all_adj.begin(), all_adj.end(), unique_adj.begin());
  std::sort(unique_adj.begin(), unique_adj.end());

  size_t unique_size;
  unique_size = std::unique(unique_adj.begin(), unique_adj.end()) -
    unique_adj.begin();
  unique_adj.resize(unique_size);

  for (size_t i = 0; i < all_adj.size(); ++i) {
    indices[i] = std::lower_bound(unique_adj.begin(), unique_adj.end(),
                                  all_adj[i]) - unique_adj.begin();
  }

  return iBase_SUCCESS;
}
