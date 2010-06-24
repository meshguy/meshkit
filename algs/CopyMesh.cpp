#include "CopyMesh.hpp"
#include "CopyVerts.hpp"
#include "../utils/ArrayManager.hpp"
#include "MBCN.h"
#include <stdlib.h>
#include <algorithm>
#include <functional>

/*
 * - call update_ce_lists to update copySets or expandSets
 */

/**\brief Get the entities copied from a set
 *
 * Given a set and a tag, return the values of the tag set on entities from that
 * set (i.e. the copied entities).
 * \param imeshImpl the iMesh instance handle
 * \param set the set containing the entities we want
 * \param local_tag the tag relating source and target entities
 * \param ents a vector which will hold our copied entities
 * \return an ITAPS error code
 */
static
int get_copied_ents(iMesh_Instance imeshImpl, iBase_EntitySetHandle set,
                    iBase_TagHandle local_tag,
                    std::vector<iBase_EntityHandle> &ents)
{
  int err;

  iBase_EntityHandle *tmp_ents = NULL;
  int tmp_ents_alloc = 0, tmp_ents_size;
  iMesh_getEntities(imeshImpl, set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                    &tmp_ents, &tmp_ents_alloc, &tmp_ents_size, &err);
  ERRORR("Failed to get ceSet entities.", iBase_FAILURE);

  ents.reserve(tmp_ents_size);

  // get copied entities
  for (int i = 0; i < tmp_ents_size; i++) {
    iBase_EntityHandle eh_tag;
    iMesh_getEHData(imeshImpl, tmp_ents[i], local_tag, &eh_tag, &err);
    if (iBase_SUCCESS == err && eh_tag)
      ents.push_back(eh_tag);
  }
  free(tmp_ents);
  return iBase_SUCCESS;
}

/**\brief Create a copy set if one doesn't exist yet
 *
 * Given a source set and a tag, create a destination (copy) set unless one
 * already exists.
 * \param imeshImpl the iMesh instance handle
 * \param local_tag the tag relating source and target sets
 * \param src the source set
 * \param dest pointer to the destination set
 * \return an ITAPS error code
 */
static
int get_dest_set(iMesh_Instance imeshImpl, iBase_TagHandle local_tag,
                 iBase_EntitySetHandle src, iBase_EntitySetHandle *dest)
{
  int err;

  if (dest == NULL) return iBase_FAILURE;

  iMesh_getEntSetEHData(imeshImpl, src, local_tag,
                        reinterpret_cast<iBase_EntityHandle*>(dest), &err);

  if (err != iBase_SUCCESS) {
    iMesh_createEntSet(imeshImpl, false, dest, &err);
    ERRORR("Failed to create copied set.", iBase_FAILURE);
    iMesh_setEntSetEHData(imeshImpl, src, local_tag,
                          reinterpret_cast<iBase_EntityHandle>(*dest), &err);
    ERRORR("Failed to tag copied set.", iBase_FAILURE);
  }

  return iBase_SUCCESS;
}

/**\brief Add copied entities/sets recursively
 *
 * Helper function for process_ce_sets. This adds any entities directly
 * contained in the current set to the dest set and then loops through the
 * contained sets and adds the children if they are also CE sets. If not, it
 * adds the entities in the child set and recurses down a level.
 * \param imeshImpl the iMesh instance handle
 * \param src the source set
 * \param current the child set to examine (start with current == src)
 * \param local_tag the tag relating source and target sets
 * \return an ITAPS error code
 */
static
int process_ce_subsets(iMesh_Instance imeshImpl, iBase_EntitySetHandle src,
                       iBase_EntitySetHandle current,
                       const std::set<iBase_EntitySetHandle> &cesets,
                       iBase_TagHandle local_tag) 
{
  int err;
  iBase_EntitySetHandle dest;

  // First, add entities directly contained in this set.
  std::vector<iBase_EntityHandle> tmp_tags;
  err = get_copied_ents(imeshImpl, current, local_tag, tmp_tags);
  ERRORR("Failed to get copied ents.", iBase_FAILURE);

  if (!tmp_tags.empty()) {
    err = get_dest_set(imeshImpl, local_tag, src, &dest);
    ERRORR("Failed to get copied set.", iBase_FAILURE);
    iMesh_addEntArrToSet(imeshImpl, &tmp_tags[0], tmp_tags.size(), dest,
                         &err);
    ERRORR("Failed to add copied entities to ce set.", iBase_FAILURE);
  }

  // Next, start looking at children.
  iBase_EntitySetHandle *children = NULL;
  int children_alloc = 0, children_size;
  iMesh_getEntSets(imeshImpl, current, 1, &children, &children_alloc,
                   &children_size, &err);
  ERRORR("Failed to get ceSet sets.", iBase_FAILURE);

  for (int i = 0; i < children_size; i++) {

    // If this child set is one of our cesets, add just the set...
    if (cesets.find(children[i]) != cesets.end()) {
      err = get_dest_set(imeshImpl, local_tag, src, &dest);
      ERRORR("Failed to get copied set.", iBase_FAILURE);

      if (src == dest) continue;

      iMesh_addEntSet(imeshImpl, children[i], dest, &err);
      ERRORR("Failed to add set to ce set.", iBase_FAILURE);
    }

    // ... otherwise, add the entities and recurse into the next level of
    // children.
    else {
      err = process_ce_subsets(imeshImpl, src, children[i], cesets, local_tag);
      ERRORR("Failed to process contained sets.", iBase_FAILURE);
    }
  }

  free(children);
  return iBase_SUCCESS;
}

// xxx - still need to do unique tags
// - for each unique tag
//   . if this set has the tag, add/append to get unique value
int process_ce_sets(iMesh_Instance imeshImpl,
                    const std::set<iBase_EntitySetHandle> &cesets,
                    iBase_TagHandle local_tag) 
{
  int err;
  std::set<iBase_EntitySetHandle>::const_iterator src;
  for (src = cesets.begin(); src != cesets.end(); ++src) {
    err = process_ce_subsets(imeshImpl, *src, *src, cesets, local_tag);
    ERRORR("Failed to process contained sets.", iBase_FAILURE);
  }
  
  return iBase_SUCCESS;
}

int tag_copy_sets(iMesh_Instance imeshImpl, iBase_TagHandle copyTag,
                  const std::set<iBase_EntitySetHandle> &copySets,
                  iBase_TagHandle tag, const char *tag_val)
{
  int err = iBase_SUCCESS;

  int tag_size;
  iMesh_getTagSizeBytes(imeshImpl, tag, &tag_size, &err);
  ERRORR("Failed to get tag size.", err);

  // allocate temp space for tag value
  std::vector<char> value;
  value.resize(tag_size);
  char *value_ptr = &value[0];
  int value_alloc = tag_size, value_size;
  
  // for each orig copy set with this tag, copy it to its copy
  std::set<iBase_EntitySetHandle>::iterator set;
  for (set = copySets.begin(); set != copySets.end(); ++set) {
    // get the tag value
    iMesh_getEntSetData(imeshImpl, *set, tag, 
                        &value_ptr, &value_alloc, &value_size, &err);
    if (iBase_TAG_NOT_FOUND == err) {
      // reset, just in case we return after this setting of err
      err = iBase_SUCCESS;
      continue;
    }
    ERRORR("Problem getting copy tag for set.", err);

    // compare to tag value if necessary
    if (tag_val && strncmp(tag_val, value_ptr, tag_size)) {
      // reset, just in case we return after this setting of err
      err = iBase_SUCCESS;
      continue;
    }      
      
    // if we got here, we should set the tag on the copy; get the copy
    iBase_EntitySetHandle copy_set;
    iMesh_getEntSetEHData(imeshImpl, *set, copyTag, 
                          reinterpret_cast<iBase_EntityHandle*>(&copy_set), 
                          &err);
    if (iBase_TAG_NOT_FOUND == err) {
      // we (probably) didn't copy anything from this set, so ignore it
      err = iBase_SUCCESS;
      continue;
    }
    ERRORR("Didn't get copied set from orig copy set.", err);

    if (copy_set != *set)
      iMesh_setEntSetData(imeshImpl, copy_set, tag, value_ptr, tag_size, &err);
    ERRORR("Failed to set data on copied set.", err);
  }
  
  return err;
}

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
  int block_alloc = *ents_allocated, block_size, num_verts;
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

CopyMesh::CopyMesh(iMesh_Instance impl) 
  : imeshImpl(impl), updatedCELists(false), copyTag(impl, "__CopyMeshTag")
{}

CopyMesh::~CopyMesh()
{
  std::vector<tag_data>::iterator i;
  for (i = copyTags.begin(); i != copyTags.end(); ++i)
    free(i->value);

  for (i = expandTags.begin(); i != expandTags.end(); ++i)
    free(i->value);
}

int CopyMesh::add_copy_expand_list(iBase_EntitySetHandle *ce_sets, int num_ce_sets,
                                   int copy_or_expand) 
{
  int err;
  for (int i = 0; i < num_ce_sets; i++) {
    if (copy_or_expand == COPY) 
      copySets.insert(ce_sets[i]);
    else if (copy_or_expand == EXPAND) 
      expandSets.insert(ce_sets[i]);
    else {
      err = iBase_INVALID_ARGUMENT;
      ERRORR("Needs to be copy or expand.", iBase_INVALID_ARGUMENT);
    }
  }

  return iBase_SUCCESS;
}

int CopyMesh::add_copy_tag(iBase_TagHandle tag_handle, const char *tag_val)
{
  assert(tag_handle != NULL);
  char *tmp = NULL;

  if (tag_val) {
    int err;
    int tag_size;
    iMesh_getTagSizeBytes(imeshImpl, tag_handle, &tag_size, &err);
    if (iBase_SUCCESS != err)
      ERRORR("Failed to get size of tag", iBase_FAILURE);
    tmp = static_cast<char*>(malloc(tag_size));
    memcpy(tmp, tag_val, tag_size);
  }

  copyTags.push_back(tag_data(tag_handle, tmp));
  return iBase_SUCCESS;
}

int CopyMesh::add_expand_tag(iBase_TagHandle tag_handle, const char *tag_val)
{
  assert(tag_handle != NULL);
  char *tmp = NULL;

  if (tag_val) {
    int err;
    int tag_size;
    iMesh_getTagSizeBytes(imeshImpl, tag_handle, &tag_size, &err);
    if (iBase_SUCCESS != err)
      ERRORR("Failed to get size of tag", iBase_FAILURE);
    tmp = static_cast<char*>(malloc(tag_size));
    memcpy(tmp, tag_val, tag_size);
  }

  expandTags.push_back(tag_data(tag_handle, tmp));
  return iBase_SUCCESS;
}



/* \brief Copy all the entities in the set
 */
int CopyMesh::copy_entities(iBase_EntitySetHandle set_handle,
                            iBase_EntityHandle **new_ents,
                            int *new_ents_allocated,
                            int *new_ents_size) 
{
  iBase_EntityHandle *ents = NULL;
  int ents_alloc = 0, ents_size;
  int err;
  iMesh_getEntities(imeshImpl, set_handle, iBase_ALL_TYPES,
                    iMesh_ALL_TOPOLOGIES, &ents, &ents_alloc, &ents_size,
                    &err);
  ERRORR("Couldn't get entities in set.", iBase_FAILURE);
  
  int result = copy_move_entities(ents, ents_size, NULL,
                                  new_ents, new_ents_allocated, new_ents_size);

  free(ents);
  return result;
}

  
/* \brief Copy all the entities
 */
int CopyMesh::copy_entities(iBase_EntityHandle *ent_handles,
                            int num_ents,
                            iBase_EntityHandle **new_ents,
                            int *new_ents_allocated,
                            int *new_ents_size) 
{
  return copy_move_entities(ent_handles, num_ents, NULL,
                            new_ents, new_ents_allocated, new_ents_size);;
}

/* \brief Copy and move all the entities
 */
int CopyMesh::copy_move_entities(iBase_EntitySetHandle set_handle,
                                 const double *dx,
                                 iBase_EntityHandle **new_ents,
                                 int *new_ents_alloc,
                                 int *new_ents_size,
                                 bool do_merge) 
{
  int err;
  iBase_EntityHandle *ents = NULL;
  int ents_alloc = 0, ents_size;
  iMesh_getEntitiesRec(imeshImpl, set_handle, 
                       iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, true,
                       &ents, &ents_alloc, &ents_size, &err);
  ERRORR("Failed to get entities from set recursively.", err);
  
  int result = copy_move_entities(ents, ents_size, dx, 
                                  new_ents, new_ents_alloc, new_ents_size,
                                  do_merge);

  free(ents);
  return result;
}


/* \brief Copy and move all the entities
 */
int CopyMesh::copy_move_entities(iBase_EntityHandle *ent_handles,
                                 int num_ents,
                                 const double *dx,
                                 iBase_EntityHandle **new_ents,
                                 int *new_ents_alloc,
                                 int *new_ents_size,
                                 bool do_merge) 
{
  double zero[3] = {0,0,0};
  if(dx == NULL)
    dx = zero;
  return copy_transform_entities(ent_handles, num_ents,
                                 CopyMoveVerts(imeshImpl, dx),
                                 new_ents, new_ents_alloc, new_ents_size,
                                 do_merge);
}

/* \brief Copy and rotate all the entities
 */
int CopyMesh::copy_rotate_entities(iBase_EntitySetHandle set_handle,
                                   const double *origin,
                                   const double *z,
                                   double theta,
                                   iBase_EntityHandle **new_ents,
                                   int *new_ents_alloc,
                                   int *new_ents_size,
                                   bool do_merge)
{
  int err;
  iBase_EntityHandle *ents = NULL;
  int ents_alloc = 0, ents_size;
  iMesh_getEntitiesRec(imeshImpl, set_handle, 
                       iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, true,
                       &ents, &ents_alloc, &ents_size, &err);
  ERRORR("Failed to get entities from set recursively.", err);
  
  int result = copy_rotate_entities(ents, ents_size, origin, z, theta, 
                                    new_ents, new_ents_alloc, new_ents_size,
                                    do_merge);

  free(ents);
  return result;
}


/* \brief Copy and rotate all the entities
 */
int CopyMesh::copy_rotate_entities(iBase_EntityHandle *ent_handles,
                                   int num_ents,
                                   const double *origin,
                                   const double *z,
                                   double theta,
                                   iBase_EntityHandle **new_ents,
                                   int *new_ents_alloc,
                                   int *new_ents_size,
                                   bool do_merge) 
{
  return copy_transform_entities(ent_handles, num_ents,
                                 CopyRotateVerts(imeshImpl, origin, z, theta),
                                 new_ents, new_ents_alloc, new_ents_size,
                                 do_merge);
}


int CopyMesh::copy_transform_entities(iBase_EntityHandle *ent_handles,
                                      int num_ents,
                                      const CopyVerts &trans,
                                      iBase_EntityHandle **new_ents,
                                      int *new_ents_allocated,
                                      int *new_ents_size,
                                      bool do_merge)
{
  int err;

  iBase_EntitySetHandle set;
  iMesh_createEntSet(imeshImpl, false, &set, &err);
  ERRORR("Couldn't create source entity set.", err);
  
  iMesh_addEntArrToSet(imeshImpl, ent_handles, num_ents, set, &err);
  ERRORR("Couldn't add entities to source entity set.", err);

  int ret = copy_transform_entities(set, trans, new_ents, new_ents_allocated,
                                    new_ents_size, do_merge);

  iMesh_destroyEntSet(imeshImpl, set, &err);
  ERRORR("Couldn't destroy source entity set.", err);

  return ret;
}

int CopyMesh::copy_transform_entities(iBase_EntitySetHandle set_handle,
                                      const CopyVerts &trans,
                                      iBase_EntityHandle **new_ents,
                                      int *new_ents_allocated,
                                      int *new_ents_size,
                                      bool do_merge)
{
  int err;
  LocalTag local_tag(imeshImpl);

  iBase_EntityHandle *ents  = NULL; int ents_alloc  = 0, ents_size;
  iBase_EntityHandle *verts = NULL; int verts_alloc = 0, verts_size;
  int *indices              = NULL; int ind_alloc   = 0, ind_size;
  int *offsets              = NULL; int off_alloc   = 0, off_size;

  iMesh_getStructure(imeshImpl, set_handle,
                     &ents,    &ents_alloc,  &ents_size,
                     &verts,   &verts_alloc, &verts_size,
                     &indices, &ind_alloc,   &ind_size,
                     &offsets, &off_alloc,   &off_size,
                     &err);
  ERRORR("Trouble getting source adjacencies.", err);

  // copy the vertices
  iBase_EntityHandle *new_verts = NULL;
  int new_verts_alloc = 0, new_verts_size;
  trans(verts, verts_size, &new_verts, &new_verts_alloc, &new_verts_size);
  ERRORR("Couldn't create new vertices.", iBase_FAILURE);
  assert(new_verts_size == verts_size);

  // set the local copy tags on vertices
  // XXX: Should this really happen? Doing so adds more entities to copy sets
  // than explicitly passed into this function. This may be a domain-specific
  // question.
  iMesh_setEHArrData(imeshImpl, verts, verts_size, local_tag, 
                     new_verts, new_verts_size, &err);
  ERRORR("Error setting local copy tag data on old vertices.", iBase_FAILURE);

  // now connect the new vertices to make the higher-dimension entities
  err = connect_the_dots(ents, ents_size, local_tag, indices, offsets,
                         new_verts);
  ERRORR("Couldn't create new entities.", iBase_FAILURE);

  // take care of copy/expand sets
  if (!updatedCELists) {
    err = update_ce_lists();
    ERRORR(" ", err);
  }

  // set the target sets for expand sets to be themselves
  std::set<iBase_EntitySetHandle>::iterator set;
  for (set = expandSets.begin(); set != expandSets.end(); ++set) {
    iMesh_setEntSetEHData(imeshImpl, *set, local_tag,
                          reinterpret_cast<iBase_EntityHandle>(*set), &err);
  }

  err = process_ce_sets(imeshImpl, copySets, local_tag);
  ERRORR("Failed to update copy/expand sets.", iBase_FAILURE);
  
  err = process_ce_sets(imeshImpl, expandSets, local_tag);
  ERRORR("Failed to update expand/expand sets.", iBase_FAILURE);

  // set the copy tag on all copied sets
  for (set = copySets.begin(); set != copySets.end(); ++set) {
    iBase_EntityHandle eh;
    iMesh_getEntSetEHData(imeshImpl, *set, local_tag, &eh, &err);
    if (iBase_SUCCESS == err) {
      iMesh_setEntSetEHData(imeshImpl, *set, copyTag, eh, &err);
      ERRORR("Failed to tag copied set with copyTag.", iBase_FAILURE);
    }
  }

  std::vector<tag_data>::iterator tag;
  for(tag = copyTags.begin(); tag != copyTags.end(); ++tag) {
    err = tag_copy_sets(imeshImpl, copyTag, copySets, tag->tag, tag->value);
    ERRORR("Failed to tag copied sets.", iBase_FAILURE);
  }

  // get all the copies
  if (new_ents) {
    iMesh_getEHArrData(imeshImpl, ents, ents_size, local_tag, 
		       new_ents, new_ents_allocated, new_ents_size, &err);
    ERRORR("Failed to get copies from local tag.", iBase_FAILURE);
  }

  free(new_verts);
  free(ents);
  free(verts);
  free(indices);
  free(offsets);

  return iBase_SUCCESS;
}

int CopyMesh::connect_the_dots(iBase_EntityHandle *ents, int ents_size,
                               iBase_TagHandle local_tag,
                               int *indices, int *offsets,
                               iBase_EntityHandle *verts)
{
  int err;

  int *topos = NULL, topos_alloc = 0, topos_size;
  iMesh_getEntArrTopo(imeshImpl, ents, ents_size,
                      &topos, &topos_alloc, &topos_size, &err);
  ERRORR("Failed to get topos of all ents.", err);
  
  // scan forward to first non-vertex
  int pos = 0;
  while (iMesh_POINT == topos[pos] && pos < topos_size) 
    pos++;
  if (pos == topos_size) return iBase_SUCCESS;

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
      ERRORR("Couldn't create new entities.", iBase_FAILURE);
    }
    else {
      // use single-entity function in this case, entity might have higher-order
      // nodes (imesh fcn doesn't have argument for # entities)
      for (int i = 0; i < num_ents; i++) {
        int status;
        iMesh_createEnt(imeshImpl, topo, &connect[i*vtx_per_ent],
                        vtx_per_ent, &new_ents[i], &status, &err);
        ERRORR("Couldn't create new entities.", iBase_FAILURE);
      }
    }

    // set the local copy tags
    iMesh_setEHArrData(imeshImpl, &ents[begin], num_ents, local_tag, 
                       &new_ents[0], num_ents, &err);
    ERRORR("Error setting local copy tag data on old ents.", iBase_FAILURE);
  }

  free(topos);
  return iBase_SUCCESS;
}

int CopyMesh::update_ce_lists() 
{
  if (updatedCELists) reset_ce_lists();

  int err;
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(imeshImpl, &root_set, &err);
  ERRORR("Trouble getting root set", iBase_FAILURE);
  
  err = update_tagged_sets(root_set, copyTags, copySets);
  ERRORR("Trouble updating tagged copy sets", iBase_FAILURE);

  err = update_tagged_sets(root_set, expandTags, expandSets);
  ERRORR("Trouble updating tagged expand sets", iBase_FAILURE);

  updatedCELists = true;
  
  return iBase_SUCCESS;
}

int CopyMesh::update_tagged_sets(iBase_EntitySetHandle from_set,
                                 const std::vector<CopyMesh::tag_data> &tags,
                                 std::set<iBase_EntitySetHandle> &tagged_sets)
{
  int err;
  std::vector<tag_data>::const_iterator tag;
  for (tag = tags.begin(); tag != tags.end(); ++tag) {
    iBase_EntitySetHandle *tmp_sets = NULL;
    int tmp_alloc = 0, tmp_size;
    iMesh_getEntSetsByTagsRec(imeshImpl, from_set, &tag->tag,
                              (tag->value ? &tag->value : NULL), 1, false,
                              &tmp_sets, &tmp_alloc, &tmp_size, &err);
    ERRORR("Couldn't get tagged sets.", iBase_FAILURE);
    tagged_sets.insert(tmp_sets, tmp_sets+tmp_size);    
    free(tmp_sets);
  }

  return iBase_SUCCESS;
}

int CopyMesh::tag_copied_sets(const char **tag_names, const char **tag_vals,
                              const int num_tags) 
{
  int err;
  
  for (int t = 0; t < num_tags; t++) {
    iBase_TagHandle tag;
    iMesh_getTagHandle(imeshImpl, tag_names[t], &tag, &err,
                       strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);

    err = tag_copy_sets(imeshImpl, copyTag, copySets, tag,
                        tag_vals ? tag_vals[t] : NULL);
    ERRORR("Failed to tag copied set.", err);
  }

  return iBase_SUCCESS;
}

int CopyMesh::tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
                              const int num_tags) 
{
  int err;
  
  for (int t = 0; t < num_tags; t++) {
    err = tag_copy_sets(imeshImpl, copyTag, copySets, tags[t],
                        tag_vals ? tag_vals[t] : NULL);
    ERRORR("Failed to tag copied set.", err);
  }
  
  return iBase_SUCCESS;
}

#ifdef TEST

#include <cfloat>
#include "MergeMesh.hpp"

iMesh_Instance impl;
iBase_EntitySetHandle root_set;

int check_num_ents(int ent_type, int expected) 
{
  int num_type, err;
  iMesh_getNumOfType(impl, root_set, ent_type, &num_type, &err);
  ERRORR("Failed to get # entities", 1);
  if (num_type != expected) {
    std::cerr << "Didn't get right # entities of dimension " << ent_type
              << "; got " << num_type
              << ", expected " << expected << std::endl;
    return 1;
  }
  return 0;
}

int parse_input(int argc, char **argv, 
                std::string &infile, std::string &outfile, 
                double *dx, 
                std::vector<std::string> &ctags, 
                std::vector<std::string> &etags, 
                std::vector<std::string> &utags,
                std::vector<char*> &cvals, 
                std::vector<char*> &evals);

int make_ents(CopyMesh *cm, iBase_EntityHandle *&ents, 
              int &ents_alloc, int &ents_size, 
              iBase_EntityHandle *&verts, 
              int &verts_alloc, int &verts_size);

int main(int argc, char **argv) 
{
  // make a mesh instance
  int err;
  iMesh_newMesh("MOAB", &impl, &err, 4);
  ERRORR("Failed to create instance.", 1);
  
  iMesh_getRootSet(impl, &root_set, &err);

  CopyMesh *cm = new CopyMesh(impl);
  MergeMesh *mm = new MergeMesh(impl);

  // some entity tag types are always copy or expand
  cm->add_expand_tag("MATERIAL_SET");
  cm->add_expand_tag("DIRICHLET_SET");
  cm->add_expand_tag("NEUMANN_SET");
  
  if (1 == argc) {
    int tag_val = 3;
    cm->add_copy_tag(std::string("GEOM_DIMENSION"), (const char*) &tag_val);

    // create vertices
    iBase_EntityHandle *ents = NULL, *verts = NULL;
    int ents_alloc = 0, ents_size, verts_alloc = 0, verts_size, 
      new_ents_alloc, new_ents_size;

    err = make_ents(cm, ents, ents_alloc, ents_size, 
                    verts, verts_alloc, verts_size);
    ERRORR("Making entities failed.", err);
    
    double dx[] = {1.0,0.0,0.0};
    iBase_EntityHandle *new_ents = ents + ents_size;
    new_ents_alloc = ents_size;
    
  
    err = cm->copy_move_entities(ents, ents_size, dx,
                                 &new_ents, &new_ents_alloc, &new_ents_size, 
                                 false);
    ERRORR("Failed to copy_move entities.", 1);

    // check # entities
    if (check_num_ents(iBase_VERTEX, 32)) return 1;
    if (check_num_ents(iBase_REGION, 6)) return 1;
  
    err = mm->merge_entities(ents, ents_size+new_ents_size, 1.0e-8);
    ERRORR("Failed to merge entities.", 1);

    // now get all vertices, put in new verts array
    iBase_EntityHandle *nverts = NULL;
    int nverts_alloc = 0, nverts_size;
    iMesh_getEntities(impl, root_set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES, 
                      &nverts, &nverts_alloc, &nverts_size, &err);
    ERRORR("Didn't get all vertices.", 1);

    // and their coords
    double *vert_coords = NULL;
    int vert_coords_alloc = 0, vert_coords_size;
    int sorder = iBase_INTERLEAVED;
    iMesh_getVtxArrCoords (impl, nverts, nverts_size, sorder, 
                           &vert_coords, &vert_coords_alloc, &vert_coords_size,
                           &err);
    ERRORR("Didn't get vtx coords.", 1);

    double xmin[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double xmax[3] = {DBL_MIN, DBL_MIN, DBL_MIN};

    for (int i = 0; i < nverts_size; i++) {
      for (int j = 0; j < 3; j++) {
        if (vert_coords[3*i+j] < xmin[j]) xmin[j] = vert_coords[3*i+j];
        if (vert_coords[3*i+j] > xmax[j]) xmax[j] = vert_coords[3*i+j];
      }
    }
    if (xmin[0] != 0.0 || xmin[1] != 0.0 || xmin[2] != 0.0 ||
        xmax[0] != 2.0 || xmax[1] != 1.0 || xmax[2] != 3.0) {
      std::cerr << "Didn't get correct min/max; output values = ("
                << xmin[0] << ", " << xmin[1] << ", " << xmin[2] << "), ("
                << xmax[0] << ", " << xmax[1] << ", " << xmax[2] << ")" << std::endl;
      return 1;
    }

    int expected = 2*verts_size - 8;
    if (nverts_size != expected) {
      std::cerr << "Didn't get right # vertices; got " << nverts_size
                << ", expected " << expected << std::endl;
    
      return 1;
    }
  
    free(verts);
    free(nverts);
    free(ents);
  }
  
  else {
    std::vector<std::string> ctags, etags, utags;
    std::vector<char*> cvals, evals;
    std::string infile, outfile;
    double dx[3] = {0.0, 0.0, 0.0};

    err = parse_input(argc, argv, 
                      infile, outfile, dx, 
                      ctags, etags, utags,
                      cvals, evals);
    if (-1 == err) return err;
    ERRORR("Couldn't parse input, exiting.", err);
    
    // read the file
    iMesh_load(impl, 0, infile.c_str(), NULL, &err, infile.length(), 0);
    ERRORR("Couldn't read mesh file.", err);
    
    // set copy/expand/unique tags
    for (unsigned int i = 0; i < ctags.size(); i++)
      cm->add_copy_tag(ctags[i], cvals[i]);
    for (unsigned int i = 0; i < etags.size(); i++)
      cm->add_expand_tag(etags[i], evals[i]);
    for (unsigned int i = 0; i < utags.size(); i++)
      cm->add_unique_tag(utags[i]);
      
    // copy
    iBase_EntityHandle* new_ents = NULL;
    int new_ents_alloc = 0;
    int new_ents_size = 0;
    err = cm->copy_move_entities(root_set, dx, &new_ents,
				 &new_ents_alloc, &new_ents_size,
                                 false);
    ERRORR("Failed to copy_move entities.", 1);
    //getting elements for merge mesh   
    iBase_EntityHandle *ents = NULL;
    int ents_alloc = 0, ents_size;
    iMesh_getEntitiesRec(impl, root_set, 
			 iBase_ALL_TYPES, iMesh_HEXAHEDRON, true,
			 &ents, &ents_alloc, &ents_size, &err);
    ERRORR("Failed to get entities from set recursively.", err);
    // merge 
    const double merge_tol =  1.0e-8;
    const int do_merge = 0;
    const int update_sets= 0;
    iBase_TagHandle merge_tag = NULL;
    err = mm->merge_entities(ents, ents_size, merge_tol,
			     do_merge, update_sets, merge_tag);
    ERRORR("Failed to merge entities.", 1);

    // export
    iMesh_save(impl, root_set, outfile.c_str(), NULL, &err, 
               outfile.length(), 0);
    ERRORR("Failed to save mesh.", 1);
    std::cout << "Wrote " << outfile << std::endl;
  }

  delete cm;
  delete mm;
  
  iMesh_dtor(impl, &err);
  ERRORR("Destructor failed.", 1);
  
  return 0;
}

int parse_input(int argc, char **argv, 
                std::string &infile, std::string &outfile, 
                double *dx, 
                std::vector<std::string> &ctags, 
                std::vector<std::string> &etags, 
                std::vector<std::string> &utags,
                std::vector<char*> &cvals, 
                std::vector<char*> &evals) 
{
  if (!strcmp(argv[1], "-h")) {
    std::cout << "Usage: " << argv[0] << " _OR_ " << std::endl
              << argv[0] << " [-c copy_tag [-cv copy_val]]* "
              << "[-e expand_tag [-ev copy_val]]* "
              << "[-u unique_tag]* " 
              << "[-x <dx> <dy> <dz>] <infile> <outfile>" << std::endl;
    std::cout << "* - repeated zero or more times" << std::endl;
    return -1;
  }
  
  int pos = 1;
  for (; pos < argc; pos++)
    {
      if (!strcmp(argv[pos],"-c")) {
	pos++;
	ctags.push_back(argv[pos]);
	if (!strcmp(argv[pos],"-cv")) {
	  pos++;
	  cvals.push_back(argv[pos]);
	}
	else cvals.push_back(NULL);
      }
      else if (!strcmp(argv[pos],"-e")) {
	pos++;
	etags.push_back(argv[pos]);
	if (!strcmp(argv[pos],"-ev")) {
	  pos++;
	  evals.push_back(argv[pos]);
	}
	else evals.push_back(NULL);
      }
      else if (!strcmp(argv[pos],"-u")) {
	pos++;
	utags.push_back(argv[pos]);
      }
      else if (!strcmp(argv[pos],"-x")) {
	pos++;
	dx[0] = atof(argv[pos++]);
	dx[1] = atof(argv[pos++]);
	dx[2] = atof(argv[pos]);
      }
      else if (argv[pos][0] != '-') {
	break;
      }
      else {
	std::cerr << "Invalid option: \"" << argv[pos] << '"' << std::endl;
      }
    }

  if (pos == argc) {
    std::cerr << "Not enough arguments, infile missing." << std::endl;
    return 1;
  }
  else
    infile = argv[pos++];

  if (pos == argc) {
    std::cerr << "Not enough arguments, outfile missing." << std::endl;
    return 1;
  }
  
  outfile = argv[pos++];
  
  return 0;
}

int make_ents(CopyMesh *cm, iBase_EntityHandle *&ents, 
              int &ents_alloc, int &ents_size, 
              iBase_EntityHandle *&verts, 
              int &verts_alloc, int &verts_size) 
{
  // make an arbitrary mesh to copy
  double vert_pos[] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
    0.0, 0.0, 2.0,
    1.0, 0.0, 2.0,
    1.0, 1.0, 2.0,
    0.0, 1.0, 2.0,
    0.0, 0.0, 3.0,
    1.0, 0.0, 3.0,
    1.0, 1.0, 3.0,
    0.0, 1.0, 3.0
  };
	 
  int connect_idx[] = {
    0, 1, 2, 3, 4, 5, 6, 7,
    4, 5, 6, 7, 8, 9, 10, 11,
    8, 9, 10, 11, 12, 13, 14, 15
  };

  int adj_table[] = {
    1, 0, 0, 1,
    1, 0, 0, 0,
    1, 0, 1, 0, 
    1, 0, 0, 1
  };
  int err;

  // set the adjacency table to request faces
  iMesh_setAdjTable(cm->impl(), adj_table, 16, &err);
  ERRORR("Failed to set adj table.", err);
	
  // create vertices
  verts = NULL; verts_alloc = 0;
  iMesh_createVtxArr(cm->impl(), 16, iBase_INTERLEAVED, vert_pos, 48,
                     &verts, &verts_alloc, &verts_size, &err);
  ERRORR("Failed to create 16 vertices.", 1);
	
  // pre-allocate ents space
  ents = (iBase_EntityHandle*) malloc(6*sizeof(iBase_EntityHandle));
  ents_alloc = 6;
	 
  // create hexes
  iBase_EntityHandle connect[sizeof(connect_idx)/sizeof(int)];
  int num_connect = sizeof(connect_idx)/sizeof(int);
  for (int i = 0; i < num_connect; i++)
    connect[i] = verts[connect_idx[i]];
	 
  int *status = NULL, status_alloc = 0, status_size;
  iMesh_createEntArr(impl,  iMesh_HEXAHEDRON, connect, 24,
                     &ents, &ents_alloc, &ents_size,
                     &status, &status_alloc, &status_size, &err);
  ERRORR("Failed to create 3 hexes.", 1);

  // create entity sets and add things to it
  iBase_EntitySetHandle esets[7];
  for (int i = 0; i < 7; i++) {
    iMesh_createEntSet(cm->impl(), false, esets+i, &err);
    ERRORR("Failed to create entity set.", err);
  }
  
  iBase_TagHandle gtag, mtag, ntag, dtag;
  iMesh_getTagHandle(cm->impl(), "GEOM_DIMENSION", &gtag, &err, 14);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "GEOM_DIMENSION", 1, iBase_INTEGER,
                    &gtag, &err, 14);
    ERRORR("Failed to create tag handle for geom dimension.", err);
  }
  
  iMesh_getTagHandle(cm->impl(), "MATERIAL_SET", &mtag, &err, 12);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "MATERIAL_SET", 1, iBase_INTEGER,
                    &mtag, &err, 12);
    ERRORR("Failed to create tag handle for material set.", err);
  }
  
  iMesh_getTagHandle(cm->impl(), "NEUMANN_SET", &ntag, &err, 11);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "NEUMANN_SET", 1, iBase_INTEGER,
                    &ntag, &err, 11);
    ERRORR("Failed to create tag handle for neumann set.", err);
  }
  
  iMesh_getTagHandle(cm->impl(), "DIRICHLET_SET", &dtag, &err, 13);
  if (iBase_TAG_NOT_FOUND == err) {
    iMesh_createTag(cm->impl(), "DIRICHLET_SET", 1, iBase_INTEGER,
                    &dtag, &err, 13);
    ERRORR("Failed to create tag handle for dirichlet set.", err);
  }
  
  // ES1: geom dim = 3, hexes 1-3, C
  iMesh_addEntArrToSet(cm->impl(), ents, 3, *esets, &err);
  ERRORR("Failed to add ents to set 1.", err);
  int dum = 3;
  iMesh_setEntSetIntData(cm->impl(),esets[0], gtag, dum, &err);
  ERRORR("Failed to set geom set tag on set 1.", err);
  
  // ES2: geom dim = 2, faces all, -
  iBase_EntityHandle *faces = NULL;
  int faces_alloc = 0, faces_size;
  iMesh_getEntities(cm->impl(), root_set, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                    &faces, &faces_alloc, &faces_size, &err);
  ERRORR("Failed to get faces.", err);
  iMesh_addEntArrToSet(cm->impl(), faces, faces_size, *(esets+1), &err);
  ERRORR("Failed to add ents to set 2.", err);
  dum = 2;
  iMesh_setEntSetIntData(cm->impl(), esets[1], gtag, dum, &err);
  ERRORR("Failed to set geom set tag on set 2.", err);
  free(faces);

  // ES3: - , skin nodes, E
  iMesh_addEntArrToSet(cm->impl(), verts, verts_size, *(esets+2), &err);
  ERRORR("Failed to add ents to set 3.", err);
  err = cm->add_copy_expand_list(esets+2, 1, CopyMesh::EXPAND);
  ERRORR("Failed to add set 3 to expand list.", err);

  // ES4: mat set, set 1, E
  iMesh_addEntSet(cm->impl(), esets[0],*(esets+3), &err);
  ERRORR("Failed to add set to set 4.", err);
  dum = 100;
  iMesh_setEntSetIntData(cm->impl(), esets[3], mtag, dum, &err);
  ERRORR("Failed to set mat set tag on set 4.", err);

  // ES5: neu set, set 2, E
  iMesh_addEntSet(cm->impl(), esets[1],*(esets+4), &err);
  ERRORR("Failed to add set to set 5.", err);
  dum = 101;
  iMesh_setEntSetIntData(cm->impl(), esets[4], ntag, dum, &err);
  ERRORR("Failed to set neu set tag on set 5.", err);

  // ES6: dir set, set 3, E
  iMesh_addEntSet(cm->impl(), esets[2],*(esets+5), &err);
  ERRORR("Failed to add set to set 6.", err);
  dum = 102;
  iMesh_setEntSetIntData(cm->impl(), esets[5], dtag, dum, &err);
  ERRORR("Failed to set dir set tag on set 6.", err);

  // ES7: -, set 4, C
  iMesh_addEntSet(cm->impl(), esets[3],*(esets+6), &err);
  ERRORR("Failed to add set to set 7.", err);
  err = cm->add_copy_expand_list(esets+6, 1, CopyMesh::COPY);
  ERRORR("Failed to add set 7 to copy list.", err);
	
  // check # entities
  if (check_num_ents(iBase_VERTEX, 16)) return 1;
  if (check_num_ents(iBase_REGION, 3)) return 1;

  return 0;
}

#endif
