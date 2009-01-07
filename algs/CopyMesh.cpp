#include "CopyMesh.hpp"
#include "MBCN.hpp"
#include <stdlib.h>
#include <algorithm>
#include <functional>

/*
 * - call update_ce_lists to update copySets or expandSets
 */

#define CHECK_SIZE(array, allocated, size, type, retval)  \
  if (0 != allocated && NULL != array && allocated < (size)) {\
    ERRORR("Allocated array not large enough to hold returned contents.", \
           iBase_MEMORY_ALLOCATION_FAILED);                             \
  }\
  if (allocated == 0 || NULL == array) {\
    array = (type*)malloc((size)*sizeof(type));\
    allocated=(size);\
    if (NULL == array) {ERRORR("Couldn't allocate array.",\
                               iBase_MEMORY_ALLOCATION_FAILED); }     \
  }

CopyMesh::CopyMesh(iMesh_Instance impl) 
        : imeshImpl(impl), updatedCELists(false)
{
    // allocate copy tag
  copyTag = 0;
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

  if (0 == copyTag && copy_or_expand == COPY) {
    iMesh_createTag(imeshImpl, "__CopyMeshTag", 1,
                    iBase_ENTITY_HANDLE, &copyTag, &err, 13);
    ERROR("Couldn't create copy mesh tag.");
  }
  
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
  int err;

    // create a tag for this local copy operation
  iBase_TagHandle local_tag;
  iMesh_createTag(imeshImpl, "local_copy", 1, iBase_ENTITY_HANDLE,
                  &local_tag, &err, 13);
  ERRORR("Failed to create local copy tag.", err);

    // create a set to hold entities to be copied
  iBase_EntitySetHandle copy_set;
  iMesh_createEntSet(imeshImpl, false, &copy_set, &err);
  ERRORR("Failed to create set.", err);

    // add entities & adj vertices
  iMesh_addEntArrToSet(imeshImpl, ent_handles, num_ents,
                       &copy_set, &err);
  ERRORR("Failed to add ents to set", err);
  iBase_EntityHandle *verts = NULL;
  int verts_alloc = 0, verts_size;
  int *offset = NULL, offset_alloc = 0, offset_size;
  int *in_set = NULL, in_set_alloc = 0, in_set_size;
  iMesh_getAdjEntities(imeshImpl, copy_set, 
                       iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, iBase_VERTEX,
                       &verts, &verts_alloc, &verts_size,
                       &offset, &offset_alloc, &offset_size,
                       &in_set, &in_set_alloc, &in_set_size,
                       &err);
  ERRORR("Failed to get adj entities in copy set", err);
  iMesh_addEntArrToSet(imeshImpl, verts, verts_size,
                       &copy_set, &err);
  ERRORR("Failed to add verts to set", err);

  err = copy_move_verts(copy_set, dx, local_tag);
  ERRORR("Failed to copy/move vertices.", iBase_FAILURE);

    // now copy entities
  err = copy_move_ents(copy_set, local_tag);
  ERRORR("Failed to copy/move entities.", iBase_FAILURE);

    // take care of copy/expand sets
  if (!updatedCELists) {
    err = update_ce_lists();
    ERRORR(" ", err);
  }
  err = process_ce_sets(copySets, COPY, local_tag);
  ERRORR("Failed to update copy/expand sets.", iBase_FAILURE);
  
  err = process_ce_sets(expandSets, EXPAND, local_tag);
  ERRORR("Failed to update expand/expand sets.", iBase_FAILURE);

  if (!copyTags.empty()) {
    err = tag_copied_sets(&copyTags[0], &copyTagVals[0], copyTags.size());
    ERRORR("Failed to tag copied sets.", iBase_FAILURE);
  }
    
  free(verts);
  free(offset);
  free(in_set);

    // destroy local tag, removing it from all entities
  iMesh_destroyTag(imeshImpl, local_tag, true, &err);
  ERRORR("Failed to force-destroy local copy tag.", iBase_FAILURE);
  
  return iBase_SUCCESS;
}

int CopyMesh::process_ce_sets(std::set<iBase_EntitySetHandle> &cesets,
                              int copy_or_expand,
                              iBase_TagHandle local_tag) 
{
  int err;
  std::vector<iBase_EntitySetHandle> copy_sets(cesets.size(), (iBase_EntitySetHandle)0);
  
    // for each ceset:
  unsigned int i;
  std::set<iBase_EntitySetHandle>::iterator sit;
  for (i = 0, sit = cesets.begin(); sit != cesets.end(); i++, sit++) {
      // - handle entities
    iBase_EntityHandle *tmp_ents = NULL;
    int tmp_ents_alloc = 0, tmp_ents_size;
    iMesh_getEntities(imeshImpl, *sit, 
                      iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                      &tmp_ents, &tmp_ents_alloc, &tmp_ents_size, &err);
    ERRORR("Failed to get ceSet entities.", iBase_FAILURE);

      // get copy tags and remove null ones
    std::vector<iBase_EntityHandle> tmp_tags(tmp_ents_size, (iBase_EntitySetHandle)0);
    iBase_EntityHandle *tmp_tags_ptr = &tmp_tags[0];
    int tmp_tags_alloc = tmp_ents_size, tmp_tags_size;
    iMesh_getEHArrData(imeshImpl, tmp_ents, tmp_ents_size, local_tag, 
                       &tmp_tags_ptr, &tmp_tags_alloc, &tmp_tags_size, &err);
    std::vector<iBase_EntityHandle>::iterator new_end, vit;
    new_end = std::remove_if(tmp_tags.begin(), tmp_tags.end(), 
                             std::bind2nd(std::equal_to<iBase_EntityHandle>(),
                                          (iBase_EntitySetHandle)0));

    // - if non-zero copied ents & copy set, make a new copy set, add copied ents
    if (new_end != tmp_tags.begin() && copy_or_expand == COPY) {
      iMesh_createEntSet(imeshImpl, false, &copy_sets[i], &err);
      ERRORR("Failed to create copied set.", iBase_FAILURE);
      iMesh_setEntSetEHData(imeshImpl, *sit, local_tag, copy_sets[i], &err);
      ERRORR("Failed to tag copied set.", iBase_FAILURE);
    }

    iBase_EntitySetHandle added_to_set = 
        (copy_or_expand == COPY ? copy_sets[i] : *sit);

    if (new_end != tmp_tags.begin()) {
      iMesh_addEntArrToSet(imeshImpl, &tmp_tags[0], new_end-tmp_tags.begin(),
                           &added_to_set, &err);
      ERRORR("Failed to add copied entities to ce set.", iBase_FAILURE);
    }
    free(tmp_ents);

      // - handle sets
    iBase_EntitySetHandle *tmp_sets = NULL;
    int tmp_sets_alloc = 0, tmp_sets_size;
    iMesh_getEntSets(imeshImpl, *sit, 1,
                     &tmp_sets, &tmp_sets_alloc, &tmp_sets_size, &err);
    ERRORR("Failed to get ceSet entities.", iBase_FAILURE);
    
      // get copy tags and save non-null ones
    tmp_tags.clear();
    iBase_EntityHandle eh_tag;
    for (int j = 0; j < tmp_sets_size; j++) {
      iMesh_getEntSetEHData(imeshImpl, tmp_sets[j], local_tag, 
                            &eh_tag, &err);
      if (iBase_SUCCESS == err && eh_tag) tmp_tags.push_back(eh_tag);
    }

      // might need to create set here
    if (!copy_sets[i] && !tmp_tags.empty() && copy_or_expand == COPY) {
      iMesh_createEntSet(imeshImpl, false, &copy_sets[i], &err);
      ERRORR("Failed to create copied set.", iBase_FAILURE);
      iMesh_setEntSetEHData(imeshImpl, *sit, local_tag, copy_sets[i], &err);
      ERRORR("Failed to tag copied set.", iBase_FAILURE);
      added_to_set = copy_sets[i];
    }

    for (vit = tmp_tags.begin(); vit != tmp_tags.end(); vit++) {
      iMesh_addEntSet(imeshImpl, *vit, &added_to_set, &err);
      ERRORR("Failed to add set to ce set.", iBase_FAILURE);
    }

    if (copy_or_expand == COPY && copyTag && copy_sets[i]) {
      iMesh_setEntSetEHData(imeshImpl, *sit, copyTag, copy_sets[i], &err);
      ERRORR("Failed to tag copied set with copyTag.", iBase_FAILURE);
    }
  }
  
    // xxx - still need to do unique tags
    // - for each unique tag
    //   . if this set has the tag, add/append to get unique value

  return iBase_SUCCESS;
}

int CopyMesh::copy_move_ents(iBase_EntitySetHandle copy_set,
                             iBase_TagHandle local_tag) 
{
  int err;
  iBase_EntityHandle *ents = NULL;
  int ents_alloc = 0, ents_size;
  iMesh_getEntities(imeshImpl, copy_set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                    &ents, &ents_alloc, &ents_size, &err);
  ERRORR("Failed to get all ents in copy set.", err);
  int *topos = NULL, topos_alloc = 0, topos_size;
  iMesh_getEntArrTopo(imeshImpl, ents, ents_size,
                      &topos, &topos_alloc, &topos_size, &err);
  ERRORR("Failed to get topos of all ents.", err);
  
    // scan forward to first non-vertex
  int pos = 0;
  while (iMesh_POINT == topos[pos] && pos < topos_size) 
    pos++;
  if (pos == topos_size) return iBase_SUCCESS;

    // get connectivity
  iBase_EntityHandle *verts = NULL;
  int verts_alloc = 0, verts_size;
  int *offset = NULL, offset_alloc = 0, offset_size;
  ents_size -= pos;
  iMesh_getEntArrAdj(imeshImpl, &ents[pos], ents_size, iBase_VERTEX,
                     &verts, &verts_alloc, &verts_size,
                     &offset, &offset_alloc, &offset_size, &err);
  ERRORR("Failed to get adj verts of all ents.", err);

    // for each run of same size & type
  iBase_EntityHandle *tmp_ents = &ents[pos];
  int *tmp_topos = &topos[pos];
  int start_idx, end_idx = 0;
  
  std::vector<iBase_EntityHandle> connect, new_ents;
  iBase_EntityHandle *connect_ptr, *new_ents_ptr;
  std::vector<int> status;
  int status_alloc, status_size, *status_ptr;
  int num_corner_verts;
  while (end_idx < ents_size) {
    
      // get next run; end_idx points to start of *next* element,
      // or offset_size if no elems left
    start_idx = end_idx++;
    int topo_start = tmp_topos[start_idx];
    int vtx_per_ent = offset[end_idx] - offset[start_idx];
    int mbcn_type;
    iMesh_MBCNType(topo_start, &mbcn_type);
    MBCN_VerticesPerEntity(mbcn_type, &num_corner_verts);
    while (end_idx < ents_size &&
           tmp_topos[end_idx] == topo_start &&
           offset[end_idx+1] - offset[end_idx] == vtx_per_ent)
      end_idx++;

      // build vector of vtx handles
    int num_ents = end_idx - start_idx;
    int totv = vtx_per_ent * num_ents;
    connect.resize(totv);
    connect_ptr = &connect[0];
    int tmp_alloc = totv, tmp_size = totv; 
    iMesh_getEHArrData(imeshImpl,
                       &verts[offset[start_idx]], totv,
                       local_tag, &connect_ptr, &tmp_alloc, &tmp_size, &err);
    ERRORR("Couldn't get connectivity for copied ents.", iBase_FAILURE);
    
      // create entities
    status.resize(num_ents); status_ptr = &status[0]; status_alloc = status.size();
    new_ents.resize(num_ents); new_ents_ptr = &new_ents[0]; tmp_alloc = num_ents;
    if (1 < num_ents && num_corner_verts == vtx_per_ent) {
      iMesh_createEntArr(imeshImpl, topo_start, 
                         &connect[0], connect.size(),
                         &new_ents_ptr, &tmp_alloc, &tmp_size,
                         &status_ptr, &status_alloc, &status_size, &err);
      ERRORR("Couldn't create new entities.", iBase_FAILURE);
    }
    else {
        // use single-entity function in this case, entity might have higher-order nodes
        // (imesh fcn doesn't have argument for # entities)
      tmp_size = 0;
      for (int i = 0; i < num_ents; i++) {
        iMesh_createEnt(imeshImpl, topo_start, 
                        &connect[i*vtx_per_ent], vtx_per_ent, 
                        new_ents_ptr+i, status_ptr, &err);
        ERRORR("Couldn't create new entities.", iBase_FAILURE);
        tmp_size++;
      }
    }

      // set the local copy tags
    iMesh_setEHArrData(imeshImpl, &tmp_ents[start_idx], num_ents, local_tag, 
                       new_ents_ptr,  tmp_size, &err);
    ERRORR("Error setting local copy tag data on old ents.", iBase_FAILURE);
  }

  return iBase_SUCCESS;
}

int CopyMesh::copy_move_verts(iBase_EntitySetHandle copy_set, 
                              const double *dx,
                              iBase_TagHandle local_tag) 
{
    // get the vertices in the order they exist in the set, w/o duplicates
  iBase_EntityHandle *verts = NULL;
  int verts_alloc = 0, verts_size;
  int err;
  iMesh_getEntities(imeshImpl, copy_set, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                    &verts, &verts_alloc, &verts_size, &err);
  ERRORR("Failed to get verts in copy set.", err);
  
  std::vector<iBase_EntityHandle> new_verts(verts_size);
  iBase_EntityHandle *new_verts_ptr = &new_verts[0];

    // get position of vertices, transform with move vector
  double *coords = NULL;
  int coords_alloc = 0, coords_size;
  int order = iBase_UNDETERMINED;
  iMesh_getVtxArrCoords(imeshImpl, verts, verts_size,
                        &order, &coords, &coords_alloc, &coords_size, &err);
  ERRORR("Failed to get vtx coords.", iBase_FAILURE);
  if (NULL != dx) {
    if (iBase_INTERLEAVED == order) {
      for (int i = 0; i < 3*verts_size; i+=3)
        coords[i] += dx[0], coords[i+1] += dx[1], coords[i+2] += dx[2];
    }
    else if (iBase_BLOCKED == order) {
      double *coordsy = coords + verts_size;
      double *coordsz = coordsy + verts_size;
      for (int i = 0; i < verts_size; i++)
        coords[i] += dx[0], coordsy[i] += dx[1], coordsz[i] += dx[2];
    }
    else {
      err = iBase_FAILURE;
      ERRORR("Couldn't determine ordering for vertex coordinates.", err);
    }
  }
  
    // copy vertices
  int tmp_size, tmp_alloc = verts_size;
  iMesh_createVtxArr(imeshImpl, verts_size, order, coords, coords_size,
                     &new_verts_ptr, &tmp_alloc, &tmp_size, &err);
  ERRORR("Couldn't create new vertices.", iBase_FAILURE);
  assert(tmp_size == verts_size);

    // set the local copy tags
  iMesh_setEHArrData(imeshImpl, verts, verts_size, local_tag, 
                     &new_verts[0],  verts_size, &err);
  ERRORR("Error setting local copy tag data on old verts.", iBase_FAILURE);

#ifndef NDEBUG
  iMesh_getEHArrData(imeshImpl, verts, verts_size, local_tag, 
                     &new_verts_ptr,  &tmp_alloc, &tmp_size, &err);
  ERRORR("Error getting local copy tag data on old verts.", iBase_FAILURE);
#endif

  free(coords);
  free(verts);
  
  return iBase_SUCCESS;
}

int CopyMesh::update_ce_lists() 
{
  if (updatedCELists) reset_ce_lists();

  int err;
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(imeshImpl, &root_set, &err);
  ERRORR("Trouble getting root set", iBase_FAILURE);
  
  if (!copyTags.empty()) {
    err = update_tagged_sets(root_set, &copyTags[0], &copyTagVals[0],
                             copyTags.size(), copySets);
    ERRORR("Trouble updating tagged copy sets", iBase_FAILURE);
  }
  
  if (!expandTags.empty()) {
    err = update_tagged_sets(root_set, &expandTags[0], &expandTagVals[0],
                             expandTags.size(), expandSets);
    ERRORR("Trouble updating tagged expand sets", iBase_FAILURE);
  }

  updatedCELists = true;
  
  return iBase_SUCCESS;
}
  
int CopyMesh::update_tagged_sets(iBase_EntitySetHandle from_set,
                                 iBase_TagHandle *tag_handles,
                                 const char **tag_vals,
                                 int num_tags,
                                 std::set<iBase_EntitySetHandle> &tagged_sets)
{
    // if xxx
  iBase_EntitySetHandle *tmp_sets;
  int tmp_alloc, tmp_size;
  int err;
  for (int i = 0; i < num_tags; i++) {
    tmp_sets = NULL;
    tmp_alloc = 0;
    iMesh_getEntSetsByTagsRec(imeshImpl, from_set, &tag_handles[i], 
                              (tag_vals[i] ? &tag_vals[i] : NULL),
                              1, 0,
                              &tmp_sets, &tmp_alloc, &tmp_size, &err);
    ERRORR("Couldn't get tagged sets.", iBase_FAILURE);
    for (int j = 0; j < tmp_size; j++) 
      tagged_sets.insert(tmp_sets[j]);
    
    free(tmp_sets);
  }

  return iBase_SUCCESS;
}
    
int CopyMesh::tag_copied_sets(const char **tag_names, const char **tag_vals,
                              const int num_tags) 
{
  int err;
  std::vector<iBase_TagHandle> tag_handles;
  
  for (int t = 0; t < num_tags; t++) {
      // get tag handle & size
    iBase_TagHandle th;
    iMesh_getTagHandle(imeshImpl, tag_names[t], &th, &err, strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);
    tag_handles.push_back(th);
  }

  return tag_copied_sets(&tag_handles[0], tag_vals, num_tags);
}

int CopyMesh::tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
                              const int num_tags) 
{
  int err = iBase_SUCCESS;
  std::vector<char> tag_space;
  char *tag_space_ptr;
  int tag_space_alloc, tag_space_size;
  assert(copyTag);
  
  for (int t = 0; t < num_tags; t++) {
      // get tag handle & size
    int tag_size;
    iMesh_getTagSizeBytes(imeshImpl, tags[t], &tag_size, &err);
    ERRORR("Failed to get tag size.", err);

      // allocate tmp space for tag value
    tag_space.resize(tag_size);
    tag_space_ptr = &tag_space[0];
    tag_space_alloc = tag_size;
    
      // for each orig copy set with this tag, copy it to its copy
    for (std::set<iBase_EntitySetHandle>::iterator sit = copySets.begin();
         sit != copySets.end(); sit++) {
        // get the tag value
      iMesh_getEntSetData(imeshImpl, *sit, tags[t], 
                          &tag_space_ptr, &tag_space_alloc, &tag_space_size, &err);
      if (iBase_TAG_NOT_FOUND == err) continue;
      ERRORR("Problem getting copy tag for set.", err);

        // compare to tag value if necessary
      if (tag_vals && tag_vals[t] && strncmp(tag_vals[t], tag_space_ptr, tag_size))
        continue;
      
        // if we got here, we should set the tag on the copy; get the copy
      iBase_EntitySetHandle copy_set;
      iMesh_getEntSetEHData(imeshImpl, *sit, copyTag, &copy_set, &err);
      ERRORR("Didn't get copied set from orig copy set.", err);
      iMesh_setEntSetData(imeshImpl, copy_set, tags[t], 
                          tag_space_ptr, tag_size, &err);
      ERRORR("Failed to set data on copied set.", err);
    }
  }
  
  return err;
}
    
#ifdef TEST

#include <cfloat>
#include "MergeMesh.hpp"

iMesh_Instance impl;
iBase_EntityHandle root_set;

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
    iBase_EntityHandle *new_ents = ents + ents_size/2;
    new_ents_alloc = ents_size/2;
    
  
    err = cm->copy_move_entities(ents, ents_size, dx,
                                 &new_ents, &new_ents_alloc, &new_ents_size, 
                                 false);
    ERRORR("Failed to copy_move entities.", 1);

      // check # entities
    if (check_num_ents(iBase_VERTEX, 32)) return 1;
    if (check_num_ents(iBase_REGION, 6)) return 1;
  
    err = mm->merge_entities(ents, ents_size+new_ents_size, 1.0e-8);
    ERRORR("Failed to merge entities.", 1);

    double xmin[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double xmax[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
    double *ncoords = NULL;
    int ncoords_alloc = 0, ncoords_size;
    int *inset = NULL, inset_alloc = 0, inset_size;
    ERRORR("Couldn't get root set.", 1);
    int sorder = iBase_INTERLEAVED;
    iMesh_getAllVtxCoords (impl, root_set,
                           &ncoords, &ncoords_alloc, &ncoords_size,
                           &inset, &inset_alloc, &inset_size, &sorder, &err);
    ERRORR("Didn't get vtx coords.", 1);

    int num_err = 0;
    for (int i = 0; i < inset_size; i++) {
      if (!inset[i]) num_err++;
      for (int j = 0; j < 3; j++) {
        if (ncoords[3*i+j] < xmin[j]) xmin[j] = ncoords[3*i+j];
        if (ncoords[3*i+j] > xmax[j]) xmax[j] = ncoords[3*i+j];
      }
    }
    if (xmin[0] != 0.0 || xmin[1] != 0.0 || xmin[2] != 0.0 ||
        xmax[0] != 2.0 || xmax[1] != 1.0 || xmax[2] != 3.0) {
      std::cerr << "Didn't get correct min/max; output values = ("
                << xmin[0] << ", " << xmin[1] << ", " << xmin[2] << "), ("
                << xmax[0] << ", " << xmax[1] << ", " << xmax[2] << ")" << std::endl;
      return 1;
    }

    iMesh_getNumOfType(impl, root_set, iBase_VERTEX, &num_err, &err);
    ERRORR("Didn't get num verts.", 1);
  
    int expected = 2*verts_size - 8;
    if (num_err != expected) {
      std::cerr << "Didn't get right # vertices; got " << num_err 
                << ", expected " << expected << std::endl;
    
      return 1;
    }
  
    free(verts);
    free(ents);
    free(ncoords);
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
    err = cm->copy_move_entities(&root_set, 1, dx, NULL, NULL, NULL,
                                 false);
    ERRORR("Failed to copy_move entities.", 1);
    
      // merge
    err = mm->merge_entities(NULL, 0, 1.0e-8);
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
  iMesh_addEntArrToSet(cm->impl(), ents, 3, esets, &err);
  ERRORR("Failed to add ents to set 1.", err);
  int dum = 3;
  iMesh_setEntSetIntData(cm->impl(), esets[0], gtag, dum, &err);
  ERRORR("Failed to set geom set tag on set 1.", err);
  
    // ES2: geom dim = 2, faces all, -
  iBase_EntityHandle *faces = NULL;
  int faces_alloc = 0, faces_size;
  iMesh_getEntities(cm->impl(), root_set, iBase_FACE, iMesh_ALL_TOPOLOGIES,
                    &faces, &faces_alloc, &faces_size, &err);
  ERRORR("Failed to get faces.", err);
  iMesh_addEntArrToSet(cm->impl(), faces, faces_size, esets+1, &err);
  ERRORR("Failed to add ents to set 2.", err);
  dum = 2;
  iMesh_setEntSetIntData(cm->impl(), esets[1], gtag, dum, &err);
  ERRORR("Failed to set geom set tag on set 2.", err);
  free(faces);

    // ES3: - , skin nodes, E
  iMesh_addEntArrToSet(cm->impl(), verts, verts_size, esets+2, &err);
  ERRORR("Failed to add ents to set 3.", err);
  err = cm->add_copy_expand_list(esets+2, 1, CopyMesh::EXPAND);
  ERRORR("Failed to add set 3 to expand list.", err);

    // ES4: mat set, set 1, E
  iMesh_addEntSet(cm->impl(), esets[0], esets+3, &err);
  ERRORR("Failed to add set to set 4.", err);
  dum = 100;
  iMesh_setEntSetIntData(cm->impl(), esets[3], mtag, dum, &err);
  ERRORR("Failed to set mat set tag on set 4.", err);

    // ES5: neu set, set 2, E
  iMesh_addEntSet(cm->impl(), esets[1], esets+4, &err);
  ERRORR("Failed to add set to set 5.", err);
  dum = 101;
  iMesh_setEntSetIntData(cm->impl(), esets[4], ntag, dum, &err);
  ERRORR("Failed to set neu set tag on set 5.", err);

    // ES6: dir set, set 3, E
  iMesh_addEntSet(cm->impl(), esets[2], esets+5, &err);
  ERRORR("Failed to add set to set 6.", err);
  dum = 102;
  iMesh_setEntSetIntData(cm->impl(), esets[5], dtag, dum, &err);
  ERRORR("Failed to set dir set tag on set 6.", err);

    // ES7: -, set 4, C
  iMesh_addEntSet(cm->impl(), esets[3], esets+6, &err);
  ERRORR("Failed to add set to set 7.", err);
  err = cm->add_copy_expand_list(esets+6, 1, CopyMesh::COPY);
  ERRORR("Failed to add set 7 to copy list.", err);
	
    // check # entities
  if (check_num_ents(iBase_VERTEX, 16)) return 1;
  if (check_num_ents(iBase_REGION, 3)) return 1;

  return 0;
}

#endif
