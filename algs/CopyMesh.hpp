#ifndef COPYMESH_HPP
#define COPYMESH_HPP

#include "iMesh_extensions.h"
#include "CopyVerts.hpp"

#include <string>
#include <vector>
#include <set>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cstring>

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

class CopyMesh 
{
public:
    /* \brief Constructor
     */
  CopyMesh(iMesh_Instance impl);
  
    /* \brief Destructor
     */
  virtual ~CopyMesh();

    /* \brief Return the imesh instance
     */
  iMesh_Instance impl() const {return imeshImpl;}

    /* \brief Return the copyTag used to indicate set copies
     */
  iBase_TagHandle copy_tag() {return copyTag;}
  
    /* \brief Set the copy or expand set list
     */
  int add_copy_expand_list(iBase_EntitySetHandle *ce_sets, int num_ce_sets,
                           int copy_or_expand);
  
    /* \brief Reset the copy and expand set lists
     */
  int reset_ce_lists();

    /* \brief Add tag indicating sets to copy w/ new entities
     */
  int add_copy_tag(const std::string &tag_name, const char *tag_val = NULL);
  
    /* \brief Add tag indicating sets to copy w/ new entities
     */
  int add_copy_tag(iBase_TagHandle tag_handle, const char *tag_val = NULL);
  
    /* \brief Add tag indicating sets to expand w/ new entities
     */
  int add_expand_tag(const std::string &tag_name, const char *tag_val = NULL);

    /* \brief Add tag indicating sets to expand w/ new entities
     */
  int add_expand_tag(iBase_TagHandle tag_handle, const char *tag_val = NULL);

    /* \brief Add tag which should have unique values
     */
  int add_unique_tag(const std::string &tag_name);

    /* \brief Add tag which should have unique values
     */
  int add_unique_tag(iBase_TagHandle tag_handle);
  
    /* \brief Return reference to copy sets 
     */
  std::set<iBase_EntitySetHandle> &copy_sets();
  
    /* \brief Return reference to expand sets 
     */
  std::set<iBase_EntitySetHandle> &expand_sets();
  
    /* \brief Return reference to unique sets 
     */
  std::set<iBase_EntitySetHandle> &unique_sets();
  
    /* \brief Copy all the entities in the set
     */
  int copy_entities(iBase_EntitySetHandle set_handle,
                    iBase_EntityHandle **new_ents = NULL,
                    int *new_ents_allocated = 0,
                    int *new_ents_size = 0);
  
    /* \brief Copy all the entities
     */
  int copy_entities(iBase_EntityHandle *ent_handles,
                    int num_ents,
                    iBase_EntityHandle **new_ents = NULL,
                    int *new_ents_allocated = 0,
                    int *new_ents_size = 0);
  
    /* \brief Copy and move all the entities
     */
/*  int copy_move_entities(iBase_EntityHandle *ent_handles,
                         int num_ents,
                         const double dx,
                         const double dy,
                         const double dz,
                         iBase_EntityHandle **new_ents = NULL,
                         int *new_ents_allocated = 0,
                         int *new_ents_size = 0);*/
  
    /* \brief Copy/move all entities in a set
     */
  int copy_move_entities(iBase_EntitySetHandle set_handle,
                         const double *dx,
                         iBase_EntityHandle **new_ents,
                         int *new_ents_alloc,
                         int *new_ents_size,
                         bool do_merge = true);
  
    /* \brief Copy and move all the entities
     */
  int copy_move_entities(iBase_EntityHandle *ent_handles,
                         int num_ents,
                         const double dx[3],
                         iBase_EntityHandle **new_ents = NULL,
                         int *new_ents_allocated = 0,
                         int *new_ents_size = 0,
                         bool do_merge = true);

  int copy_transform_entities(iBase_EntitySetHandle set_handle,
                              const CopyVerts &trans,
                              iBase_EntityHandle **new_ents = NULL,
                              int *new_ents_allocated = 0,
                              int *new_ents_size = 0,
                              bool do_merge = true);

  int copy_transform_entities(iBase_EntityHandle *ent_handles,
                              int num_ents,
                              const CopyVerts &trans,
                              iBase_EntityHandle **new_ents = NULL,
                              int *new_ents_allocated = 0,
                              int *new_ents_size = 0,
                              bool do_merge = true);

  int update_ce_lists();
  
  int process_ce_sets(std::set<iBase_EntitySetHandle> &cesets,
                      int copy_or_expand,
                      iBase_TagHandle local_tag);

    /* \brief Tag copied sets with indicated tag from original set
     */
  int tag_copied_sets(const char **tag_names, const char **tag_vals,
                      const int num_tags);
  
    /* \brief Tag copied sets with indicated tag from original set
     */
  int tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
                      const int num_tags);
  
  enum {COPY = 0, EXPAND, UNIQUE};
  
private:
  
    //- copy entities, from connect indices in ind, types in topos, and
    //- # verts in offset
  int copy_move_ents(iBase_EntitySetHandle copy_set,
                     iBase_TagHandle local_tag);
  
    //- copy/move vertices, putting results in new_ents
  int copy_transform_verts(iBase_EntitySetHandle copy_set,
                           const CopyVerts &copier,
                           iBase_TagHandle local_tag);

    //- get the copy/expand sets based on copy/expand tags
  int get_copy_expand_sets(iBase_EntitySetHandle *&copy_sets,
                           int &num_copy_sets,
                           iBase_EntitySetHandle *&expand_sets,
                           int &num_expand_sets);
  
    //- get the sets tagged with the given vector of tags/values
  int get_tagged_sets(iBase_EntitySetHandle from_set,
                      iBase_TagHandle *tag_handles,
                      const char **tag_vals,
                      int num_tags,
                      iBase_EntitySetHandle *&tagged_sets,
                      int &num_tagged_sets);
  
  int update_tagged_sets(iBase_EntitySetHandle from_set,
                         iBase_TagHandle *tag_handles,
                         const char **tag_vals,
                         int num_tags,
                         std::set<iBase_EntitySetHandle> &tagged_sets);

    //- interface instance
  iMesh_Instance imeshImpl;

    //- tags indicating which sets should be expanded to include new ents
  std::vector<iBase_TagHandle> expandTags;

    //- tag values indicating which sets should be expanded to include new ents
  std::vector<const char*> expandTagVals;

    //- expand sets, input or found from searching; use set to allow multiple additions
  std::set<iBase_EntitySetHandle> expandSets;

    //- tags indicating which sets should be copied with new entities
  std::vector<iBase_TagHandle> copyTags;
  
    //- tag values indicating which sets should be copied to include new ents
  std::vector<const char*> copyTagVals;

    //- copy sets, input or found from searching; use set to allow multiple additions
  std::set<iBase_EntitySetHandle> copySets;

    //- tags indicating which sets should be copied with new entities
  std::vector<iBase_TagHandle> uniqueTags;
  
    //- tag values indicating which sets should be copied to include new ents
  std::vector<const char*> uniqueTagVals;

    //- unique sets, input or found from searching
  std::set<iBase_EntitySetHandle> uniqueSets;

    //- flag denoting whether copy/expand sets lists have been updated
  bool updatedCELists;
  
    //- tag storing copy-to tag
  iBase_TagHandle copyTag;
  
};

inline CopyMesh::~CopyMesh() 
{
  if (0 != copyTag) {
    int err;
    iMesh_destroyTag(imeshImpl, copyTag, true, &err);
    if (iBase_SUCCESS != err) ERROR("Failed to destroy copyTag.");
  }
}

inline int CopyMesh::add_copy_tag(const std::string &tag_name, const char *tag_val) 
{
  iBase_TagHandle tag_handle = 0;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err, tag_name.length());
  if (iBase_SUCCESS != err) {
    std::string tmp_str("Failed to get handle for tag ");
    tmp_str += tag_name;
    ERRORR(tmp_str, iBase_FAILURE);
  }

  return add_copy_tag(tag_handle, tag_val);
}

  /* \brief Add tag indicating sets to copy w/ new entities
   */
inline int CopyMesh::add_copy_tag(iBase_TagHandle tag_handle, const char *tag_val)
{
  assert(0 != tag_handle);
  copyTags.push_back(tag_handle);
  int err = iBase_SUCCESS;

  if (tag_val) {
    int tag_size;
    iMesh_getTagSizeBytes(imeshImpl, tag_handle, &tag_size, &err);
    if (iBase_SUCCESS != err) {
      std::string tmp_str("Failed to get size of tag");
      ERRORR(tmp_str, iBase_FAILURE);
    }
    char *tmp_mem = (char*) malloc(tag_size);
    memcpy(tmp_mem, tag_val, tag_size);
    copyTagVals.push_back(tmp_mem);
  }
  else
    copyTagVals.push_back(NULL);

  if (!copyTag) {
    iMesh_createTag(imeshImpl, "__CopyMeshTag", 1,
                    iBase_ENTITY_HANDLE, &copyTag, &err, 13);
    ERROR("Couldn't create copy mesh tag.");
  }
  
  return err;
}

inline int CopyMesh::add_unique_tag(const std::string &tag_name) 
{
  iBase_TagHandle tag_handle = 0;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err, tag_name.length());
  if (iBase_SUCCESS != err) {
    std::string tmp_str("Failed to get handle for tag ");
    tmp_str += tag_name;
    ERRORR(tmp_str, iBase_FAILURE);
  }

  return add_unique_tag(tag_handle);
}

  /* \brief Add tag indicating sets to unique w/ new entities
   */
inline int CopyMesh::add_unique_tag(iBase_TagHandle tag_handle)
{
  assert(0 != tag_handle);
  uniqueTags.push_back(tag_handle);
  return iBase_SUCCESS;
}

  /* \brief Add tag indicating sets to expand w/ new entities
   */
inline int CopyMesh::add_expand_tag(const std::string &tag_name, const char *tag_val)
{
  iBase_TagHandle tag_handle = 0;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err, tag_name.length());
  if (iBase_SUCCESS != err) {
    std::string tmp_str("Failed to get handle for tag ");
    tmp_str += tag_name;
    ERRORR(tmp_str, iBase_FAILURE);
  }

  return add_expand_tag(tag_handle, tag_val);
}

  /* \brief Add tag indicating sets to expand w/ new entities
   */
inline int CopyMesh::add_expand_tag(iBase_TagHandle tag_handle, const char *tag_val)
{
  assert(0 != tag_handle);
  expandTags.push_back(tag_handle);
  int err = iBase_SUCCESS;

  if (tag_val) {
    int tag_size;
    iMesh_getTagSizeBytes(imeshImpl, tag_handle, &tag_size, &err);
    if (iBase_SUCCESS != err) {
      std::string tmp_str("Failed to get size of tag");
      ERRORR(tmp_str, iBase_FAILURE);
    }
    char *tmp_mem = (char*) malloc(tag_size);
    memcpy(tmp_mem, tag_val, tag_size);
    expandTagVals.push_back(tmp_mem);
  }
  else
    expandTagVals.push_back(NULL);

  return (int) err;
}

inline int CopyMesh::reset_ce_lists() 
{
  copySets.clear();
  expandSets.clear();
  updatedCELists = false;
  return iBase_SUCCESS;
}

  /* \brief Return reference to copy sets 
   */
inline std::set<iBase_EntitySetHandle> &CopyMesh::copy_sets() {return copySets;}
  
  /* \brief Return reference to expand sets 
   */
inline std::set<iBase_EntitySetHandle> &CopyMesh::expand_sets() {return expandSets;}
  
  /* \brief Return reference to unique sets 
   */
inline std::set<iBase_EntitySetHandle> &CopyMesh::unique_sets() {return uniqueSets;}

#endif
