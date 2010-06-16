#ifndef COPYMESH_HPP
#define COPYMESH_HPP

#include "iMesh_extensions.h"
#include "CopyVerts.hpp"
#include "LocalTag.hpp"

#include <string>
#include <vector>
#include <set>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

int process_ce_sets(iMesh_Instance imeshImpl,
                    std::set<iBase_EntitySetHandle> &cesets,
                    iBase_TagHandle local_tag);

int tag_copy_sets(iMesh_Instance imeshImpl, iBase_TagHandle copyTag,
                  const std::set<iBase_EntitySetHandle> &copySets,
                  iBase_TagHandle tag, const char *tag_val);

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
                        int *err);

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
  
  /* \brief Copy/move all entities in a set
   */
  int copy_move_entities(iBase_EntitySetHandle set_handle,
                         const double *dx,
                         iBase_EntityHandle **new_ents = NULL,
                         int *new_ents_alloc = 0,
                         int *new_ents_size = 0,
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

  int copy_rotate_entities(iBase_EntitySetHandle set_handle,
                           const double *origin,
                           const double *z,
                           double theta,
                           iBase_EntityHandle **new_ents = NULL,
                           int *new_ents_allocated = 0,
                           int *new_ents_size = 0,
                           bool do_merge = true);

  int copy_rotate_entities(iBase_EntityHandle *ent_handles,
                           int num_ents,
                           const double *origin,
                           const double *z,
                           double theta,
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
  struct tag_data
  {
    tag_data(iBase_TagHandle tag, char *value)
      : tag(tag), value(value)
    {}

    iBase_TagHandle tag;
    char *value;
  };

  int connect_the_dots(iBase_EntityHandle *ents, int size,
                       iBase_TagHandle local_tag, int *indices, int *offsets,
                       iBase_EntityHandle *verts);

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
                         const std::vector<tag_data> &tags,
                         std::set<iBase_EntitySetHandle> &tagged_sets);

  iMesh_Instance imeshImpl;

  std::vector<tag_data> expandTags;
  std::set<iBase_EntitySetHandle> expandSets;

  std::vector<tag_data> copyTags;
  std::set<iBase_EntitySetHandle> copySets;

  std::vector<iBase_TagHandle> uniqueTags;
  std::set<iBase_EntitySetHandle> uniqueSets;

  // flag denoting whether copy/expand sets lists have been updated
  bool updatedCELists;
  
  // tag storing copy-to tag
  LocalTag copyTag;
};

inline int CopyMesh::add_copy_tag(const std::string &tag_name,
                                  const char *tag_val) 
{
  iBase_TagHandle tag_handle;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err,
                     tag_name.size());
  if (iBase_SUCCESS != err)
    ERRORR("Failed to get handle for tag "+tag_name, iBase_FAILURE);

  return add_copy_tag(tag_handle, tag_val);
}

/* \brief Add tag indicating sets to copy w/ new entities
 */
inline int CopyMesh::add_copy_tag(iBase_TagHandle tag_handle,
                                  const char *tag_val)
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

/* \brief Add tag indicating sets to expand w/ new entities
 */
inline int CopyMesh::add_expand_tag(const std::string &tag_name,
                                    const char *tag_val)
{
  iBase_TagHandle tag_handle;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err,
                     tag_name.size());
  if (iBase_SUCCESS != err)
    ERRORR("Failed to get handle for tag "+tag_name, iBase_FAILURE);

  return add_expand_tag(tag_handle, tag_val);
}

/* \brief Add tag indicating sets to expand w/ new entities
 */
inline int CopyMesh::add_expand_tag(iBase_TagHandle tag_handle,
                                    const char *tag_val)
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

inline int CopyMesh::add_unique_tag(const std::string &tag_name) 
{
  iBase_TagHandle tag_handle;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err,
                     tag_name.size());
  if (iBase_SUCCESS != err)
    ERRORR("Failed to get handle for tag "+tag_name, iBase_FAILURE);

  return add_unique_tag(tag_handle);
}

/* \brief Add tag indicating sets to unique w/ new entities
 */
inline int CopyMesh::add_unique_tag(iBase_TagHandle tag_handle)
{
  assert(tag_handle != NULL);
  uniqueTags.push_back(tag_handle);
  return iBase_SUCCESS;
}

inline int CopyMesh::reset_ce_lists()
{
  std::vector<tag_data>::iterator i;
  for (i = copyTags.begin(); i != copyTags.end(); ++i)
    free(i->value);

  for (i = expandTags.begin(); i != expandTags.end(); ++i)
    free(i->value);

  copySets.clear();
  expandSets.clear();
  updatedCELists = false;
  return iBase_SUCCESS;
}

/* \brief Return reference to copy sets 
 */
inline std::set<iBase_EntitySetHandle> &CopyMesh::copy_sets()
{
  return copySets;
}
  
/* \brief Return reference to expand sets 
 */
inline std::set<iBase_EntitySetHandle> &CopyMesh::expand_sets()
{
  return expandSets;
}
  
/* \brief Return reference to unique sets 
 */
inline std::set<iBase_EntitySetHandle> &CopyMesh::unique_sets()
{
  return uniqueSets;
}

#endif
