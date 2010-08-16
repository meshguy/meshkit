#ifndef COPYMESH_HPP
#define COPYMESH_HPP

#include "iMesh_extensions.h"
#include "Transform.hpp"
#include "LocalTag.hpp"
#include "CESets.hpp"
#include "MKException.hpp"

#include <cassert>
#include <string>
#include <vector>
#include <set>

class CopyMesh 
{
public:
  /* \brief Constructor
   *
   * Create a new CopyMesh instance
   * \param impl the iMesh instance handle for the mesh
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

  void update_sets();

  /* \brief Reset the copy and expand set lists
   */
  void reset_sets();

  /* \brief Add tag which should have unique values
   */
  void add_unique_tag(const std::string &tag_name);

  /* \brief Add tag which should have unique values
   */
  void add_unique_tag(iBase_TagHandle tag_handle);
  
  /* \brief Return reference to copy sets 
   */
  CESets &copy_sets();
  
  /* \brief Return reference to expand sets 
   */
  CESets &expand_sets();
  
  /* \brief Return reference to unique sets 
   */
  std::set<iBase_EntitySetHandle> &unique_sets();
  
  void copy(iBase_EntitySetHandle set_handle,
            const copy::Transform &trans = copy::Identity(),
            iBase_EntityHandle **new_ents = NULL,
            int *new_ents_allocated = 0,
            int *new_ents_size = 0,
            bool do_merge = true);

  void copy(iBase_EntityHandle *ent_handles,
            int num_ents,
            const copy::Transform &trans = copy::Identity(),
            iBase_EntityHandle **new_ents = NULL,
            int *new_ents_allocated = 0,
            int *new_ents_size = 0, 
            bool do_merge = true);

  /* \brief Tag copied sets with indicated tag from original set
   */
  void tag_copied_sets(const char **tag_names, const char **tag_vals,
                       const int num_tags);
  
  /* \brief Tag copied sets with indicated tag from original set
   */
  void tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
                       const int num_tags);
  
private:
  //- get the copy/expand sets based on copy/expand tags
  void get_copy_expand_sets(iBase_EntitySetHandle *&copy_sets,
                            int &num_copy_sets,
                            iBase_EntitySetHandle *&expand_sets,
                            int &num_expand_sets);
  
  //- get the sets tagged with the given vector of tags/values
  void get_tagged_sets(iBase_EntitySetHandle from_set,
                       iBase_TagHandle *tag_handles,
                       const char **tag_vals,
                       int num_tags,
                       iBase_EntitySetHandle *&tagged_sets,
                       int &num_tagged_sets);
  
  iMesh_Instance imeshImpl;

  // tag storing copy-to tag
  LocalTag copyTag;

  CESets copySets;
  CESets expandSets;

  std::vector<iBase_TagHandle> uniqueTags;
  std::set<iBase_EntitySetHandle> uniqueSets;
};

inline void CopyMesh::add_unique_tag(const std::string &tag_name) 
{
  iBase_TagHandle tag_handle;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err,
                     tag_name.size());
  check_error(imeshImpl, err);

  add_unique_tag(tag_handle);
}

inline void CopyMesh::add_unique_tag(iBase_TagHandle tag_handle)
{
  assert(tag_handle != NULL);
  uniqueTags.push_back(tag_handle);
}

inline void CopyMesh::reset_sets()
{
  copySets.clear();
  expandSets.clear();
}

inline CESets &CopyMesh::copy_sets()
{
  return copySets;
}
  
inline CESets &CopyMesh::expand_sets()
{
  return expandSets;
}
  
inline std::set<iBase_EntitySetHandle> &CopyMesh::unique_sets()
{
  return uniqueSets;
}

#endif
