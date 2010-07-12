#ifndef COPYMESH_HPP
#define COPYMESH_HPP

#include "iMesh_extensions.h"
#include "CopyVerts.hpp"
#include "LocalTag.hpp"
#include "CESets.hpp"
#include "MKException.hpp"

#include <cassert>
#include <string>
#include <vector>
#include <set>

/**\brief Get the entities and unique adjacent vertices of a set
 *
 * Return an array of all entities in a set, sorted by topology, and all
 * unique adjacent vertices. The adjacent vertices can be accessed with an
 * offset array indexing into an index buffer.
 * \param imeshImpl the iMesh instance handle
 * \param set the set of entities from which to query
 * \param ents pointer to an array of entity handles in the set
 * \param unique_adj pointer to an array of unique vertices adjacent to |ents|
 * \param indices index buffer into |unique_adj|
 * \param offsets offset array indicating start and end of |indices| for each
                  entity in |ents|
 * \param err pointer to an ITAPS error code 
 */
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
  
  /* \brief Copy all the entities in the set
   */
  void copy_entities(iBase_EntitySetHandle set_handle,
                     iBase_EntityHandle **new_ents = NULL,
                     int *new_ents_allocated = 0,
                     int *new_ents_size = 0);
  
  /* \brief Copy all the entities
   */
  void copy_entities(iBase_EntityHandle *ent_handles,
                     int num_ents,
                     iBase_EntityHandle **new_ents = NULL,
                     int *new_ents_allocated = 0,
                     int *new_ents_size = 0);
  
  /* \brief Copy/move all entities in a set
   */
  void copy_move_entities(iBase_EntitySetHandle set_handle,
                          const double *dx,
                          iBase_EntityHandle **new_ents = NULL,
                          int *new_ents_alloc = 0,
                          int *new_ents_size = 0,
                          bool do_merge = true);
  
  /* \brief Copy and move all the entities
   */
  void copy_move_entities(iBase_EntityHandle *ent_handles,
                          int num_ents,
                          const double dx[3],
                          iBase_EntityHandle **new_ents = NULL,
                          int *new_ents_allocated = 0,
                          int *new_ents_size = 0,
                          bool do_merge = true);

  void copy_rotate_entities(iBase_EntitySetHandle set_handle,
                            const double *origin,
                            const double *z,
                            double theta,
                            iBase_EntityHandle **new_ents = NULL,
                            int *new_ents_allocated = 0,
                            int *new_ents_size = 0,
                            bool do_merge = true);

  void copy_rotate_entities(iBase_EntityHandle *ent_handles,
                            int num_ents,
                            const double *origin,
                            const double *z,
                            double theta,
                            iBase_EntityHandle **new_ents = NULL,
                            int *new_ents_allocated = 0,
                            int *new_ents_size = 0,
                            bool do_merge = true);

  void copy_transform_entities(iBase_EntitySetHandle set_handle,
                               const CopyVerts &trans,
                               iBase_EntityHandle **new_ents = NULL,
                               int *new_ents_allocated = 0,
                               int *new_ents_size = 0,
                               bool do_merge = true);

  void copy_transform_entities(iBase_EntityHandle *ent_handles,
                               int num_ents,
                               const CopyVerts &trans,
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
  void connect_the_dots(iBase_EntityHandle *ents, int size,
                        iBase_TagHandle local_tag, int *indices, int *offsets,
                        iBase_EntityHandle *verts);

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
