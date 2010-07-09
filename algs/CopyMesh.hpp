#ifndef COPYMESH_HPP
#define COPYMESH_HPP

#include "iMesh_extensions.h"
#include "CopyVerts.hpp"
#include "LocalTag.hpp"
#include "MKException.hpp"

#include <cassert>
#include <string>
#include <vector>
#include <set>

/**\brief Add newly-created entities/sets to a collection of sets
 *
 * Given a collection of copy, expand, or extrude source sets and a tag, create 
 * a destination (copy) set unless one already exists. Fill this set with any
 * new entities/sets created from those in the source set.
 * \param imeshImpl the iMesh instance handle
 * \param cesets a collection of source sets
 * \param local_tag the tag relating source and target entities/sets
 */
void process_ce_sets(iMesh_Instance imeshImpl,
                     const std::set<iBase_EntitySetHandle> &cesets,
                     iBase_TagHandle local_tag);

/**\brief Tag a collection of copied sets
 *
 * Given a collection of source sets and a tag |copyTag| relating sources to
 * destinations, apply a tag |tag| to the destination sets if the tag exists on
 * the corresponding source.
 * \param imeshImpl the iMesh instance handle
 * \param copyTag the tag relating sources and destinations
 * \param cesets a collection of source sets
 * \param tag the tag to set on the destinations
 * \param tag_val if non-NULL, only set |tag| on the destination if the source's
 *                tag matches this value
 */
void tag_copy_sets(iMesh_Instance imeshImpl, iBase_TagHandle copyTag,
                   const std::set<iBase_EntitySetHandle> &copySets,
                   iBase_TagHandle tag, const char *tag_val);

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
  
  /* \brief Set the copy or expand set list
   */
  void add_copy_expand_list(iBase_EntitySetHandle *ce_sets, int num_ce_sets,
                            int copy_or_expand);
  
  /* \brief Reset the copy and expand set lists
   */
  void reset_ce_lists();

  /* \brief Add tag indicating sets to copy w/ new entities
   */
  void add_copy_tag(const std::string &tag_name, const char *tag_val = NULL);
  
  /* \brief Add tag indicating sets to copy w/ new entities
   */
  void add_copy_tag(iBase_TagHandle tag_handle, const char *tag_val = NULL);
  
  /* \brief Add tag indicating sets to expand w/ new entities
   */
  void add_expand_tag(const std::string &tag_name, const char *tag_val = NULL);

  /* \brief Add tag indicating sets to expand w/ new entities
   */
  void add_expand_tag(iBase_TagHandle tag_handle, const char *tag_val = NULL);

  /* \brief Add tag which should have unique values
   */
  void add_unique_tag(const std::string &tag_name);

  /* \brief Add tag which should have unique values
   */
  void add_unique_tag(iBase_TagHandle tag_handle);
  
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

  void update_ce_lists();
  
  /* \brief Tag copied sets with indicated tag from original set
   */
  void tag_copied_sets(const char **tag_names, const char **tag_vals,
                       const int num_tags);
  
  /* \brief Tag copied sets with indicated tag from original set
   */
  void tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
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
  
  void update_tagged_sets(iBase_EntitySetHandle from_set,
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

inline void CopyMesh::add_copy_tag(const std::string &tag_name,
                                   const char *tag_val) 
{
  iBase_TagHandle tag_handle;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err,
                     tag_name.size());
  check_error(imeshImpl, err);

  add_copy_tag(tag_handle, tag_val);
}

inline void CopyMesh::add_expand_tag(const std::string &tag_name,
                                    const char *tag_val)
{
  iBase_TagHandle tag_handle;
  int err;
  iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err,
                     tag_name.size());
  check_error(imeshImpl, err);

  add_expand_tag(tag_handle, tag_val);
}

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

inline void CopyMesh::reset_ce_lists()
{
  std::vector<tag_data>::iterator i;
  for (i = copyTags.begin(); i != copyTags.end(); ++i)
    free(i->value);

  for (i = expandTags.begin(); i != expandTags.end(); ++i)
    free(i->value);

  copySets.clear();
  expandSets.clear();
  updatedCELists = false;
}

inline std::set<iBase_EntitySetHandle> &CopyMesh::copy_sets()
{
  return copySets;
}
  
inline std::set<iBase_EntitySetHandle> &CopyMesh::expand_sets()
{
  return expandSets;
}
  
inline std::set<iBase_EntitySetHandle> &CopyMesh::unique_sets()
{
  return uniqueSets;
}

#endif
