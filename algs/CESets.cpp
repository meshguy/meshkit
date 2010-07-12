#include "CESets.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iMesh_extensions.h>

CESets::~CESets()
{
  std::vector<tag_data>::iterator i;
  for (i = tags_.begin(); i != tags_.end(); ++i)
    free(i->value);
}

void CESets::add_tag(iBase_TagHandle tag_handle, const char *value)
{
  assert(tag_handle != NULL);
  char *tmp = NULL;

  if (value) {
    int err;
    int tag_size;
    iMesh_getTagSizeBytes(impl_, tag_handle, &tag_size, &err);
    check_error(impl_, err);

    tmp = static_cast<char*>(malloc(tag_size));
    memcpy(tmp, value, tag_size);
  }

  tags_.push_back(tag_data(tag_handle, tmp));
}

void CESets::add_tag(const std::string &tag_name, const char *value)
{
  iBase_TagHandle tag_handle;
  int err;
  iMesh_getTagHandle(impl_, tag_name.c_str(), &tag_handle, &err,
                     tag_name.size());
  check_error(impl_, err);

  add_tag(tag_handle, value);
}

void CESets::update_tagged_sets()
{
  int err;
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(impl_, &root_set, &err);
  check_error(impl_, err);

  std::vector<tag_data>::const_iterator tag;
  for (tag = tags_.begin(); tag != tags_.end(); ++tag) {
    SimpleArray<iBase_EntitySetHandle> tmp_sets;
    iMesh_getEntSetsByTagsRec(impl_, root_set, &tag->tag,
                              (tag->value ? &tag->value : NULL), 1, false,
                              ARRAY_INOUT(tmp_sets), &err);
    check_error(impl_, err);
    sets_.insert(tmp_sets.begin(), tmp_sets.end());
  }
}

/**\brief Get the entities copied from a set
 *
 * Given a set and a tag, return the values of the tag set on entities from that
 * set (i.e. the copied entities).
 * \param imeshImpl the iMesh instance handle
 * \param set the set containing the entities we want
 * \param local_tag the tag relating source and target entities
 * \param ents a vector which will hold our copied entities
 */
static
void get_copied_ents(iMesh_Instance imeshImpl, iBase_EntitySetHandle set,
                     iBase_TagHandle local_tag,
                     std::vector<iBase_EntityHandle> &ents)
{
  int err;

  SimpleArray<iBase_EntityHandle> tmp_ents;
  iMesh_getEntities(imeshImpl, set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                    ARRAY_INOUT(tmp_ents), &err);
  check_error(imeshImpl, err);

  ents.reserve(tmp_ents.size());

  // get copied entities
  for (int i = 0; i < tmp_ents.size(); i++) {
    iBase_EntityHandle eh_tag;
    iMesh_getEHData(imeshImpl, tmp_ents[i], local_tag, &eh_tag, &err);
    if (iBase_SUCCESS == err && eh_tag)
      ents.push_back(eh_tag);
  }
}

/**\brief Create a copy set if one doesn't exist yet
 *
 * Given a source set and a tag, create a destination (copy) set unless one
 * already exists.
 * \param imeshImpl the iMesh instance handle
 * \param local_tag the tag relating source and target sets
 * \param src the source set
 * \param dest pointer to the destination set
 */
static
void get_dest_set(iMesh_Instance imeshImpl, iBase_TagHandle local_tag,
                  iBase_EntitySetHandle src, iBase_EntitySetHandle *dest)
{
  int err;

  if (dest == NULL) return;

  iMesh_getEntSetEHData(imeshImpl, src, local_tag,
                        reinterpret_cast<iBase_EntityHandle*>(dest), &err);

  if (err != iBase_SUCCESS) {
    iMesh_createEntSet(imeshImpl, false, dest, &err);
    check_error(imeshImpl, err);
    iMesh_setEntSetEHData(imeshImpl, src, local_tag,
                          reinterpret_cast<iBase_EntityHandle>(*dest), &err);
    check_error(imeshImpl, err);
  }
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
 */
static
void process_ce_subsets(iMesh_Instance imeshImpl, iBase_EntitySetHandle src,
                        iBase_EntitySetHandle current,
                        const std::set<iBase_EntitySetHandle> &cesets,
                        iBase_TagHandle local_tag) 
{
  int err;
  iBase_EntitySetHandle dest;

  // First, add entities directly contained in this set.
  std::vector<iBase_EntityHandle> tmp_tags;
  get_copied_ents(imeshImpl, current, local_tag, tmp_tags);

  if (!tmp_tags.empty()) {
    get_dest_set(imeshImpl, local_tag, src, &dest);
    iMesh_addEntArrToSet(imeshImpl, &tmp_tags[0], tmp_tags.size(), dest,
                         &err);
    check_error(imeshImpl, err);
  }

  // Next, start looking at children.
  SimpleArray<iBase_EntitySetHandle> children;
  iMesh_getEntSets(imeshImpl, current, 1, ARRAY_INOUT(children), &err);
  check_error(imeshImpl, err);

  for (int i = 0; i < children.size(); i++) {

    // If this child set is one of our cesets, add just the set...
    if (cesets.find(children[i]) != cesets.end()) {
      get_dest_set(imeshImpl, local_tag, src, &dest);
      if (src == dest) continue;

      iMesh_addEntSet(imeshImpl, children[i], dest, &err);
      check_error(imeshImpl, err);
    }

    // ... otherwise, add the entities and recurse into the next level of
    // children.
    else {
      process_ce_subsets(imeshImpl, src, children[i], cesets, local_tag);
    }
  }
}

// xxx - still need to do unique tags
// - for each unique tag
//   . if this set has the tag, add/append to get unique value
void process_ce_sets(iMesh_Instance imeshImpl,
                     const std::set<iBase_EntitySetHandle> &cesets,
                     iBase_TagHandle local_tag) 
{
  std::set<iBase_EntitySetHandle>::const_iterator src;
  for (src = cesets.begin(); src != cesets.end(); ++src)
    process_ce_subsets(imeshImpl, *src, *src, cesets, local_tag);
}

void tag_copy_sets(iMesh_Instance imeshImpl, iBase_TagHandle copyTag,
                  const std::set<iBase_EntitySetHandle> &copySets,
                  iBase_TagHandle tag, const char *tag_val)
{
  int err;

  int tag_size;
  iMesh_getTagSizeBytes(imeshImpl, tag, &tag_size, &err);
  check_error(imeshImpl, err);

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
    if (err == iBase_TAG_NOT_FOUND)
      continue;
    check_error(imeshImpl, err);

    // compare to tag value if necessary
    if (tag_val && strncmp(tag_val, value_ptr, tag_size))
      continue;
      
    // if we got here, we should set the tag on the copy; get the copy
    iBase_EntitySetHandle copy_set;
    iMesh_getEntSetEHData(imeshImpl, *set, copyTag, 
                          reinterpret_cast<iBase_EntityHandle*>(&copy_set), 
                          &err);
    if (err == iBase_TAG_NOT_FOUND)
      continue;
    check_error(imeshImpl, err);

    if (copy_set != *set) {
      iMesh_setEntSetData(imeshImpl, copy_set, tag, value_ptr, tag_size, &err);
      check_error(imeshImpl, err);
    }
  }
}
