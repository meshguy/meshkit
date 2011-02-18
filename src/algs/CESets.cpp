#include "meshkit/CESets.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iMesh_extensions.h>

#include <meshkit/Error.hpp>
#include "SimpleArray.hpp"

namespace MeshKit {

CESets::~CESets()
{
  std::vector<tag_data>::iterator i;
  for (i = tags_.begin(); i != tags_.end(); ++i)
    free(i->value);
}

void CESets::add_tag(iMesh::TagHandle tag_handle, const char *value)
{
  assert(tag_handle != NULL);
  char *tmp = NULL;

  if (value) {
    int tag_size;
    IBERRCHK(mesh_->getTagSizeBytes(tag_handle, tag_size), *mesh_);

    tmp = static_cast<char*>(malloc(tag_size));
    memcpy(tmp, value, tag_size);
  }

  tags_.push_back(tag_data(tag_handle, tmp));
}

void CESets::add_tag(const std::string &tag_name, const char *value)
{
  iMesh::TagHandle tag_handle;
  IBERRCHK(mesh_->getTagHandle(tag_name.c_str(), tag_handle), *mesh_);
  add_tag(tag_handle, value);
}

void CESets::update_tagged_sets()
{
  int err;
  iMesh::EntitySetHandle root_set = mesh_->getRootSet();

  std::vector<tag_data>::const_iterator tag;
  for (tag = tags_.begin(); tag != tags_.end(); ++tag) {
    SimpleArray<iMesh::EntitySetHandle> tmp_sets;
    iMesh_getEntSetsByTagsRec(mesh_->instance(), root_set, &tag->tag,
                              (tag->value ? &tag->value : NULL), 1, false,
                              ARRAY_INOUT(tmp_sets), &err);
    IBERRCHK(err, *mesh_);
    sets_.insert(tmp_sets.begin(), tmp_sets.end());
  }
}

void link_expand_sets(const CESets &ce_sets, iMesh::TagHandle local_tag)
{
  CESets::set_iterator set;
  for (set = ce_sets.sbegin(); set != ce_sets.send(); ++set) {
    IBERRCHK(ce_sets.imesh_instance()->setEntSetEHData(*set, local_tag,
      reinterpret_cast<iMesh::EntityHandle>(*set)), *ce_sets.imesh_instance());
  }
}

/**\brief Get the entities copied from a set
 *
 * Given a set and a tag, return the values of the tag set on entities from that
 * set (i.e. the copied entities).
 * \param mesh the iMesh instance handle
 * \param set the set containing the entities we want
 * \param local_tag the tag relating source and target entities
 * \param ents a vector which will hold our copied entities
 */
static
void get_copied_ents(iMesh *mesh, iMesh::EntitySetHandle set,
                     iMesh::TagHandle local_tag,
                     std::vector<iMesh::EntityHandle> &ents)
{
  std::vector<iMesh::EntityHandle> tmp_ents;
  IBERRCHK(mesh->getEntities(set, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, 
                             tmp_ents), *mesh);

  ents.reserve(tmp_ents.size());

  // get copied entities
  for (size_t i = 0; i < tmp_ents.size(); i++) {
    iMesh::EntityHandle eh_tag;
    iMesh::Error err = mesh->getEHData(tmp_ents[i], local_tag, eh_tag);
    if (err == iBase_SUCCESS && eh_tag)
      ents.push_back(eh_tag);
  }
}

/**\brief Create a copy set if one doesn't exist yet
 *
 * Given a source set and a tag, create a destination (copy) set unless one
 * already exists.
 * \param mesh the iMesh instance handle
 * \param local_tag the tag relating source and target sets
 * \param src the source set
 * \param dest the destination set
 */
static
void get_dest_set(iMesh *mesh, iMesh::TagHandle local_tag,
                  iMesh::EntitySetHandle src, iMesh::EntitySetHandle &dest)
{
  iMesh::Error err;

  if (dest == NULL) return;

  err = mesh->getEntSetEHData(src, local_tag,
                              reinterpret_cast<iMesh::EntityHandle&>(dest));

  if (err != iBase_SUCCESS) {
    IBERRCHK(mesh->createEntSet(false, dest), *mesh);
    IBERRCHK(mesh->setEntSetEHData(src, local_tag,
                                   reinterpret_cast<iMesh::EntityHandle>(dest)),
             *mesh);
  }
}

/**\brief Add copied entities/sets recursively
 *
 * Helper function for process_ce_sets. This adds any entities directly
 * contained in the current set to the dest set and then loops through the
 * contained sets and adds the children if they are also CE sets. If not, it
 * adds the entities in the child set and recurses down a level.
 * \param mesh the iMesh instance handle
 * \param src the source set
 * \param current the child set to examine (start with current == src)
 * \param cesets the copy or expand sets to operate on
 * \param local_tag the tag relating source and target sets
 */
static
void process_ce_subsets(iMesh *mesh, iMesh::EntitySetHandle src,
                        iMesh::EntitySetHandle current,
                        const std::set<iMesh::EntitySetHandle> &cesets,
                        iMesh::TagHandle local_tag) 
{
  iMesh::EntitySetHandle dest;

  // First, add entities directly contained in this set.
  std::vector<iMesh::EntityHandle> tmp_tags;
  get_copied_ents(mesh, current, local_tag, tmp_tags);

  if (!tmp_tags.empty()) {
    get_dest_set(mesh, local_tag, src, dest);
    IBERRCHK(mesh->addEntArrToSet(&tmp_tags[0], tmp_tags.size(), dest), *mesh);
  }

  // Next, start looking at children.
  std::vector<iMesh::EntitySetHandle> children;
  IBERRCHK(mesh->getEntSets(current, 0, children), *mesh);

  for (size_t i = 0; i < children.size(); i++) {

    // If this child set is one of our cesets, add just the set...
    if (cesets.find(children[i]) != cesets.end()) {
      get_dest_set(mesh, local_tag, src, dest);
      if (src == dest) continue;

      IBERRCHK(mesh->addEntSet(children[i], dest), *mesh);
    }

    // ... otherwise, add the entities and recurse into the next level of
    // children.
    else {
      process_ce_subsets(mesh, src, children[i], cesets, local_tag);
    }
  }
}

// xxx - still need to do unique tags
// - for each unique tag
//   . if this set has the tag, add/append to get unique value
void process_ce_sets(iMesh *mesh,
                     const std::set<iMesh::EntitySetHandle> &cesets,
                     iMesh::TagHandle local_tag) 
{
  std::set<iMesh::EntitySetHandle>::const_iterator src;
  for (src = cesets.begin(); src != cesets.end(); ++src)
    process_ce_subsets(mesh, *src, *src, cesets, local_tag);
}

void tag_copy_sets(iMesh *mesh, iMesh::TagHandle copyTag,
                   const std::set<iMesh::EntitySetHandle> &copySets,
                   iMesh::TagHandle tag, const char *tag_val)
{
  iMesh::Error err;

  int tag_size;
  IBERRCHK(mesh->getTagSizeBytes(tag, tag_size), *mesh);

  // allocate temp space for tag value
  std::vector<char> value;
  value.resize(tag_size);
  char *value_ptr = &value[0];
  int value_alloc = tag_size, value_size;
  
  // for each orig copy set with this tag, copy it to its copy
  std::set<iMesh::EntitySetHandle>::iterator set;
  for (set = copySets.begin(); set != copySets.end(); ++set) {
    // get the tag value
    err = mesh->getEntSetData(*set, tag, &value[0]);
    if (err == iBase_TAG_NOT_FOUND)
      continue;
    IBERRCHK(err, *mesh);

    // compare to tag value if necessary
    if (tag_val && strncmp(tag_val, value_ptr, tag_size))
      continue;
      
    // if we got here, we should set the tag on the copy; get the copy
    iMesh::EntitySetHandle copy_set;
    err = mesh->getEntSetEHData(
      *set, copyTag, reinterpret_cast<iMesh::EntityHandle&>(copy_set));
    if (err == iBase_TAG_NOT_FOUND)
      continue;
    IBERRCHK(err, *mesh);

    if (copy_set != *set)
      IBERRCHK(mesh->setEntSetData(copy_set, tag, value_ptr), *mesh);
  }
}

void tag_copy_sets(const CESets &ce_sets, iMesh::TagHandle local_tag,
                   iMesh::TagHandle copy_tag)
{
  // set the copy tag on all copied sets
  for (CESets::const_set_iterator set = ce_sets.sbegin(); set != ce_sets.send();
       ++set) {
    iMesh::Error err;
    iMesh::EntityHandle eh;
    err = ce_sets.imesh_instance()->getEntSetEHData(*set, local_tag, eh);
    if (err == iBase_SUCCESS) {
      IBERRCHK(ce_sets.imesh_instance()->setEntSetEHData(*set, copy_tag, eh),
               *ce_sets.imesh_instance());
    }
  }

  // tag the newly-created sets
  for (CESets::const_tag_iterator tag = ce_sets.tbegin();
       tag != ce_sets.tend(); ++tag) {
    tag_copy_sets(ce_sets.imesh_instance(), copy_tag, ce_sets.sets(),
                  tag->tag, tag->value);
  }
}

} // namespace MeshKit
