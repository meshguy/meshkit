#ifndef METADATA_HPP
#define METADATA_HPP

#include <set>
#include <string>
#include <vector>

#include <iBase.h>
#include <iMesh.h>

namespace MeshKit {

class CESets
{
public:
  struct tag_data
  {
    tag_data(iBase_TagHandle tag, char *value)
      : tag(tag), value(value)
    {}

    iBase_TagHandle tag;
    char *value;
  };

  typedef std::vector<tag_data> tag_type;
  typedef std::set<iBase_EntitySetHandle> set_type;
  typedef tag_type::iterator       tag_iterator;
  typedef tag_type::const_iterator const_tag_iterator;
  typedef set_type::iterator       set_iterator;
  typedef set_type::const_iterator const_set_iterator;

  CESets(iMesh_Instance impl) : impl_(impl)
  {}

  ~CESets();

  iMesh_Instance impl() const { return impl_; }

  void add_set(iBase_EntitySetHandle set)
  {
    sets_.insert(set);
  }

  template <typename T>
  void add_sets(T begin, T end)
  {
    sets_.insert(begin, end);
  }

  void add_tag(iBase_TagHandle tag_handle, const char *value = NULL);
  void add_tag(const std::string &tag_name, const char *value = NULL);
  void update_tagged_sets();

  void clear()
  {
    sets_.clear();
  }

  tag_type & tags()             { return tags_; }
  const tag_type & tags() const { return tags_; }
  set_type & sets()             { return sets_; }
  const set_type & sets() const { return sets_; }

  tag_iterator tbegin()             { return tags_.begin(); }
  const_tag_iterator tbegin() const { return tags_.begin(); }
  tag_iterator tend()               { return tags_.end(); }
  const_tag_iterator tend() const   { return tags_.end(); }
  set_iterator sbegin()             { return sets_.begin(); }
  const_set_iterator sbegin() const { return sets_.begin(); }
  set_iterator send()               { return sets_.end(); }
  const_set_iterator send() const   { return sets_.end(); }

private:
  iMesh_Instance impl_;
  tag_type tags_;
  set_type sets_;
};

/**\brief Set the target sets for expand sets to be themselves
 */
void link_expand_sets(const CESets &ce_sets, iBase_TagHandle local_tag);

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

void tag_copy_sets(const CESets &ce_sets, iBase_TagHandle local_tag,
                   iBase_TagHandle copy_tag);

} // namespace MeshKit
#endif
