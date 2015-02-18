#ifndef MESHKIT_CD_SETS_HPP
#define MESHKIT_CD_SETS_HPP

#include "meshkitalgs_export.h"

#include <set>
#include <string>
#include <vector>

#include "meshkit/MKCore.hpp"

#include <iBase.h>
#include <iMesh.h>

namespace MeshKit {

class MESHKITALGS_EXPORT CESets
{
public:
  struct tag_data
  {
    tag_data(iMesh::TagHandle tag, char *value)
      : tag(tag), value(value)
    {}

    iMesh::TagHandle tag;
    char *value;
  };

  typedef std::vector<tag_data> tag_type;
  typedef std::set<iMesh::EntitySetHandle> set_type;
  typedef tag_type::iterator       tag_iterator;
  typedef tag_type::const_iterator const_tag_iterator;
  typedef set_type::iterator       set_iterator;
  typedef set_type::const_iterator const_set_iterator;

  CESets(MKCore *mkcore) : mesh_(mkcore->imesh_instance())
  {}

  ~CESets();

  iMesh * imesh_instance() const { return mesh_; }

  void add_set(iMesh::EntitySetHandle set)
  {
    sets_.insert(set);
  }

  template <typename T>
  void add_sets(T begin, T end)
  {
    sets_.insert(begin, end);
  }

  void add_tag(iMesh::TagHandle tag_handle, const char *value = NULL);
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
  iMesh *mesh_;
  tag_type tags_;
  set_type sets_;
};

/**\brief Set the target sets for expand sets to be themselves
 */
void MESHKITALGS_EXPORT link_expand_sets(const CESets &ce_sets, iMesh::TagHandle local_tag);

/**\brief Add newly-created entities/sets to a collection of sets
 *
 * Given a collection of copy, expand, or extrude source sets and a tag, create 
 * a destination (copy) set unless one already exists. Fill this set with any
 * new entities/sets created from those in the source set.
 * \param mesh the iMesh instance handle
 * \param cesets a collection of source sets
 * \param local_tag the tag relating source and target entities/sets
 */
void MESHKITALGS_EXPORT process_ce_sets(iMesh *mesh,
                     const std::set<iMesh::EntitySetHandle> &cesets,
                     iMesh::TagHandle local_tag);

/**\brief Tag a collection of copied sets
 *
 * Given a collection of source sets and a tag |copyTag| relating sources to
 * destinations, apply a tag |tag| to the destination sets if the tag exists on
 * the corresponding source.
 * \param mesh the iMesh instance handle
 * \param copyTag the tag relating sources and destinations
 * \param copySets a collection of source sets
 * \param tag the tag to set on the destinations
 * \param tag_val if non-NULL, only set |tag| on the destination if the source's
 *                tag matches this value
 */
void MESHKITALGS_EXPORT tag_copy_sets(iMesh *mesh, iMesh::TagHandle copyTag,
                   const std::set<iMesh::EntitySetHandle> &copySets,
                   iMesh::TagHandle tag, const char *tag_val);

/**\brief Tag a collection of copied sets
 *
 * Given a CESets instance and a tag |local_tag| relating sources to
 * destinations, move the contents of the local tag to |copy_tag| and 
 * apply the tags in CESets to the destination sets if the tag exists on the
 * corresponding source.
 * \param ce_sets the CESets to operate on
 * \param local_tag the local tag relating copied sets with their originals
 * \param copy_tag the copy tag to receive local_tag's data
 */
void MESHKITALGS_EXPORT tag_copy_sets(const CESets &ce_sets, iMesh::TagHandle local_tag,
                   iMesh::TagHandle copy_tag);

} // namespace MeshKit
#endif
