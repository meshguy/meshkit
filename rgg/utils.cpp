#include "utils.hpp"
#include "SimpleArray.hpp"
#include <cstring>

static int get_ce_sets(CESets &ce_sets, iBase_EntitySetHandle orig_set,
                       const char **tag_names, const char **tag_vals,
                       int num_tags)
{
  for (int i = 0; i < num_tags; i++) {
    int err;
    iBase_TagHandle tmp_tag;
    iMesh_getTagHandle(ce_sets.impl(), tag_names[i], &tmp_tag, &err,
                       strlen(tag_names[i]));
    ERRORR("Failed to get tag handle.", iBase_FAILURE);
    SimpleArray<iBase_EntitySetHandle> tmp_sets;
    iMesh_getEntSetsByTagsRec(ce_sets.impl(), orig_set, &tmp_tag, &tag_vals[i],
                              1, 0, ARRAY_INOUT(tmp_sets), &err);
    ERRORR("Failure getting sets by tags.", err);
  
    if (tmp_sets.size() != 0)
      ce_sets.add_sets(tmp_sets.begin(), tmp_sets.end());
  }

  return iBase_SUCCESS;
}

int get_copy_sets(CopyMesh *cm, iBase_EntitySetHandle orig_set,
                  const char **tag_names, const char **tag_vals,
                  int num_tags)
{
  return get_ce_sets(cm->copy_sets(), orig_set, tag_names, tag_vals, num_tags);
}

int get_expand_sets(CopyMesh *cm, iBase_EntitySetHandle orig_set,
                    const char **tag_names, const char **tag_vals,
                    int num_tags)
{
  return get_ce_sets(cm->expand_sets(), orig_set, tag_names, tag_vals,
                     num_tags);
}

int extend_expand_sets(CopyMesh *cm) 
{
    // check expand sets for any contained sets which aren't already copy sets, 
    // and add them to the list
  std::set<iBase_EntitySetHandle>::iterator sit;
  for (sit = cm->expand_sets().sets().begin();
       sit != cm->expand_sets().sets().end(); ++sit) {
    int err;
    SimpleArray<iBase_EntitySetHandle> sets;
    iMesh_getEntSets(cm->impl(), *sit, 1, ARRAY_INOUT(sets), &err);
    ERRORR("Failed to get contained sets.", err);

    if (sets.size() != 0)
      cm->copy_sets().add_sets(sets.begin(), sets.end());
  }
  
  return iBase_SUCCESS;
}
