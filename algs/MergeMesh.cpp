#include "MergeMesh.hpp"

#ifdef MOAB
#include "MBSkinner.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"
#endif

int MergeMesh::merge_entities(iBase_EntityHandle *elems,
                              int elems_size,
                              const double merge_tol,
                              const int do_merge,
                              const int update_sets,
                              iBase_TagHandle merge_tag)  
{
  mergeTol = merge_tol;
  mergeTolSq = merge_tol*merge_tol;
  
#ifdef MOAB  
  MBRange tmp_elems;
  MBEntityHandle *tmp_elems_ptr = (MBEntityHandle*) elems;
  std::copy(tmp_elems_ptr, tmp_elems_ptr+elems_size, mb_range_inserter(tmp_elems));
  MBErrorCode result = merge_entities(tmp_elems, do_merge, update_sets, (MBTag) merge_tag);
  return (MB_SUCCESS == result ? iBase_SUCCESS : iBase_FAILURE);
#else
  return iBase_NOT_SUPPORTED;
#endif

}

int MergeMesh::perform_merge(iBase_TagHandle merge_tag) 
{
#ifdef MOAB
    // put into a range
  return perform_merge((MBTag) merge_tag);
#else
  return iBase_NOT_SUPPORTED;
#endif
}



#ifdef MOAB
MBErrorCode MergeMesh::merge_entities(MBRange &elems,
                                      const int do_merge,
                                      const int update_sets,
                                      MBTag merge_tag) 
{
    // get the skin of the entities
  MBSkinner skinner(mbImpl);
  MBRange skin_range;
  MBErrorCode result = skinner.find_skin(elems, 0, skin_range);
  if (MB_SUCCESS != result) return result;

    // create a tag to mark merged-to entity; reuse tree_root
  MBEntityHandle tree_root = 0;
  if (0 == merge_tag) {
    result = mbImpl->tag_create("__merge_tag", sizeof(MBEntityHandle), 
                                MB_TAG_DENSE, MB_TYPE_HANDLE,
                                mbMergeTag, &tree_root);
    if (MB_SUCCESS != result) return result;
  }
  else mbMergeTag = merge_tag;
  
    // build a kd tree with the vertices
  MBAdaptiveKDTree kd(mbImpl, true);
  result = kd.build_tree(skin_range, tree_root);
  if (MB_SUCCESS != result) return result;

    // find matching vertices, mark them
  result = find_merged_to(tree_root, mbMergeTag);
  if (MB_SUCCESS != result) return result;
  
    // merge them if requested
  if (do_merge) {
    result = perform_merge(mbMergeTag);
    if (MB_SUCCESS != result) return result;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MergeMesh::perform_merge(MBTag merge_tag) 
{
  if (mbImpl->type_from_handle(*deadEnts.rbegin()) != MBVERTEX) 
    return MB_FAILURE;
  
  std::vector<MBEntityHandle> merge_tag_val(deadEnts.size());
  MBErrorCode result = mbImpl->tag_get_data(merge_tag, deadEnts, &merge_tag_val[0]);
  if (MB_SUCCESS != result) return result;
  
  MBRange::iterator rit;
  unsigned int i;
  for (rit = deadEnts.begin(), i = 0; rit != deadEnts.end(); rit++, i++) {
    assert(merge_tag_val[i]);
    result = mbImpl->merge_entities(merge_tag_val[i], *rit, false, false);
    if (MB_SUCCESS != result) return result;
  }
  
  return mbImpl->delete_entities(deadEnts);
}

MBErrorCode MergeMesh::find_merged_to(MBEntityHandle &tree_root, MBTag merge_tag) 
{
  MBAdaptiveKDTreeIter iter;
  MBAdaptiveKDTree tree(mbImpl);
  
  MBErrorCode result = tree.get_tree_iterator(tree_root, iter);
  if (MB_SUCCESS != result) return result;
  
    // evaluate vertices in this leaf
  MBRange leaf_range, leaf_range2, leaves_checked;
  std::vector<double> coords;
  std::vector<MBEntityHandle> merge_tag_val, leaves_out;

  while (result != MB_ENTITY_NOT_FOUND) {
    leaves_checked.insert(iter.handle());
    result = mbImpl->get_entities_by_handle(iter.handle(), leaf_range);
    if (MB_SUCCESS != result) return result;
    coords.resize(3*leaf_range.size());
    merge_tag_val.resize(leaf_range.size());
    result = mbImpl->get_coords(leaf_range, &coords[0]);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->tag_get_data(merge_tag, leaf_range, &merge_tag_val[0]);
    if (MB_SUCCESS != result) return result;
    MBRange::iterator rit;
    unsigned int i;
    bool inleaf_merged, outleaf_merged = false;
    unsigned int lr_size = leaf_range.size();
    
    for (i = 0, rit = leaf_range.begin(); i != lr_size; rit++, i++) {
      if (0 != merge_tag_val[i]) continue;
      MBCartVect from(&coords[3*i]);
      inleaf_merged = false;

        // check close-by leaves too
      leaves_out.clear();
      result = tree.leaves_within_distance(tree_root, from.array(), mergeTol,
                                           leaves_out);
      leaf_range2.clear();
      for (std::vector<MBEntityHandle>::iterator vit = leaves_out.begin();
           vit != leaves_out.end(); vit++) {
        if (leaves_checked.find(*vit) == leaves_checked.end()) {
          result = mbImpl->get_entities_by_handle(*vit, leaf_range2, MBInterface::UNION);
          if (MB_SUCCESS != result) return result;
        }
      }
      if (!leaf_range2.empty()) {
        coords.resize(3*(lr_size+leaf_range2.size()));
        merge_tag_val.resize(lr_size+leaf_range2.size());
        result = mbImpl->get_coords(leaf_range2, &coords[3*lr_size]);
        if (MB_SUCCESS != result) return result;
        result = mbImpl->tag_get_data(merge_tag, leaf_range2, &merge_tag_val[lr_size]);
        if (MB_SUCCESS != result) return result;
        outleaf_merged = false;
      }

        // check other verts in this leaf
      for (unsigned int j = i+1; j < merge_tag_val.size(); j++) {
        MBEntityHandle to_ent = j >= lr_size ? leaf_range2[j-lr_size] : 
            leaf_range[j];
        
        if (*rit == to_ent) continue;
        
        if ((from - MBCartVect(&coords[3*j])).length_squared() < mergeTolSq) {
          merge_tag_val[j] = *rit;
          if (j < lr_size) inleaf_merged = true;
          else outleaf_merged = true;
          deadEnts.insert(to_ent);
        }
      }

      if (outleaf_merged) {
        result = mbImpl->tag_set_data(merge_tag, leaf_range2, &merge_tag_val[leaf_range.size()]);
        if (MB_SUCCESS != result) return result;
      }
    }
    
    if (inleaf_merged) {
      result = mbImpl->tag_set_data(merge_tag, leaf_range, &merge_tag_val[0]);
      if (MB_SUCCESS != result) return result;
    }

    result = iter.step();
  }

  return MB_SUCCESS;
}

#endif // ifdef MOAB
