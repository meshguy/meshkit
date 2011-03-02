#define HAVE_MOAB 1

#include "meshkit/MKCore.hpp"
#include "meshkit/MergeMesh.hpp"

#ifdef HAVE_MOAB

//KDTree required
#include "MBSkinner.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBRange.hpp"
#include "MBCartVect.hpp"

#include <algorithm>
#include <vector>

namespace MeshKit {

MergeMesh::MergeMesh(MKCore *mkcore, const MEntVector &me_vec)
  : MeshScheme(mkcore, me_vec),

    mbImpl( (reinterpret_cast<MBiMesh *>(this->mk_core()->imesh_instance()) )->mbImpl),
    mergeTol(1e-3),
    do_merge(true),
    update_sets(false),
    merge_tag(0)
{}

MergeMesh::~MergeMesh() 
{
  if (merge_tag) 
    mbImpl->tag_delete(merge_tag); //not the behaviour we want when merge_tag is passed in, fixme
}
  
// static registration of this  mesh scheme                                                                 
moab::EntityType MergeMesh_tps[] = { moab::MBVERTEX,
				    moab::MBEDGE,
				    moab::MBTRI,
				    moab::MBHEX,
				    moab::MBMAXTYPE};

const moab::EntityType* MergeMesh::output_types(){ 
  return MergeMesh_tps; 
}

const char* MergeMesh::name() {
  return "MergeMesh";
}
  
void MergeMesh::setup_this(){
  //no-op
}

void MergeMesh::execute_this(){
  ModelEnt *me = mentSelection.begin()->first;
  moab::EntityHandle set=me->mesh_handle(); // this is the initial set                                       

  

}


MBErrorCode MergeMesh::merge_entities(MBRange &elems)
{
  //elems is a bunch of 3d elements in 3d, or a bunch of 2d elements in 2d
  //  if (result != MB_SUCCESS)
  //throw Error(MK_FAILURE, "Something");

  // get the skin of the entities
  MBSkinner skinner(mbImpl);
  MBRange skin_range;
  MBErrorCode result = skinner.find_skin(elems, 0, skin_range);
  if (MB_SUCCESS != result) 
    return result;

  // create a tag to mark merged-to entity; reuse tree_root
  MBEntityHandle tree_root = 0;
  if (0 == merge_tag) {
    result = mbImpl->tag_create("__merge_tag", sizeof(MBEntityHandle), 
                                MB_TAG_DENSE, MB_TYPE_HANDLE,
                                merge_tag, &tree_root);
    if (MB_SUCCESS != result) 
      return result;
  }
  
  // build a kd tree with the vertices
  MBAdaptiveKDTree kd(mbImpl, true);
  result = kd.build_tree(skin_range, tree_root);
  if (MB_SUCCESS != result) 
    return result;

  // find matching vertices, mark them
  result = find_merged_to(tree_root);
  if (MB_SUCCESS != result) 
    return result;
  
  // merge them if requested
  if (do_merge) {
    result = perform_merge();
    if (MB_SUCCESS != result) 
      return result;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MergeMesh::perform_merge()
{
  if (deadEnts.size()==0){
    std::cout << "\nWarning: Geometries don't have a common face; Nothing to merge" << std::endl;
    return MB_SUCCESS; //nothing to merge carry on with the program
  }
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

MBErrorCode MergeMesh::find_merged_to(MBEntityHandle &tree_root)
{
  MBAdaptiveKDTreeIter iter;
  MBAdaptiveKDTree tree(mbImpl);

  double mergeTolSq = mergeTol*mergeTol;

  // evaluate vertices in this leaf
  MBRange leaf_range, leaf_range2;
  std::vector<MBEntityHandle> sorted_leaves;
  std::vector<double> coords;
  std::vector<MBEntityHandle> merge_tag_val, leaves_out;
  
  MBErrorCode result = tree.get_tree_iterator(tree_root, iter);
  if (MB_SUCCESS != result) return result;
  while (result == MB_SUCCESS) {
    sorted_leaves.push_back( iter.handle() );
    result = iter.step();
  }
  if (result != MB_ENTITY_NOT_FOUND)
    return result;
  std::sort( sorted_leaves.begin(), sorted_leaves.end() );
  
  std::vector<MBEntityHandle>::iterator it;
  for (it = sorted_leaves.begin(); it != sorted_leaves.end(); ++it) {

    leaf_range.clear();
    result = mbImpl->get_entities_by_handle(*it, leaf_range);
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
        if (*vit > *it) { // if we haven't visited this leaf yet in the outer loop
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
          if (j < lr_size){
	    inleaf_merged = true;}
          else{
	    outleaf_merged = true;}

          deadEnts.insert(to_ent);
        }

      }
      if (outleaf_merged) {
	result = mbImpl->tag_set_data(merge_tag, leaf_range2, &merge_tag_val[leaf_range.size()]);
        if (MB_SUCCESS != result) return result;
	outleaf_merged = false;
      }
      if (inleaf_merged) {
	result = mbImpl->tag_set_data(merge_tag, leaf_range, &merge_tag_val[0]);
	if (MB_SUCCESS != result) return result;
      }

    }
  }
  return MB_SUCCESS;
}

}//namespace MeshKit 
#endif // ifdef MOAB
