#include "meshkit/MergeMesh.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

#include "moab/Skinner.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"

namespace MeshKit
{
  // static registration of this  mesh scheme
  moab::EntityType MergeMesh_tps[] = { moab::MBVERTEX,
				       moab::MBEDGE,
				       moab::MBTRI,
				       moab::MBHEX,
				       moab::MBMAXTYPE};

  const moab::EntityType* MergeMesh::output_types()
  { return MergeMesh_tps; }
  
  MergeMesh::MergeMesh(MKCore *mkcore, const MEntVector &me_vec)
    : MeshScheme(mkcore, me_vec),
      imeshImpl(mkcore->imesh_instance()), mbImpl (mkcore->moab_instance())
  {
    mergeTol = 1e-4;
    updateSets = 0;
    doMerge = 1;
    mergeTag = NULL;
  }

  MergeMesh::~MergeMesh()
  {
    if (mbMergeTag) mbImpl->tag_delete(mbMergeTag);
  }

  bool MergeMesh::add_modelent(ModelEnt *model_ent)
  {
    return MeshOp::add_modelent(model_ent);
  }

  void MergeMesh::set_merge_params(double merge_tol, int updatesets,
				   int domerge, iBase_TagHandle merge_tag)
  {
    mergeTol = merge_tol;
    updateSets = updatesets;
    doMerge = domerge;
    mergeTag = merge_tag;
  }

  void MergeMesh::setup_this()
  {
    std::cout << "IN SETUPTHIS MERGEMESH " << std::endl;

  }

  void MergeMesh::execute_this()
  {
    std::cout << "IN EXTHIS MERGEMESH " << std::endl;

    std::vector<iBase_EntityHandle> ents(mentSelection.size());
    int dim = 3;
    for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {
      ModelEnt *me = mit->first;
      dim  = me->dimension();
    }

    if(dim == 3){
      imeshImpl->getEntities(0, iBase_REGION, iMesh_ALL_TOPOLOGIES,
			     ents);
    }
    else if (dim == 2){
      imeshImpl->getEntities(0, iBase_FACE, iMesh_ALL_TOPOLOGIES,
			     ents);
    }
    merge_entities(&ents[0], ents.size(), mergeTol, doMerge, updateSets,
		   mergeTag);
  }

  void MergeMesh::merge_entities(iBase_EntityHandle *elems,
				 int elems_size,
				 const double merge_tol,
				 const int do_merge,
				 const int update_sets,
				 iBase_TagHandle merge_tag)  
  {
    mergeTol = merge_tol;
    mergeTolSq = merge_tol*merge_tol;
  
    moab::Range tmp_elems;
    tmp_elems.insert((moab::EntityHandle*)elems, (moab::EntityHandle*)elems + elems_size);
    moab::ErrorCode result = merge_entities(tmp_elems, do_merge, update_sets,
					(moab::Tag)merge_tag);
    if (moab::MB_SUCCESS != result)
      exit(0);
  }

  void MergeMesh::perform_merge(iBase_TagHandle merge_tag) 
  {
    // put into a range
    moab::ErrorCode result = perform_merge((moab::Tag) merge_tag);
    if (moab::MB_SUCCESS != result)
      exit(0);
  }

  moab::ErrorCode MergeMesh::merge_entities(moab::Range &elems,
					const int do_merge,
					const int update_sets,
					moab::Tag merge_tag) 
  {
    // get the skin of the entities
    moab::Skinner skinner(mbImpl);
    moab::Range skin_range;
    moab::ErrorCode result = skinner.find_skin(0, elems, 0, skin_range);
    if (moab::MB_SUCCESS != result) return result;

    // create a tag to mark merged-to entity; reuse tree_root
    moab::EntityHandle tree_root = 0;
    if (0 == merge_tag) {
      //result = mbImpl->tag_create("__merge_tag", sizeof(moab::EntityHandle), 
      //                            MB_TAG_DENSE, MB_TYPE_HANDLE,
      //                            mbMergeTag, &tree_root);
      result = mbImpl->tag_get_handle("__merge_tag", 1, moab::MB_TYPE_HANDLE,
                                      mbMergeTag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT,
				      &tree_root);
      if (moab::MB_SUCCESS != result) return result;
    }
    else mbMergeTag = merge_tag;
  
    // build a kd tree with the vertices
    moab::AdaptiveKDTree kd(mbImpl);
    result = kd.build_tree(skin_range, &tree_root);
    if (moab::MB_SUCCESS != result) return result;

    // find matching vertices, mark them
    result = find_merged_to(kd, tree_root, mbMergeTag);
    if (moab::MB_SUCCESS != result) return result;
  
    // merge them if requested
    if (do_merge) {
      result = perform_merge(mbMergeTag);
      if (moab::MB_SUCCESS != result) return result;
    }
  
    return moab::MB_SUCCESS;
  }

  moab::ErrorCode MergeMesh::perform_merge(moab::Tag merge_tag) 
  {
    if (deadEnts.size()==0){
      std::cout << "\nWarning: Geometries don't have a common face; Nothing to merge" << std::endl;
      return moab::MB_SUCCESS; //nothing to merge carry on with the program
    }
    if (mbImpl->type_from_handle(*deadEnts.rbegin()) != moab::MBVERTEX)
      return moab::MB_FAILURE;
  
    std::vector<moab::EntityHandle> merge_tag_val(deadEnts.size());
    moab::ErrorCode result = mbImpl->tag_get_data(merge_tag, deadEnts, &merge_tag_val[0]);
    if (moab::MB_SUCCESS != result) return result;
  
    moab::Range::iterator rit;
    unsigned int i;
    for (rit = deadEnts.begin(), i = 0; rit != deadEnts.end(); rit++, i++) {
      assert(merge_tag_val[i]);
      result = mbImpl->merge_entities(merge_tag_val[i], *rit, false, false);
      if (moab::MB_SUCCESS != result) return result;
    }
  
    return mbImpl->delete_entities(deadEnts);
  }

  moab::ErrorCode MergeMesh::find_merged_to(moab::AdaptiveKDTree & tree, moab::EntityHandle &tree_root, moab::Tag merge_tag)
  {
    moab::AdaptiveKDTreeIter iter;
    //moab::AdaptiveKDTree tree(mbImpl);
  
    // evaluate vertices in this leaf
    moab::Range leaf_range, leaf_range2;
    std::vector<moab::EntityHandle> sorted_leaves;
    std::vector<double> coords;
    std::vector<moab::EntityHandle> merge_tag_val, leaves_out;
  
    moab::ErrorCode result = tree.get_tree_iterator(tree_root, iter);
    if (moab::MB_SUCCESS != result) return result;
    while (result == moab::MB_SUCCESS) {
      sorted_leaves.push_back( iter.handle() );
      result = iter.step();
    }
    if (result != moab::MB_ENTITY_NOT_FOUND)
      return result;
    std::sort( sorted_leaves.begin(), sorted_leaves.end() );
  
    std::vector<moab::EntityHandle>::iterator it;
    for (it = sorted_leaves.begin(); it != sorted_leaves.end(); ++it) {

      leaf_range.clear();
      result = mbImpl->get_entities_by_handle(*it, leaf_range);
      if (moab::MB_SUCCESS != result) return result;
      coords.resize(3*leaf_range.size());
      merge_tag_val.resize(leaf_range.size());
      result = mbImpl->get_coords(leaf_range, &coords[0]);
      if (moab::MB_SUCCESS != result) return result;
      result = mbImpl->tag_get_data(merge_tag, leaf_range, &merge_tag_val[0]);
      if (moab::MB_SUCCESS != result) return result;
      moab::Range::iterator rit;
      unsigned int i;
      bool inleaf_merged, outleaf_merged = false;
      unsigned int lr_size = leaf_range.size();
    
      for (i = 0, rit = leaf_range.begin(); i != lr_size; rit++, i++) {
	if (0 != merge_tag_val[i]) continue;
	moab::CartVect from(&coords[3*i]);
	inleaf_merged = false;

	// check close-by leaves too
	leaves_out.clear();
	result = tree.distance_search( from.array(), mergeTol,
					     leaves_out);
	leaf_range2.clear();
	for (std::vector<moab::EntityHandle>::iterator vit = leaves_out.begin();
	     vit != leaves_out.end(); vit++) {
	  if (*vit > *it) { // if we haven't visited this leaf yet in the outer loop
	    result = mbImpl->get_entities_by_handle(*vit, leaf_range2, moab::Interface::UNION);
	    if (moab::MB_SUCCESS != result) return result;
	  }
	}
	if (!leaf_range2.empty()) {
	  coords.resize(3*(lr_size+leaf_range2.size()));
	  merge_tag_val.resize(lr_size+leaf_range2.size());
	  result = mbImpl->get_coords(leaf_range2, &coords[3*lr_size]);
	  if (moab::MB_SUCCESS != result) return result;
	  result = mbImpl->tag_get_data(merge_tag, leaf_range2, &merge_tag_val[lr_size]);
	  if (moab::MB_SUCCESS != result) return result;
	  outleaf_merged = false;
	}

	// check other verts in this leaf
	for (unsigned int j = i+1; j < merge_tag_val.size(); j++) {
	  moab::EntityHandle to_ent = j >= lr_size ? leaf_range2[j-lr_size] : 
	    leaf_range[j];
        
	  if (*rit == to_ent) continue;
        
          if ((from - moab::CartVect(&coords[3*j])).length_squared() < mergeTolSq) {
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
	  if (moab::MB_SUCCESS != result) return result;
	  outleaf_merged = false;
	}
	if (inleaf_merged) {
	  result = mbImpl->tag_set_data(merge_tag, leaf_range, &merge_tag_val[0]);
	  if (moab::MB_SUCCESS != result) return result;
	}

      }
    }
    return moab::MB_SUCCESS;
  }
}
