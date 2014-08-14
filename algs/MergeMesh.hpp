#ifndef MESHKIT_MERGEMESH_HPP
#define MESHKIT_MERGEMESH_HPP

#include "iMesh.h"

#include "MBiMesh.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp"

class MergeMesh 
{
public:
    /* \brief Constructor
     */
  MergeMesh(iMesh_Instance impl);
  
    /* \brief Destructor
     */
  virtual ~MergeMesh();

    /* \brief Merge vertices in elements passed in
     */
  void merge_entities(iBase_EntityHandle *elems,
                      int elems_size,
                      const double merge_tol,
                      const int do_merge = true,
                      const int update_sets = false,
                      iBase_TagHandle merge_tag = 0);
  
    /* \brief Perform the actual merge between entities
     */
  void perform_merge(iBase_TagHandle merged_to);
private:
  iMesh_Instance imeshImpl;

  double mergeTol, mergeTolSq;

  iBase_TagHandle mergeTag;
  

    //- given a kdtree, set tag on vertices in leaf nodes with vertices
    //- to which they should be merged
  MBErrorCode find_merged_to(MBAdaptiveKDTree & tree, MBEntityHandle &tree_root, MBTag merge_tag);
  
  MBErrorCode merge_entities(MBRange &elems,
                             const int do_merge,
                             const int update_sets,
                             MBTag merge_tag);
  
    //- perform the actual merge
  MBErrorCode perform_merge(MBTag merged_to);

  MBInterface *mbImpl;

    //- the tag pointing to the entity to which an entity will be merged
  MBTag mbMergeTag;

    //- entities which will go away after the merge
  MBRange deadEnts;
  
};

inline MergeMesh::MergeMesh(iMesh_Instance impl) 
        : imeshImpl(impl)
{
  mbImpl = (reinterpret_cast<MBiMesh*> (impl))->mbImpl;
}

inline MergeMesh::~MergeMesh() 
{
  if (mbMergeTag) mbImpl->tag_delete(mbMergeTag);
}

#endif

