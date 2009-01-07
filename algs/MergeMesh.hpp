#ifndef MERGEMESH_HPP
#define MERGEMESH_HPP

#include "iMesh.h"

#ifdef MOAB
#include "MBInterface.hpp"
#include "MBRange.hpp"
#endif

#include <string>
#include <vector>
#include <assert.h>
#include <iostream>

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

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
  int merge_entities(iBase_EntityHandle *elems,
                     int elems_size,
                     const double merge_tol,
                     const int do_merge = true,
                     const int update_sets = false,
                     iBase_TagHandle merge_tag = 0);
  
    /* \brief Perform the actual merge between entities
     */
  int perform_merge(iBase_TagHandle merged_to);
  
  
private:

  iMesh_Instance imeshImpl;

  double mergeTol, mergeTolSq;

  iBase_TagHandle mergeTag;
  
#ifdef MOAB

    //- given a kdtree, set tag on vertices in leaf nodes with vertices
    //- to which they should be merged
  MBErrorCode find_merged_to(MBEntityHandle &tree_root, MBTag merged_to);
  
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
#endif
  
};

inline MergeMesh::MergeMesh(iMesh_Instance impl) 
        : imeshImpl(impl)
{
#ifdef MOAB
  mbImpl = (MBInterface*) impl;
#endif
}

inline MergeMesh::~MergeMesh() 
{
#ifdef MOAB
  if (mbMergeTag) mbImpl->tag_delete(mbMergeTag);
#endif
}

#endif

