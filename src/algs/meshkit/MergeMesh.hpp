#ifndef MESHKIT_MERGEMESH_HPP
#define MESHKIT_MERGEMESH_HPP
#ifdef HAVE_MOAB

//todo:
//- nice case for when we don't have MOAB

#include "iMesh.h"

#include "meshkit/MeshScheme.hpp"

#include "MBiMesh.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"


namespace MeshKit {


class MergeMesh : public MeshScheme
{
private:

  MBInterface *mbImpl;
  double mergeTol;
  bool do_merge;
  bool update_sets;
  MBTag merge_tag; // tag pointing to the entity to which an entity will be merged

  MBRange deadEnts; // entities which will go away after the merge

  //- given a kdtree, set tag on vertices in leaf nodes with vertices
  //- to which they should be merged
  MBErrorCode find_merged_to(MBEntityHandle &tree_root);  
  MBErrorCode merge_entities(MBRange &elems);
  
  //- perform the actual merge
  MBErrorCode perform_merge();

public:
  MergeMesh(MKCore *mkcore, const MEntVector &me_vec);
  virtual ~MergeMesh();

  /* \brief Get class name */
  static const char* name();
  const moab::EntityType* output_types();
  void setup_this();
  void execute_this();

  void set_merge_tol(double v){mergeTol = v;};
  void set_do_merge(bool v){do_merge = v;};
  void set_update_sets(bool v){update_sets = v;};
  void set_merge_tag(MBTag v){merge_tag = v;};
};


}//namespace MeshKit
#endif //HAVE_MOAB
#endif //MESHKIT_MERGEMESH_HPP

