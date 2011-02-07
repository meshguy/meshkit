// Author : Rajeev Jain
// Feb 01, 2011
// CopyMesh is an algorithm is used to copy move mesh. 
// It takes care of the meta-data propagation on the final mesh.

#ifndef COPYMESH_HPP
#define COPYMESH_HPP

#include <vector>
#include <map>
#include <sys/resource.h>

#include "meshkit/Types.h"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"
#include "iMesh.hh"

#include "iMesh_extensions.h"
#include "utils/Transform.hpp"
#include "LocalTag.hpp"
#include "utils/CESets.hpp"
#include "MKException.hpp"


namespace MeshKit {

  class MKCore;
    

  /** \class CopyMesh CopyMesh.hpp "meshkit/CopyMesh.hpp"
   * \brief A simple class for meshing copy moving meshes
   *
   * INPUT: ModelEnts representing meshes
   * MESH TYPE(S): ALL TYPES
   * OUTPUT: Copied mesh along with the existing mesh
   * DEPENDENCIES: (none)
   * CONSTRAINTS: (none)
   *
   * This class performs the trivial task of copy moving meshes.  There can be 
   * multiple instances of this class, and therefore it is pointed to and managed by ....
   *
   * Each instance of this class stores all the ModelEnt's representing the mesh data,
   * and after execution after meshing new entities are created and tag propagation happens.
   */
  class CopyMesh : public MeshScheme
  {
  public:

    //! Bare constructor
    CopyMesh(MKCore *mkcore, const MEntVector &me_vec);
  
    //! Destructor
    virtual ~CopyMesh();

    /** \brief Factory function registered with MeshOpFactory
     * \param mkcore MKCore instance for the factory
     * \param me_vec ModelEnts to which this scheme will be applied
     */

    /* \brief Return the imesh instance
     */
    iMesh_Instance impl() const {return imeshImpl;}

    static MeshOp *factory(MKCore *mkcore, const MEntVector &me_vec);
  
    /** \brief Function returning whether this scheme can mesh the specified entity
     * 
     * Used by MeshOpFactory to find scheme for an entity.
     * \param me ModelEnt being queried
     * \return If true, this scheme can mesh the specified ModelEnt
     */
    static bool can_mesh(ModelEnt *me);
  
    /** \brief Return the mesh entity types operated on by this scheme
     * \param tps Entity types generated by this scheme
     */
    virtual void mesh_types(std::vector<moab::EntityType> &tps);
  
    /** \brief Re-implemented here so we can check topological dimension of model_ent
     * \param model_ent ModelEnt being added
     */
    virtual bool add_modelent(ModelEnt *model_ent);

    //! Setup is a no-op, but must be provided since it's pure virtual
    virtual void setup_this();

    //! The only setup/execute function we need, since meshing vertices is trivial
    virtual void execute_this();

    // CopyMesh function local //////////////////
    /* \brief Return the copyTag used to indicate set copies
     */
    iBase_TagHandle copy_tag() {return copyTag;}

    void update_sets();

    /* \brief Reset the copy and expand set lists
     */
    void reset_sets();

    /* \brief Add tag which should have unique values
     */
    void add_unique_tag(const std::string &tag_name);

    /* \brief Add tag which should have unique values
     */
    void add_unique_tag(iBase_TagHandle tag_handle);

    /* \brief Return reference to copy sets
     */
    CESets &copy_sets();

    /* \brief Return reference to expand sets
     */
    CESets &expand_sets();

    /* \brief Return reference to unique sets
     */
    std::set<iBase_EntitySetHandle> &unique_sets();
  
    void copy(iBase_EntitySetHandle set_handle,
	      const copy::Transform &trans = copy::Identity(),
	      iBase_EntityHandle **new_ents = NULL,
	      int *new_ents_allocated = 0,
	      int *new_ents_size = 0,
	      bool do_merge = true);

    void copy(iBase_EntityHandle *ent_handles,
	      int num_ents,
	      const copy::Transform &trans = copy::Identity(),
	      iBase_EntityHandle **new_ents = NULL,
	      int *new_ents_allocated = 0,
	      int *new_ents_size = 0,
	      bool do_merge = true);

    /* \brief Tag copied sets with indicated tag from original set
     */
    void tag_copied_sets(const char **tag_names, const char **tag_vals,
			 const int num_tags);

    /* \brief Tag copied sets with indicated tag from original set
     */
    void tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
			 const int num_tags);

    /* \brief Set the location to copy move the mesh
     */
    void set_location(const double []);
  
  private:
    //- get the copy/expand sets based on copy/expand tags
    void get_copy_expand_sets(iBase_EntitySetHandle *&copy_sets,
			      int &num_copy_sets,
			      iBase_EntitySetHandle *&expand_sets,
			      int &num_expand_sets);

    //- get the sets tagged with the given vector of tags/values
    void get_tagged_sets(iBase_EntitySetHandle from_set,
			 iBase_TagHandle *tag_handles,
			 const char **tag_vals,
			 int num_tags,
			 iBase_EntitySetHandle *&tagged_sets,
			 int &num_tagged_sets);

    // mesh instance
    iMesh_Instance imeshImpl;

    // tag storing copy-to tag
    LocalTag copyTag;

    CESets copySets;
    CESets expandSets;

    std::vector<iBase_TagHandle> uniqueTags;
    std::set<iBase_EntitySetHandle> uniqueSets;
    double m_x[3];
  };

  inline void CopyMesh::add_unique_tag(const std::string &tag_name)
  {
    iBase_TagHandle tag_handle;
    int err;
    iMesh_getTagHandle(imeshImpl, tag_name.c_str(), &tag_handle, &err,
		       tag_name.size());
    check_error(imeshImpl, err);

    add_unique_tag(tag_handle);
  }

  inline void CopyMesh::add_unique_tag(iBase_TagHandle tag_handle)
  {
    assert(tag_handle != NULL);
    uniqueTags.push_back(tag_handle);
  }

  inline void CopyMesh::reset_sets()
  {
    copySets.clear();
    expandSets.clear();
  }

  inline CESets &CopyMesh::copy_sets()
  {
    return copySets;
  }
  
  inline CESets &CopyMesh::expand_sets()
  {
    return expandSets;
  }
  
  inline std::set<iBase_EntitySetHandle> &CopyMesh::unique_sets()
  {
    return uniqueSets;
  }

};
#endif

  