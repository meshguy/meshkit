#ifndef MESHKIT_MERGE_MESH_HPP
#define MESHKIT_MERGE_MESH_HPP

#include "MKException.hpp"

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"

#include <cassert>
#include <string.h>
#include <vector>
#include <set>

#include "iMesh.h"

#ifdef HAVE_MOAB
#include "MBiMesh.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp"
#endif

namespace MeshKit {

  class MKCore;


  /** \class MergeMesh MergeMesh.hpp "meshkit/MergeMesh.hpp"
   * \brief A simple class for merging meshes 
   *
   * INPUT: ModelEnts representing meshes
   * MESH TYPE(S): ALL TYPES
   * OUTPUT: Copied mesh along with the existing mesh
   * DEPENDENCIES: (none)
   * CONSTRAINTS: (none)
   *
   * This class performs the trivial task of merging meshes.  There can be
   * multiple instances of this class, and therefore it is pointed to and managed by ....
   *
   * Each instance of this class stores all the ModelEnt's representing the mesh data,
   * and after execution after meshing new entities are created and tag propagation happens.
   */


class MergeMesh : public MeshScheme
{
public:
  /* \brief Constructor
   *
   * Create a new MergeMesh instance
   * \param impl the iMesh instance handle for the mesh
   */


  MergeMesh(MKCore *mkcore, const MEntVector &me_vec);
  
  /* \brief Destructor
   */
  virtual ~MergeMesh();


    /**\brief Get class name */
    static const char* name();

    /**\brief Function returning whether this scheme can mesh entities of t
     *        the specified dimension.
     *\param dim entity dimension
     */
    static bool can_mesh(iBase_EntityType dim);

    /** \brief Function returning whether this scheme can mesh the specified entity
     *
     * Used by MeshOpFactory to find scheme for an entity.
     * \param me ModelEnt being queried
     * \return If true, this scheme can mesh the specified ModelEnt
     */
    static bool can_mesh(ModelEnt *me);

    /**\brief Get list of mesh entity types that can be generated.
     *\return array terminated with \c moab::MBMAXTYPE
     */
    static const moab::EntityType* output_types();

    /** \brief Return the mesh entity types operated on by this scheme
     * \return array terminated with \c moab::MBMAXTYPE
     */
    virtual const moab::EntityType* mesh_types_arr() const;

    /** \brief Re-implemented here so we can check topological dimension of model_ent
     * \param model_ent ModelEnt being added
     */
    virtual bool add_modelent(ModelEnt *model_ent);

    //! Setup is a no-op, but must be provided since it's pure virtual
    virtual void setup_this();

    //! The only setup/execute function we need, since meshing vertices is trivial
    virtual void execute_this();

  /* \brief Return the imesh instance
   */
  iMesh* impl() const {return imeshImpl;}



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
  void set_merge_params(double mergeTol = 1e-4, int update_sets = 0,
			  int do_merge = 1, iBase_TagHandle mergeTag = NULL);
private:
  
  iMesh *imeshImpl;
  double mergeTol, mergeTolSq;
  iBase_TagHandle mergeTag;
  int updateSets;
  int doMerge;

#ifdef HAVE_MOAB

    //- given a kdtree, set tag on vertices in leaf nodes with vertices
    //- to which they should be merged
  MBErrorCode find_merged_to(MBAdaptiveKDTree & tree, MBEntityHandle &tree_root, MBTag merged_to);

  
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

  inline const char* MergeMesh::name()
  {
    return "MergeMesh";
  }

  inline bool MergeMesh::can_mesh(iBase_EntityType)
  {
    return false;
  }

  inline bool MergeMesh::can_mesh(ModelEnt *)
  {
    return true;
  }

  inline const moab::EntityType* MergeMesh::mesh_types_arr() const
  {
    return output_types();
  }

} // namespace

#endif
