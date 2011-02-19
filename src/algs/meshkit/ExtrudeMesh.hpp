#ifndef EXTRUDEMESH_HPP
#define EXTRUDEMESH_HPP

#include <vector>
#include <set>

#include "meshkit/Types.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MeshScheme.hpp"
#include "meshkit/ModelEnt.hpp"

#include "meshkit/CESets.hpp"
#include "meshkit/LocalTag.hpp"
#include "meshkit/Transform.hpp"

#include "meshkit/iMesh.hpp"


namespace MeshKit {

  class MKCore;


  /** \class ExtrudeMesh ExtrudeMesh.hpp "meshkit/ExtrudeMesh.hpp"
   * \brief A simple class for extruding meshes
   *
   * INPUT: ModelEnts representing meshes
   * MESH TYPE(S): ALL TYPES
   * OUTPUT: Copied mesh along with the existing mesh
   * DEPENDENCIES: (none)
   * CONSTRAINTS: (none)
   *
   * This class performs the trivial task of extruding meshes.  There can be
   * multiple instances of this class, and therefore it is pointed to and managed by ....
   *
   * Each instance of this class stores all the ModelEnt's representing the mesh data,
   * and after execution after meshing new entities are created and tag propagation happens.
   */
  class ExtrudeMesh : public MeshScheme
  {
  public:
    //! Bare constructor
    ExtrudeMesh(MKCore *mkcore, const MEntVector &me_vec);

    //! Destructor
    virtual ~ExtrudeMesh();

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

    // ExtrudeMesh function local //////////////////

    /* \brief Return the tag used to indicate extruded sets
     */
    iBase_TagHandle extrude_tag();

    /* \brief Return the tag used to indicate copied sets
     */
    iBase_TagHandle copy_tag();

    /* \brief Set the transform function for this operation
     */
    void set_transform(const Extrude::Transform &transform);

    /* \brief Turn copying of faces on (true) or off (false)
     */
    void copy_faces(bool copy);

    void update_sets();

    /* \brief Reset the copy and expand set lists
     */
    void reset_sets();

    /* \brief Return reference to extrude sets
     */
    CESets &extrude_sets();

    /* \brief Return reference to copy sets
     */
    CESets &copy_sets();

    /* \brief Return reference to expand sets
     */
    CESets &expand_sets();
  private:
    void do_extrude(iMesh::EntitySetHandle set_handle);

    void get_normals(iBase_EntityHandle *verts, int *indices,
                     int *offsets, int size, const Vector<3> &dv,
                     std::vector<int> &normals);

    void connect_up_dots(iBase_EntityHandle *src, int size,
                         iBase_TagHandle local_tag, int *norms, int *inds,
                         int *offs, iBase_EntityHandle *pre,
                         iBase_EntityHandle *post);
    void connect_up_dots(iBase_EntityHandle *src, int size,
                         iBase_TagHandle local_tag,
                         int *pre_norms,  int *pre_inds,  int *pre_offs,
                         iBase_EntityHandle *pre,
                         int *post_norms, int *post_inds, int *post_offs,
                         iBase_EntityHandle *post);

    iMesh *mesh;                   // mesh instance
    LocalTag extrudeTag;           // tag storing extrude-to tag
    LocalTag copyTag;              // tag storing copy-to tag
    Extrude::Transform *transform; // transform function for extrusion
    bool copyFaces;

    CESets extrudeSets;
    CESets copySets;
    CESets expandSets;
  };

  inline const char* ExtrudeMesh::name()
  {
    return "ExtrudeMesh";
  }

  inline bool ExtrudeMesh::can_mesh(iBase_EntityType)
  {
    // Given just a dimension, ExtrudeMesh can't do anything since it doesn't
    // know what to extrude.
    return false;
  }

  inline bool ExtrudeMesh::can_mesh(ModelEnt *me)
  {
    return me->dimension() < 3;
  }

  inline const moab::EntityType* ExtrudeMesh::mesh_types_arr() const
  {
    return output_types();
  }

  inline void ExtrudeMesh::set_transform(const Extrude::Transform &trans)
  {
    delete transform;
    transform = trans.clone();
  }

  inline void ExtrudeMesh::copy_faces(bool copy)
  {
    copyFaces = copy;
  }

  inline void ExtrudeMesh::reset_sets()
  {
    extrudeSets.clear();
    copySets.clear();
    expandSets.clear();
  }

  inline CESets &ExtrudeMesh::copy_sets()
  {
    return copySets;
  }

  inline CESets &ExtrudeMesh::expand_sets()
  {
    return expandSets;
  }

  inline void
  ExtrudeMesh::connect_up_dots(iBase_EntityHandle *src, int size,
                               iBase_TagHandle local_tag, int *norms, int *inds,
                               int *offs, iBase_EntityHandle *pre,
                               iBase_EntityHandle *post)
  {
    connect_up_dots(src, size, local_tag,
                    norms, inds, offs, pre,
                    norms, inds, offs, post);
  }

} // namespace MeshKit
#endif


