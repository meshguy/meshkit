#include "meshkit/CopyMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"

#include "CopyUtils.hpp"
#include "meshkit/LocalSet.hpp"
#include "SimpleArray.hpp"

namespace MeshKit 
{
 
  // static registration of this  mesh scheme
  moab::EntityType CopyMesh_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBHEX};
  iBase_EntityType CopyMesh_mtp = iBase_REGION;
  static RegisterMeshOp<CopyMesh,true> INIT("CopyMesh", &CopyMesh_mtp, 1, 
                                             CopyMesh_tps,
                                             sizeof(CopyMesh_tps)/sizeof(CopyMesh_tps[0]) );


  CopyMesh::CopyMesh(MKCore *mkcore, const MEntVector &me_vec)
    : MeshScheme(mkcore, me_vec),
      imeshImpl(mkcore->imesh_instance()->instance()),
      copyTag(mkcore, "__CopyMeshTag"),
      copySets(mkcore->imesh_instance()->instance()),
      expandSets(mkcore->imesh_instance()->instance())
  {}

  CopyMesh::~CopyMesh()
  {}


  bool CopyMesh::can_mesh(ModelEnt *me)
  {
    if (me->dimension() == 3 || me->dimension() == 2) return true;
    else return false;
  }
  
  void CopyMesh::mesh_types(std::vector<moab::EntityType> &tps)
  {
    tps.push_back(moab::MBVERTEX);
    tps.push_back(moab::MBEDGE);
    tps.push_back(moab::MBTRI);
    tps.push_back(moab::MBHEX);
  }
    
  bool CopyMesh::add_modelent(ModelEnt *model_ent)
  {
    // make sure this represents geometric volumes
    if (3 != model_ent->dimension() || 2 != model_ent->dimension())
      throw Error(MK_WRONG_DIMENSION, "Wrong dimension entity added to CopyMesh.");

    return MeshOp::add_modelent(model_ent);
  }
  /*
    void CopyMesh::setup()
    {
    setup_this();
    }
  */
  void CopyMesh::setup_this()
  {
  
  }
  /*
    void CopyMesh::execute()
    {
    execute_this();
    }
  */
  void CopyMesh::execute_this()
  {
    SimpleArray<iBase_EntityHandle> orig_ents(mentSelection.size());

    int i = 0;
    for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {
      ModelEnt *me = mit->first;
      orig_ents[i++] = reinterpret_cast<iBase_EntityHandle> (me->mesh_handle());
    }

    LocalSet set(this->mk_core());
  
    IBERRCHK(imeshImpl.addEntArrToSet(ARRAY_IN(orig_ents), set), "FIXME");
    do_copy(set);
  }

  void CopyMesh::do_copy(iBase_EntitySetHandle set_handle,
                         const Copy::Transform &trans,
                         iBase_EntityHandle **new_ents,
                         int *new_ents_allocated,
                         int *new_ents_size,
                         bool do_merge)
  {
    int err;
    LocalTag local_tag(this->mk_core());

    SimpleArray<iBase_EntityHandle> ents;
    SimpleArray<iBase_EntityHandle> verts;
    SimpleArray<int> indices;
    SimpleArray<int> offsets;

    iMesh_getStructure(imeshImpl.instance(), set_handle, ARRAY_INOUT(ents),
                       ARRAY_INOUT(verts), ARRAY_INOUT(indices),
                       ARRAY_INOUT(offsets), &err);
    IBERRCHK(err, "FIXME");

    // copy the vertices
    SimpleArray<iBase_EntityHandle> new_verts;
    trans(imeshImpl, ARRAY_IN(verts), ARRAY_INOUT(new_verts));
    assert(new_verts.size() == verts.size());

    // set the local copy tags on vertices
    // XXX: Should this really happen? Doing so adds more entities to copy sets
    // than explicitly passed into this function. This may be a domain-specific
    // question.
    iMesh_setEHArrData(imeshImpl.instance(), ARRAY_IN(verts), local_tag,
                       ARRAY_IN(new_verts), &err);
    IBERRCHK(err, "FIXME");

    // now connect the new vertices to make the higher-dimension entities
    connect_the_dots(imeshImpl.instance(),
                     ARRAY_IN(ents), local_tag, &indices[0],
                     &offsets[0], &new_verts[0]);

    // take care of copy/expand sets
    update_sets();

    link_expand_sets(expandSets, local_tag);

    process_ce_sets(imeshImpl.instance(), copySets.sets(), local_tag);
    process_ce_sets(imeshImpl.instance(), expandSets.sets(), local_tag);

    tag_copy_sets(copySets, local_tag, copyTag);

    // get all the copies
    if (new_ents) {
      iMesh_getEHArrData(imeshImpl.instance(), ARRAY_IN(ents), local_tag,
                         new_ents, new_ents_allocated, new_ents_size, &err);
      IBERRCHK(err, "FIXME");
    }
  }

  void CopyMesh::update_sets()
  {
    copySets.update_tagged_sets();
    expandSets.update_tagged_sets();
  }

  void CopyMesh::tag_copied_sets(const char **tag_names, const char **tag_vals,
                                 const int num_tags)
  {
    int err;
  
    for (int t = 0; t < num_tags; t++) {
      iBase_TagHandle tag;
      iMesh_getTagHandle(imeshImpl.instance(), tag_names[t], &tag, &err,
                         strlen(tag_names[t]));
      IBERRCHK(err, "FIXME");

      tag_copy_sets(imeshImpl.instance(), copyTag, copySets.sets(), tag,
                    tag_vals ? tag_vals[t] : NULL);
    }
  }

  void CopyMesh::tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
                                 const int num_tags)
  {
    for (int t = 0; t < num_tags; t++)
      tag_copy_sets(imeshImpl.instance(), copyTag, copySets.sets(), tags[t],
                    tag_vals ? tag_vals[t] : NULL);
  }

} // namespace MeshKit
