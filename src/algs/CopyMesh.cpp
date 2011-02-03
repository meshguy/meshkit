#include "meshkit/CopyMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"

#include "utils/CopyUtils.hpp"
#include "LocalSet.hpp"
#include "SimpleArray.hpp"

namespace MeshKit 
{
  //// static registration of this  mesh scheme
  //  moab::EntityType CopyMesh_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBHEX};
  //  iBase_EntityType CopyMesh_mtp = iBase_REGION;
  //
  //  static int success = MKCore::register_meshop("CopyMesh", &CopyMesh_mtp, 1, CopyMesh_tps, 3,
  //					       CopyMesh::factory, MeshOp::canmesh_edge);

  MeshOp *CopyMesh::factory(MKCore *mkcore, const MEntVector &me_vec)
  {
    return new CopyMesh(mkcore, me_vec);
  }


  CopyMesh::CopyMesh(MKCore *mkcore, const MEntVector &me_vec)
    : MeshScheme(mkcore, me_vec),
      imeshImpl(reinterpret_cast<iMesh_Instance> (mkcore->mb_imesh())),
      copyTag(reinterpret_cast<iMesh_Instance> (mkcore->mb_imesh()), "__CopyMeshTag"),
      copySets(reinterpret_cast<iMesh_Instance> (mkcore->mb_imesh())),
      expandSets(reinterpret_cast<iMesh_Instance> (mkcore->mb_imesh()))
  {
    m_x[0] = 0.0;
    m_x[1] = 0.0;
    m_x[2] = 0.0;
  }

  CopyMesh::~CopyMesh()
  {

  }


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
    int err = 0, i  = 0;
    int orig_ents_alloc = 0, orig_ents_size = 0;
    iBase_EntityHandle *orig_ents = NULL;

    iBase_EntityHandle *new_ents;
    int new_ents_alloc, new_ents_size;
    new_ents = NULL;
    new_ents_alloc = 0;
    new_ents_size = 0;

    for (MEntSelection::iterator mit = mentSelection.begin(); mit != mentSelection.end(); mit++) {
      ModelEnt *me = mit->first;
      orig_ents[orig_ents_size] = reinterpret_cast<iBase_EntityHandle> (me->mesh_handle());
      orig_ents_size++;
    }
    std::cout << orig_ents_size << " value of orig_ents_size ? shouldn't be zero!" << std::endl;
    copy(orig_ents , i, copy::Translate(m_x),
	 &new_ents, &new_ents_alloc, &new_ents_size, false);
  }

  void CopyMesh::set_location(const double x[])
  {
    m_x[0] = x[0];
    m_x[1] = x[1];
    m_x[2] = x[2];
  }

  void CopyMesh::copy(iBase_EntityHandle *ent_handles,
		      int num_ents,
		      const copy::Transform &trans,
		      iBase_EntityHandle **new_ents,
		      int *new_ents_allocated,
		      int *new_ents_size,
		      bool do_merge)
  {
    int err;

    LocalSet set(imeshImpl);
  
    iMesh_addEntArrToSet(imeshImpl, ent_handles, num_ents, set, &err);
    check_error(imeshImpl, err);

    copy(set, trans, new_ents, new_ents_allocated, new_ents_size, do_merge);
  }

  void CopyMesh::copy(iBase_EntitySetHandle set_handle,
		      const copy::Transform &trans,
		      iBase_EntityHandle **new_ents,
		      int *new_ents_allocated,
		      int *new_ents_size,
		      bool do_merge)
  {
    int err;
    LocalTag local_tag(imeshImpl);

    SimpleArray<iBase_EntityHandle> ents;
    SimpleArray<iBase_EntityHandle> verts;
    SimpleArray<int> indices;
    SimpleArray<int> offsets;

    iMesh_getStructure(imeshImpl, set_handle, ARRAY_INOUT(ents),
		       ARRAY_INOUT(verts), ARRAY_INOUT(indices),
		       ARRAY_INOUT(offsets), &err);
    check_error(imeshImpl, err);

    // copy the vertices
    SimpleArray<iBase_EntityHandle> new_verts;
    trans(imeshImpl, ARRAY_IN(verts), ARRAY_INOUT(new_verts));
    assert(new_verts.size() == verts.size());

    // set the local copy tags on vertices
    // XXX: Should this really happen? Doing so adds more entities to copy sets
    // than explicitly passed into this function. This may be a domain-specific
    // question.
    iMesh_setEHArrData(imeshImpl, ARRAY_IN(verts), local_tag,
		       ARRAY_IN(new_verts), &err);
    check_error(imeshImpl, err);

    // now connect the new vertices to make the higher-dimension entities
    connect_the_dots(imeshImpl, ARRAY_IN(ents), local_tag, &indices[0],
		     &offsets[0], &new_verts[0]);

    // take care of copy/expand sets
    update_sets();

    link_expand_sets(expandSets, local_tag);

    process_ce_sets(imeshImpl, copySets.sets(), local_tag);
    process_ce_sets(imeshImpl, expandSets.sets(), local_tag);

    tag_copy_sets(copySets, local_tag, copyTag);

    // get all the copies
    if (new_ents) {
      iMesh_getEHArrData(imeshImpl, ARRAY_IN(ents), local_tag,
			 new_ents, new_ents_allocated, new_ents_size, &err);
      check_error(imeshImpl, err);
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
      iMesh_getTagHandle(imeshImpl, tag_names[t], &tag, &err,
			 strlen(tag_names[t]));
      check_error(imeshImpl, err);

      tag_copy_sets(imeshImpl, copyTag, copySets.sets(), tag,
		    tag_vals ? tag_vals[t] : NULL);
    }
  }

  void CopyMesh::tag_copied_sets(iBase_TagHandle *tags, const char **tag_vals,
				 const int num_tags)
  {
    for (int t = 0; t < num_tags; t++)
      tag_copy_sets(imeshImpl, copyTag, copySets.sets(), tags[t],
		    tag_vals ? tag_vals[t] : NULL);
  }

} // namespace MeshKit
