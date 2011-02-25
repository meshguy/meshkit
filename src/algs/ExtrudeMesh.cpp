#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/LocalSet.hpp"

#include "CopyUtils.hpp"

#include "iMesh_extensions.h"

namespace MeshKit
{

  // static registration of this  mesh scheme
  moab::EntityType ExtrudeMesh_tps[] = { moab::MBVERTEX,
                                         moab::MBEDGE,
                                         moab::MBTRI,
                                         moab::MBHEX,
                                         moab::MBMAXTYPE};

  const moab::EntityType* ExtrudeMesh::output_types()
    { return ExtrudeMesh_tps; }

  ExtrudeMesh::ExtrudeMesh(MKCore *mkcore, const MEntVector &me_vec)
    : MeshScheme(mkcore, me_vec),
      mesh(mkcore->imesh_instance()),
      extrudeTag(mkcore, "__ExtrudeMeshTag"),
      copyTag(mkcore, "__CopyMeshTag"),
      transform(0),
      copyFaces(false),
      extrudeSets(mkcore),
      copySets(mkcore),
      expandSets(mkcore)
  {}

  ExtrudeMesh::~ExtrudeMesh()
  {}


  bool ExtrudeMesh::add_modelent(ModelEnt *model_ent)
  {
    return MeshOp::add_modelent(model_ent);
  }

  void ExtrudeMesh::setup_this()
  {}

  void ExtrudeMesh::execute_this()
  {
    std::vector<iMesh::EntityHandle> orig_ents(mentSelection.size());

    int i = 0;
    for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {
      ModelEnt *me = mit->first;
      orig_ents[i++] = reinterpret_cast<iBase_EntityHandle> (me->mesh_handle());
    }

    LocalSet set(this->mk_core());

    IBERRCHK(mesh->addEntArrToSet(&orig_ents[0], orig_ents.size(), set), *mesh);
    do_extrude(set);
  }

  void ExtrudeMesh::update_sets()
  {
    copySets.update_tagged_sets();
    expandSets.update_tagged_sets();
  }

  void ExtrudeMesh::do_extrude(iBase_EntitySetHandle src)
  {
    assert(transform && transform->steps() > 0);

    update_sets();

    std::vector<iBase_EntityHandle> ents;
    std::vector<iBase_EntityHandle> verts;
    std::vector<int> indices;
    std::vector<int> offsets;

    IBERRCHK(iMesh_getStructure(mesh->instance(), src, ents, verts,
                                indices, offsets), *mesh);

    if (ents.size() == 0) return;

    std::vector<iBase_EntityHandle> curr;
    std::vector<iBase_EntityHandle> next;
    std::vector<int> normals;

    LocalTag local_tag(this->mk_core());

    curr.resize(verts.size());
    next.resize(verts.size());
    transform->transform(1, mesh, verts, next);

    // Get the offset between vertices between steps
    Vector<3> xa, xb, dx;
    IBERRCHK(mesh->getVtxCoord(next[0],  xa[0], xa[1], xa[2]), *mesh);
    IBERRCHK(mesh->getVtxCoord(verts[0], xb[0], xb[1], xb[2]), *mesh);
    dx = xa-xb;

    get_normals(verts, indices, offsets, dx, normals);

    // Make the first set of volumes
    connect_up_dots(&ents[0], ents.size(), local_tag, &normals[0], &indices[0],
                    &offsets[0], &verts[0], &next[0]);

    // Now do the rest
    for (int i=2; i<=transform->steps(); i++) {
      std::swap(curr, next);
      transform->transform(i, mesh, verts, next);
      connect_up_dots(&ents[0], ents.size(), local_tag, &normals[0],
                      &indices[0], &offsets[0], &curr[0], &next[0]);
    }

    tag_copy_sets(extrudeSets, local_tag, extrudeTag);

    if (copyFaces) {
      // set the local copy tags on vertices
      // XXX: Should this really happen? Doing so adds more entities to copysets
      // than explicitly passed into this function. This may be a domain-
      // specific question.
      IBERRCHK(mesh->setEHArrData(&verts[0], verts.size(), local_tag,
                                  &next[0]), *mesh);

      connect_the_dots(mesh, local_tag, ents, indices, offsets, next);

      link_expand_sets(expandSets, local_tag);

      process_ce_sets(mesh, copySets.sets(), local_tag);
      process_ce_sets(mesh, expandSets.sets(), local_tag);

      tag_copy_sets(copySets, local_tag, copyTag);
    }
  }

  // calculate the normals for each face (1 = towards v, -1 = away from v)
  // TODO: this can fail with non-convex faces
  void ExtrudeMesh::get_normals(const std::vector<iBase_EntityHandle> &verts,
                                const std::vector<int> &indices,
                                const std::vector<int> &offsets,
                                const Vector<3> &dv, std::vector<int> &normals)
  {
    size_t size = offsets.size() - 1;
    normals.resize(size);

    for(size_t i=0; i<size; i++) {
      Vector<3> a, b;
      iBase_EntityHandle curr_verts[3];

      if(offsets[i+1] - offsets[i] > 2) { // face
        for(int j=0; j<3; j++)
          curr_verts[j] = verts[indices[ offsets[i]+j ]];

        std::vector< Vector<3> > coords(3);
        IBERRCHK(mesh->getVtxArrCoords(curr_verts, 3, iBase_INTERLEAVED,
                                       vec2ptr(coords)), *mesh);

        a = coords[1] - coords[0];
        b = coords[2] - coords[1];
        normals[i] = (vector_product(a, b) % dv) > 0 ? 1:-1;
      }
      else if(offsets[i+1] - offsets[i] == 2) { // line
        normals[i] = 1; // TODO: figure out a way of distinguishing swapped
                        // lines
      }
      else // vertex
        normals[i] = 1;
    }
  }

  void ExtrudeMesh::connect_up_dots(
    iBase_EntityHandle *src, int size, iBase_TagHandle local_tag,
    int *pre_norms,  int *pre_inds,  int *pre_offs,  iBase_EntityHandle *pre,
    int *post_norms, int *post_inds, int *post_offs, iBase_EntityHandle *post)
  {
    for(int i=0; i<size; i++) {
      int count = pre_offs[i+1] - pre_offs[i];

      // If the normal is facing in the wrong direction (away from the
      // translation) we add the vertices in reverse order. Otherwise, we go
      // in the usual order. If count is 2, then we are creating quads and so
      // need to swap the order of the post set of verts.

      int dx = pre_norms [i];
      int dy = post_norms[i] * (count == 2 ? -1:1);
      int x  = (dx == 1) ? pre_offs [i] : pre_offs [i+1]-1;
      int y  = (dy == 1) ? post_offs[i] : post_offs[i+1]-1;

      iBase_EntityHandle *nodes = new iBase_EntityHandle[count*2];
      for(int j=0; j<count; j++) {
        nodes[j]       = pre [ pre_inds [x + dx*j] ];
        nodes[j+count] = post[ post_inds[y + dy*j] ];
      }

      iBase_EntityHandle out;

      iMesh::Error err;
      if(count == 4)      // quad
        err = mesh->createEnt(iMesh_HEXAHEDRON, nodes, 8, out);
      else if(count == 3) // tri
        err = mesh->createEnt(iMesh_PRISM, nodes, 6, out);
      else if(count == 2) // line
        err = mesh->createEnt(iMesh_QUADRILATERAL, nodes, 4, out);
      else if(count == 1) // vertex
        err = mesh->createEnt(iMesh_LINE_SEGMENT, nodes, 2, out);
      else
        throw Error(iBase_FAILURE, "Couldn't extrude face; unusual shape.");

      IBERRCHK(err, *mesh);
      delete[] nodes;

      IBERRCHK(mesh->setEHData(src[i], local_tag, out), *mesh);
    }

    process_ce_sets(mesh, extrudeSets.sets(), local_tag);
  }
} // namespace MeshKit
