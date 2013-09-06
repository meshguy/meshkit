#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/CopyMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/LocalSet.hpp"

#include "iMesh_extensions.h"

#define CHKERR(err) do {                        \
  if ((err) != iBase_SUCCESS)                   \
    return iBase_ErrorType(err);                \
  } while(false)

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

  int ExtrudeMesh::getStructure(iMesh_Instance instance,
                             iBase_EntitySetHandle set,
                             std::vector<iBase_EntityHandle> &ents,
                             std::vector<iBase_EntityHandle> &unique_adj,
                             std::vector<int> &indices,
                             std::vector<int> &offsets)
  {
      // 1) Get source entities, making sure verts are first
      int num;
      int err;
      iMesh_getNumOfTypeRec(instance, set, iBase_ALL_TYPES, true, &num, &err);
      CHKERR(err);

      ents.resize(num);
      offsets.resize(num+1);

      iBase_EntityHandle *block = &ents[0];
      int block_alloc = ents.size(), block_size, num_verts = 0;
      for (int t = iMesh_POINT; t < iMesh_ALL_TOPOLOGIES && block_alloc; ++t) {
          iMesh_getEntitiesRec(instance, set, iBase_ALL_TYPES, t, true,
                               &block, &block_alloc, &block_size, &err);
          CHKERR(err);

          block_alloc -= block_size;
          block += block_size;
          if (t == iMesh_POINT)
              num_verts = block_size;
      }

      // 2) Get verts adjacent to all source entitites (verts are adj to themselves)
      std::vector<iBase_EntityHandle> all_adj(ents.begin(), ents.begin()+num_verts);

      // first, fill the vertex-vertex adjacencies
      for (int i = 0; i < num_verts; ++i)
          offsets[i] = i;

      iBase_EntityHandle *tmp_adj = NULL;
      int tmp_adj_alloc = 0, tmp_adj_size;
      int *tmp_off = &offsets[num_verts];
      int tmp_off_alloc = offsets.size() - num_verts, tmp_off_size;
      iMesh_getEntArrAdj(instance, &ents[num_verts], ents.size()-num_verts,
                         iBase_VERTEX, &tmp_adj, &tmp_adj_alloc, &tmp_adj_size,
                         &tmp_off, &tmp_off_alloc, &tmp_off_size, &err);
      CHKERR(err);

      // shift all the offsets to account for vertices
      for(unsigned int i = num_verts; i < offsets.size(); ++i)
          offsets[i] += num_verts;

      all_adj.reserve(all_adj.size() + tmp_adj_size);
      all_adj.insert(all_adj.end(), tmp_adj, tmp_adj+tmp_adj_size);
      free(tmp_adj);

      // 3) Get unique adjacent vertices and offsets
      unique_adj.resize(all_adj.size());
      indices.resize(all_adj.size());
      std::copy(all_adj.begin(), all_adj.end(), unique_adj.begin());
      std::sort(unique_adj.begin(), unique_adj.end());

      size_t unique_size;
      unique_size = std::unique(unique_adj.begin(), unique_adj.end()) -
              unique_adj.begin();
      unique_adj.resize(unique_size);

      for (size_t i = 0; i < all_adj.size(); ++i) {
          indices[i] = std::lower_bound(unique_adj.begin(), unique_adj.end(),
                                        all_adj[i]) - unique_adj.begin();
      }

      return 0;
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

//    IBERRCHK(iMesh_getStructure(mesh->instance(), src, ents, verts,
//                                indices, offsets), *mesh);
    getStructure(mesh->instance(), src, ents, verts,
                                    indices, offsets);

    if (ents.size() == 0) return;

    std::vector<iBase_EntityHandle> curr;
    std::vector<iBase_EntityHandle> next;
    std::vector<int> normals;

    LocalTag local_extrude_tag(this->mk_core());
    LocalTag local_copy_tag(this->mk_core());

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
    connect_up_dots(&ents[0], ents.size(), local_extrude_tag, &normals[0],
                    &indices[0], &offsets[0], &verts[0], &next[0]);

    // Now do the rest
    for (int i=2; i<=transform->steps(); i++) {
      std::swap(curr, next);
      transform->transform(i, mesh, verts, next);
      connect_up_dots(&ents[0], ents.size(), local_extrude_tag, &normals[0],
                      &indices[0], &offsets[0], &curr[0], &next[0]);
    }

    tag_copy_sets(extrudeSets, local_extrude_tag, extrudeTag);

    if (copyFaces) {
      // set the local copy tags on vertices
      // XXX: Should this really happen? Doing so adds more entities to copysets
      // than explicitly passed into this function. This may be a domain-
      // specific question.
      IBERRCHK(mesh->setEHArrData(&verts[0], verts.size(), local_copy_tag,
                                  &next[0]), *mesh);

      connect_the_dots(mesh, local_copy_tag, ents, indices, offsets, next);

      link_expand_sets(expandSets, local_copy_tag);

      process_ce_sets(mesh, copySets.sets(), local_copy_tag);
      process_ce_sets(mesh, expandSets.sets(), local_copy_tag);

      tag_copy_sets(copySets, local_copy_tag, copyTag);
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
