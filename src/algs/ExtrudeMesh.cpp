#include "meshkit/ExtrudeMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/LocalSet.hpp"

#include "CopyUtils.hpp"
#include "SimpleArray.hpp"

#include "iMesh_extensions.h"
#include "vec_utils.hpp"

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
      copyFaces(false),
      extrudeTag(mkcore, "__ExtrudeMeshTag"),
      copyTag(mkcore, "__CopyMeshTag"),
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
    SimpleArray<iMesh::EntityHandle> orig_ents(mentSelection.size());

    int i = 0;
    for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {
      ModelEnt *me = mit->first;
      orig_ents[i++] = reinterpret_cast<iBase_EntityHandle> (me->mesh_handle());
    }

    LocalSet set(this->mk_core());

    IBERRCHK(mesh->addEntArrToSet(ARRAY_IN(orig_ents), set), *mesh);
    do_extrude(set);
  }

  void ExtrudeMesh::update_sets()
  {
    copySets.update_tagged_sets();
    expandSets.update_tagged_sets();
  }

  void ExtrudeMesh::do_extrude(iBase_EntitySetHandle src)
  {
    assert(steps > 0);

    update_sets();

    int err;

    SimpleArray<iBase_EntityHandle> ents;
    SimpleArray<iBase_EntityHandle> adj;
    SimpleArray<int> indices;
    SimpleArray<int> offsets;

    iMesh_getStructure(mesh->instance(), src, ARRAY_INOUT(ents),
                       ARRAY_INOUT(adj), ARRAY_INOUT(indices),
                       ARRAY_INOUT(offsets), &err);
    IBERRCHK(err, *mesh);

    if (adj.size() == 0) return;

    SimpleArray<iBase_EntityHandle> curr;
    SimpleArray<iBase_EntityHandle> next;
    std::vector<int> normals;

    LocalTag local_tag(this->mk_core());

    curr.resize(adj.size());
    next.resize(adj.size());
    transform(1, *mesh, ARRAY_IN(adj), ARRAY_INOUT(next));

    // Get the offset between vertices between steps
    Vector<3> xa, xb, dx;
    IBERRCHK(mesh->getVtxCoord(next[0], xa[0], xa[1], xa[2]), *mesh);
    IBERRCHK(mesh->getVtxCoord(curr[0], xb[0], xb[1], xb[2]), *mesh);
    dx = xa-xb;

    get_normals(&adj[0], &indices[0], &offsets[0], ents.size(), dx, normals);

    // Make the first set of volumes
    connect_up_dots(ARRAY_IN(ents), local_tag, &normals[0], &indices[0],
                    &offsets[0], &adj[0], &next[0]);

    // Now do the rest
    for (int i=2; i<=steps; i++) { // XXX!!
      std::swap(curr, next);
      transform(i, *mesh, ARRAY_IN(adj), ARRAY_INOUT(next));
      connect_up_dots(ARRAY_IN(ents), local_tag, &normals[0], &indices[0],
                      &offsets[0], &curr[0], &next[0]);
    }

    tag_copy_sets(extrudeSets, local_tag, extrudeTag);

    if (copyFaces) {
      connect_the_dots(mesh->instance(), ARRAY_IN(ents), local_tag, &indices[0],
                       &offsets[0], &next[0]);

      link_expand_sets(expandSets, local_tag);

      process_ce_sets(mesh, copySets.sets(), local_tag);
      process_ce_sets(mesh, expandSets.sets(), local_tag);

      tag_copy_sets(copySets, local_tag, copyTag);
    }
  }

  // calculate the normals for each face (1 = towards v, -1 = away from v)
  // TODO: this can fail with non-convex faces
  void ExtrudeMesh::get_normals(iBase_EntityHandle *verts, int *indices,
                                int *offsets, int size, const Vector<3> &dv,
                                std::vector<int> &normals)
  {
    int err;
    normals.resize(size);

    for(int i=0; i<size; i++) {
      double res[3], a[3], b[3];
      iBase_EntityHandle curr_verts[3];

      if(offsets[i+1] - offsets[i] > 2) { // face
        for(int j=0; j<3; j++)
          curr_verts[j] = verts[indices[ offsets[i]+j ]];

        SimpleArray<double> coords;
        iMesh_getVtxArrCoords(mesh->instance(), curr_verts, 3, iBase_INTERLEAVED,
                              ARRAY_INOUT(coords), &err);
        IBERRCHK(err, *mesh);

        for(int j=0; j<3; j++) {
          a[j] = coords[1*3 + j] - coords[0*3 + j];
          b[j] = coords[2*3 + j] - coords[1*3 + j];
        }
        normals[i] = (dot( cross(res,a,b),dv.data() ) > 0) ? 1:-1;
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
