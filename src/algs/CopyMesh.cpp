#include "meshkit/CopyMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "meshkit/LocalSet.hpp"

#include "iMesh_extensions.h"
#include "MBCN.h"

#define CHKERR(err) do {                        \
    if ((err) != iBase_SUCCESS)                   \
    return iBase_ErrorType(err);                \
    } while(false)


namespace MeshKit
{

// static registration of this  mesh scheme
moab::EntityType CopyMesh_tps[] = { moab::MBVERTEX,
                                    moab::MBEDGE,
                                    moab::MBTRI,
                                    moab::MBHEX,
                                    moab::MBMAXTYPE};
const moab::EntityType* CopyMesh::output_types()
{ return CopyMesh_tps; }

CopyMesh::CopyMesh(MKCore *mkcore, const MEntVector &me_vec)
    : MeshScheme(mkcore, me_vec),
      mesh(mkcore->imesh_instance()),
      copyTag(mkcore, "__CopyMeshTag"),
      transform(new Copy::Identity()),
      copySets(mkcore),
      expandSets(mkcore)
{
  flag_process_ce_set = false;
}

CopyMesh::~CopyMesh()
{
    delete transform;
}


bool CopyMesh::add_modelent(ModelEnt *model_ent)
{
    return MeshOp::add_modelent(model_ent);
}

void CopyMesh::setup_this()
{
}

void CopyMesh::execute_this()
{
    iBase_EntitySetHandle t = NULL;
    for (MEntSelection::iterator mit = mentSelection.begin();
         mit != mentSelection.end(); mit++) {
        ModelEnt *me = mit->first;
        //orig_ents[i++] = reinterpret_cast<iBase_EntityHandle> (me->mesh_handle());
        t = reinterpret_cast<iBase_EntitySetHandle> (me->mesh_handle());
    }

    LocalSet set(this->mk_core());
    std::vector<iBase_EntityHandle> original;
    original.clear();
    iBase_EntityType type_t;
    IBERRCHK(mesh->getEntType((iBase_EntityHandle)t, type_t), *mesh);
    // check if this is a set(4)
    if((int) type_t !=4){
        std::cout << "Invalid model ent set, bailing out.." << std::endl;
        exit(0);
    }

    IBERRCHK(mesh->getEntities(t, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, original), *mesh);
    IBERRCHK(mesh->addEntArrToSet(&original[0], original.size(), set), *mesh);
    do_copy(set);
}

int CopyMesh::getStructure(iMesh_Instance instance,
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

void CopyMesh::do_copy(iMesh::EntitySetHandle set_handle)
{
    assert(transform);

    LocalTag local_tag(this->mk_core());

    std::vector<iBase_EntityHandle> ents;
    std::vector<iBase_EntityHandle> verts;
    std::vector<int> indices;
    std::vector<int> offsets;
    //    IBERRCHK(iMesh_getStructure(mesh->instance(), set_handle, ents, verts,
    //                                indices, offsets), *mesh);
    getStructure(mesh->instance(), set_handle, ents, verts,
                 indices, offsets);
    // copy the vertices
    std::vector<iBase_EntityHandle> new_verts;
    transform->transform(mesh, verts, new_verts);
    assert(new_verts.size() == verts.size());

    // set the local copy tags on vertices
    // XXX: Should this really happen? Doing so adds more entities to copy sets
    // than explicitly passed into this function. This may be a domain-specific
    // question.
    IBERRCHK(mesh->setEHArrData(&verts[0], verts.size(), local_tag,
             &new_verts[0]), *mesh);

    // now connect the new vertices to make the higher-dimension entities
    connect_the_dots(mesh, local_tag, ents, indices, offsets, new_verts);

    // take care of copy/expand sets
    update_sets();

    link_expand_sets(expandSets, local_tag);

    if(flag_process_ce_set == true)
      process_ce_sets(mesh, copySets.sets(), local_tag);
    process_ce_sets(mesh, expandSets.sets(), local_tag);

    tag_copy_sets(copySets, local_tag, copyTag);
}

void CopyMesh::update_sets()
{
    copySets.update_tagged_sets();
    expandSets.update_tagged_sets();
}

void CopyMesh::tag_copied_sets(const char **tag_names, const char **tag_vals,
                               const int num_tags)
{
    for (int t = 0; t < num_tags; t++) {
        iMesh::TagHandle tag;
        IBERRCHK(mesh->getTagHandle(tag_names[t], tag), *mesh);

        tag_copy_sets(mesh, copyTag, copySets.sets(), tag,
                      tag_vals ? tag_vals[t] : NULL);
    }
}

void CopyMesh::tag_copied_sets(iMesh::TagHandle *tags, const char **tag_vals,
                               const int num_tags)
{
    for (int t = 0; t < num_tags; t++)
        tag_copy_sets(mesh, copyTag, copySets.sets(), tags[t],
                      tag_vals ? tag_vals[t] : NULL);
}

void connect_the_dots(iMesh *mesh, iBase_TagHandle local_tag,
                      const std::vector<iBase_EntityHandle> &ents,
                      const std::vector<int> &indices,
                      const std::vector<int> &offsets,
                      const std::vector<iBase_EntityHandle> &verts)
{
    std::vector<iMesh_EntityTopology> topos(ents.size());
    IBERRCHK(mesh->getEntArrTopo(&ents[0], ents.size(), &topos[0]), *mesh);

    // scan forward to first non-vertex
    size_t pos = 0;
    while (pos < topos.size() && iMesh_POINT == topos[pos])
        pos++;
    if (pos == topos.size()) return;

    // for each run of same size & type
    std::vector<iBase_EntityHandle> connect, new_ents;
    size_t begin, end = pos;
    while (end < ents.size()) {
        // get next run; end points to start of *next* element,
        // or ents_size if no elems left
        begin = end++;

        iMesh_EntityTopology topo = topos[begin];
        int vtx_per_ent = offsets[end] - offsets[begin];
        while (end < ents.size() &&
               topos[end] == topo &&
               offsets[end+1] - offsets[end] == vtx_per_ent)
            end++;
        size_t num_ents = end - begin;

        int mbcn_type;
        int num_corner_verts;
        iMesh_MBCNType(topo, &mbcn_type);
        MBCN_VerticesPerEntity(mbcn_type, &num_corner_verts);

        // build vector of vtx handles
        connect.resize(vtx_per_ent * num_ents);
        for (size_t i = 0; i < connect.size(); i++)
            connect[i] = verts[indices[offsets[begin] + i]];

        // create entities
        new_ents.resize(num_ents);

        if (num_corner_verts == vtx_per_ent) {
            IBERRCHK(mesh->createEntArr(topo, &connect[0], connect.size(),
                     &new_ents[0]), *mesh);
        }
        else {
            // use single-entity function in this case, entity might have higher-
            // order nodes (imesh fcn doesn't have argument for # entities)
            for (size_t i = 0; i < num_ents; i++) {
                IBERRCHK(mesh->createEnt(topo, &connect[i*vtx_per_ent],
                         vtx_per_ent, new_ents[i]), *mesh);
            }
        }

        // set the local copy tags
        IBERRCHK(mesh->setEHArrData(&ents[begin], num_ents, local_tag,
                                    &new_ents[0]), *mesh);
    }
}

} // namespace MeshKit
