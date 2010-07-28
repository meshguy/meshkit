#include "ExtrudeMesh.hpp"

#include "CopyUtils.hpp"
#include "LocalSet.hpp"
#include "SimpleArray.hpp"
#include "vec_utils.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

static double * vtx_diff(double *res, iMesh_Instance mesh, iBase_EntityHandle a,
                         iBase_EntityHandle b)
{
  int err;
  double xa[3],xb[3];

  iMesh_getVtxCoord(mesh,a,xa+0,xa+1,xa+2,&err);
  check_error(mesh, err);
  iMesh_getVtxCoord(mesh,b,xb+0,xb+1,xb+2,&err);
  check_error(mesh, err);

  for(int i=0; i<3; i++)
    res[i] = xa[i]-xb[i];
  return res;
}


ExtrudeMesh::ExtrudeMesh(iMesh_Instance mesh)
  : impl_(mesh),
    copy_tag_(mesh, "__CopyMeshTag"), extrude_tag_(mesh, "__ExtrudeMeshTag"),
    copy_sets_(mesh), expand_sets_(mesh), extrude_sets_(mesh)
{}

ExtrudeMesh::~ExtrudeMesh()
{}

void ExtrudeMesh::update_sets()
{
  copy_sets_.update_tagged_sets();
  expand_sets_.update_tagged_sets();
  extrude_sets_.update_tagged_sets();
}

void ExtrudeMesh::reset_sets()
{
  copy_sets_.clear();
  expand_sets_.clear();
  extrude_sets_.clear();
}

void ExtrudeMesh::translate(iBase_EntityHandle *src, int size, int steps,
                            const double *dx, bool copy_faces)
{
  extrude(src, size, extrude::Translate(dx, steps), copy_faces);
}

void ExtrudeMesh::translate(iBase_EntitySetHandle src, int steps,
                            const double *dx, bool copy_faces)
{
  extrude(src, extrude::Translate(dx, steps), copy_faces);
}

void ExtrudeMesh::translate(iBase_EntityHandle *src, iBase_EntityHandle *dest,
                            int size, int steps)
{
  int err;
  LocalSet src_set(impl_), dest_set(impl_);
  
  iMesh_addEntArrToSet(impl_, src, size, src_set, &err);
  check_error(impl_, err);
  iMesh_addEntArrToSet(impl_, dest, size, dest_set, &err);
  check_error(impl_, err);

  translate(src_set, dest_set, steps);
}

void ExtrudeMesh::translate(iBase_EntitySetHandle src, 
                            iBase_EntitySetHandle dest, int steps)
{
  int err;

  // Deduce the per-step displacement vector "dx"
  // Note: we assume that src and dest are the same shape, etc.
  double dx[3];
  double coords[2][3];

  iBase_EntitySetHandle ends[] = { src, dest };
  for(int i=0; i<2; i++) {
    iMesh_EntityIterator iter;
    iMesh_initEntIter(impl_, ends[i], iBase_FACE, iMesh_ALL_TOPOLOGIES, &iter,
                      &err);
    check_error(impl_, err);

    iBase_EntityHandle face;
    int has_data;
    iMesh_getNextEntIter(impl_, iter, &face, &has_data, &err);
    check_error(impl_, err);
    if(!has_data) {
      iMesh_endEntIter(impl_, iter, &err);
      throw MKException(iBase_FAILURE, "Couldn't get a face");
    }

    iMesh_endEntIter(impl_, iter, &err);
    check_error(impl_, err);

    SimpleArray<iBase_EntityHandle> verts;
    iMesh_getEntAdj(impl_, face, iBase_VERTEX, ARRAY_INOUT(verts), &err);
    check_error(impl_, err);

    iMesh_getVtxCoord(impl_, verts[0], coords[i]+0, coords[i]+1, coords[i]+2,
                      &err);
    check_error(impl_, err);
  }

  for(int i=0; i<3; i++)
    dx[i] = (coords[1][i]-coords[0][i]);

  extrude(src, dest, extrude::Translate(dx, steps));
}

void ExtrudeMesh::rotate(iBase_EntityHandle *src, int size, int steps,
                         const double *origin, const double *z, double angle,
                         bool copy_faces)
{
  extrude(src, size, extrude::Rotate(origin, z, angle, steps), copy_faces);
}

void ExtrudeMesh::rotate(iBase_EntitySetHandle src, int steps,
                         const double *origin, const double *z, double angle,
                         bool copy_faces)
{
  extrude(src, extrude::Rotate(origin, z, angle, steps), copy_faces);
}

void ExtrudeMesh::rotate(iBase_EntityHandle *src, iBase_EntityHandle *dest,
                         int size, int steps, const double *origin,
                         const double *z, double angle)
{
  extrude(src, dest, size, extrude::Rotate(origin, z, angle, steps));
}

void ExtrudeMesh::rotate(iBase_EntitySetHandle src, iBase_EntitySetHandle dest,
                         int steps, const double *origin, const double *z,
                         double angle)
{
  extrude(src, dest, extrude::Rotate(origin, z, angle, steps));
}

void ExtrudeMesh::extrude(iBase_EntityHandle *src, int size, 
                          const extrude::Transform &trans, bool copy_faces)
{
  int err;
  LocalSet set(impl_);
  
  iMesh_addEntArrToSet(impl_, src, size, set, &err);
  check_error(impl_, err);

  extrude(set, trans, copy_faces);
}

void ExtrudeMesh::extrude(iBase_EntityHandle *src, iBase_EntityHandle *dest,
                          int size, const extrude::Transform &trans)
{
  int err;
  LocalSet src_set(impl_), dest_set(impl_);
  
  iMesh_addEntArrToSet(impl_, src, size, src_set, &err);
  check_error(impl_, err);
  iMesh_addEntArrToSet(impl_, dest, size, dest_set, &err);
  check_error(impl_, err);

  extrude(src_set, dest_set, trans);
}

void ExtrudeMesh::extrude(iBase_EntitySetHandle src,
                          const extrude::Transform &trans, bool copy_faces)
{
  assert(trans.steps() > 0);

  update_sets();

  int err;

  SimpleArray<iBase_EntityHandle> ents;
  SimpleArray<iBase_EntityHandle> adj;
  SimpleArray<int> indices;
  SimpleArray<int> offsets;

  iMesh_getStructure(impl_, src, ARRAY_INOUT(ents), ARRAY_INOUT(adj),
                     ARRAY_INOUT(indices), ARRAY_INOUT(offsets), &err);
  check_error(impl_, err);

  double dx[3];
  SimpleArray<iBase_EntityHandle> curr;
  SimpleArray<iBase_EntityHandle> next;
  std::vector<int> normals;

  LocalTag local_tag(impl_);

  curr.resize(adj.size());
  next.resize(adj.size());
  trans(1, impl_, ARRAY_IN(adj), ARRAY_INOUT(next));

  vtx_diff(dx, impl_, next[0], adj[0]);
  normals = get_normals(&adj[0], &indices[0], &offsets[0], ents.size(), dx);

  // Make the first set of volumes
  connect_up_dots(ARRAY_IN(ents), local_tag, &normals[0], &indices[0],
                  &offsets[0], &adj[0], &next[0]);

  // Now do the rest
  for (int i=2; i<=trans.steps(); i++) {
    std::swap(curr, next);
    trans(i, impl_, ARRAY_IN(adj), ARRAY_INOUT(next));
    connect_up_dots(ARRAY_IN(ents), local_tag, &normals[0], &indices[0],
                    &offsets[0], &curr[0], &next[0]);
  }

  if (copy_faces) {
    connect_the_dots(impl_, ARRAY_IN(ents), local_tag, &indices[0], &offsets[0],
                     &next[0]);
  }

  tag_all_sets(local_tag);
}

void ExtrudeMesh::extrude(iBase_EntitySetHandle src, iBase_EntitySetHandle dest,
                          const extrude::Transform &trans)
{
  assert(trans.steps() > 0);

  update_sets();

  int err;

  SimpleArray<iBase_EntityHandle> ents;
  SimpleArray<iBase_EntityHandle> adj;
  SimpleArray<int> indices;
  SimpleArray<int> offsets;

  iMesh_getStructure(impl_, src, ARRAY_INOUT(ents), ARRAY_INOUT(adj),
                     ARRAY_INOUT(indices), ARRAY_INOUT(offsets), &err);
  check_error(impl_, err);

  SimpleArray<iBase_EntityHandle> ents2;
  SimpleArray<iBase_EntityHandle> adj2;
  SimpleArray<int> indices2;
  SimpleArray<int> offsets2;

  iMesh_getStructure(impl_, dest, ARRAY_INOUT(ents2), ARRAY_INOUT(adj2),
                     ARRAY_INOUT(indices2), ARRAY_INOUT(offsets2), &err);
  check_error(impl_, err);

  double dx[3];
  SimpleArray<iBase_EntityHandle> curr;
  SimpleArray<iBase_EntityHandle> next;
  std::vector<int> normals;

  LocalTag local_tag(impl_);

  if (trans.steps() == 0) {
    next.resize(adj.size());
    std::copy(adj.begin(), adj.end(), next.begin());
  }
  else {
    curr.resize(adj.size());
    next.resize(adj.size());
    trans(1, impl_, ARRAY_IN(adj), ARRAY_INOUT(next));

    vtx_diff(dx, impl_, next[0], adj[0]);
    normals = get_normals(&adj[0], &indices[0], &offsets[0], ents.size(), dx);

    // Make the first set of volumes
    connect_up_dots(ARRAY_IN(ents), local_tag, &normals[0], &indices[0],
                    &offsets[0], &adj[0], &next[0]);

    // Make the inner volumes
    for(int i=2; i<trans.steps(); i++) {
      std::swap(curr, next);
      trans(i, impl_, ARRAY_IN(adj), ARRAY_INOUT(next));
      connect_up_dots(ARRAY_IN(ents), local_tag, &normals[0], &indices[0],
                      &offsets[0], &curr[0], &next[0]);
    }
  }

  // Connect to the destination set
  vtx_diff(dx, impl_, adj2[indices2[ offsets2[0] ]],
                      next[indices [ offsets [0] ]]);
  if(normals.empty())
    normals = get_normals(&adj[0], &indices[0], &offsets[0], ents.size(), dx);
  std::vector<int> normals2 = get_normals(&adj2[0], &indices2[0],
                                          &offsets2[0], ents.size(), dx);

  connect_up_dots(ARRAY_IN(ents), local_tag,
                  &normals[0],  &indices[0],  &offsets[0],  &next[0],
                  &normals2[0], &indices2[0], &offsets2[0], &adj2[0]);

  tag_all_sets(local_tag);
}

// calculate the normals for each face (1 = towards v, -1 = away from v)
// TODO: this can fail with non-convex faces
std::vector<int> ExtrudeMesh::get_normals(iBase_EntityHandle *verts,
                                          int *indices, int *offsets, int size,
                                          double *dv)
{
  int err;
  std::vector<int> normals(size);

  for(int i=0; i<size; i++) {
    double res[3], a[3], b[3];
    iBase_EntityHandle curr_verts[3];

    if(offsets[i+1] - offsets[i] > 2) { // face
      for(int j=0; j<3; j++)
        curr_verts[j] = verts[indices[ offsets[i]+j ]];

      SimpleArray<double> coords;
      iMesh_getVtxArrCoords(impl_, curr_verts, 3, iBase_INTERLEAVED,
                            ARRAY_INOUT(coords), &err);
      check_error(impl_, err);

      for(int j=0; j<3; j++) {
        a[j] = coords[1*3 + j] - coords[0*3 + j];
        b[j] = coords[2*3 + j] - coords[1*3 + j];
      }
      normals[i] = (dot( cross(res,a,b),dv ) > 0) ? 1:-1;
    }
    else if(offsets[i+1] - offsets[i] == 2) { // line
      normals[i] = 1; // TODO: figure out a way of distinguishing swapped
                      // lines
    }
    else // vertex
      normals[i] = 1;
  }

  return normals;
}

void ExtrudeMesh::connect_up_dots(
  iBase_EntityHandle *src, int size, iBase_TagHandle local_tag,
  int *pre_norms,  int *pre_inds,  int *pre_offs,  iBase_EntityHandle *pre,
  int *post_norms, int *post_inds, int *post_offs, iBase_EntityHandle *post)
{
  int err;

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

    int status;
    iBase_EntityHandle out;

    if(count == 4)      // quad
      iMesh_createEnt(impl_, iMesh_HEXAHEDRON, nodes, 8, &out, &status, &err);
    else if(count == 3) // tri
      iMesh_createEnt(impl_, iMesh_PRISM, nodes, 6, &out, &status, &err);
    else if(count == 2) // line
      iMesh_createEnt(impl_, iMesh_QUADRILATERAL, nodes, 4, &out, &status,
                      &err);
    else if(count == 1) // vertex
      iMesh_createEnt(impl_, iMesh_LINE_SEGMENT, nodes, 2, &out, &status, &err);
    else
      throw MKException(iBase_FAILURE, "Couldn't extrude face; unusual shape.");

    check_error(impl_, err);

    iMesh_setEHData(impl_, src[i], local_tag, out, &err);
    check_error(impl_, err);
    delete[] nodes;
  }

  process_ce_sets(impl_, extrude_sets_.sets(), local_tag);
}

void ExtrudeMesh::tag_all_sets(iBase_TagHandle local_tag)
{
  int err;
  std::set<iBase_EntitySetHandle>::iterator set;
  std::vector<CESets::tag_data>::iterator tag;

  // set the target sets for expand sets to be themselves
  for (set = expand_sets_.sets().begin(); set != expand_sets_.sets().end();
       ++set) {
    iMesh_setEntSetEHData(impl_, *set, local_tag,
                          reinterpret_cast<iBase_EntityHandle>(*set), &err);
    check_error(impl_, err);
  }

  process_ce_sets(impl_, copy_sets_.sets(), local_tag);
  process_ce_sets(impl_, expand_sets_.sets(), local_tag);

  // set the copy tag on all copied sets
  for (set = copy_sets_.sets().begin(); set != copy_sets_.sets().end(); ++set) {
    iBase_EntityHandle eh;
    iMesh_getEntSetEHData(impl_, *set, local_tag, &eh, &err);
    if (err == iBase_SUCCESS) {
      iMesh_setEntSetEHData(impl_, *set, copy_tag_, eh, &err);
      check_error(impl_, err);
    }
  }

  // set the extrude tag on all extruded sets
  for(set = extrude_sets_.sets().begin(); set != extrude_sets_.sets().end();
      ++set) {
    iBase_EntityHandle eh;
    iMesh_getEntSetEHData(impl_, *set, local_tag, &eh, &err);
    if (err == iBase_SUCCESS) {
      iMesh_setEntSetEHData(impl_, *set, extrude_tag_, eh, &err);
      check_error(impl_, err);
    }
  }

  // tag the newly-created sets
  for(tag = copy_sets_.tags().begin(); tag != copy_sets_.tags().end(); ++tag) {
    tag_copy_sets(impl_, copy_tag_, copy_sets_.sets(), tag->tag,
                  tag->value);
  }

  for(tag = extrude_sets_.tags().begin(); tag != extrude_sets_.tags().end();
      ++tag) {
    tag_copy_sets(impl_, extrude_tag_, extrude_sets_.sets(), tag->tag,
                  tag->value);
  }
}


#ifdef TEST
#include <cmath>
#include <iostream>

#define CHECK_ERR(str) do {                       \
    if (err != iBase_SUCCESS) {                   \
      std::cerr << (str) << std::endl;            \
      return 1;                                   \
    }                                             \
  } while(false)

#define CHECK_THROW(expr) do {                    \
    try {                                         \
      expr;                                       \
    }                                             \
    catch(const std::exception &e) {              \
      std::cerr << e.what() << std::endl;         \
      return 1;                                   \
    }                                             \
  } while(false)

#define RUN_TEST(function) do {                   \
    std::cerr << "  " << #function << "... ";     \
    if(function() == 0) {                         \
      std::cerr << "passed" << std::endl;         \
    }                                             \
    else {                                        \
      std::cerr << "failed" << std::endl;         \
      result = 1;                                 \
    }                                             \
  } while(false)

int test_translate_ents()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  CHECK_ERR("Couldn't create mesh.");

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  CHECK_ERR("Couldn't get root set.");

  ExtrudeMesh *ext = new ExtrudeMesh(mesh);

  double coords[] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,
    2, 1, 0,
    9, 9, 9,
  };

  iBase_EntityHandle *verts = NULL;
  int alloc=0, size;
  iMesh_createVtxArr(mesh, 6, iBase_INTERLEAVED, coords, 3*6, &verts, &alloc,
                     &size, &err);
  CHECK_ERR("Couldn't create vertex array");

  iBase_EntityHandle quad;
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts, 4, &quad, &status, &err);
  CHECK_ERR("Couldn't create entity");

  iBase_EntityHandle line;
  iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, verts, 2, &line, &status, &err);
  CHECK_ERR("Couldn't create entity");

  iBase_EntityHandle tri;
  iBase_EntityHandle tri_verts[] = { verts[1], verts[2], verts[4] };
  iMesh_createEnt(mesh, iMesh_TRIANGLE, tri_verts, 3, &tri, &status, &err);
  CHECK_ERR("Couldn't create entity");

  free(verts);

  iBase_EntityHandle faces[] = {quad, tri, line};
  double v[] = { 0, 0, 5 };
  int steps = 50;
  CHECK_THROW( ext->translate(faces, 3, steps, v) );

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  CHECK_ERR("Couldn't get number of vertices.");
  if(count != 5*(steps+1)+1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  if(count != 2+1*steps)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  if(count != steps*2)
    return 1;

#ifdef TESTSAVE
  // VisIt doesn't like 1-d objects in pseudocolor volume graphs
  iMesh_deleteEnt(mesh, line, &err);
  assert(err == 0);

  const char *file = "test1.vtk";
  iMesh_save(mesh, root_set, file, "", &err, strlen(file), 0);
  assert(err == 0);
#endif

  delete ext;
  iMesh_dtor(mesh, &err);
  CHECK_ERR("Couldn't destroy mesh.");

  return 0;
}

int test_translate_ents_to_dest()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  CHECK_ERR("Couldn't create mesh.");

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  CHECK_ERR("Couldn't get root set.");

  ExtrudeMesh *ext = new ExtrudeMesh(mesh);

  double coords[] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,

    2, 1, 0,
    9, 9, 9,
    0, 0, 1,
    1, 0, 1,

    1, 1, 1,
    0, 1, 1,
    2, 1, 1,
  };
  iBase_EntityHandle *verts = NULL;
  int alloc = 0, size;

  iMesh_createVtxArr(mesh, 11, iBase_INTERLEAVED, coords, 3*11, &verts, &alloc,
                     &size, &err);
  CHECK_ERR("Couldn't create vertex array");

  iBase_EntityHandle quad[2];
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts+0, 4, quad+0, &status, &err);
  CHECK_ERR("Couldn't create entity");

  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts+6, 4, quad+1, &status, &err);
  CHECK_ERR("Couldn't create entity");

  iBase_EntityHandle tri[2];
  iBase_EntityHandle tri_verts[] = { verts[1], verts[2], verts[4],
                                     verts[7], verts[8], verts[10] };
  iMesh_createEnt(mesh, iMesh_TRIANGLE, tri_verts+0, 3, tri+0, &status, &err);
  CHECK_ERR("Couldn't create entity");

  iMesh_createEnt(mesh, iMesh_TRIANGLE, tri_verts+3, 3, tri+1, &status, &err);
  CHECK_ERR("Couldn't create entity");

  free(verts);

  iBase_EntityHandle pre[]  = {quad[0], tri[0]};
  iBase_EntityHandle post[] = {quad[1], tri[1]};
  int steps = 5;
  CHECK_THROW( ext->translate(pre, post, 2, steps) );

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  CHECK_ERR("Couldn't get number of vertices.");
  if(count != 5*(steps+1)+1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  if(count != 4)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  if(count != steps*2)
    return 1;

#ifdef TESTSAVE
  const char *file = "test2.vtk";
  iMesh_save(mesh, root_set, file, "", &err, strlen(file), 0);
  assert(err == 0);
  #endif

  delete ext;
  iMesh_dtor(mesh, &err);
  CHECK_ERR("Couldn't destroy mesh.");

  return 0;
}

int test_translate_set()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  CHECK_ERR("Couldn't create mesh.");

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  CHECK_ERR("Couldn't get root set.");

  ExtrudeMesh *ext = new ExtrudeMesh(mesh);

  double coords[] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,
  };
  iBase_EntityHandle *verts = NULL;
  int alloc = 0, size;

  iMesh_createVtxArr(mesh, 4, iBase_INTERLEAVED, coords, 3*4, &verts, &alloc,
                     &size, &err);
  CHECK_ERR("Couldn't create vertex array");

  iBase_EntityHandle quad;
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts, 4, &quad, &status, &err);
  CHECK_ERR("Couldn't create entity");

  free(verts);

  iBase_EntityHandle faces[] = { quad };
  double v[] = { 0, 0, 5 };
  int steps = 5;
  CHECK_THROW( ext->translate(faces, 1, steps, v) );

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  CHECK_ERR("Couldn't get number of vertices.");
  if(count != 4*(steps+1))
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  if(count != 1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  if(count != steps)
    return 1;

#ifdef TESTSAVE
  const char *file = "test3.vtk";
  iMesh_save(mesh, root_set, file, "", &err, strlen(file), 0);
  assert(err == 0);
#endif

  delete ext;
  iMesh_dtor(mesh, &err);
  CHECK_ERR("Couldn't destroy mesh.");

  return 0;
}

int test_rotate_set()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  CHECK_ERR("Couldn't create mesh.");

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  CHECK_ERR("Couldn't get root set.");

  ExtrudeMesh *ext = new ExtrudeMesh(mesh);

  double coords[] = {
    0, 0, 0,
    0, 1, 0,
    0, 1, 1,
    0, 0, 1,
  };
  iBase_EntityHandle *verts = NULL;
  int alloc = 0, size;

  iMesh_createVtxArr(mesh, 4, iBase_INTERLEAVED, coords, 3*4, &verts, &alloc,
                     &size, &err);
  CHECK_ERR("Couldn't create vertex array");

  iBase_EntityHandle quad;
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts, 4, &quad, &status, &err);
  CHECK_ERR("Couldn't create entity");

  free(verts);

  iBase_EntityHandle faces[] = { quad };

  iBase_EntitySetHandle set;
  iMesh_createEntSet(mesh, false, &set, &err);
  CHECK_ERR("Couldn't create entity set");
  
  iMesh_addEntArrToSet(mesh, faces, 1, set, &err);
  CHECK_ERR("Couldn't add entities to set");

  iBase_TagHandle tag;
  iMesh_createTag(mesh, "my_tag", 1, iBase_BYTES, &tag, &err, 6);
  CHECK_ERR("Couldn't create tag");
  iMesh_setEntSetData(mesh, set, tag, "x", 1, &err);
  CHECK_ERR("Couldn't set tag data on set.");
  CHECK_THROW( ext->extrude_sets().add_tag(tag, "x") );

  int steps = 200;
  double origin[] = { 0, -3, 0 };
  double z[] = { 1, 1, 1 };
  double angle = 2*3.14159;
  CHECK_THROW( ext->rotate(set, steps, origin, z, angle) );

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  CHECK_ERR("Couldn't get number of vertices.");
  if(count != 4*(steps+1))
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  if(count != 1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  if(count != steps)
    return 1;

  iMesh_getNumEntSets(mesh, root_set, 0, &count, &err);
  CHECK_ERR("Couldn't get number of entity sets.");
  if(count != 2)
    return 1;

#ifdef TESTSAVE
  const char *file = "test4.vtk";
  iMesh_save(mesh, root_set, file, "", &err, strlen(file), 0);
  assert(err == 0);
#endif

  delete ext;
  iMesh_dtor(mesh, &err);
  CHECK_ERR("Couldn't destroy mesh.");

  return 0;
}

int main()
{
  int result = 0;

  RUN_TEST(test_translate_ents);
  RUN_TEST(test_translate_ents_to_dest);
  RUN_TEST(test_translate_set);
  RUN_TEST(test_rotate_set);

  return result;
}

#endif
