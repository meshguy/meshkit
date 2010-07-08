#include "ExtrudeMesh.hpp"

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
  ERROR("Couldn't get vertex coords");
  iMesh_getVtxCoord(mesh,b,xb+0,xb+1,xb+2,&err);
  ERROR("Couldn't get vertex coords");

  for(int i=0; i<3; i++)
    res[i] = xa[i]-xb[i];
  return res;
}


ExtrudeMesh::ExtrudeMesh(iMesh_Instance mesh)
  : impl_(mesh), updated_set_(false),
    copy_tag_   (mesh, "__CopyMeshTag"),
    extrude_tag_(mesh, "__ExtrudeMeshTag")
{}

ExtrudeMesh::~ExtrudeMesh()
{
  std::vector<tag_data>::iterator i;
  for(i = extrude_tags_.begin(); i != extrude_tags_.end(); ++i)
    free(i->value);
}

int ExtrudeMesh::add_copy_tag(const std::string &tag_name,
                              const char *tag_val)
{
  iBase_TagHandle tag_handle = 0;
  int err;
  iMesh_getTagHandle(impl_, tag_name.c_str(), &tag_handle, &err,
                     tag_name.length());
  if(err != iBase_SUCCESS)
    ERRORR("Failed to get handle for tag "+tag_name, iBase_FAILURE);

  return add_copy_tag(tag_handle, tag_val);
}

int ExtrudeMesh::add_copy_tag(iBase_TagHandle tag_handle,
                              const char *tag_val)
{
  char *tmp = NULL;

  if(tag_val) {
    int err;
    int tag_size;
    iMesh_getTagSizeBytes(impl_, tag_handle, &tag_size, &err);
    if (err != iBase_SUCCESS)
      ERRORR("Failed to get size of tag", iBase_FAILURE);
    tmp = static_cast<char*>(malloc(tag_size));
    memcpy(tmp, tag_val, tag_size);
  }

  copy_tags_.push_back(tag_data(tag_handle, tmp));
  return iBase_SUCCESS;
}

int ExtrudeMesh::add_expand_tag(const std::string &tag_name,
                                const char *tag_val)
{
  iBase_TagHandle tag_handle = 0;
  int err;
  iMesh_getTagHandle(impl_, tag_name.c_str(), &tag_handle, &err,
                     tag_name.length());
  if(err != iBase_SUCCESS)
    ERRORR("Failed to get handle for tag "+tag_name, iBase_FAILURE);

  return add_expand_tag(tag_handle, tag_val);
}

int ExtrudeMesh::add_expand_tag(iBase_TagHandle tag_handle,
                                const char *tag_val)
{
  char *tmp = NULL;

  if(tag_val) {
    int err;
    int tag_size;
    iMesh_getTagSizeBytes(impl_, tag_handle, &tag_size, &err);
    if (err != iBase_SUCCESS)
      ERRORR("Failed to get size of tag", iBase_FAILURE);
    tmp = static_cast<char*>(malloc(tag_size));
    memcpy(tmp, tag_val, tag_size);
  }

  expand_tags_.push_back(tag_data(tag_handle, tmp));
  return iBase_SUCCESS;
}

int ExtrudeMesh::add_extrude_tag(const std::string &tag_name,
                                 const char *tag_val)
{
  iBase_TagHandle tag_handle = 0;
  int err;
  iMesh_getTagHandle(impl_, tag_name.c_str(), &tag_handle, &err,
                     tag_name.length());
  if(err != iBase_SUCCESS)
    ERRORR("Failed to get handle for tag "+tag_name, iBase_FAILURE);

  return add_extrude_tag(tag_handle, tag_val);
}

int ExtrudeMesh::add_extrude_tag(iBase_TagHandle tag_handle,
                                 const char *tag_val)
{
  char *tmp = NULL;

  if(tag_val) {
    int err;
    int tag_size;
    iMesh_getTagSizeBytes(impl_, tag_handle, &tag_size, &err);
    if (err != iBase_SUCCESS)
      ERRORR("Failed to get size of tag", iBase_FAILURE);
    tmp = static_cast<char*>(malloc(tag_size));
    memcpy(tmp, tag_val, tag_size);
  }

  extrude_tags_.push_back(tag_data(tag_handle, tmp));
  return iBase_SUCCESS;
}

int ExtrudeMesh::update_sets()
{
  if(updated_set_)
    reset_sets();

  int err;
  iBase_EntitySetHandle root;
  iMesh_getRootSet(impl_, &root, &err);
  ERRORR("Trouble getting root set", iBase_FAILURE);
  
  for(std::vector<tag_data>::iterator i=copy_tags_.begin();
      i!=copy_tags_.end(); ++i) {
    SimpleArray<iBase_EntitySetHandle> tmp_sets;
    iMesh_getEntSetsByTagsRec(impl_, root, &i->tag,
                              (i->value ? &i->value:NULL), 1, 0,
                              ARRAY_INOUT(tmp_sets), &err);
    ERRORR("Couldn't get tagged sets.", iBase_FAILURE);

    copy_sets_.insert(tmp_sets.begin(), tmp_sets.end());
  }

  for(std::vector<tag_data>::iterator i=expand_tags_.begin();
      i!=expand_tags_.end(); ++i) {
    SimpleArray<iBase_EntitySetHandle> tmp_sets;
    iMesh_getEntSetsByTagsRec(impl_, root, &i->tag,
                              (i->value ? &i->value:NULL), 1, 0,
                              ARRAY_INOUT(tmp_sets), &err);
    ERRORR("Couldn't get tagged sets.", iBase_FAILURE);

    expand_sets_.insert(tmp_sets.begin(), tmp_sets.end());
  }

  for(std::vector<tag_data>::iterator i=extrude_tags_.begin();
      i!=extrude_tags_.end(); ++i) {
    SimpleArray<iBase_EntitySetHandle> tmp_sets;
    iMesh_getEntSetsByTagsRec(impl_, root, &i->tag,
                              (i->value ? &i->value:NULL), 1, 0,
                              ARRAY_INOUT(tmp_sets), &err);
    ERRORR("Couldn't get tagged sets.", iBase_FAILURE);

    extrude_sets_.insert(tmp_sets.begin(), tmp_sets.end());
  }

  updated_set_ = true;
  return iBase_SUCCESS;
}

int ExtrudeMesh::reset_sets()
{
  copy_sets_.clear();
  expand_sets_.clear();
  extrude_sets_.clear();
  updated_set_ = false;

  return iBase_SUCCESS;
}

int ExtrudeMesh::translate(iBase_EntityHandle *src, int size, int steps,
                           const double *dx, bool copy_faces)
{
  return extrude(src, size, steps, CopyMoveVerts(impl_, dx, steps), copy_faces);
}

int ExtrudeMesh::translate(iBase_EntitySetHandle src, int steps,
                           const double *dx, bool copy_faces)
{
  return extrude(src, steps, CopyMoveVerts(impl_, dx, steps) ,copy_faces);
}

int ExtrudeMesh::rotate(iBase_EntityHandle *src, int size, int steps,
                        const double *origin, const double *z, double angle,
                        bool copy_faces)
{
  return extrude(src, size, steps,
                 CopyRotateVerts(impl_, origin, z, angle, steps), copy_faces);
}

int ExtrudeMesh::rotate(iBase_EntitySetHandle src, int steps,
                        const double *origin, const double *z, double angle,
                        bool copy_faces)
{
  return extrude(src, steps, CopyRotateVerts(impl_, origin, z, angle),
                 copy_faces);
}

int ExtrudeMesh::extrude(iBase_EntityHandle *src, int size, int steps,
                         const CopyVerts &trans, bool copy_faces)
{
  int err;
  iBase_EntitySetHandle set;
  iMesh_createEntSet(impl_, false, &set, &err);
  ERRORR("Couldn't create source entity set.", err);
  
  iMesh_addEntArrToSet(impl_, src, size, set, &err);
  ERRORR("Couldn't add entities to source entity set.", err);

  int ret = extrude(set, steps, trans, copy_faces);

  iMesh_destroyEntSet(impl_, set, &err);
  ERRORR("Couldn't destroy source entity set.", err);

  return ret;
}

int ExtrudeMesh::extrude(iBase_EntitySetHandle src, int steps,
                         const CopyVerts &trans, bool copy_faces)
{
  assert(steps > 0);

  if(!updated_set_)
    update_sets();

  int err;

  SimpleArray<iBase_EntityHandle> ents;
  SimpleArray<iBase_EntityHandle> adj;
  SimpleArray<int> indices;
  SimpleArray<int> offsets;

  iMesh_getStructure(impl_, src, ARRAY_INOUT(ents), ARRAY_INOUT(adj),
                     ARRAY_INOUT(indices), ARRAY_INOUT(offsets), &err);
  ERRORR("Trouble getting source adjacencies.", err);

  double dx[3];
  int *normals = 0;

  LocalTag local_copy_tag(impl_);
  LocalTag local_extrude_tag(impl_);

  SimpleArray<iBase_EntityHandle> curr(adj.size());
  SimpleArray<iBase_EntityHandle> next(adj.size());
  trans(1, ARRAY_IN(adj), ARRAY_INOUT(next));

  vtx_diff(dx, impl_, next[0], adj[0]);
  normals = get_normals(&adj[0], &indices[0], &offsets[0], ents.size(), dx);

  // Make the first set of volumes
  connect_higher_dots(ARRAY_IN(ents), local_extrude_tag, &normals[0],
                      &indices[0], &offsets[0], &adj[0], &next[0]);

  // Make the inner volumes
  for(int i=2; i<=steps; i++) {
    std::swap(curr, next);
    trans(i, ARRAY_IN(adj), ARRAY_INOUT(next));
    connect_higher_dots(ARRAY_IN(ents), local_extrude_tag, &normals[0],
                        &indices[0], &offsets[0], &curr[0], &next[0]);
  }

  if(copy_faces)
    connect_the_dots(impl_, ARRAY_IN(ents), local_copy_tag,
                     &indices[0], &offsets[0], &next[0]);

  free(normals);

  err = process_ce_sets(impl_, copy_sets_, local_copy_tag);
  ERRORR("Failed to update copy/expand sets.", iBase_FAILURE);
  
  err = process_ce_sets(impl_, expand_sets_, local_extrude_tag);
  ERRORR("Failed to update expand/expand sets.", iBase_FAILURE);

  std::set<iBase_EntitySetHandle>::iterator set;
  // set the copy tag on all copied sets
  for(set = copy_sets_.begin(); set != copy_sets_.end(); ++set) {
    iBase_EntityHandle eh;
    iMesh_getEntSetEHData(impl_, *set, local_copy_tag, &eh, &err);
    if(err == iBase_SUCCESS) {
      iMesh_setEntSetEHData(impl_, *set, copy_tag_, eh, &err);
      ERRORR("Failed to tag extruded set with extrude tag.", iBase_FAILURE);
    }
  }

  // set the extrude tag on all extruded sets
  for(set = extrude_sets_.begin(); set != extrude_sets_.end(); ++set) {
    iBase_EntityHandle eh;
    iMesh_getEntSetEHData(impl_, *set, local_extrude_tag, &eh, &err);
    if(err == iBase_SUCCESS) {
      iMesh_setEntSetEHData(impl_, *set, extrude_tag_, eh, &err);
      ERRORR("Failed to tag extruded set with extrude tag.", iBase_FAILURE);
    }
  }

  std::vector<tag_data>::iterator tag;
  for(tag = copy_tags_.begin(); tag != copy_tags_.end(); ++tag) {
    err = tag_copy_sets(impl_, copy_tag_, copy_sets_, tag->tag,
                        tag->value);
    ERRORR("Failed to tag copied sets.", iBase_FAILURE);
  }

  for(tag = extrude_tags_.begin(); tag != extrude_tags_.end(); ++tag) {
    err = tag_copy_sets(impl_, extrude_tag_, extrude_sets_, tag->tag,
                        tag->value);
    ERRORR("Failed to tag copied sets.", iBase_FAILURE);
  }

  return iBase_SUCCESS;
}

// calculate the normals for each face (1 = towards v, -1 = away from v)
// TODO: this can fail with non-convex faces
int * ExtrudeMesh::get_normals(iBase_EntityHandle *verts, int *indices,
                               int *offsets, int size, double *dv)
{
  int err;
  int *normals = (int*)malloc(size*sizeof(int));

  for(int i=0; i<size; i++) {
    double res[3], a[3], b[3];
    iBase_EntityHandle curr_verts[3];

    if(offsets[i+1] - offsets[i] > 2) { // face
      for(int j=0; j<3; j++)
        curr_verts[j] = verts[indices[ offsets[i]+j ]];

      SimpleArray<double> coords;
      iMesh_getVtxArrCoords(impl_, curr_verts, 3, iBase_INTERLEAVED,
                            ARRAY_INOUT(coords), &err);
      ERROR("Couldn't get vertex coordinates.");

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

void ExtrudeMesh::connect_higher_dots(iBase_EntityHandle *src, int size,
                                      iBase_TagHandle local_extrude_tag,
                                      int *normals, int *indices, int *offsets,
                                      iBase_EntityHandle *pre,
                                      iBase_EntityHandle *post)
{
  int err;

  for(int i=0; i<size; i++) {
    int count = offsets[i+1] - offsets[i];

    // If the normal is facing in the wrong direction (away from the
    // translation) we add the vertices in reverse order. Otherwise, we go
    // in the usual order. If count is 2, then we are creating quads and so
    // need to swap the order of the post set of verts.

    int dx = normals[i];
    int dy = normals[i] * (count == 2 ? -1:1);
    int x  = (dx == 1) ? offsets[i] : offsets[i+1]-1;
    int y  = (dy == 1) ? offsets[i] : offsets[i+1]-1;

    iBase_EntityHandle *nodes = new iBase_EntityHandle[count*2];
    for(int j=0; j<count; j++) {
      nodes[j]       = pre [ indices[x + dx*j] ];
      nodes[j+count] = post[ indices[y + dy*j] ];
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
      std::cerr << "Couldn't extrude face; unusual shape." << std::endl;

    ERROR("Couldn't create extruded face.");

    iMesh_setEHData(impl_, src[i], local_extrude_tag, out, &err);
    ERROR("Couldn't set local tag data.");
    delete[] nodes;
  }

  process_ce_sets(impl_, extrude_sets_, local_extrude_tag);
}

#ifdef TEST
#include <cmath>

int test1()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

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
  ERRORR("Couldn't create vertex array", 1);

  iBase_EntityHandle quad;
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts, 4, &quad, &status, &err);
  ERRORR("Couldn't create entity", 1);

  iBase_EntityHandle line;
  iMesh_createEnt(mesh, iMesh_LINE_SEGMENT, verts, 2, &line, &status, &err);
  ERRORR("Couldn't create entity", 1);

  iBase_EntityHandle tri;
  iBase_EntityHandle tri_verts[] = { verts[1], verts[2], verts[4] };
  iMesh_createEnt(mesh, iMesh_TRIANGLE, tri_verts, 3, &tri, &status, &err);
  ERRORR("Couldn't create entity", 1);

  free(verts);

  iBase_EntityHandle faces[] = {quad, tri, line};
  double v[] = { 0, 0, 5 };
  int steps = 50;
  err = ext->translate(faces, 3, steps, v);
  ERRORR("Couldn't extrude mesh", 1);

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  ERRORR("Couldn't get number of vertices.", 1);
  if(count != 5*(steps+1)+1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  ERRORR("Couldn't get number of faces.", 1);
  if(count != 2+1*steps)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  ERRORR("Couldn't get number of regions.", 1);
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
  ERRORR("Couldn't destroy mesh.", 1);

  return 0;
}

int test2()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

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
  ERRORR("Couldn't create vertex array", 1);

  iBase_EntityHandle quad[2];
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts+0, 4, quad+0, &status, &err);
  ERRORR("Couldn't create entity", 1);

  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts+6, 4, quad+1, &status, &err);
  ERRORR("Couldn't create entity", 1);

  iBase_EntityHandle tri[2];
  iBase_EntityHandle tri_verts[] = { verts[1], verts[2], verts[4],
                                     verts[7], verts[8], verts[10] };
  iMesh_createEnt(mesh, iMesh_TRIANGLE, tri_verts+0, 3, tri+0, &status, &err);
  ERRORR("Couldn't create entity",1);

  iMesh_createEnt(mesh, iMesh_TRIANGLE, tri_verts+3, 3, tri+1, &status, &err);
  ERRORR("Couldn't create entity",1);

  free(verts);

  iBase_EntityHandle pre[]  = {quad[0], tri[0]};
  iBase_EntityHandle post[] = {quad[1], tri[1]};
  int steps = 5;
//  err = ext->translate(pre, post, 2, steps);
  ERRORR("Couldn't extrude mesh", 1);

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  ERRORR("Couldn't get number of vertices.", 1);
  if(count != 5*(steps+1)+1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  ERRORR("Couldn't get number of faces.", 1);
  if(count != 4)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  ERRORR("Couldn't get number of regions.", 1);
  if(count != steps*2)
    return 1;

#ifdef TESTSAVE
  const char *file = "test2.vtk";
  iMesh_save(mesh, root_set, file, "", &err, strlen(file), 0);
  assert(err == 0);
  #endif

  delete ext;
  iMesh_dtor(mesh, &err);
  ERRORR("Couldn't destroy mesh.", 1);

  return 0;
}

int test3()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

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
  ERRORR("Couldn't create vertex array",1);

  iBase_EntityHandle quad;
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts, 4, &quad, &status, &err);
  ERRORR("Couldn't create entity", 1);

  free(verts);

  iBase_EntityHandle faces[] = { quad };
  double v[] = { 0, 0, 5 };
  int steps = 5;
  err = ext->translate(faces, 1, steps, v);
  ERRORR("Couldn't extrude mesh", 1);

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  ERRORR("Couldn't get number of vertices.", 1);
  if(count != 4*(steps+1))
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  ERRORR("Couldn't get number of faces.", 1);
  if(count != 1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  ERRORR("Couldn't get number of regions.", 1);
  if(count != steps)
    return 1;

#ifdef TESTSAVE
  const char *file = "test3.vtk";
  iMesh_save(mesh, root_set, file, "", &err, strlen(file), 0);
  assert(err == 0);
#endif

  delete ext;
  iMesh_dtor(mesh, &err);
  ERRORR("Couldn't destroy mesh.", 1);

  return 0;
}

int test4()
{
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

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
  ERRORR("Couldn't create vertex array", 1);

  iBase_EntityHandle quad;
  int status;
  iMesh_createEnt(mesh, iMesh_QUADRILATERAL, verts, 4, &quad, &status, &err);
  ERRORR("Couldn't create entity", 1);

  free(verts);

  iBase_EntityHandle faces[] = { quad };

  iBase_EntitySetHandle set;
  iMesh_createEntSet(mesh, false, &set, &err);
  ERRORR("Couldn't create entity set", 1);
  
  iMesh_addEntArrToSet(mesh, faces, 1, set, &err);
  ERRORR("Couldn't add entities to set", 1);

  iBase_TagHandle tag;
  iMesh_createTag(mesh, "my_tag", 1, iBase_BYTES, &tag, &err, 6);
  iMesh_setEntSetData(mesh, set, tag, "x", 1, &err);
  ext->add_extrude_tag(tag, "x");

  int steps = 200;
  double origin[] = { 0, -3, 0 };
  double z[] = { 1, 1, 1 };
  double angle = 2*3.14159;
  err = ext->rotate(set, steps, origin, z, angle);
  ERRORR("Couldn't extrude mesh", 1);

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  ERRORR("Couldn't get number of vertices.", 1);
  if(count != 4*(steps+1))
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  ERRORR("Couldn't get number of faces.", 1);
  if(count != 1)
    return 1;

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  ERRORR("Couldn't get number of regions.", 1);
  if(count != steps)
    return 1;

  iMesh_getNumEntSets(mesh, root_set, 0, &count, &err);
  ERRORR("Couldn't get number of entity sets.", 1);
  if(count != 2)
    return 1;

#ifdef TESTSAVE
  const char *file = "test4.vtk";
  iMesh_save(mesh, root_set, file, "", &err, strlen(file), 0);
  assert(err == 0);
#endif

  delete ext;
  iMesh_dtor(mesh, &err);
  ERRORR("Couldn't destroy mesh.", 1);

  return 0;
}

int main()
{
  if(test1())
    return 1;
//  if(test2())
//    return 1;
  if(test3())
    return 1;
  if(test4())
    return 1;

  return 0;
}

#endif
