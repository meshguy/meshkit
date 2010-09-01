#include "ExtrudeMesh.hpp"
#include "TestFramework.hpp"
#include <cmath>

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
  CHECK_THROW( ext->extrude(faces, 3, extrude::Translate(v, steps)) );

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  CHECK_ERR("Couldn't get number of vertices.");
  CHECK_TRUE(count == 5*(steps+1)+1, "Invalid number of vertices.");

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  CHECK_TRUE(count == 2+1*steps, "Invalid number of faces.");

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  CHECK_TRUE(count == steps*2, "Invalid number of regions.");

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
  CHECK_TRUE(count == 5*(steps+1)+1, "Invalid number of vertices.");

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  CHECK_TRUE(count == 4, "Invalid number of faces.");

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  CHECK_TRUE(count == steps*2, "Invalid number of regions.");

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
  CHECK_THROW( ext->extrude(faces, 1, extrude::Translate(v, steps)) );

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  CHECK_ERR("Couldn't get number of vertices.");
  CHECK_TRUE(count == 4*(steps+1), "Invalid number of vertices.");

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  CHECK_TRUE(count == 1, "Invalid number of faces.");

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  CHECK_TRUE(count == steps, "Invalid number of regions.");

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
  CHECK_THROW( ext->extrude(set, extrude::Rotate(origin, z, angle, steps)) );

  int count;
  iMesh_getNumOfType(mesh, 0, iBase_VERTEX, &count, &err);
  CHECK_ERR("Couldn't get number of vertices.");
  CHECK_TRUE(count == 4*(steps+1), "Invalid number of vertices.");

  iMesh_getNumOfType(mesh, 0, iBase_FACE, &count, &err);
  CHECK_ERR("Couldn't get number of faces.");
  CHECK_TRUE(count == 1, "Invalid number of faces.");

  iMesh_getNumOfType(mesh, 0, iBase_REGION, &count, &err);
  CHECK_ERR("Couldn't get number of regions.");
  CHECK_TRUE(count == steps, "Invalid number of regions.");

  iMesh_getNumEntSets(mesh, root_set, 0, &count, &err);
  CHECK_ERR("Couldn't get number of entity sets.");
  CHECK_TRUE(count == 2, "Invalid number of entity sets.");

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