#include "CutCellMesh.hpp"
#include "iGeom.h"

#include <iostream>

#define DEFAULT_TEST_FILE_GEOM "sphere.stl"
#define DEFAULT_TEST_FILE_MESH "sphere.stl"
#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

int load_and_cutcell_test(const char *filename_geom, const char *filename_mesh);

int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename_geom = 0;
  const char *filename_mesh = 0;

  if (argc == 3) {
    filename_geom = argv[1];
    filename_mesh = argv[2];
  }
  else {
    printf("Usage: %s <geom_filename> <mesh_filename>\n", argv[0]);
    if (argc != 1) return 1;
    printf("  No file specified.  Defaulting to: %s %s\n", DEFAULT_TEST_FILE_GEOM,
	   DEFAULT_TEST_FILE_MESH);
    filename_geom = DEFAULT_TEST_FILE_GEOM;
    filename_mesh = DEFAULT_TEST_FILE_MESH;
  }

  if (load_and_cutcell_test(filename_geom, filename_mesh)) return 1;

  return 0;
}

int load_and_cutcell_test(const char *filename_geom, const char *filename_mesh)
{/*
  // initialize the Geom
  int err;
  iGeom_Instance geom;
  iGeom_newGeom("", &geom, &err, 0);

  // load geom file
  iGeom_load(geom, filename_geom, NULL, &err, strlen(filename_geom), 0);
  ERRORR("can not load a geometry.", 1);
 */
  // initialize the Mesh
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

  // load a mesh file, just surface triangle now
  // should be changed
  iMesh_load(mesh, root_set, filename_mesh, NULL, &err, strlen(filename_mesh), 0);
  ERRORR("Couldn't load mesh file.", 1);

  // make cut cell mesher
  CutCellMesh *ccm = new CutCellMesh(mesh, root_set);

  // get initial cartesian cell division
  err = ccm->do_mesh();
  ERRORR("Couldn't cut-cell mesh.", 1);

  std::cerr << "Cut-cell mesh is succesfully produced." << std::endl;

  delete ccm;
  iMesh_dtor(mesh, &err);
  ERRORR("Couldn't destroy mesh.", 1);

  return 0;
}
