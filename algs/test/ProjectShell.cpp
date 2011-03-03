#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ProjectShell.hpp"

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#define DEFAULT_TEST_FILE STRINGIFY(SRCDIR) "/3k-tri-cube.h5m"
#define DEFAULT_OUT_FILE "ProjShell_3k-tri-cube.h5m"

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

int projectOnDirection(const char *filename_mesh, double direction[3], const char * newFile);

int main(int argc, char* argv[] )
{
  // check command line arg
  const char *filename_mesh = 0;
  double direction[3] = {0., 0., 1.};
  const char *newFile = 0;
  if (argc == 6)
    {
      filename_mesh = argv[1];
      direction[0] = atof(argv[2]);
      direction[1] = atof(argv[3]);
      direction[2] = atof(argv[4]);
      // do project mesh
      // new projectShell mesh
      newFile = argv[5];
    }
  else {
    printf("Usage: %s <mesh_filename> <dir_x> <dir_y> <dir_z> <newFile>\n", argv[0]);
    if (argc != 1) return 1;
    printf("No file specified.  Defaulting to: %s %g %g %g %s\n",
	   DEFAULT_TEST_FILE, direction[0], direction[1], direction[2], DEFAULT_OUT_FILE);
    filename_mesh = DEFAULT_TEST_FILE;
    newFile = DEFAULT_OUT_FILE;
  }
  
  if (projectOnDirection(filename_mesh, direction, newFile)  ) return 1;
  
  return 0;
}

int projectOnDirection(const char *filename_mesh, double direction[3], const char *newFile )
{
  // initialize the Mesh
  int err;
  iMesh_Instance mesh;
  iMesh_newMesh("", &mesh, &err, 0);
  ERRORR("Couldn't create mesh.", 1);

  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

  // read geometry and establish facet mesh
  clock_t start_time = clock();
  iMesh_load(mesh, root_set, filename_mesh, NULL, &err, strlen(filename_mesh), 0);
  ERRORR("Couldn't load mesh file.", 1);
  clock_t load_time = clock();


  // make project shell mesher
  ProjectShell *psm = new ProjectShell(mesh, root_set, direction);

  
  err = psm->project();
  ERRORR("Couldn't project shell mesh.", 1);
  // write the output mesh 
 
  iMesh_Instance newMesh;
  iMesh_newMesh("", &newMesh, &err, 0);
  ERRORR("Couldn't create new mesh.", 1);

  iBase_EntitySetHandle new_root_set;
  iMesh_getRootSet(newMesh, &new_root_set, &err);
  
  ERRORR("Couldn't get new root set.", 1);

  err = psm->writeNewMesh(newMesh);
  ERRORR ("Couldn't create new mesh.", 1);
 
  iMesh_save(newMesh, new_root_set, newFile, NULL, &err, strlen(newFile), 0);

  clock_t mesh_time = clock();

  std::cout << "Project Shell mesh is succesfully produced." << std::endl;
  std::cout << "Time including loading is "
	    << (double) (mesh_time - start_time)/CLOCKS_PER_SEC
	    << " secs, Time excluding loading is "
	    << (double) (mesh_time - load_time)/CLOCKS_PER_SEC
	    << std::endl;

  delete psm;
  iMesh_dtor(mesh, &err);
  ERRORR("Couldn't destroy mesh.", 1);

  return 0;
}
