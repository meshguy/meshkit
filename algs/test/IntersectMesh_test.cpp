/*
 * IntersectMesh_test.cpp
 *
 *  This test is for intersection of 2d meshes, on top of each other, for
 *   a domain decomposition algorithm; it leverages the method from ProjectShell algorithm.
 *
 *  inputs are 2 meshes, vtk format, output will be another mesh, m3, with the common part
 *    intersected
 *
 *  Created on: Aug 19, 2010
 *      Author: iulian
 */

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "IntersectMesh.hpp"
#include "TestFramework.hpp"

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

int main(int argc, char* argv[]) {
   // check command line arg
   const char *filename_mesh1 = 0;
   const char *filename_mesh2 = 0;
   const char *newFile = 0;
   if (argc == 4) {
      filename_mesh1 = argv[1];
      filename_mesh2 = argv[2];
      newFile = argv[3];
   } else {
      printf("Usage: %s <mesh_filename1> <mesh_filename2>  <newFile>\n",
            argv[0]);
      if (argc != 1)
         return 1;
      const char * f1 = STRINGIFY(SRCDIR) "/m1.vtk";
      const char * f2 = STRINGIFY(SRCDIR) "/m2.vtk";
      const char* f3 = "intx.vtk";
      printf("No files specified.  Defaulting to: %s  %s  %s\n", f1, f2, f3);
      filename_mesh1 = f1;
      filename_mesh2 = f2;
      newFile = f3;
   }

   // initialize the Meshes
   int err;
   iMesh_Instance mesh1;
   iMesh_newMesh("", &mesh1, &err, 0);
   CHECK_ERR("Couldn't create mesh 1.");

   iBase_EntitySetHandle root_set1;
   iMesh_getRootSet(mesh1, &root_set1, &err);
   CHECK_ERR("Couldn't get root set 1.");

   iMesh_Instance mesh2;
   iMesh_newMesh("", &mesh2, &err, 0);
   CHECK_ERR("Couldn't create mesh 2.");

   iBase_EntitySetHandle root_set2;
   iMesh_getRootSet(mesh2, &root_set2, &err);
   CHECK_ERR("Couldn't get root set 2.");

   // read meshes
   clock_t start_time = clock();
   iMesh_load(mesh1, root_set1, filename_mesh1, NULL, &err, strlen(
         filename_mesh1), 0);
   CHECK_ERR("Couldn't load mesh file 1.");

   iMesh_load(mesh2, root_set2, filename_mesh2, NULL, &err, strlen(
         filename_mesh2), 0);
   CHECK_ERR("Couldn't load mesh file 2.");

   clock_t load_time = clock();

   // make intersect mesher
   IntersectMesh *psm = new IntersectMesh(mesh1, root_set1, mesh2, root_set2);

   iMesh_Instance newMesh;
   iMesh_newMesh("", &newMesh, &err, 0);
   CHECK_ERR("Couldn't create new mesh.");

   iBase_EntitySetHandle new_root_set;
   iMesh_getRootSet(newMesh, &new_root_set, &err);

   CHECK_ERR("Couldn't get new root set.");

   err = psm->compute(newMesh, new_root_set);
   CHECK_ERR("Couldn't intersect meshes.");
   // write the output mesh

   iMesh_save(newMesh, new_root_set, newFile, NULL, &err, strlen(newFile), 0);
   CHECK_ERR("Couldn't save new mesh");
   clock_t mesh_time = clock();

   std::cout << "intersected mesh is succesfully produced." << std::endl;
   std::cout << "Time including loading is " << (double) (mesh_time
         - start_time) / CLOCKS_PER_SEC << " secs, Time excluding loading is "
         << (double) (mesh_time - load_time) / CLOCKS_PER_SEC << std::endl;

   delete psm;
   iMesh_dtor(mesh1, &err);
   CHECK_ERR("Couldn't destroy mesh 1.");
   iMesh_dtor(mesh2, &err);

   CHECK_ERR("Couldn't destroy mesh 2.");
   iMesh_dtor(newMesh, &err);

   return 0;
}
