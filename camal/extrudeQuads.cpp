/**
 * \file extrudeQuads.cpp
 *
 * \brief from one mesh file, with quads, on bottom, generate 
 *  hexa meshes, by extrusion, in the z direction
 *
 * compute a thickness between those 2 sheets (top loaded as geo file)
 *
 *
 */

#include "iMesh.h"
#include "iGeom.h"
#include "iRel.h"

#include "stdlib.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"
#include "moab/Skinner.hpp"
#include "MBiMesh.hpp"
#include <vector>

#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

//bool debug_surf_eval = false;


// usually a NEUMANN set
bool create_a_set_with_tag_value(iMesh_Instance meshIface, iBase_EntityHandle * ents, int num_ents, char * tag_name,
      int size_name_tag, int value)
{
   // will first create a new entity set in the mesh; will add the
   // entities;
   // will associate the tag with the given value (create the tag if not existent yet)
   // NEUMANN_SET
   int result = iBase_SUCCESS;

   bool isList = false;
   iBase_EntitySetHandle mesh_set;

   iMesh_createEntSet(meshIface, isList, &mesh_set, &result);
   if (iBase_SUCCESS != result) return false;

   iBase_TagHandle tag_handle;
   iMesh_getTagHandle(meshIface, tag_name, &tag_handle, &result, size_name_tag);
   assert (0 != tag_handle);
   iMesh_setEntSetIntData(meshIface, mesh_set, tag_handle, value, &result);
   if (iBase_SUCCESS != result) return false;
   // add the entities to the set
   //
   iMesh_addEntArrToSet(meshIface, ents, num_ents,
                          mesh_set, &result);
   if (iBase_SUCCESS != result) return false;

   return true;

}
int main(int argc, char *argv[])
{
  // Check command line arg
  std::string bottom_mesh;
  std::string top_filename;
  std::string out_mesh;
  int nbLayers = 0;
  double grading = 1.5;
  std::vector<double> grades;
  bool smooth = true;

  if (argc < 4) {
    std::cout << "Usage: " << argv[0]
        << " <bottom_mesh> <top_filename> <out_mesh> [-n layers] [-g grading ] <ratios> "

    << std::endl;

    return 0;
  } else {
    bottom_mesh = argv[1];
    top_filename = argv[2];
    out_mesh = argv[3];
    int argno = 4;
    while (argno < argc) {
      if (!strcmp(argv[argno], "-n")) {
        argno++;
        nbLayers = atoi(argv[argno]);
        argno++;
      } else if (!strcmp(argv[argno], "-g")) {
        argno++;
        grading = atof(argv[argno]);
        argno++;
      } else {
        grades.push_back(atof(argv[argno]));
        argno++;
      }
    }
  }
  int nbGrades = (int) grades.size();
  if (nbGrades > 0) {
    // add all grades, then divide by total to get the new accumulated ratios
    double total = 0.;
    int i;
    for (i = 0; i < nbGrades; i++) {
      total += grades[i];
      //grades[i] = total; // accumulated so far
    }
    for (i = 0; i < nbGrades; i++) {
      grades[i] /= total;
    }
    nbLayers = nbGrades;
  } else {
    nbGrades = nbLayers;
    for (int i = 0; i < nbLayers; i++)
      grades.push_back(1. / nbLayers);
  }

  clock_t start_time = clock();
  int err = 0;
  // read initial mesh (quad mesh surface bottom)
  iMesh_Instance mesh1;
  iMesh_newMesh(0, &mesh1, &err, 0);
  assert(iBase_SUCCESS == err);

  iGeom_Instance geom2;
  iGeom_newGeom(0, &geom2, &err, 0);
  assert(iBase_SUCCESS == err);

  // read bottom mesh
  // we use Smooth MOAB:
  char * options = NULL; // "SMOOTH;";
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh1, &root_set, &err);
  ERRORR("Couldn't get root set.", 1);

  // read  mesh

  iMesh_load(mesh1, root_set, bottom_mesh.c_str(), options, &err,
      bottom_mesh.size(), 0);
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR : can not load a mesh from " << bottom_mesh
        << std::endl;
    return 1;
  }
  clock_t load_time1 = clock();

  char * opts2 = "SMOOTH;";
  if (smooth)
    iGeom_load(geom2, top_filename.c_str(), opts2, &err, top_filename.size(), 8);
  else
    iGeom_load(geom2, top_filename.c_str(), 0, &err, top_filename.size(), 0);

  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR : can not load a geometry from " << top_filename
        << std::endl;
    return 1;
  }

  clock_t load_time2 = clock();// load the smooth file

  double direction[3] = { 0., 0., 1. }; // normalized
  // first get all nodes from bottom mesh
  // get the quads and the vertices in one shot
  iBase_EntityHandle *quads = NULL;
  int quads_alloc = 0;
  iBase_EntityHandle *vert_adj = NULL;
  int vert_adj_alloc = 0, vert_adj_size, numQuads;
  int * offsets = NULL, offsets_alloc = 0, indices_size;
  int * indices = NULL, indices_alloc = 0, offsets_size;
  iMesh_getAdjEntIndices(mesh1, root_set, iBase_FACE, iMesh_QUADRILATERAL,
      iBase_VERTEX, &quads, &quads_alloc, &numQuads, &vert_adj,
      &vert_adj_alloc, &vert_adj_size, &indices, &indices_alloc, &indices_size,
      &offsets, &offsets_alloc, &offsets_size, &err);

  ERRORR("failed to get quads.", 1);

  // create one set with desired tag (1)
  create_a_set_with_tag_value(mesh1, quads, numQuads, "NEUMANN_SET", 11, 1);
  /* get the coordinates in one array */

  int vert_coords_alloc = 0;
  double * xyz = 0; // not allocated
  int vertex_coord_size = 0;

  iMesh_getVtxArrCoords(mesh1, vert_adj, vert_adj_size, iBase_INTERLEAVED,
      &xyz, &vert_coords_alloc, &vertex_coord_size, &err);
  ERRORR("failed to get vertex coordinates of entities.", 1);

  // then, go to shoot rays to decide the sweeping thickness
  int & numNodes = vert_adj_size; // to not change too much
  double * dArr = new double[numNodes];
  // then, go to shoot rays
  int j = 0;
  int numRaysIntersected = 0;
  // double factorFloating = (1.-937./1026.);
  for (j = 0; j < numNodes; j++) {
    //dArr[j] = xyz[j * 3 + 2];
    dArr[j] = 0; // no intersection situation, marked by 0
    // for a point, see if it is inside the polygon, with winding number

    iBase_EntityHandle * intersect_entity_handles = NULL;
    int intersect_entity_handles_allocated = 0, intersect_entity_handles_size =
        0;
    double * intersect_coords = NULL;
    int intersect_coords_allocated = 0, intersect_coords_size = 0;
    double * param_coords = NULL;
    int param_coords_allocated = 0, param_coords_size = 0;
    double pos[3] = { xyz[3 * j], xyz[3 * j + 1], xyz[3 * j + 2] };
    iGeom_getPntRayIntsct(geom2, pos[0], pos[1], pos[2], direction[0],
        direction[1], direction[2], &intersect_entity_handles,
        &intersect_entity_handles_allocated, &intersect_entity_handles_size,
        iBase_INTERLEAVED, &intersect_coords, &intersect_coords_allocated,
        &intersect_coords_size, &param_coords, &param_coords_allocated,
        &param_coords_size, &err);
    // get the first coordinate
    if (err != 0 || intersect_entity_handles_size == 0 || param_coords == NULL)
      continue;
    double zTop = intersect_coords[2]; // the z of the top, intersection computation
    numRaysIntersected++;
    // consider only the first intersection point
    dArr[j] = param_coords[0]; // the first intersection only
    // this is the thickness
    // intersect point has z = intersect_coords[0]
    free(intersect_entity_handles);
    free(intersect_coords);
    free(param_coords);

  }
  // to xyz coordinates, add some thickness, and (thick/layers), to create
  // a next layer, and the hexas
  iBase_EntityHandle * layer1 = vert_adj;
  // create in a loop layer 2 , and vertices in layer 2, and hexas
  // between layer 1 and layer2
  for (int i = 0; i < nbLayers; i++) {
    for (int j = 0; j < numNodes; j++) {
      for (int k = 0; k < 3; k++)
        xyz[3 * j + k] += direction[k] * dArr[j] * grades[i];
    }
    // now create new vertices at this position, layer 2
    iBase_EntityHandle * newVerts = NULL; // no vertices yet
    int size1, size2;
    iMesh_createVtxArr(mesh1,
    /*in*/numNodes,
    /*in*/iBase_INTERLEAVED,
    /*in*/xyz,
    /*in*/numNodes * 3,
    /*inout*/&newVerts,
    /*inout*/&size1,
    /*inout*/&size2,
    /*out*/&err);
    ERRORR("failed to create verts on new layer", 1);// also, careful with the grounding line...

    // the vertices are identified as the index in vert_adj
    // indices are the quads
    // i is index in quads
    //  offsets[i], offsets[i+1] are offsets in indices
    // vertices of quad [i] have indices  indices[offsets[i]+j],
    //                    j=0:(offsets[ i+1]-offsets[i])-1
    // start copy
    int numHexa = numQuads;
    long int * adjacency = new long int[8 * numHexa];
    iBase_EntityHandle * conn = (iBase_EntityHandle *) adjacency;
    for (int L = 0; L < numHexa; L++) {

      for (int k = 0; k < 4; k++) {
        conn[8 * L + k] = layer1[indices[offsets[L] + k]];
        conn[8 * L + k + 4] = newVerts[indices[offsets[L] + k]];
      }
    }
    iBase_EntitySetHandle orig_set;
    iMesh_createEntSet(mesh1, 0, &orig_set, &err);
    int n = numHexa;
    int junk1 = n, junk2 = n, junk3 = n, junk4 = n;
    int * stat = new int[numHexa];
    int* ptr2 = stat;
    int ierr;
    iBase_EntityHandle * new_entity_handles = NULL;
    iMesh_createEntArr(mesh1, iMesh_HEXAHEDRON, conn, 8 * n,
        &new_entity_handles, &junk1, &junk2, &ptr2, &junk3, &junk4, &err);
    // end copy
    // at the end, layer 1 becomes layer 2
    layer1 = newVerts;
    // if last layer, create top quads too
    if (i == nbLayers - 1) {
      // last layer, create top quads too
      for (int L = 0; L < numQuads; L++) {

        for (int k = 0; k < 4; k++) {
          //conn[8*L+k]   = layer1[ indices[offsets[L]+k ] ];
          conn[4 * L + k] = newVerts[indices[offsets[L] + k]];
        }
      }

      iMesh_createEntArr(mesh1, iMesh_QUADRILATERAL, conn, 4 * n,
          &new_entity_handles, &junk1, &junk2, &ptr2, &junk3, &junk4, &err);
      // create another set
      create_a_set_with_tag_value(mesh1, new_entity_handles, numQuads, "NEUMANN_SET", 11, 2); // top
    }
  }

  clock_t compute_time = clock();

  // now get MOAB , to use skin on the hexas, to get the rest of quads
  moab::Interface * mb = reinterpret_cast<MBiMesh*> (mesh1)->mbImpl;
  // get all hexas, and get the skin
  moab::Range hexas, iniQuads;
  moab::ErrorCode rval = mb->get_entities_by_type(0, MBHEX, hexas);
  rval = mb->get_entities_by_type(0, MBQUAD, iniQuads);
  moab::Skinner skinner(mb);
  moab::Range allQuads;
  rval = skinner.find_skin(hexas,2, allQuads);

  moab::Range latQuads = subtract(allQuads, iniQuads);

  moab::EntityHandle latSet;
  mb->create_meshset(moab::MESHSET_SET, latSet);
  mb->add_entities(latSet, latQuads);

  // add range and tag
  moab::Tag ntag;
  rval = mb->tag_get_handle( "NEUMANN_SET", ntag );
  int val = 3;
  rval = mb->tag_set_data( ntag, &latSet, 1, (void*)&val );


  assert(iBase_SUCCESS == err);
  iMesh_save(mesh1, root_set, out_mesh.c_str(), 0, &err, out_mesh.length(), 0);
  if (iBase_SUCCESS != err) {
    std::cerr << "ERROR saving mesh to " << out_mesh << std::endl;
    return 1;
  }
  clock_t out_time = clock();
  std::cout << "Total time is " << (double) (out_time - start_time)
      / CLOCKS_PER_SEC << " s\n  load bottom : " << (double) (load_time1
      - start_time) / CLOCKS_PER_SEC << " s\n  load top : "
      << (double) (load_time2 - load_time1) / CLOCKS_PER_SEC
      << " s\n  compute time : " << (double) (compute_time - load_time2)
      / CLOCKS_PER_SEC << " s\n  write time : " << (double) (out_time
      - compute_time) / CLOCKS_PER_SEC << std::endl;
  std::cout << "num rays intersected : " << numRaysIntersected << "\n";

  free(xyz);
  delete[] dArr;
  return 0;
}

