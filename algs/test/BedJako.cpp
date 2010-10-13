#include <iostream>

#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include "math.h"
#include "iMesh.h"
#include "vec_utils.hpp"
#include <stdlib.h>
#include <cstring>

// triangle stuff
#define REAL double
#define ANSI_DECLARATORS 1
#define VOID int
extern "C" {
#include "triangle.h"
}

#define DEFAULT_TEST_FILE "gridData.txt"
#define DEFAULT_BED_FILE "output.smf"

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}


int main(int argc, char* argv[]) {
   int err;
   // check command line arg
   const char *filename_data = 0;

   char * output_bed_file = NULL;
   char * output_top_file = NULL;

   double fretting = 0;
   if (argc == 4) {
      filename_data = argv[1];
      output_bed_file = argv[2];
      fretting = atof(argv[3]);

   } else {
      printf("Usage: %s <vertex_filename> <output> <fretting>\n", argv[0]);
      if (argc != 1)
         return 1;
      printf("No file specified.  Defaulting to: %s %s 2.\n", DEFAULT_TEST_FILE,
            DEFAULT_BED_FILE );
      filename_data = DEFAULT_TEST_FILE;
      output_bed_file = DEFAULT_BED_FILE;
      fretting = 2.; //
   }

   // read vertices

   clock_t start_time = clock();

   // read the file with the vertex data
   std::ifstream datafile(filename_data, std::ifstream::in);
   if (!datafile) {
      std::cout << "can't read file\n";
      return 1;
   }
   //
   char temp[100];
   datafile.getline(temp, 100);// first line has no data in it


   std::vector<double> lats, longits, beds;
   while (!datafile.eof()) {
      datafile.getline(temp, 100);
      int id = 0;
      double latit, longit, bed, thickness, surfdata;
      sscanf(temp, "%lf %lf %lf", &latit, &longit, &bed);

      lats.push_back(latit);
      longits.push_back(longit);
      beds.push_back(bed);
   }
   datafile.close();
   int numNodes = lats.size();
   std::cout << "num nodes " << numNodes << "\n";

   /*// write the poly file

   std::ofstream myfile;
   myfile.open("JakoBed.poly");
   myfile << numNodes << " 2 0 0 \n";
   for (int i = 0; i < numNodes; i++) {
      // we just need a reasonable connectivity file
      myfile << std::setprecision(11) << (i + 1) << " " << lats[i] << " "
            << longits[i] << " " << beds[i] << " \n";
   }
   myfile << "0 \n";
   myfile << "0 \n";

   myfile.close();*/

   struct triangulateio in, out;
   /* Define input points. */

   in.numberofpoints = numNodes;
   in.numberofpointattributes = 0;
   in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
   for (int j = 0; j < in.numberofpoints; j++) {
      in.pointlist[2 * j] = lats[j];
      in.pointlist[2 * j + 1] = longits[j];
   }
   in.pointattributelist = (REAL *) NULL;

   in.pointmarkerlist = (int *) NULL;

   in.numberofsegments = 0;
   in.numberofholes = 0;
   in.numberofregions = 0;
   in.regionlist = (REAL *) NULL;

   out.pointlist = (REAL *) NULL;

   out.pointattributelist = (REAL *) NULL;
   out.trianglelist = (int *) NULL;
   out.triangleattributelist = (REAL *) NULL;

   out.numberofpoints = 0;
   out.numberofpointattributes = 0;

   out.pointmarkerlist = (int *) NULL;

   out.numberofsegments = 0;
   out.segmentlist = (int*)NULL;
   out.segmentmarkerlist = (int*) NULL;
   out.numberofholes = 0;
   out.numberofregions = 0;
   out.regionlist = (REAL *) NULL;

   // triangulate

   //std::string opts("pc");
   //triangulate(opts.c_str(), &in, &out, (struct triangulateio *) NULL);// no voronoi
   char * opts="pc";
   triangulate(opts, &in, &out, (struct triangulateio *) NULL);// no voronoi

   // run triangle

   // here run triangle in the background
   //system(" triangle -pc JakoBed ");
   iMesh_Instance mesh;
   iMesh_newMesh("", &mesh, &err, 0);
   ERRORR("Couldn't create mesh.", 1);

   iBase_EntitySetHandle root_set;
   iMesh_getRootSet(mesh, &root_set, &err);
   ERRORR("Couldn't get root set.", 1);

   //
   iBase_EntityHandle * newVerts = NULL; // no vertices yet
   // iBase_INTERLEAVED
   int size1, size2;
   // first, only xyz = lat , long, 0
   double * xyz = new double[3 * numNodes];

   for (int i = 0; i < numNodes; i++) {
      xyz[3 * i] = lats[i];
      xyz[3 * i + 1] = longits[i];
      xyz[3 * i + 2] = beds[i];// 0;// no height first
   }
   iMesh_createVtxArr(mesh,
   /*in*/numNodes,
   /*in*/iBase_INTERLEAVED,
   /*in*/xyz,
   /*in*/numNodes * 3,
   /*inout*/&newVerts,
   /*inout*/&size1,
   /*inout*/&size2,
   /*out*/&err);
   if (err != 0) {
      std::cout << "can't create vertices\n";
      exit(1);
   }
   //
   // read the new file with the element connectivity
   /*std::ifstream elefile("JakoBed.1.ele", std::ifstream::in);
   //
   elefile.getline(temp, 100);

   long int numElements;
   char * pEnd;
   numElements = strtol(temp, &pEnd, 10);*/
   int numTriangles = out.numberoftriangles;
   std::cout<< " numTriangles " << numTriangles << std::endl;
   // adjacency information: stored in an array 3*numElements
   long int * adjacency = new long int[3 * numTriangles];
   iBase_EntityHandle * conn = (iBase_EntityHandle *) adjacency;

   int numGoodTriangles = 0;

   for (int L = 0; L < numTriangles; L++) {
      /*elefile.getline(temp, 100);
      int id = strtol(temp, &pEnd, 10);*/
      int k = 0;
      int v[3];
      for (k = 0; k < 3; k++) {
         /*int indexInV = strtol(pEnd, &pEnd, 10) - 1; // it is 0 based
         v[k] = indexInV; // conn[3*L+k] = newVerts[indexInV];*/
         int newV = out.trianglelist[3*L+k]-1;
         v[k] = newV;
      }

      int good = 1;
      for (k = 0; k < 3; k++) {
         int k1 = (k + 1) % 3;
         double lenEdge = dist2(&xyz[3 * v[k]], &xyz[3 * v[k1]]);

         if (lenEdge > fretting) {
            good = 0;
            break;
         }
      }
      if (good) {
         for (k = 0; k < 3; k++) {
            conn[3 * numGoodTriangles + k] = newVerts[v[k]];
         }
         numGoodTriangles++;
      }
   }
   std::cout << "initial triangles: " << numTriangles << " after trim:"
         << numGoodTriangles << std::endl;

   // size of
   //  then entity set, then elements (triangles)

   iBase_EntitySetHandle orig_set;
   iMesh_createEntSet(mesh, 0, &orig_set, &err);
   int n = numGoodTriangles;
   int junk1 = n, junk2 = n, junk3 = n, junk4 = n;
   int * stat = new int[numGoodTriangles];
   int* ptr2 = stat;

   iBase_EntityHandle * new_entity_handles = NULL;
   iMesh_createEntArr(mesh, iMesh_TRIANGLE, conn, 3 * numGoodTriangles,
         &new_entity_handles, &junk1, &junk2, &ptr2, &junk3, &junk4, &err);
   ERRORR("Couldn't create bed mesh.", 1);

   iMesh_addEntArrToSet(mesh, new_entity_handles, numGoodTriangles, orig_set,
         &err);
   ERRORR("Couldn't add to set.", 1);
   // write this mesh as an smf file, then run it with 100000 faces reduction
   //char * bedSmfFile = "Bed.smf";
   std::cout << "Saving in file " << output_bed_file << std::endl;
   iMesh_save(mesh, root_set, output_bed_file, NULL, &err, strlen(
         output_bed_file), 0);

   delete[] xyz;
   free(newVerts);
   free (in.pointlist);
   // we should also free the triangle generated stuff
   //
   return 0;
}
