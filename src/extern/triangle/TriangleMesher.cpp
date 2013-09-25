/*
 * TriangleMesher.cpp
 *
 *  Created on: Sep 23, 2011
 *      Author: iulian
 */

#include "meshkit/TriangleMesher.hpp"
#include "meshkit/ModelEnt.hpp"

// triangle stuff
#define REAL double
#define ANSI_DECLARATORS 1
#define VOID int
extern "C" {
#include "triangle.h"
}

namespace MeshKit {

//Entity Type initialization for triangle meshing
moab::EntityType TriangleMesher_tps[] = { moab::MBVERTEX, moab::MBEDGE,
    moab::MBTRI, moab::MBMAXTYPE };
const moab::EntityType* TriangleMesher::output_types()
{
  return TriangleMesher_tps;
}

TriangleMesher::TriangleMesher(MKCore *mk_core, const MEntVector &me_vec) :
  MeshScheme(mk_core, me_vec)
{
}

TriangleMesher::~TriangleMesher()
{
  // TODO Auto-generated destructor stub
}

void TriangleMesher::setup_this()
{

  // nothing to do here, maybe later?
}

//---------------------------------------------------------------------------//
// execute function: Generate the mesh for decimation
void TriangleMesher::execute_this()
{
  ModelEnt * me = mentSelection.begin()->first;
  moab::Range & outR = mentSelection.begin()->second;
  moab::EntityHandle mset = me->mesh_handle(); // this is the initial set

  double frt2 = _fretting*_fretting;
  // the vertices and triangles are collected from initial moab set
  // the result will be put in the same set, and the range will store result
  struct triangulateio in, out;
  /* Define input points. */
  // first, collect the number of points
  // we may have to collect the edges and triangles in the initial set
  // this will be later on, now, concentrate on the nodes, as input
  // eventually, edges and triangles can be input for a planar complex linear graph
  //  for Triangle
  moab::Range verts;
  moab::ErrorCode rval = mk_core()->moab_instance()->get_entities_by_type(mset,
      moab::MBVERTEX, verts);
  MBERRCHK(rval, mk_core()->moab_instance());
  // if there are other elements, edges and triangles, they should be
  // probably retrieved here too, with their vertices (connectivity)
  // as it is now, get only the vertices

  unsigned int numNodes = verts.size();

  // allocate for node position
  std::vector<double> coords;
  coords.resize(3 * numNodes);
  rval = mk_core()->moab_instance()->get_coords(verts, &(coords[0]));

  MBERRCHK(rval, mk_core()->moab_instance());
  // now populate triangulateio data structure with desired data
  in.numberofpoints = numNodes;
  in.numberofpointattributes = 0;

  // depending on direction, get the right coordinates
  int i1 = (_dir % 3), i2 = (_dir + 1) % 3;// corresponding to direction 3 would be x and y
  // direction 1 would be y and z, direction 2 would be z and x
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
  for (int j = 0; j < in.numberofpoints; j++) {
    in.pointlist[2 * j] = coords[3 * j + i1];
    in.pointlist[2 * j + 1] = coords[3 * j + i2];
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
  out.segmentlist = (int*) NULL;
  out.segmentmarkerlist = (int*) NULL;
  out.numberofholes = 0;
  out.numberofregions = 0;
  out.regionlist = (REAL *) NULL;

  //char * options = const_cast<char *>(_opts.c_str());
  triangulate(_opts, &in, &out, (struct triangulateio *) NULL);// no voronoi

  // now, we can assume all nodes are original, or get new ones;
  // prefer to get new ones, although it is a waste to recreate them again
  // it would also cover cases when new nodes are created (by some refinement option)


  int numTriangles = out.numberoftriangles;
  std::cout << " numTriangles " << numTriangles << std::endl;
  // adjacency information: stored in an array 3*numElements
  //moab::EntityHandle * conn = new moab::EntityHandle [3 * numTriangles];
  //iBase_EntityHandle * conn = (iBase_EntityHandle *) adjacency;

  int numGoodTriangles = 0;

  for (int L = 0; L < numTriangles; L++) {
    /*elefile.getline(temp, 100);
     int id = strtol(temp, &pEnd, 10);*/
    int k = 0;
    int v[3];
    for (k = 0; k < 3; k++) {
      /*int indexInV = strtol(pEnd, &pEnd, 10) - 1; // it is 0 based
       v[k] = indexInV; // conn[3*L+k] = newVerts[indexInV];*/
      int newV = out.trianglelist[3 * L + k] - 1;
      v[k] = newV;
    }

    int good = 1;
    for (k = 0; k < 3; k++) {
      int k1 = (k + 1) % 3;
      double len2=0;
      for (int j=0; j<3; j++ )
      {
        double df = (coords[3*v[k]+j] - coords[3*v[k1]+j]);
        len2+= df*df;
      }
      if (len2 > frt2) {
        good = 0;
        std::cout<< "length2: " << len2 <<"\n";
        break;
      }
    }
    if (good) {
      moab::EntityHandle conn[3];

      for (k = 0; k < 3; k++) {
        conn[k] = verts[v[k]];
      }
      moab::EntityHandle triangle;
      rval = mk_core()->moab_instance()->create_element(moab::MBTRI,
          conn, 3, triangle);
      MBERRCHK(rval, mk_core()->moab_instance());
      numGoodTriangles++;
      outR.insert(triangle);


    }
  }
  // add it to the initial set (which serve as output, too!!)
  rval = mk_core()->moab_instance()->add_entities(mset, outR);
  MBERRCHK(rval, mk_core()->moab_instance());

  std::cout << "initial triangles: " << numTriangles << " after trim:"
      << numGoodTriangles << std::endl;

  // size of
  //  then entity set, then elements (triangles)

  // create triangles with this connectivity
  // because of fretting, create a triangle one by one


}

}
