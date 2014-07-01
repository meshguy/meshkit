#ifndef DIJKSHORT_H
#define DIJKSHORT_H

#include "meshkit/Mesh.hpp"
#include "meshkit/basic_math.hpp"

#ifdef USE_HASHMAP
#include <tr1/unordered_map>
#endif

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////

class  DijkstraShortestPath {
public:
     DijkstraShortestPath( Mesh *m, int mtype = 0) {
          vsrc = 0;
          vdst = 0;
          mesh = m;
          filter = NULL;
          distance_measure_method = mtype;
          complete_path_known = 0;
     }

     const NodeSequence& getPath(Vertex *vs, Vertex *vd = NULL ) {
          assert( !vs->isRemoved() );
          filter = NULL;
          if( vs != vsrc ) complete_path_known = 0;
          vsrc = vs;
          vdst = vd;
          if(!complete_path_known) fastmarching();
          traceback();
          return nodepath;
     }

     const NodeSequence& getPath(Vertex *vs, MeshFilter *f) {
          assert( vs );
          assert( !vs->isRemoved() );
          filter = f;
          if( vs != vsrc ) complete_path_known = 0;
          vsrc = vs;
          vdst = NULL;
          if(!complete_path_known) fastmarching();
          traceback();
          return nodepath;
     }

private:
     // Input parameters ...
     Mesh *mesh;
     Vertex *vsrc, *vdst;
     MeshFilter  *filter;

     // Output data ...
     NodeSequence  nodepath;

     // Local data ...
     int   distance_measure_method;
     bool  complete_path_known;
     struct LVertex {

          LVertex() {
               distance = 0.0;
               previous = NULL;
               vertex   = NULL;
          }

          size_t getID() const {
               return vertex->getID();
          }

          bool operator > ( const LVertex &rhs) const {
               return distance > rhs.distance;
          }

          bool operator < ( const LVertex &rhs) const {
               return distance < rhs.distance;
          }

          double   distance;   // Shortest distance from the source to this point
          Vertex   *vertex;    // Current Vertex
          Vertex   *previous;  // Previous vertex
     };

#ifdef USE_HASHMAP
     std::tr1::unordered_map<Vertex*, LVertex> vmap;
     std::tr1::unordered_map<Vertex*,LVertex>::const_iterator miter;
#else
     std::map<Vertex*, LVertex> vmap;
     std::map<Vertex*, LVertex>::const_iterator miter;
#endif

     std::priority_queue<LVertex, vector<LVertex>, greater<LVertex> > vertexQ;

     int   atomicOp( LVertex &node);
     void  fastmarching();  // Fast Marching Style algorithm O(nlogn) with heap
     void  traceback();

     double getCost(const LVertex &vi, const LVertex &vj ) const {
          double val = vi.distance;
          if(distance_measure_method == TOPOLOGICAL_DISTANCE ) {
               val += 1.0;
          } else {
               const Point3D &p1 = vi.vertex->getXYZCoords();
               const Point3D &p2 = vj.vertex->getXYZCoords();
               double d  = Math::length2(p1,p2);
               assert( d > 0.0);
               val += d;
          }
          return val;
     }
};

namespace Jaal {
int dijkstra_shortest_path_test();
};

////////////////////////////////////////////////////////////////////////////////

#endif
