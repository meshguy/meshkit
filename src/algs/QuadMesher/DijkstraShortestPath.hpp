#ifndef DIJKSHORT_H
#define DIJKSHORT_H

#include "Mesh.hpp"

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////

class  DijkstraShortestPath
{
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

        double   distance;
        Vertex   *previous, *vertex;
    };

    vector<LVertex> vnodes;

    std::priority_queue<LVertex, vector<LVertex>, greater<LVertex> > vertexQ;

    void  initialize();
    int   atomicOp( LVertex &node);
    void  fastmarching();  // Fast Marching Style algorithm O(nlogn) with heap
    void  traceback();

    double getCost(const LVertex &vi, const LVertex &vj ) const {
        double val = vi.distance;
        if(distance_measure_method == TOPOLOGICAL_DISTANCE ) {
            val += 1;
        } else {
            Point3D p1 = vi.vertex->getXYZCoords();
            Point3D p2 = vj.vertex->getXYZCoords();
            val += Math::length2(p1,p2);;
        }
        return val;
    }
};

namespace Jaal
{
int dijkstra_shortest_path_test();
};

////////////////////////////////////////////////////////////////////////////////

#endif
