#ifndef EDGE_MESHER_H
#define EDGE_MESHER_H

#include "itaps.h"

using namespace std;

double get_edge_length(const Vertex *v0, const Vertex *v1);

////////////////////////////////////////////////////////////////////////////////

class EdgeMesher
{

public:

    int discretize(GeomMesh &geomesh, iBase_EntityHandle gEdge, int nSegments);

protected:
    GeomMesh geomesh;

    iBase_EntitySetHandle discretize_close_edge(iBase_EntityHandle gEdge, int nSegments);
    iBase_EntitySetHandle discretize_open_edge(iBase_EntityHandle gEdge, int nSegments);

    vector<Vertex*> getInsertedNodes() const
    {
        return insertedNodes;
    }

    mutable vector<Vertex*> insertedNodes;
};

////////////////////////////////////////////////////////////////////////////////

#endif
