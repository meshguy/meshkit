#ifndef VOL_MESHER_H
#define VOL_MESHER_H

#include "itaps.h"

#ifdef NETGEN
namespace nglib
{
#include <nglib.h>
}
using namespace nglib;
#endif

#include <tetgen.h>

using namespace std;

struct VolumeMesh
{

    size_t getSize(int dim) const
    {
        if (dim == 0) return nodes.size();
        if (dim == 3) return cells.size();
        return 0;
    }
    vector<Cell> cells;
    vector<Vertex*> nodes;
};

struct VolumeMesher
{
    int discretize(GeomMesh &geomesh, iBase_EntityHandle gCell);

    virtual VolumeMesh getMesh(const vector<Face> &face) const = 0;

    vector<Vertex*> getInsertedNodes() const
    {
        return insertedNodes;
    }

protected:
    mutable vector<Vertex*> insertedNodes;
    int checkSurfMesh();
};

#ifdef NETGEN
struct NetGenVolumeMesher : public VolumeMesher
{
    static VolumeMesher * newObject();

    NetGenVolumeMesher()
    {
        Ng_Init();
    }

    ~NetGenVolumeMesher()
    {
        Ng_Exit();
    }

    VolumeMesh getMesh(const vector<Face> &edge) const;
};
#endif

struct TetGenVolumeMesher : public VolumeMesher
{
    static VolumeMesher * newObject();

    TetGenVolumeMesher()
    {
      //exactinit();
    }

    ~TetGenVolumeMesher()
    {
    }

    VolumeMesh getMesh(const vector<Face> &faces) const;
};

struct VolumeMesherFactory
{
    static VolumeMesher * getProduct(const string & s);
};

void example_volume_mesher(const string &filename);

#endif
