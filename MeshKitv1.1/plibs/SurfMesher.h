#ifndef SURF_MESHER_H
#define SURF_MESHER_H

#include "itaps.h"

#ifdef NETGEN
#include "ITAP_NetGen_SurfMesh.h"

extern "C"
{
// using instead of triangle.h because it conflicts with OCC
#include <TriangleMeshGenerator.h>  
}

namespace nglib
{
#include <nglib.h>
}

using namespace nglib;
#endif

using namespace std;

double get_edge_length(const Vertex *v0, const Vertex *v1);

////////////////////////////////////////////////////////////////////////////////

struct SurfaceMesh
{

    size_t getSize(int dim) const
    {
        if (dim == 0) return nodes.size();
        if (dim == 2) return faces.size();
        return 0;
    }

    vector<Face> faces;
    vector<Vertex*> nodes;

    int laplacian_uv_smoothing(vector<Vertex*> &restricted);
    int laplacian_xyz_smoothing(vector<Vertex*> &restricted);
};

////////////////////////////////////////////////////////////////////////////////

class SurfaceMesher
{
public:

    int discretize(GeomMesh &geomesh, iBase_EntityHandle gFace);

    virtual SurfaceMesh getMesh(const vector<Edge> &edge) const = 0;

protected:
    iGeom_Instance geometry;
    iBase_EntityHandle faceHandle;
    iBase_TagHandle    geom_id_tag;

    vector<Vertex*> getInsertedNodes() const
    {
        return insertedNodes;
    }

    mutable vector<Vertex*> insertedNodes;

    int saveAs(const vector<Edge> &segments, const string &filename);
    int saveAs(SurfaceMesh &uvmesh, const string &filename);


};

////////////////////////////////////////////////////////////////////////////////
#ifdef TRIANGLE
struct TriangleSurfaceMesher : public SurfaceMesher
{
    static SurfaceMesher * newObject();

    TriangleSurfaceMesher()
    {
        exactinit();
    }


    ~TriangleSurfaceMesher()
    { }

    SurfaceMesh getMesh(const vector<Edge> &edge) const;
};
#endif

////////////////////////////////////////////////////////////////////////////////
#ifdef NETGEN
struct NetGenSurfaceMesher : public SurfaceMesher
{
    static SurfaceMesher * newObject();

    ~NetGenSurfaceMesher()
    {
    }

    SurfaceMesh getMesh(const vector<Edge> &edge) const;
};
#endif
////////////////////////////////////////////////////////////////////////////////

struct TetGenSurfaceMesher : public SurfaceMesher
{
    static SurfaceMesher * newObject();

    ~TetGenSurfaceMesher()
    {
    }

    SurfaceMesh getMesh(const vector<Edge> &edge) const
    {
    }
};

////////////////////////////////////////////////////////////////////////////////

struct SurfaceMesherFactory
{
    static SurfaceMesher * getProduct(const string & s);
};

////////////////////////////////////////////////////////////////////////////////

struct FlipEdge
{

    FlipEdge()
    {
        faces[0] = NULL;
        faces[1] = NULL;

        connect[0] = NULL;
        connect[1] = NULL;

        oppositeNodes[0] = NULL;
        oppositeNodes[1] = NULL;
    }
    Face * faces[2];
    Vertex * connect[2];
    Vertex * oppositeNodes[2];
    bool isFlipAllowed() const;
};

////////////////////////////////////////////////////////////////////////////////

class EdgeFlipping
{
public:

    EdgeFlipping(SurfaceMesh &m)
    {
        mesh = m;
    }

    void execute();

private:
    SurfaceMesh mesh;

    vector<Face*> faces;
    std::map<Vertex*, set<Face*> > vFaceMap;

    void build_relations();
    int commit(const FlipEdge &f);

    Vertex* opposite_vertex(Face *f, Vertex *v0, Vertex *v1);
    FlipEdge build_edge(Vertex *v0, Vertex *v1);

    void saveAs(const string &s);
};

////////////////////////////////////////////////////////////////////////////////
#endif
