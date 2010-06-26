#ifndef GEDGE_H
#define GEDGE_H

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include <boost/array.hpp>
#include <vector>

#include <ANN/ANN.h>

#include "SimpleArray.hpp"

using namespace std;

typedef boost::array<double, 3 > Vec3D;
typedef boost::array<double, 2 > Point2D;
typedef boost::array<double, 3 > Point3D;

inline double square_length(const Point3D &p1, const Point3D &p2) 
{
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double dz = p2[2] - p1[2];
    return dx * dx + dy * dy + dz*dz;
}

inline Point3D linear_interpolation01(const Point3D &p1, const Point3D &p2, double t) 
{
    assert( t >= 0.0 && t <= 1.0);

    Point3D  p3d;
    p3d[0] = (1-t)*p1[0] + t*p2[0];
    p3d[1] = (1-t)*p1[1] + t*p2[1];
    p3d[2] = (1-t)*p1[2] + t*p2[2];
    return p3d;
    
}

struct HO_Points 
{
     HO_Points() { nx = 0; ny = 0; nz = 0; numHPoints = 0; nodeHandles = 0; }
    int nx, ny, nz;
    int numHPoints;
    iBase_EntityHandle* nodeHandles;
};

class GEdge 
{
public:

    static int get_quad_edge_number(const iBase_EntityHandle *facenodes,
            const iBase_EntityHandle *edgenodes, int &side_no, int &orientation);

    static int get_hex_edge_number(const iBase_EntityHandle *cellnodes,
            const iBase_EntityHandle *edgenodes, int &side_no, int &orientation);

    GEdge(iGeom_Instance g, iBase_EntityHandle h,
            iMesh_Instance m = 0, iRel_Instance a = 0, iRel_RelationHandle r = 0);

    ~GEdge();

    iBase_EntityHandle getVertex(int i) const;

    Point3D getXYZCoords(double u) const;
    Point3D getClosestPoint(const Point3D &p) const;

    Point2D getURange() const;
    double getUCoord(const Point3D &p) const;
    double getUCoord(const Point3D &p, double nearto) const;
    int getUCoord(double ustart, double dist, double &uguess) const;

    Vec3D getNormal(double u) const;
    Vec3D getFirstDer(double u) const;
    Vec3D getSecondDer(double u) const;

    double getLength() const;
    double getLength(double ustart, double uend) const;

    bool isPeriodic() const;

    // Discretize the curve using equal parametric distance, The argument persistent = 1 will update the
    // Coordinates in the data structures also.
    bool is_u_uniformly_discretized();
    void uniform_u_discretization(int n, vector<Point3D> &xyz, vector<double> &u, bool persistent = 0);

    // Discretize the curve using equal xyz distance
    bool is_xyz_uniformly_discretized();
    void uniform_xyz_discretization(int n, vector<Point3D> &xyz, vector<double> &u, bool persistent = 0);

    void projectHigherOrderNodes(const vector<double> &gllnodes);

private:

    iGeom_Instance geometry;
    iMesh_Instance mesh;
    iRel_Instance assoc;
    iRel_RelationHandle rel;

    iBase_EntityHandle gEdgeHandle;

    vector<double> gllnodes, arclength_ratio;
    iBase_TagHandle horder_tag; // Higher order elements tag on mesh entities

    bool periodic;
    double umin, umax;
    double xmin, xmax, ymin, ymax, zmin, zmax;

    // Stuff for K-D tree
    int numNodes;
    double *uCoord;
    int numNeighs;
    ANNkd_tree *kdtree;
    ANNpointArray kdnodes;

    mutable vector<ANNidx> annIdx;
    mutable vector<ANNdist> anndist;

    void refine_search(const Point3D &qPoint, double u0, double u1, double &urefine, double &mindist2) const;

    void projectHigherOrderNodes(iBase_EntityHandle mEdgeHandle);

    void build_kdtree();
    void delete_kdtree();
};

#endif
