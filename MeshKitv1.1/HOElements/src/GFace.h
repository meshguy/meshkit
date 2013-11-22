#ifndef GFACE_H
#define GFACE_H

#include <boost/array.hpp>
#include <boost/foreach.hpp>

#include <vector>

#include <ANN/ANN.h>

#include <meshkit/iGeom.hpp>
#include <meshkit/iMesh.hpp>
#include <meshkit/iRel.hpp>

#include "GEdge.h"

#include "TFIMap.h"

using namespace std;

typedef boost::array<double, 3 > Vec3D;
typedef boost::array<double, 2 > Point2D;
typedef boost::array<double, 3 > Point3D;

class TFIMap;

class GFace {
public:

     static int get_hex_face_number(const vector<iBase_EntityHandle> &cellnodes,
                                    const vector<iBase_EntityHandle> &edgenodes,
                                    int &side_no, int &orientation);

     GFace(iBase_EntityHandle h, iGeom *g, iMesh *m = 0, iRel *r = 0, iRel::PairHandle *a = 0);

     ~GFace();

     // iGeom at present don't support edge loops or ordering of vertices in
     // EntAdj function, so this is just an helper function

     iBase_EntityHandle getVertex(int i) const;

     int getUVCoords(const Point2D &uvstart, double dist, Point2D &uvguess) const;
     int getUVCoords(const Point3D &xyz, Point2D &puv) const;
     int getUVCoords(const Point3D &xyz, const Point2D &nearto, Point2D &puv) const;

     int getXYZCoords(const Point2D &uv, Point3D &xyz) const;
     int getClosestPoint(const Point3D &xyz, Point3D &psurf) const;

     int getNormal(const Point2D &uv, Vec3D &v) const;

     int getFirstDer(const Point2D &uv, Vec3D &du, Vec3D &dv) const;
     int getSecondDer(const Point2D &uv, Vec3D &v) const;

     double getGeodesicLength(const Point2D &u0, const Point2D &u1) const;

     double getArea() const;

     bool isPeriodic(int i) const;
     bool hasSeam() const;
     void saveAs(const string &s) const;

     void projectEdgeHigherOrderNodes(const vector<double> &gnodes);
     void projectFaceHigherOrderNodes(const vector<double> &gnodes);

private:

     iGeom *geom;
     iMesh *mesh;
     iRel  *rel;
     iRel::PairHandle *relPair;

     iBase_EntityHandle gFaceHandle;

     vector<double> gllnodes, arclength_ratio;
     iBase_TagHandle horder_tag; // Higher order elements tag on mesh entities

     double umin, umax, vmin, vmax;
     double xmin, xmax, ymin, ymax, zmin, zmax;
     double xlength, ylength, zlength, maxlength;

     bool periodic[2];

     // Stuff for K-D tree
     double *uvCoords;
     int numNeighs;
     ANNkd_tree *kdtree;
     ANNpointArray kdnodes;
     mutable vector<ANNidx> annIdx;
     mutable vector<ANNdist> anndist;

     void build_kdtree();
     void delete_kdtree();

     void projectEdgeHigherOrderNodes(iBase_EntityHandle medgehandle);
     void projectFaceHigherOrderNodes(iBase_EntityHandle mfacehandle);
};

// This is an extension to iGeom as the current iGeom doesn't return
// entities in some particular order.
void iGeom_getEntAdj_Ext(iGeom *geom,
                         iBase_EntityHandle entity_handle,
                         int to_dimension,
                         std::vector<iBase_EntityHandle> &entities);

#endif
