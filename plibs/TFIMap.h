#ifndef TFI_MAP_H
#define TFI_MAP_H

#include "GFace.h"

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>
#include <MBCN.hpp>

#include <set>
#include <map>
#include <list>

#include "SpectralElements.h"

class GFace;

using namespace std;

class TFIMap {
public:
    static const int PARAMETRIC_TFI = 0;

    static const int PHYSICAL_TFI = 1;

    static double linear_interpolation(double r, double x0, double x1);

    static double bilinear_interpolation(double r, double s, double *valCorners);

    static double bilinear_interpolation(double r, double s,
            double x00, double x10, double x11, double x01);

    static double trilinear_interpolation(double r, double s, double t, double *valCorners);

    static double trilinear_interpolation(double r, double s, double t,
            double x000, double x100, double x110, double x010,
            double x001, double x101, double x111, double x011);

    static double transfinite_blend(double r, double s,
            double x00, double x10, double x11, double x01,
            double xr0, double x1s, double xr1, double x0s);

    static double transfinite_blend(double r, double s, double t,
            double x000, double xr00, double x100,
            double x0s0, double xrs0, double x1s0,
            double x010, double xr10, double x110,
            double x00t, double xr0t, double x10t,
            double x0st, double x1st,
            double x01t, double xr1t, double x11t,
            double x001, double xr01, double x101,
            double x0s1, double xrs1, double x1s1,
            double x011, double xr11, double x111);

    static void blend_from_corners(vector<double> &x,
            const vector<double> &glxnodes);

    static void blend_from_corners(vector<double> &x,
            const vector<double> &glxnodes,
            const vector<double>& glynodes);

    static void blend_from_corners(vector<double> &x,
            const vector<double> &glxnodes,
            const vector<double> &glynodes,
            const vector<double> &glznodes);

    static void blend_from_edges(vector<double> &x,
            const vector<double> &glxnodes,
            const vector<double> &glynodes);

    static void blend_from_edges(vector<double> &x,
            const vector<double> &glxnodes,
            const vector<double> &glynodes,
            const vector<double> &glznodes);

    static void blend_from_faces(vector<double> &x,
            const vector<double> &glxnodes,
            const vector<double> &glynodes,
            const vector<double> &glznodes);

    TFIMap(iMesh_Instance m, iGeom_Instance g = 0, iRel_Instance a = 0, iRel_RelationHandle r = 0) 
    {
        setMesh(m);
        hoelement_order = 2;
        tfi_type = PHYSICAL_TFI;

        if( g )  {
           setGeometry(g, a, r);
           tfi_type = PARAMETRIC_TFI;
        }
    }

    void setGeometry(iGeom_Instance g, iRel_Instance a, iRel_RelationHandle r) 
    {
        geometry = g;
        assoc = a;
        rel = r;
    }

    void setHigherElementOrder(int n) {
        hoelement_order = n;
    }

    void setMesh(iMesh_Instance m) 
    {
        mesh = m;
        int err;
        const char *tagname = "HO_POINTS";
        int namelen = strlen(tagname);
        iMesh_getTagHandle(mesh, tagname, &horder_tag, &err, namelen);
    }

    void setType(int t) {
        tfi_type = t;
    }

    int getTFI2D(iBase_EntityHandle gfacehandle, const vector<double> &gnodes);
    int getTFI3D(iBase_EntityHandle gcellhandle, const vector<double> &gnodes);

    void projectHigherOrderNodes2D( iBase_EntityHandle mFaceHandle, const vector<double> &gnodes);
    void projectHigherOrderNodes3D( iBase_EntityHandle mCellHandle, const vector<double> &gnodes);

    void saveAs( const string &s);
private:
    int tfi_type;

    iMesh_Instance mesh;
    iGeom_Instance geometry;
    iRel_Instance assoc;
    iRel_RelationHandle rel;

    iBase_TagHandle horder_tag; // Higher order elements tag on mesh entities
    iBase_EntityHandle gFaceHandle;

    struct CoEdge {

        CoEdge() {
            direction = 0;
            imaginary = 0;
        }
        Point3D getXYZCoords(double u); // U is relative ( 0-1);
        Point2D getUVCoords(double u);  // U is relative ( 0-1);
        int direction;
        bool imaginary; // wtheher the edge is fictitious
        iGeom_Instance geometry;
        GFace *geomface;
        iBase_EntityHandle edgeHandle;
    };

    CoEdge coedge[4];
    GFace *geomface; // An helper class for geometric face.

    int hoelement_order;
    int nx, ny, nz;

    vector<double> gllnodes; // Gauss-Labotto nodal points.
    vector<double> arclength_ratio;
    vector<iBase_EntityHandle> structured_nodes;
    vector<double> uCorner, vCorner;

    int search_edge(set<iBase_EntityHandle> &aset,
            iBase_EntityHandle n0,
            iBase_EntityHandle n1, iBase_EntityHandle &edgeHandle);

    void project_face(iBase_EntityHandle edgeHandle);

    void project_edge_u(iBase_EntityHandle edgeHandle, int dir, double u0, double uN, double s);
    void project_edge_v(iBase_EntityHandle edgeHandle, int dir, double v0, double vN, double r);

    int init_coedges2();
    int init_coedges4();
    int init_coedges();

    //
    // iMesh may give faces in an unstructured way, even though the mesh
    // be topologically structured. TFI require mesh to be structured and
    // if the underlying mesh is unstructured and we want TFI over it,
    // then other tricks are applied and this case is separately 
    // considered. The following method will fill the "structured_nodes"
    // data member in the structured order.
    //
    int reconstruct_structured_mesh2d();
    int reconstruct_structured_mesh3d();

    void physicalTFI2D();
    void parametricTFI2D();
};


#endif
