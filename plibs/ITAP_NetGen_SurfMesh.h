#ifndef ITAPGEOMETRY
#define ITAPGEOMETRY

#ifdef NETGEN
#include <meshing.hpp>
#include <linalg.hpp>

using namespace netgen;
#endif

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include "SimpleArray.h"
#include "SearchUV.h"

#define PARAMETERSPACE -1
#define PLANESPACE      1

class SingularMatrixException
{
};

class UVBoundsException
{
};

class ITAP_Surface
{
public:

    ITAP_Surface(const iGeom_Instance &g, const iBase_EntityHandle &fh, int aprojecttype)
    {
        geometry = g;
        faceHandle = fh;
        projecttype = aprojecttype;

        int err;
        iGeom_getEntUVRange(geometry, faceHandle, &umin, &vmin, &umax, &vmax, &err);
        assert( !err );

        umin -= fabs(umax - umin) / 100.0;
        vmin -= fabs(vmax - vmin) / 100.0;
        umax += fabs(umax - umin) / 100.0;
        vmax += fabs(vmax - vmin) / 100.0;
    };

    ~ITAP_Surface()
    { }; 

    void setProjectionType( int p ) { projecttype = p; }
    int  getProjectionType()  const { return projecttype; }

    void Project(Point < 3 > & p, PointGeomInfo & gi);

    int  Project(const Point < 3 > &pOff, Point < 3 > &pOn,
                 Point < 2 > &pUV, double tolerance);

    void GetNormalVector(const Point < 3 > & p,
                         const PointGeomInfo & geominfo,
                         Vec < 3 > & n) const;

    /**
      Defines tangential plane in ap1.
      The local x-coordinate axis point to the direction of ap2 */
    void DefineTangentialPlane(const Point < 3 > & ap1,
                               const PointGeomInfo & geominfo1,
                               const Point < 3 > & ap2,
                               const PointGeomInfo & geominfo2);


    /// Transforms 3d point p3d to local coordinates pplane
    void ToPlane(const Point < 3 > & p3d, const PointGeomInfo & geominfo,
                 Point < 2 > & pplane, double h, int & zone) const;

    /// Transforms point pplane in local coordinates to 3d point
    void FromPlane(const Point < 2 > & pplane,
                   Point < 3 > & p3d,
                   PointGeomInfo & gi,
                   double h);
protected:
    Point < 3 > p1, p2;

    /// in plane, directed p1->p2
    Vec < 3 > ex;

    /// in plane
    Vec < 3 > ey;

    /// outer normal direction
    Vec < 3 > ez;

    /// normal vector in p2
    Vec < 3 > n2;

    /// average normal vector
    Vec < 3 > nmid;

    // for transformation to parameter space
    Point < 2 > psp1, psp2;
    Vec < 2 > psex, psey;

    Mat < 2, 2 > Amat, Amatinv;

    // UV Bounds
    double umin, umax, vmin, vmax;

private:
    iGeom_Instance geometry;
    iBase_EntityHandle faceHandle;
    int projecttype;

    int project(Point < 3 > & p, double tol = 1.0E-12) const;
};

////////////////////////////////////////////////////////////////////////////////

class ITAP_NetGen_SurfaceMesh : public Meshing2
{
    ITAP_Surface surface;

    Mesh  mesh;
    iGeom_Instance geometry;

public:
    ITAP_NetGen_SurfaceMesh( const iGeom_Instance &geometry, 
                             const iBase_EntityHandle &fh,
                             const Box <3> &bbox, int aprojecttype);

    int GetProjectionType()
    {
        return surface.getProjectionType();
    }

//    void execute();

protected:


    virtual void DefineTransformation(const Point3d & p1, const Point3d & p2,
                                      const PointGeomInfo * geominfo1,
                                      const PointGeomInfo * geominfo2);
    virtual void TransformToPlain(const Point3d & locpoint,
                                  const MultiPointGeomInfo & geominfo,
                                  Point2d & plainpoint, double h, int & zone);

    virtual int TransformFromPlain(Point2d & plainpoint,
                                   Point3d & locpoint, PointGeomInfo & gi,
                                   double h);

    virtual double CalcLocalH(const Point3d & p, double gh) const;
};

int NetGen_SurfaceMeshing(iGeom_Instance &g, Mesh &m, int projecttype);

///////////////////////////////////////////////////////////////////////////////

class ITAP_SurfaceMesh_Optimization : public MeshOptimize2d
{
    iGeom_Instance geometry;

public:
    ITAP_SurfaceMesh_Optimization(const iGeom_Instance &ageometry);

    virtual void ProjectPoint(int surfind, Point < 3 > & p) const;
    virtual void ProjectPoint2(int surfind, int surfind2, Point < 3 > & p) const;
    virtual int ProjectPointGI(int surfind, Point < 3 > & p, PointGeomInfo & gi) const;
    virtual void GetNormalVector(int surfind, const Point < 3 > & p, Vec < 3 > & n) const;
    virtual void GetNormalVector(int surfind, const Point < 3 > & p, PointGeomInfo & gi,
                                 Vec < 3 > & n) const;
    virtual int CalcPointGeomInfo(int surfind, PointGeomInfo& gi, const Point < 3 > & p3) const;
};

///////////////////////////////////////////////////////////////////////////////

class ITAP_Surface_Refinement : public Refinement
{
    const iGeom_Instance geometry;

public:

    ITAP_Surface_Refinement(const iGeom_Instance &ageometry);

    virtual ~ITAP_Surface_Refinement();

    virtual void PointBetween(const Point < 3 > & p1, const Point < 3 > & p2,
                              double secpoint, int surfi,
                              const PointGeomInfo & gi1,
                              const PointGeomInfo & gi2,
                              Point < 3 > & newp, PointGeomInfo & newgi);

    virtual void PointBetween(const Point < 3 > & p1, const Point < 3 > & p2,
                              double secpoint, int surfi1, int surfi2,
                              const EdgePointGeomInfo & ap1,
                              const EdgePointGeomInfo & ap2,
                              Point < 3 > & newp, EdgePointGeomInfo & newgi);

    virtual void ProjectToSurface(Point < 3 > & p, int surfi);

    virtual void ProjectToSurface(Point < 3 > & p, int surfi, PointGeomInfo & gi);
};

#endif
