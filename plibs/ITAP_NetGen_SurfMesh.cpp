#include <meshing.hpp>

#include "ITAP_NetGen_SurfMesh.h"

inline double 
Det3 (double a00, double a01, double a02,
      double a10, double a11, double a12,
      double a20, double a21, double a22)
{
  return a00*a11*a22 + a01*a12*a20 + a10*a21*a02 - a20*a11*a02 - a10*a01*a22 - a21*a12*a00;
}


Vec<3> cross_product( const Vec<3> &a, const Vec<3> &b)
{
      Vec<3> result;

      result(0) = a[1] * b[2] - b[1] * a[2];
      result(1) = -a[0] * b[2] + b[0] * a[2];
      result(2) = a[0] * b[1] - b[0] * a[1];

      return result;
}

Vec<3> cross_product( const SimpleArray<double> &a, const SimpleArray<double> &b)
{
      Vec< 3 > result;

      result(0) = a[1] * b[2] - b[1] * a[2];
      result(1) = -a[0] * b[2] + b[0] * a[2];
      result(2) = a[0] * b[1] - b[0] * a[1];

      return result;
}

void ITAP_Surface::GetNormalVector(const Point < 3 > & p,
                                  const PointGeomInfo & geominfo,
                                  Vec < 3 > & n) const
{
    int err;
    // Not Quite sure why both geominfo and Point p is passed.

    double u = geominfo.u;
    double v = geominfo.v;
    // iGeom_getEntUVtoXYZ( geometry, faceHandle, u, v, &x, &y, &z, &err);

    double nx, ny, nz;
    iGeom_getEntNrmlUV(geometry, faceHandle, u, v, &nx, &ny, &nz, &err);
    assert( !err );

    double mag = nx * nx + ny * ny + nz*nz;

    if (fabs(mag) > 1.0E-10)
    {
        mag = 1.0/sqrt(mag);
        nx *= mag;
        ny *= mag;
        nz *= mag;
    }

    n(0) = nx;
    n(1) = ny;
    n(2) = nz;
}

//////////////////////////////////////////////////////////////////////////////

void ITAP_Surface::DefineTangentialPlane(const Point < 3 > & ap1,
                                        const PointGeomInfo & geominfo1,
                                        const Point < 3 > & ap2,
                                        const PointGeomInfo & geominfo2)
{
    int err;

    if (projecttype == PLANESPACE)
    {
        p1 = ap1;
        p2 = ap2;

        GetNormalVector(p1, geominfo1, ez);

        ex = p2 - p1;
        ex -= (ex * ez) * ez;
        ex.Normalize();
        ey = cross_product(ez, ex);

        GetNormalVector(p2, geominfo2, n2);

        nmid = 0.5 * (n2 + ez);

        ez = nmid;
        ez.Normalize();

        ex = (p2 - p1).Normalize();
        ez -= (ez * ex) * ex;
        ez.Normalize();
        ey = cross_product(ez, ex);
        nmid = ez;
    }
    else
    {
        if ((geominfo1.u < umin) ||
                (geominfo1.u > umax) ||
                (geominfo2.u < umin) ||
                (geominfo2.u > umax) ||
                (geominfo1.v < vmin) ||
                (geominfo1.v > vmax) ||
                (geominfo2.v < vmin) ||
                (geominfo2.v > vmax)) throw UVBoundsException();

        p1 = ap1;
        p2 = ap2;
        psp1 = Point < 2 > (geominfo1.u, geominfo1.v);
        psp2 = Point < 2 > (geominfo2.u, geominfo2.v);

        Vec < 3 > n;
        GetNormalVector(p1, geominfo1, n);

        double u = geominfo1.u;
        double v = geominfo1.v;

        SimpleArray<double> du, dv;
        iGeom_getEnt1stDrvt(geometry, faceHandle, u, v,
                            ARRAY_INOUT(du), ARRAY_INOUT(dv), &err);

        DenseMatrix D1(3, 2), D1T(2, 3), DDTinv(2, 2);
        D1(0, 0) = du[0];
        D1(1, 0) = du[1];
        D1(2, 0) = du[2];
        D1(0, 1) = dv[0];
        D1(1, 1) = dv[1];
        D1(2, 1) = dv[2];

        Transpose(D1, D1T);
        DenseMatrix D1TD1(3, 3);

        D1TD1 = D1T*D1;
        if (D1TD1.Det() == 0) throw SingularMatrixException();

        CalcInverse(D1TD1, DDTinv);
        DenseMatrix Y(3, 2);
        Vec < 3 > y1 = (ap2 - ap1).Normalize();
        Vec < 3 > y2 = cross_product(n, y1).Normalize();
        for (int i = 0; i < 3; i++)
        {
            Y(i, 0) = y1(i);
            Y(i, 1) = y2(i);
        }

        DenseMatrix A(2, 2);
        A = DDTinv * D1T * Y;
        DenseMatrix Ainv(2, 2);

        if (A.Det() == 0) throw SingularMatrixException();

        CalcInverse(A, Ainv);

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                Amat(i, j) = A(i, j);
                Amatinv(i, j) = Ainv(i, j);
            }

        Vec < 2 > temp = Amatinv * (psp2 - psp1);

        double r = temp.Length();
        double alpha = -atan2(temp(1), temp(0));
        DenseMatrix R(2, 2);
        R(0, 0) = cos(alpha);
        R(1, 0) = -sin(alpha);
        R(0, 1) = sin(alpha);
        R(1, 1) = cos(alpha);

        A = A*R;

        if (A.Det() == 0) throw SingularMatrixException();

        CalcInverse(A, Ainv);

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                Amat(i, j) = A(i, j);
                Amatinv(i, j) = Ainv(i, j);
            }

        temp = Amatinv * (psp2 - psp1);
    }
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_Surface::ToPlane(const Point < 3 > & p3d,
                          const PointGeomInfo & geominfo,
                          Point < 2 > & pplane,
                          double h, int & zone) const
{
    if (projecttype == PLANESPACE)
    {
        Vec < 3 > p1p, n;
        GetNormalVector(p3d, geominfo, n);

        p1p = p3d - p1;
        pplane(0) = (p1p * ex) / h;
        pplane(1) = (p1p * ey) / h;

        if (n * nmid < 0)
            zone = -1;
        else
            zone = 0;

    }
    else
    {
        pplane = Point < 2 > (geominfo.u, geominfo.v);
        pplane = Point < 2 > (1 / h * (Amatinv * (pplane - psp1)));
        //      pplane = Point<2> (h * (Amatinv * (pplane-psp1)));
        //      pplane = Point<2> (1/h * ((pplane-psp1)));

        zone = 0;
    };
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_Surface::FromPlane(const Point < 2 > & pplane,
                            Point < 3 > & p3d,
                            PointGeomInfo & gi,
                            double h)
{
    int err;
    double x, y, z;

    if (projecttype == PLANESPACE)
    {
        p3d = p1 + (h * pplane(0)) * ex + (h * pplane(1)) * ey;
        Project(p3d, gi);
    }
    else
    {
        //      Point<2> pspnew = Point<2>(1/h * (Amat * Vec<2>(pplane)) + Vec<2>(psp1));
        Point < 2 > pspnew = Point < 2 > (h * (Amat * Vec < 2 > (pplane)) + Vec < 2 > (psp1));
        //      Point<2> pspnew = Point<2>(h * (Vec<2>(pplane)) + Vec<2>(psp1));
        gi.u = pspnew(0);
        gi.v = pspnew(1);
//      gi.trignum = 1;
        iGeom_getEntUVtoXYZ( geometry, faceHandle, gi.u, gi.v, &x, &y, &z, &err);
        assert( !err );
        p3d = Point <3>(x, y, z);
    }
}


/////////////////////////////////////////////////////////////////////////////

int ITAP_Surface::Project(const Point < 3 > &pOff, Point < 3 > &pOn,
                         Point < 2 > &pUV, double tolerance)
{
    cout << " Project " << endl;
    int err;
    double x[3], xold[3];
    double dx, dy, dz;

    iGeom_getEntClosestPt(geometry, faceHandle, pOff[0], pOff[1], pOff[2],
                          &x[0], &x[1], &x[2], &err);

    dx = pOff[0] - x[0];
    dy = pOff[1] - x[1];
    dz = pOff[2] - x[2];

    double sqr_tol = tolerance*tolerance;

    if (dx * dx + dy * dy + dz * dz < sqr_tol) return 0;

    xold[0] = x[0];
    xold[1] = x[1];
    xold[2] = x[2];

    double u, v;
    SearchUV searchUV;
    searchUV.getUV(geometry, faceHandle, x[0], x[1], x[2], u, v);

    SimpleArray<double> du, dv;
    iGeom_getEnt1stDrvt(geometry, faceHandle, u, v, ARRAY_INOUT(du), ARRAY_INOUT(dv), &err);

    int count = 0;
    double det, lambda, mu;

    Vec<3> n;

    while (1)
    {
        count++;

        n = cross_product(du, dv);

        det = Det3(n[0], du[0], dv[0],
                   n[1], du[1], dv[1],
                   n[2], du[2], dv[2]);

        if (det < 1e-15) return false;

        lambda = Det3(n[0], pOff[0] - x[0], dv[0],
                      n[1], pOff[1] - x[1], dv[1],
                      n[2], pOff[2] - x[2], dv[2]) / det;

        mu = Det3(n[0], du[0], pOff[0] - x[0],
                  n[1], du[1], pOff[1] - x[1],
                  n[2], du[2], pOff[2] - x[2]) / det;

        u += lambda;
        v += mu;

        iGeom_getEntUVtoXYZ(geometry, faceHandle, u, v, &x[0], &x[1], &x[2], &err);

        dx = xold[0] - x[0];
        dy = xold[1] - x[1];
        dz = xold[2] - x[2];

        if (dx * dx + dy * dy + dz * dz < sqr_tol || count >= 100) break;

        xold[0] = x[0];
        xold[1] = x[1];
        xold[2] = x[2];

        du.clear();
        dv.clear();
        iGeom_getEnt1stDrvt(geometry, faceHandle, u, v, ARRAY_INOUT(du), ARRAY_INOUT(dv), &err);
    }

    if (count >= 100) return 2;

    pOn(0) = x[0];
    pOn(1) = x[1];
    pOn(2) = x[2];

    pUV(0) = u;
    pUV(1) = v;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void ITAP_Surface::Project(Point < 3 > & p, PointGeomInfo & gi)
{
/*
    gp_Pnt pnt(p(0), p(1), p(2));

    double u, v;
    Handle(ShapeAnalysis_Surface) su = new ShapeAnalysis_Surface(occface);
    gp_Pnt2d suval = su->ValueOfUV(pnt, BRep_Tool::Tolerance(topods_face));
    suval.Coord(u, v);
    pnt = occface->Value(u, v);

    gi.u = u;
    gi.v = v;

    gi.trignum = 1;

    p = Point < 3 > (pnt.X(), pnt.Y(), pnt.Z());
*/
    cout << " CSV: Exit " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////
//
ITAP_NetGen_SurfaceMesh::ITAP_NetGen_SurfaceMesh(const iGeom_Instance &g,
                                      const iBase_EntityHandle &fh,
                                     const Box < 3 > & abb, int aprojecttype)
: Meshing2(Box <3>(abb.PMin(), abb.PMax())), surface(g, fh, aprojecttype)
{ }


/////////////////////////////////////////////////////////////////////////////

void ITAP_NetGen_SurfaceMesh::DefineTransformation(const Point3d & p1, const Point3d & p2,
                                             const PointGeomInfo * geominfo1,
                                             const PointGeomInfo * geominfo2)
{
    surface.DefineTangentialPlane(p1, *geominfo1, p2, *geominfo2);
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_NetGen_SurfaceMesh::TransformToPlain(const Point3d & locpoint,
                                         const MultiPointGeomInfo & geominfo,
                                         Point2d & planepoint,
                                         double h, int & zone)
{
    Point < 2 > hp;
    surface.ToPlane(locpoint, geominfo.GetPGI(1), hp, h, zone);
    planepoint.X() = hp(0);
    planepoint.Y() = hp(1);
}


/////////////////////////////////////////////////////////////////////////////

int ITAP_NetGen_SurfaceMesh::TransformFromPlain(Point2d & planepoint,
                                          Point3d & locpoint,
                                          PointGeomInfo & gi,
                                          double h)
{
    Point < 3 > hp;
    Point < 2 > hp2(planepoint.X(), planepoint.Y());
    surface.FromPlane(hp2, hp, gi, h);
    locpoint = hp;
    return 0;
}

/////////////////////////////////////////////////////////////////////////////

double ITAP_NetGen_SurfaceMesh::CalcLocalH(const Point3d & p, double gh) const
{
    return gh;
}

/////////////////////////////////////////////////////////////////////////////

ITAP_SurfaceMesh_Optimization::ITAP_SurfaceMesh_Optimization(const iGeom_Instance & ageometry)
: MeshOptimize2d(), geometry(ageometry)
{
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_SurfaceMesh_Optimization::ProjectPoint(int surfind, Point < 3 > & p) const
{
/*
    geometry.Project(surfind, p);
*/
   cout << " CSV: Exit " << endl;
   exit(0);
}

/////////////////////////////////////////////////////////////////////////////

int ITAP_SurfaceMesh_Optimization::ProjectPointGI(int surfind, Point<3> & p, PointGeomInfo & gi) const
{
/*
    double u = gi.u;
    double v = gi.v;

    Point < 3 > hp = p;
    if (geometry.FastProject(surfind, hp, u, v))
    {
        p = hp;
        return 1;
    }
    ProjectPoint(surfind, p);
    return CalcPointGeomInfo(surfind, gi, p);
*/
    cout << " CSV: Exit " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_SurfaceMesh_Optimization::ProjectPoint2(int surfind, int surfind2,
                                               Point < 3 > & p) const
{
/*
    TopExp_Explorer exp0, exp1;
    bool done = false;
    Handle(Geom_Curve) c;

    for (exp0.Init(geometry.fmap(surfind), TopAbs_EDGE); !done && exp0.More(); exp0.Next())
        for (exp1.Init(geometry.fmap(surfind2), TopAbs_EDGE); !done && exp1.More(); exp1.Next())
        {
            if (TopoDS::Edge(exp0.Current()).IsSame(TopoDS::Edge(exp1.Current())))
            {
                done = true;
                double s0, s1;
                c = BRep_Tool::Curve(TopoDS::Edge(exp0.Current()), s0, s1);
            }
        }

    gp_Pnt pnt(p(0), p(1), p(2));
    GeomAPI_ProjectPointOnCurve proj(pnt, c);
    pnt = proj.NearestPoint();
    p(0) = pnt.X();
    p(1) = pnt.Y();
    p(2) = pnt.Z();
*/
    cout << " CSV: Exit " << endl;
    exit(0);

}

/////////////////////////////////////////////////////////////////////////////

void ITAP_SurfaceMesh_Optimization::
GetNormalVector(int surfind, const Point < 3 > & p, PointGeomInfo & geominfo, Vec < 3 > & n) const
{
/*
    gp_Pnt pnt;
    gp_Vec du, dv;

    Handle(Geom_Surface) occface;
    occface = BRep_Tool::Surface(TopoDS::Face(geometry.fmap(surfind)));

    occface->D1(geominfo.u, geominfo.v, pnt, du, dv);

    n = cross_product(Vec < 3 > (du.X(), du.Y(), du.Z()),
                      Vec < 3 > (dv.X(), dv.Y(), dv.Z()));
    n.Normalize();

    if (geometry.fmap(surfind).Orientation() == TopAbs_REVERSED) n = -1 * n;
*/
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_SurfaceMesh_Optimization::
GetNormalVector(int surfind, const Point < 3 > & p, Vec < 3 > & n) const
{
/*
    Standard_Real u, v;

    gp_Pnt pnt(p(0), p(1), p(2));

    Handle(Geom_Surface) occface;
    occface = BRep_Tool::Surface(TopoDS::Face(geometry.fmap(surfind)));

    Handle(ShapeAnalysis_Surface) su = new ShapeAnalysis_Surface(occface);
    gp_Pnt2d suval = su->ValueOfUV(pnt, BRep_Tool::Tolerance(TopoDS::Face(geometry.fmap(surfind))));
    suval.Coord(u, v);
    pnt = occface->Value(u, v);

    gp_Vec du, dv;
    occface->D1(u, v, pnt, du, dv);

    n = cross_product(Vec3d(du.X(), du.Y(), du.Z()),
                      Vec3d(dv.X(), dv.Y(), dv.Z()));
    n.Normalize();

    if (geometry.fmap(surfind).Orientation() == TopAbs_REVERSED) n = -1 * n;
*/
    cout << " CSV Exit " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////

int ITAP_SurfaceMesh_Optimization::
CalcPointGeomInfo(int surfind, PointGeomInfo& gi, const Point < 3 > & p) const
{
/*
    Standard_Real u, v;

    gp_Pnt pnt(p(0), p(1), p(2));

    Handle(Geom_Surface) occface;
    occface = BRep_Tool::Surface(TopoDS::Face(geometry.fmap(surfind)));

    Handle(ShapeAnalysis_Surface) su = new ShapeAnalysis_Surface(occface);
    gp_Pnt2d suval = su->ValueOfUV(pnt, BRep_Tool::Tolerance(TopoDS::Face(geometry.fmap(surfind))));
    suval.Coord(u, v);

    gi.u = u;
    gi.v = v;
*/
    cout << " CSV Exit " << endl;
    exit(0);
    return 1;
}

/////////////////////////////////////////////////////////////////////////////

ITAP_Surface_Refinement::ITAP_Surface_Refinement(const iGeom_Instance & ageometry)
: Refinement(), geometry(ageometry)
{
}
/////////////////////////////////////////////////////////////////////////////

ITAP_Surface_Refinement::~ITAP_Surface_Refinement()
{
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_Surface_Refinement::
PointBetween(const Point < 3 > & p1, const Point < 3 > & p2, double secpoint,
             int surfi,
             const PointGeomInfo & gi1,
             const PointGeomInfo & gi2,
             Point < 3 > & newp, PointGeomInfo & newgi)
{
/*
    Point < 3 > hnewp;
    hnewp = p1 + secpoint * (p2 - p1);

    if (surfi > 0)
    {

        double u = gi1.u + secpoint * (gi2.u - gi1.u);
        double v = gi1.v + secpoint * (gi2.v - gi1.v);

        if (!geometry.FastProject(surfi, hnewp, u, v))
        {
            geometry.Project(surfi, hnewp);
        }

        newgi.trignum = 1;
        newgi.u = u;
        newgi.v = v;
    }

    newp = hnewp;
*/
    cout << " CSV: Exit " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_Surface_Refinement::
PointBetween(const Point < 3 > & p1, const Point < 3 > & p2, double secpoint,
             int surfi1, int surfi2,
             const EdgePointGeomInfo & ap1,
             const EdgePointGeomInfo & ap2,
             Point < 3 > & newp, EdgePointGeomInfo & newgi)
{
/*
    double s0, s1;

    Point < 3 > hnewp = p1 + secpoint * (p2 - p1);
    gp_Pnt pnt(hnewp(0), hnewp(1), hnewp(2));
    GeomAPI_ProjectPointOnCurve proj(pnt, BRep_Tool::Curve(TopoDS::Edge(geometry.emap(ap1.edgenr)), s0, s1));
    pnt = proj.NearestPoint();
    hnewp = Point < 3 > (pnt.X(), pnt.Y(), pnt.Z());
    newp = hnewp;
    newgi = ap1;
*/
    cout <<  " CSV: Exit " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_Surface_Refinement::ProjectToSurface(Point < 3 > & p, int surfi)
{
/*
    if (surfi > 0) geometry.Project(surfi, p);
*/
    cout << " CSV: Exit " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////

void ITAP_Surface_Refinement::ProjectToSurface(Point < 3 > & p, int surfi, PointGeomInfo & gi)
{
/*
    if (surfi > 0)
        if (!geometry.FastProject(surfi, p, gi.u, gi.v))
        {
            cout << "Fast projection to surface fails! Using ITAP projection" << endl;
            geometry.Project(surfi, p);
        }
*/
    cout << " CSV: Exit " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////

