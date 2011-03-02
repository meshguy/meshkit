#ifndef GFXGEOM_3D_INCLUDED // -*- C++ -*-
#define GFXGEOM_3D_INCLUDED

#include "Vec3.h"
#include "Vec4.h"
#include "Array.h"

class Bounds
{
public:

    Vec3 min, max;
    Vec3 center;
    double radius;
    unsigned int points;

    Bounds() { reset(); }

    void reset();
    void addPoint(const Vec3&);
    void complete();
};

class Plane
{
    //
    // A plane is defined by the equation:  n*p + d = 0
    Vec3 n;
    double d;

public:

    Plane() : n(0,0,1) { d=0; } // -- this will define the XY plane
    Plane(const Vec3& p, const Vec3& q, const Vec3& r) { calcFrom(p,q,r); }
    Plane(const array<Vec3>& verts) { calcFrom(verts); }
    Plane(const Plane& p) { n=p.n; d=p.d; }

    void calcFrom(const Vec3& p, const Vec3& q, const Vec3& r);
    void calcFrom(const array<Vec3>&);

    bool isValid() const { return n[X]!=0.0 || n[Y]!=0.0 || n[Z]!= 0.0; }
    void markInvalid() { n[X] = n[Y] = n[Z] = 0.0; }

    double distTo(const Vec3& p) const { return n*p + d; }
    const Vec3& normal() const { return n; }

    void coeffs(double *a, double *b, double *c, double *dd) const {
        *a=n[X]; *b=n[Y]; *c=n[Z]; *dd=d;
    }
    Vec4 coeffs() const { return Vec4(n,d); }
};

// GFXGEOM_3D_INCLUDED
#endif
