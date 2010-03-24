#ifndef GFXMATH_VEC4_INCLUDED // -*- C++ -*-
#define GFXMATH_VEC4_INCLUDED

/************************************************************************

  4D Vector class.

  $Id: Vec4.h,v 1.10 1997/03/17 22:52:27 garland Exp $

 ************************************************************************/

class Vec4 {
private:
    double elt[4];

protected:
    inline void copy(const Vec4& v);

public:
    //
    // Standard constructors
    //
    Vec4(double x=0, double y=0, double z=0, double w=0) {
	elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=w;
    }
#ifdef GFXMATH_VEC3_INCLUDED
    Vec4(const Vec3& v,double w) {elt[0]=v[0];elt[1]=v[1];elt[2]=v[2];elt[3]=w;}
#endif
    Vec4(const Vec4& v) { copy(v); }
    Vec4(const double *v) { elt[0]=v[0]; elt[1]=v[1]; elt[2]=v[2]; elt[3]=v[3]; }

    //
    // Access methods
    //
#ifdef SAFETY
    double& operator()(int i)       { assert(i>=0 && i<4); return elt[i]; }
    double  operator()(int i) const { assert(i>=0 && i<4); return elt[i]; }
#else
    double& operator()(int i)       { return elt[i]; }
    double  operator()(int i) const { return elt[i]; }
#endif
    double& operator[](int i)             { return elt[i]; }
    const double& operator[](int i) const { return elt[i]; }

    double *raw()             { return elt; }
    const double *raw() const { return elt; }

    //
    // Comparison methods
    //
    inline bool operator==(const Vec4&) const;
    inline bool operator!=(const Vec4&) const;

    //
    // Assignment and in-place arithmetic methods
    //
    inline void set(double x, double y, double z, double w){
	elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=w;
    }
    inline Vec4& operator=(const Vec4& v);
    inline Vec4& operator+=(const Vec4& v);
    inline Vec4& operator-=(const Vec4& v);
    inline Vec4& operator*=(double s);
    inline Vec4& operator/=(double s);

    //
    // Binary arithmetic methods
    //
    inline Vec4 operator+(const Vec4& v) const;
    inline Vec4 operator-(const Vec4& v) const;
    inline Vec4 operator-() const;

    inline Vec4 operator*(double s) const;
    inline Vec4 operator/(double s) const;
    inline double operator*(const Vec4& v) const;
};



////////////////////////////////////////////////////////////////////////
//
// Method definitions
//

inline void Vec4::copy(const Vec4& v)
{
    elt[0]=v.elt[0]; elt[1]=v.elt[1]; elt[2]=v.elt[2]; elt[3]=v.elt[3];
}

inline bool Vec4::operator==(const Vec4& v) const
{
    double dx=elt[X]-v[X],  dy=elt[Y]-v[Y],  dz=elt[Z]-v[Z],  dw=elt[W]-v[W];
    return (dx*dx + dy*dy + dz*dz + dw*dw) < FEQ_EPS2;
}

inline bool Vec4::operator!=(const Vec4& v) const
{
    double dx=elt[X]-v[X],  dy=elt[Y]-v[Y],  dz=elt[Z]-v[Z],  dw=elt[W]-v[W];
    return (dx*dx + dy*dy + dz*dz + dw*dw) > FEQ_EPS2;
}

inline Vec4& Vec4::operator=(const Vec4& v)
{
    copy(v);
    return *this;
}

inline Vec4& Vec4::operator+=(const Vec4& v)
{
    elt[0] += v[0];   elt[1] += v[1];   elt[2] += v[2];   elt[3] += v[3];
    return *this;
}

inline Vec4& Vec4::operator-=(const Vec4& v)
{
    elt[0] -= v[0];   elt[1] -= v[1];   elt[2] -= v[2];   elt[3] -= v[3];
    return *this;
}

inline Vec4& Vec4::operator*=(double s)
{
    elt[0] *= s;   elt[1] *= s;   elt[2] *= s;  elt[3] *= s;
    return *this;
}

inline Vec4& Vec4::operator/=(double s)
{
    elt[0] /= s;   elt[1] /= s;   elt[2] /= s;   elt[3] /= s;
    return *this;
}

inline Vec4 Vec4::operator+(const Vec4& v) const
{
    return Vec4(elt[0]+v[0], elt[1]+v[1], elt[2]+v[2], elt[3]+v[3]);
}

inline Vec4 Vec4::operator-(const Vec4& v) const
{
    return Vec4(elt[0]-v[0], elt[1]-v[1], elt[2]-v[2], elt[3]-v[3]);
}

inline Vec4 Vec4::operator-() const
{
    return Vec4(-elt[0], -elt[1], -elt[2], -elt[3]);
}

inline Vec4 Vec4::operator*(double s) const
{
    return Vec4(elt[0]*s, elt[1]*s, elt[2]*s, elt[3]*s);
}

inline Vec4 Vec4::operator/(double s) const
{
    return Vec4(elt[0]/s, elt[1]/s, elt[2]/s, elt[3]/s);
}

inline double Vec4::operator*(const Vec4& v) const
{
    return elt[0]*v[0] + elt[1]*v[1] + elt[2]*v[2] + elt[3]*v[3];
}

// Make scalar multiplication commutative
inline Vec4 operator*(double s, const Vec4& v) { return v*s; }



////////////////////////////////////////////////////////////////////////
//
// Primitive function definitions
//

//
// Code adapted from VecLib4d.c in Graphics Gems V
inline Vec4 cross(const Vec4& a, const Vec4& b, const Vec4& c)
{
    Vec4 result;

    double d1 = (b[Z] * c[W]) - (b[W] * c[Z]);
    double d2 = (b[Y] * c[W]) - (b[W] * c[Y]);
    double d3 = (b[Y] * c[Z]) - (b[Z] * c[Y]);
    double d4 = (b[X] * c[W]) - (b[W] * c[X]);
    double d5 = (b[X] * c[Z]) - (b[Z] * c[X]);
    double d6 = (b[X] * c[Y]) - (b[Y] * c[X]);

    result[X] = - a[Y] * d1 + a[Z] * d2 - a[W] * d3;
    result[Y] =   a[X] * d1 - a[Z] * d4 + a[W] * d5;
    result[Z] = - a[X] * d2 + a[Y] * d4 - a[W] * d6;
    result[W] =   a[X] * d3 - a[Y] * d5 + a[Z] * d6;

    return result;
}

inline double norm(const Vec4& v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
}

inline double norm2(const Vec4& v)
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
}

inline double length(const Vec4& v) { return norm(v); }

inline double unitize(Vec4& v)
{
    double l=norm2(v);
    if( l!=1.0 && l!=0.0 )
    {
	l = sqrt(l);
	v /= l;
    }
    return l;
}



////////////////////////////////////////////////////////////////////////
//
// Misc. function definitions
//

inline std::ostream& operator<<(std::ostream& out, const Vec4& v)
{
    return
	out << "[" << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << "]";
}


#ifdef GFXGL_INCLUDED
inline void glV(const Vec4& v) { glVertex(v[X], v[Y], v[Z], v[W]); }
inline void glC(const Vec4& v) { glColor(v[X], v[Y], v[Z], v[W]); }
#endif

// GFXMATH_VEC4_INCLUDED
#endif
