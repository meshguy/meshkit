#include "CopyVerts.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

CopyVerts::CopyVerts(iMesh_Instance impl) : impl_(impl)
{}

void CopyVerts::operator ()(int n,iBase_EntityHandle *src,int src_size,
                            iBase_EntityHandle **dest,int *dest_alloc,
                            int *dest_size) const
{
    int err;

    double *coords=0;
    int coords_alloc=0,coords_size=0;
    iMesh_getVtxArrCoords(impl_,src,src_size,iBase_INTERLEAVED,
                          &coords,&coords_alloc,&coords_size,&err);
    assert(!err);

    for(int i=0; i<src_size; i++)
        transform(n,i,coords + i*3);

    iMesh_createVtxArr(impl_,src_size,iBase_INTERLEAVED,coords,coords_size,
                       dest,dest_alloc,dest_size,&err);
    assert(!err);
    assert(*dest_size == src_size); // Sanity check

    free(coords);
}

CopyMoveVerts::CopyMoveVerts(iMesh_Instance impl,const double *dv)
    : CopyVerts(impl)
{
    memcpy(dv_,dv,sizeof(dv_));
}

void CopyMoveVerts::transform(int n,int,double *coords) const
{
    coords[0] += dv_[0]*n;
    coords[1] += dv_[1]*n;
    coords[2] += dv_[2]*n;
}

// TODO: these should be in their own file maybe?
static double * cross(double *res,const double *a,const double *b)
{
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
    return res;
}

static double dot(const double *a,const double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double dist(const double *a)
{
    return sqrt(dot(a,a));
}

static double * normalize(double *res,const double *a)
{
    double d = dist(a);
    for(int i=0; i<3; i++)
        res[i] = a[i]/d;
}

CopyRotateVerts::CopyRotateVerts(iMesh_Instance impl,const double *origin,
                                 const double *z,double theta)
    : CopyVerts(impl),theta_(theta)
{
    memcpy(origin_,origin,sizeof(origin_));
    normalize(z_,z);
}


void CopyRotateVerts::transform(int n,int,double *coords) const
{
    // uses Rodrigues' rotation formula
    double tmp[3];
    for(int i=0; i<3; i++)
        tmp[i] = coords[i]-origin_[i];

    double x[3];
    cross(x,z_,tmp);

    double a = cos(theta_*n);
    double b = sin(theta_*n);
    double c = dot(tmp,z_)*(1-cos(theta_*n));

    for(int i=0; i<3; i++)
        coords[i] = a*tmp[i] + b*x[i] + c*z_[i] + origin_[i];
}
