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

CopyMoveVerts::CopyMoveVerts(iMesh_Instance impl,double *dv) : CopyVerts(impl)
{
    memcpy(dv_,dv,sizeof(dv_));
}

void CopyMoveVerts::transform(int n,int i,double *coords) const
{
    coords[0] += dv_[0]*n;
    coords[1] += dv_[1]*n;
    coords[2] += dv_[2]*n;
}

CopyRotateVerts::CopyRotateVerts(iMesh_Instance impl,double *origin,
                                 double *angles)
    : CopyVerts(impl)
{
    memcpy(origin_,origin,sizeof(origin_));
    memcpy(angles_,angles,sizeof(angles_));
}

void CopyRotateVerts::transform(int n,int i,double *coords) const
{
    double s[3],c[3];
    for(int i=0; i<3; i++)
    {
        s[i] = sin(angles_[i]*n);
        c[i] = cos(angles_[i]*n);
    }

    double tmp[3] = { coords[0]-origin_[0],
                      coords[1]-origin_[1],
                      coords[2]-origin_[2] };

    coords[0] =
        ( c[2]*c[0] - c[1]*s[0]*s[2]) * tmp[0] +
        ( c[2]*s[0] + c[1]*c[0]*s[2]) * tmp[1] +
        ( s[2]*s[1])                  * tmp[2] +
        origin_[0];

    coords[1] =
        (-s[2]*c[0] - c[1]*s[0]*c[2]) * tmp[0] +
        (-s[2]*s[0] + c[1]*c[0]*c[2]) * tmp[1] +
        ( c[2]*s[1])                  * tmp[2] +
        origin_[1];

    coords[2] = 
        ( s[1]*s[0]) * tmp[0] +
        (-s[1]*c[0]) * tmp[1] +
        ( c[1])      * tmp[2] +
        origin_[2];
}
