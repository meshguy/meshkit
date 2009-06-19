#include "CopyVerts.hpp"

#include <cassert>
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
