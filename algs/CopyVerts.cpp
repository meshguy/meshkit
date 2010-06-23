#include "CopyVerts.hpp"

#include "vec_utils.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

CopyVerts::CopyVerts(iMesh_Instance impl,int max)
  : impl_(impl),max_(max)
{}

void CopyVerts::operator ()(int n, iBase_EntityHandle *src, int src_size,
                            iBase_EntityHandle **dest, int *dest_alloc,
                            int *dest_size) const
{
  int err;

  double *coords = NULL;
  int coords_alloc = 0, coords_size;
  iMesh_getVtxArrCoords(impl_, src, src_size, iBase_INTERLEAVED, &coords,
                        &coords_alloc, &coords_size, &err);
  assert(!err);

  for(int i=0; i<src_size; i++)
    transform(n,i,coords + i*3);

  iMesh_createVtxArr(impl_, src_size, iBase_INTERLEAVED, coords, coords_size,
                     dest, dest_alloc, dest_size, &err);
  assert(!err);
  assert(*dest_size == src_size); // Sanity check

  free(coords);
}

CopyMoveVerts::CopyMoveVerts(iMesh_Instance impl, const double *dv, int max)
  : CopyVerts(impl, max)
{
  dv_[0] = dv[0]/max;
  dv_[1] = dv[1]/max;
  dv_[2] = dv[2]/max;
}

void CopyMoveVerts::transform(int n, int, double *coords) const
{
  coords[0] += dv_[0]*n;
  coords[1] += dv_[1]*n;
  coords[2] += dv_[2]*n;
}

CopyRotateVerts::CopyRotateVerts(iMesh_Instance impl, const double *origin,
                                 const double *z, double theta, int max)
  : CopyVerts(impl, max), theta_(theta/max)
{
  memcpy(origin_, origin, sizeof(origin_));
  normalize(z_, z);
}


void CopyRotateVerts::transform(int n, int, double *coords) const
{
  // uses Rodrigues' rotation formula
  double tmp[3];
  for(int i=0; i<3; i++)
    tmp[i] = coords[i]-origin_[i];

  double x[3];
  cross(x,z_,tmp);

  double a = cos(theta_*n);
  double b = sin(theta_*n);
  double c = dot(tmp, z_)*(1-cos(theta_*n));

  for(int i=0; i<3; i++)
    coords[i] = a*tmp[i] + b*x[i] + c*z_[i] + origin_[i];
}
