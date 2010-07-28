#include "Transform.hpp"

#include "vec_utils.hpp"
#include "MKException.hpp"
#include "SimpleArray.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

void copy::Transform::operator ()(iMesh_Instance impl, iBase_EntityHandle *src,
                                  int src_size, iBase_EntityHandle **dest,
                                  int *dest_alloc, int *dest_size) const
{
  int err;

  SimpleArray<double> coords;
  iMesh_getVtxArrCoords(impl, src, src_size, iBase_INTERLEAVED,
                        ARRAY_INOUT(coords), &err);
  check_error(impl, err);

  for (int i=0; i<src_size; i++)
    transform(&coords[i*3]);

  iMesh_createVtxArr(impl, src_size, iBase_INTERLEAVED, ARRAY_IN(coords),
                     dest, dest_alloc, dest_size, &err);
  check_error(impl, err);

  assert(*dest_size == src_size); // Sanity check
}

extrude::Transform::Transform(int steps) : steps_(steps)
{}

void extrude::Transform::operator ()(int step, iMesh_Instance impl,
                                     iBase_EntityHandle *src, int src_size,
                                     iBase_EntityHandle **dest, int *dest_alloc,
                                     int *dest_size) const
{
  int err;

  SimpleArray<double> coords;
  iMesh_getVtxArrCoords(impl, src, src_size, iBase_INTERLEAVED,
                        ARRAY_INOUT(coords), &err);
  check_error(impl, err);

  for (int i=0; i<src_size; i++)
    transform(step, &coords[i*3]);

  iMesh_createVtxArr(impl, src_size, iBase_INTERLEAVED, ARRAY_IN(coords),
                     dest, dest_alloc, dest_size, &err);
  check_error(impl, err);

  assert(*dest_size == src_size); // Sanity check
}

copy::Translate::Translate(const double *dv)
{
  memcpy(dv_, dv, sizeof(dv_));
}

void copy::Translate::transform(double *coords) const
{
  coords[0] += dv_[0];
  coords[1] += dv_[1];
  coords[2] += dv_[2];
}

copy::Rotate::Rotate(const double *origin, const double *z, double dtheta)
  : dtheta_(dtheta)
{
  memcpy(origin_, origin, sizeof(origin_));
  normalize(z_, z);
}

void copy::Rotate::transform(double *coords) const
{
  // uses Rodrigues' rotation formula
  double tmp[3];
  for(int i=0; i<3; i++)
    tmp[i] = coords[i]-origin_[i];

  double x[3];
  cross(x,z_,tmp);

  double a = cos(dtheta_);
  double b = sin(dtheta_);
  double c = dot(tmp, z_)*(1-cos(dtheta_));

  for(int i=0; i<3; i++)
    coords[i] = a*tmp[i] + b*x[i] + c*z_[i] + origin_[i];
}



extrude::Translate::Translate(const double *dv, int steps) : Transform(steps)
{
  dv_[0] = dv[0] / steps;
  dv_[1] = dv[1] / steps;
  dv_[2] = dv[2] / steps;
}

void extrude::Translate::transform(int step, double *coords) const
{
  coords[0] += dv_[0]*step;
  coords[1] += dv_[1]*step;
  coords[2] += dv_[2]*step;
}

extrude::Rotate::Rotate(const double *origin, const double *z, double dtheta,
                     int steps) : Transform(steps), dtheta_(dtheta/steps)
{
  memcpy(origin_, origin, sizeof(origin_));
  normalize(z_, z);
}

void extrude::Rotate::transform(int step, double *coords) const
{
  // uses Rodrigues' rotation formula
  double tmp[3];
  for(int i=0; i<3; i++)
    tmp[i] = coords[i]-origin_[i];

  double x[3];
  cross(x,z_,tmp);

  double a = cos(dtheta_*step);
  double b = sin(dtheta_*step);
  double c = dot(tmp, z_)*(1-cos(dtheta_*step));

  for(int i=0; i<3; i++)
    coords[i] = a*tmp[i] + b*x[i] + c*z_[i] + origin_[i];
}
