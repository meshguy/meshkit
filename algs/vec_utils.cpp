#include "vec_utils.hpp"
#include <cmath>

double * cross(double *res, const double *a, const double *b)
{
  res[0] = a[1]*b[2] - a[2]*b[1];
  res[1] = a[2]*b[0] - a[0]*b[2];
  res[2] = a[0]*b[1] - a[1]*b[0];
  return res;
}

double dot(const double *a, const double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double dist(const double *a)
{
  return std::sqrt(dot(a,a));
}

double * normalize(double *res, const double *a)
{
  double d = dist(a);
  for(int i=0; i<3; i++)
    res[i] = a[i]/d;
  return res;
}
