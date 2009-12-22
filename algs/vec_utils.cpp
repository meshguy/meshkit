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
  if (0.==d)// no division with 0
     d=1.;
  for(int i=0; i<3; i++)
    res[i] = a[i]/d;
  return res;
}

// distance in 2D
double dist2(const double *a, const double *b)
{
  return std::sqrt((b[0]-a[0])*(b[0]-a[0])+(b[1]-a[1])*(b[1]-a[1]));
}
// area of a triangle (with sign)
double area2D(const double * a, const double * b, const double *c) // could be negative or 0!
{
  double area = (b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]);
  return area/2;
}
