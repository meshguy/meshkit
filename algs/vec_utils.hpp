#ifndef VEC_UTILS_H
#define VEC_UTILS_H

double * cross(double *res, const double *a, const double *b);
double dot(const double *a, const double *b);
double dist(const double *a);
double * normalize(double *res, const double *a);
double dist2(const double *a, const double * b);
double area2D(const double * a, const double * b, const double *c); // could be negative or 0!

#endif
