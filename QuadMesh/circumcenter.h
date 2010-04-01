#ifndef CCENTER_H
#define CCENTER_H

#include "geomPredicates.h"
#include "basic_math.h"

void TetCircumCenter( double *a, double *b, double *c, double *d, double *r, double *p);
void TriCircumCenter2D( double *a, double *b, double *c, double *r, double *p);
void TriCircumCenter3D( double *a, double *b, double *c, double *r, double *p);
void TriCircumCenter3D( double *a, double *b, double *c, double *r);

namespace UnitTest
{
   int test_tet_circumcenter();
   int test_tri2d_circumcenter();
   int test_tri3d_circumcenter();
}
   
#endif
