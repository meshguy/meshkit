/*
 * SizingFunctionVar.cpp
 *
 *  Created on: Jun 24, 2011
 *      Author: iulian
 */

#include "meshkit/SizingFunctionVar.hpp"

namespace MeshKit {

SizingFunctionVar::SizingFunctionVar(MKCore *mkcore, int num_int, double int_size):
    SizingFunction(mkcore, num_int, int_size)
{
  // TODO Auto-generated constructor stub

}

SizingFunctionVar::~SizingFunctionVar()
{
  // TODO Auto-generated destructor stub
}

void SizingFunctionVar::set_linear_coeff(double * fixedPoint, double * coeff)
{
  fixed[0] = fixedPoint[0];
  fixed[1] = fixedPoint[1];
  fixed[2] = fixedPoint[2];
  a = coeff[0]; b= coeff[1]; c=coeff[2]; d = coeff[3];
}

// another version, easier to pythonize
// this is just because I don't understand swig
void SizingFunctionVar::set_coeff(double x0, double y0, double z0, double a1,
      double b1, double c1, double d1)
{
  fixed[0] = x0;
  fixed[1] = y0;
  fixed[2] = z0;
  a = a1; b= b1; c=c1; d = d1;
}

double SizingFunctionVar::size(double *xyz ) const
{
  double sz = d;
  if (xyz)
  {
    sz = sz+(xyz[0]-fixed[0])*a + (xyz[1]-fixed[1])*b + (xyz[2]-fixed[2])*c ;
  }
  return sz;
}
}
