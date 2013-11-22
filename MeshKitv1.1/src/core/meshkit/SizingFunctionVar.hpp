/*
 * SizingFunctionVar.hpp
 *
 *  Created on: Jun 24, 2011
 *      Author: iulian
 */

#ifndef SIZINGFUNCTIONVAR_HPP_
#define SIZINGFUNCTIONVAR_HPP_

#include "meshkit/SizingFunction.hpp"

namespace MeshKit {

class SizingFunctionVar: public MeshKit::SizingFunction
{
public:
  SizingFunctionVar(MKCore *mkcore, int num_int = -1, double int_size = -1.0);
  virtual ~SizingFunctionVar();

  // the function would be something like
  //  a(x-x0)+b(y-y0)+c(z-z0)+d
  // at fixed point, the mesh size returned would be d
  // these should be either read from a file, or
  // passed as arguments to the test / usage scenario
  void set_linear_coeff(double * fixedPoint, double * coeff);

  // another version, easier to pythonize
  void set_coeff(double x0, double y0, double c0, double a,
      double b, double c, double d);

  virtual double size(double *xyz = NULL) const;
  // will be used by edge meshers to decide are we or not variables
  // use a different edge mesh strategy then (different enum)
  virtual bool variable() {return true;}

private:
  double a, b, c, d;
  double fixed[3];

};

}

#endif /* SIZINGFUNCTIONVAR_HPP_ */
