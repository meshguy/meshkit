#include "meshkit/Transform.hpp"

#include <cassert>
#include <cmath>

namespace MeshKit {

Vector<3> rodrigues(const Vector<3> &pt, const Vector<3> &z, double dtheta) {
  Vector<3> x = vector_product(z, pt);

  double a = cos(dtheta);
  double b = sin(dtheta);
  double c = (pt % z)*(1-a);

  return a*pt + b*x + c*z;
}

} // namespace MeshKit
