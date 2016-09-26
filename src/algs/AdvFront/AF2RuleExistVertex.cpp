#include "meshkit/AF2RuleExistVertex.hpp"

AF2RuleExistVertex::AF2RuleExistVertex(double refXCoord, double refYCoord,
    double coeffAlpha, double coeffBravo, double coeffCharlie)
{
  x = refXCoord;
  y = refYCoord;
  a = coeffAlpha;
  b = coeffBravo;
  c = coeffCharlie;
}

double AF2RuleExistVertex::getX() const
{
  return x;
}

double AF2RuleExistVertex::getY() const
{
  return y;
}

bool AF2RuleExistVertex::isMatching(
    AF2Point2D const& matchPoint, double maxDeviation) const
{
  double dx = matchPoint.getX() - x;
  double dy = matchPoint.getY() - y;
  double deviation = a*dx*dx + b*dx*dy + c*dy*dy;
  return (deviation < maxDeviation);
}
