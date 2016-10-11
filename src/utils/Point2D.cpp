#include "meshkit/Point2D.hpp"

Point2D::Point2D() : x(0), y(0)
{
  // the initializers of x and y are all that is needed
}

Point2D::Point2D(double xVal, double yVal) : x(xVal), y(yVal)
{
  // the initializers of x and y are all that is needed
}

double Point2D::getX() const
{
  return x;
}

double Point2D::getY() const
{
  return y;
}

bool operator==(const Point2D & onePnt, const Point2D & otherPnt)
{
  return onePnt.getX() == otherPnt.getX() &&
      onePnt.getY() == otherPnt.getY();
}

bool operator!=(const Point2D & onePnt, const Point2D & otherPnt)
{
  return onePnt.getX() != otherPnt.getX() ||
      onePnt.getY() != otherPnt.getY();
}

std::ostream& operator<<(std::ostream& outStream, const Point2D & point)
{
  outStream << "(" << point.getX() << ", " << point.getY() << ")";
  return outStream;
}
