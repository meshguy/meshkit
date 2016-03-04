#include "meshkit/AF2Edge2D.hpp"

AF2Edge2D::AF2Edge2D(const AF2Point2D* start, const AF2Point2D* end) :
    startPnt(start), endPnt(end)
{
  // no work to do beyond the member initializers
}

const AF2Point2D* AF2Edge2D::getStart() const
{
  return startPnt;
}

const AF2Point2D* AF2Edge2D::getEnd() const
{
  return endPnt;
}
