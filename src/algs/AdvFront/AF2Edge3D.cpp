#include "meshkit/AF2Edge3D.hpp"

AF2Edge3D::AF2Edge3D(const AF2Point3D* start, const AF2Point3D* end) :
    startPnt(start), endPnt(end)
{
  // no work to do beyond the member initializers
}

const AF2Point3D* AF2Edge3D::getStart() const
{
  return startPnt;
}

const AF2Point3D* AF2Edge3D::getEnd() const
{
  return endPnt;
}
