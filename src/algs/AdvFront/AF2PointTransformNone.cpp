#include "meshkit/AF2PointTransformNone.hpp"

AF2Point2D AF2PointTransformNone::transformPoint(
    AF2Point2D const & point, AF2VertexBinding & vBinding) const
{
  return point;
}

AF2PointTransformNone* AF2PointTransformNone::clone() const
{
  return new AF2PointTransformNone();
}
