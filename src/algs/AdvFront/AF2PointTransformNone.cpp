#include "meshkit/AF2PointTransformNone.hpp"

MeshKit::Vector<2> AF2PointTransformNone::transformPoint(
    MeshKit::Vector<2> const & point, AF2VertexBinding & vBinding) const
{
  return point;
}

AF2PointTransformNone* AF2PointTransformNone::clone() const
{
  return new AF2PointTransformNone();
}
