#include "meshkit/AF2RuleNewVertex.hpp"

AF2RuleNewVertex::AF2RuleNewVertex(AF2Point2D const & rfrncPoint,
        const AF2PointTransform* const & pntTrnsfrm)
{
  referencePoint = rfrncPoint;
  pointTransform = pntTrnsfrm->clone();
}

AF2RuleNewVertex::~AF2RuleNewVertex()
{
  delete pointTransform;
}

AF2RuleNewVertex::AF2RuleNewVertex(const AF2RuleNewVertex & toCopy)
{
  referencePoint = toCopy.referencePoint;
  pointTransform = toCopy.pointTransform->clone();
}

AF2RuleNewVertex& AF2RuleNewVertex::operator=(
    const AF2RuleNewVertex & rhs)
{
  // directly copy members that are not using a pointer
  referencePoint = rhs.referencePoint;

  // copy constructor functionality,
  // but to other parts of memory, not yet to this
  AF2PointTransform* otherPointTransform;
  otherPointTransform = rhs.pointTransform->clone();

  // destructor functionality
  delete pointTransform;

  // transfer ownership from other parts of memory to this object
  pointTransform = otherPointTransform;
  otherPointTransform = NULL; // not necessary, but to be explicit

  // return this
  return *this;
}

AF2Point2D AF2RuleNewVertex::getLocation(
    AF2Binding const & vertexBinding) const
{
  return pointTransform->transformPoint(referencePoint, vertexBinding);
}
