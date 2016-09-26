#include "meshkit/AF2RuleExistEdge.hpp"

AF2RuleExistEdge::AF2RuleExistEdge(const AF2RuleExistVertex* firstVtx,
    const AF2RuleExistVertex* secondVtx,
    double coeffAlpha, double coeffBravo, double coeffCharlie)
{
  firstVertexPtr = firstVtx;
  secondVertexPtr = secondVtx;
  vecDiffX = secondVtx->getX() - firstVtx->getX();
  vecDiffY = secondVtx->getY() - firstVtx->getY();
  a = coeffAlpha;
  b = coeffBravo;
  c = coeffCharlie;
}

const AF2RuleExistVertex* AF2RuleExistEdge::getStart() const
{
  return firstVertexPtr;
}

const AF2RuleExistVertex* AF2RuleExistEdge::getEnd() const
{
  return secondVertexPtr;
}

bool AF2RuleExistEdge::isMatching(AF2Point2D const & startPnt,
    AF2Point2D const & endPnt, double maxDeviation) const
{
  double matchDiffX = endPnt.getX() - startPnt.getX();
  double matchDiffY = endPnt.getY() - startPnt.getY();
  double dx = matchDiffX - vecDiffX;
  double dy = matchDiffY - vecDiffY;
  double deviation = a*dx*dx + b*dx*dy + c*dy*dy;
  return (deviation < maxDeviation);
}
