#include "meshkit/AF2RuleNewTriangle.hpp"

// C++
#include <stdexcept>

AF2RuleNewTriangle::AF2RuleNewTriangle(
    int firstIndex, int secondIndex, int thirdIndex)
{
  triVtxIndices[0] = firstIndex;
  triVtxIndices[1] = secondIndex;
  triVtxIndices[2] = thirdIndex;
}

int AF2RuleNewTriangle::getNumVertices() const
{
  return 3;
}

int AF2RuleNewTriangle::getVertexIndex(int vtxNum) const
{
  if (vtxNum < 0 || vtxNum > 2)
  {
    throw std::range_error(
        "AF2RuleNewTriangle vertex number must be between 0 and 2 inclusive.");
  }
  return triVtxIndices[vtxNum];
}
