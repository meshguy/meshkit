#include "meshkit/AF2RuleNewTriangle.hpp"

// MeshKit
#include "meshkit/Error.hpp"

AF2RuleNewTriangle::AF2RuleNewTriangle(
    unsigned int firstIndex, unsigned int secondIndex, unsigned int thirdIndex)
{
  triVtxIndices[0] = firstIndex;
  triVtxIndices[1] = secondIndex;
  triVtxIndices[2] = thirdIndex;
}

unsigned int AF2RuleNewTriangle::getNumVertices() const
{
  return 3u;
}

unsigned int AF2RuleNewTriangle::getVertexIndex(unsigned int vtxNum) const
{
  if (vtxNum > 2)
  {
    MeshKit::Error badArg(MeshKit::MK_BAD_INPUT);
    badArg.set_string(
        "AF2RuleNewTriangle vertex number must not be greater than 2.");
    throw badArg;
  }
  return triVtxIndices[vtxNum];
}
