#include "meshkit/AF2RuleNewEdge.hpp"

AF2RuleNewEdge::AF2RuleNewEdge(int firstIndex, int secondIndex)
{
  a = firstIndex;
  b = secondIndex;
}

int AF2RuleNewEdge::getFirstIndex() const
{
  return a;
}

int AF2RuleNewEdge::getSecondIndex() const
{
  return b;
}
