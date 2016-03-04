#include "meshkit/AF2VertexBinding.hpp"
#include "meshkit/Error.hpp"

const AF2Point2D AF2VertexBinding::NOT_BOUND(0, 0);

AF2VertexBinding::AF2VertexBinding(
    std::list<const AF2RuleExistVertex*> const & verticesToBind)
{
  validCache = false;
  typedef std::list<const AF2RuleExistVertex*>::const_iterator ItrType;
  for (ItrType itr = verticesToBind.begin();
      itr != verticesToBind.end(); ++itr)
  {
    bindingMap[*itr] = &NOT_BOUND;
  }
}

void AF2VertexBinding::bindVertex(const AF2RuleExistVertex* vertexPtr,
    const AF2Point2D* pointPtr)
{
  std::map<const AF2RuleExistVertex*, const AF2Point2D*>::iterator itr;
  itr = bindingMap.find(vertexPtr);
  if (itr == bindingMap.end())
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("Vertex does not exist as key in binding map.");
    throw badArg;
  }
  validCache = false;
  itr->second = pointPtr;
}

const AF2Point2D* AF2VertexBinding::getBoundValue(
    const AF2RuleExistVertex* vertexPtr) const
{
  std::map<const AF2RuleExistVertex*, const AF2Point2D*>::const_iterator itr;
  itr = bindingMap.find(vertexPtr);
  if (itr == bindingMap.end())
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("Vertex does not exist as key in binding map.");
    throw badArg;
  }
  return itr->second;
}
