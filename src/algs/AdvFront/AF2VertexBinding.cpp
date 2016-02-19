#include "meshkit/AF2VertexBinding.hpp"
#include "meshkit/Error.hpp"

AF2VertexBinding::AF2VertexBinding(
    std::list<const AF2RuleExistVertex*> const & verticesToBind)
{
  MeshKit::Vector<2> zeroPnt;
  zeroPnt.zero();
  validCache = false;
  typedef std::list<const AF2RuleExistVertex*>::const_iterator ItrType;
  for (ItrType itr = verticesToBind.begin();
      itr != verticesToBind.end(); ++itr)
  {
    bindingMap[*itr] = zeroPnt;
  }
}

void AF2VertexBinding::bindVertex(const AF2RuleExistVertex* vertexPtr,
    MeshKit::Vector<2> const & coordinates)
{
  std::map<const AF2RuleExistVertex*, MeshKit::Vector<2> >::iterator itr;
  itr = bindingMap.find(vertexPtr);
  if (itr == bindingMap.end())
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("Vertex does not exist as key in binding map.");
    throw badArg;
  }
  validCache = false;
  itr->second = coordinates;
}
