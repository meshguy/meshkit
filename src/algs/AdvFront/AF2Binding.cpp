#include "meshkit/AF2Binding.hpp"
#include "meshkit/Error.hpp"

typedef std::map<const AF2RuleExistEdge*,
    const AF2Edge2D*>::iterator EdgeBindMapItr;
typedef std::map<const AF2RuleExistEdge*,
    const AF2Edge2D*>::const_iterator EdgeBindMapConstItr;

AF2Binding::VtxBindRec::VtxBindRec() :
    pointPtr(NULL), numBoundEdges(0), xplctBound(false)
{
  // do nothing beyond the initializers
}

void AF2Binding::bind(
    const AF2RuleExistEdge* ruleEdgePtr, const AF2Edge2D* ngbhdEdgePtr)
{
  EdgeBindMapItr itr = edgeBindMap.find(ruleEdgePtr);
  if (itr == edgeBindMap.end())
  {
    // the rule edge is not bound to anything
    if (edgesInUse.find(ngbhdEdgePtr) == edgesInUse.end() &&
        isConsistent(ruleEdgePtr->getStart(), ngbhdEdgePtr->getStart()) &&
        isConsistent(ruleEdgePtr->getEnd(), ngbhdEdgePtr->getEnd()))
    {
      // nothing else is bound to the neighborhood edge
      // and binding the endpoints of the edge is consistent, so . . .
      // . . . bind the endpoints of the edge,
      bind(ruleEdgePtr->getStart(), ngbhdEdgePtr->getStart(), false);
      bind(ruleEdgePtr->getEnd(), ngbhdEdgePtr->getEnd(), false);
      // . . . bind the edge itself,
      edgeBindMap[ruleEdgePtr] = ngbhdEdgePtr;
      // . . . and mark the edge as being in use
      edgesInUse.insert(ngbhdEdgePtr);
    }
    else
    {
      // something else is bound to the neighborhood edge already or there
      // is an inconsistency in binding one or more of the endpoints
      throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "Binding the rule edge to the neighborhood edge is inconsistent with existing bindings.");
    }
  }
  else
  {
    if (itr->second == ngbhdEdgePtr)
    {
      // the rule edge is already bound to the neighborhood edge
      // this should imply that the endpoints of the rule edge are implicitly
      // bound to the endpoints of the neighborhood edge

      // There is nothing to do in this case.
    }
    else
    {
      // the rule edge is bound to some other neighborhood edge
      throw MeshKit::Error(MeshKit::MK_BAD_INPUT,
          "The rule edge is already bound to some other neighborhood edge.");
    }
  }
}

void AF2Binding::bind(
    const AF2RuleExistVertex* ruleVertexPtr, const AF2Point2D* ngbhdVertexPtr)
{
  bind(ruleVertexPtr, ngbhdVertexPtr, true);
}

void AF2Binding::bind(const AF2RuleExistVertex* ruleVertexPtr,
    const AF2Point2D* ngbhdVertexPtr, bool isExplicit)
{
  typedef std::map<const AF2RuleExistVertex*,
      AF2Binding::VtxBindRec>::iterator VertexBindMapItr;
  VertexBindMapItr itr = vertexBindMap.find(ruleVertexPtr);
  if (itr == vertexBindMap.end())
  {
    // the rule vertex was not previously bound
    if (verticesInUse.find(ngbhdVertexPtr) == verticesInUse.end())
    {
      // nothing else was bound to the neighborhood vertex
      // insert the binding from the rule vertex to the neighborhood vertex
      VtxBindRec* bindRecordPtr = &(vertexBindMap[ruleVertexPtr]);
      bindRecordPtr->pointPtr = ngbhdVertexPtr;
      if (isExplicit)
      {
        bindRecordPtr->xplctBound = true;
      }
      else
      {
        bindRecordPtr->numBoundEdges = 1;
      }
      // mark the neighborhood vertex as being in use
      verticesInUse.insert(ngbhdVertexPtr);
    }
    else
    {
      // something else is bound to the neighborhood vertex already
      throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "Some other rule vertex is already bound to the neighborhood vertex.");
    }
  }
  else
  {
    VtxBindRec* bindRecordPtr = &(itr->second);
    if (bindRecordPtr->pointPtr == ngbhdVertexPtr)
    {
      // the rule vertex is already bound to the correct neighborhood vertex
      if (isExplicit)
      {
        bindRecordPtr->xplctBound = true;
      }
      else
      {
        ++(bindRecordPtr->numBoundEdges);
      }
    }
    else
    {
      // the rule vertex is already bound to some other neighborhood vertex
      throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "The rule vertex is already bound to some other neighborhood vertex.");
    }
  }
}

const AF2Point2D* AF2Binding::getBoundValue(
    const AF2RuleExistVertex* ruleVertexPtr) const
{
  typedef std::map<const AF2RuleExistVertex*,
      AF2Binding::VtxBindRec>::const_iterator VertexBindMapConstItr;
  VertexBindMapConstItr itr = vertexBindMap.find(ruleVertexPtr);
  if (itr == vertexBindMap.end())
  {
    // the rule vertex is not bound to anything
    return NULL;
  }
  else
  {
    return itr->second.pointPtr;
  }
}

bool AF2Binding::isConsistent(
    const AF2RuleExistEdge* ruleEdgePtr, const AF2Edge2D* ngbhdEdgePtr) const
{
  EdgeBindMapConstItr itr = edgeBindMap.find(ruleEdgePtr);
  if (itr == edgeBindMap.end())
  {
    // the rule edge is not bound to anything
    if (edgesInUse.find(ngbhdEdgePtr) == edgesInUse.end())
    {
      // nothing else is bound to the neighborhood edge
      // check consistency of the endpoints of the edge
      return isConsistent(ruleEdgePtr->getStart(), ngbhdEdgePtr->getStart())
          && isConsistent(ruleEdgePtr->getEnd(), ngbhdEdgePtr->getEnd());
    }
    else
    {
      // something else is bound to the neighborhood edge already
      return false;
    }
  }
  else
  {
    if (itr->second == ngbhdEdgePtr)
    {
      // the rule edge is already bound to the neighborhood edge
      // this should imply that the endpoints of the rule edge are
      // bound to the endpoints of the neighborhood edge
      return true;
    }
    else
    {
      // the rule edge is bound to some other neighborhood edge
      return false;
    }
  }
}

bool AF2Binding::isConsistent(const AF2RuleExistVertex* ruleVertexPtr,
    const AF2Point2D* ngbhdVertexPtr) const
{
  typedef std::map<const AF2RuleExistVertex*,
      AF2Binding::VtxBindRec>::const_iterator VertexBindMapConstItr;
  VertexBindMapConstItr itr = vertexBindMap.find(ruleVertexPtr);
  if (itr == vertexBindMap.end())
  {
    // the rule vertex is not bound to anything
    if (verticesInUse.find(ngbhdVertexPtr) == verticesInUse.end())
    {
      // nothing else is bound to the neighborhood vertex
      return true;
    }
    else
    {
      // something else is bound to the neighborhood vertex already
      return false;
    }
  }
  else
  {
    const VtxBindRec* bindRecordPtr = &(itr->second);
    if (bindRecordPtr->pointPtr == ngbhdVertexPtr)
    {
      // the rule vertex is already bound to the neighborhood vertex
      return true;
    }
    else
    {
      // the rule vertex is bound to some other neighborhood vertex
      return false;
    }
  }
}

void AF2Binding::release(const AF2RuleExistEdge* ruleEdgePtr)
{
  EdgeBindMapItr itr = edgeBindMap.find(ruleEdgePtr);
  if (itr == edgeBindMap.end())
  {
    // the rule edge is not bound to anything
    throw MeshKit::Error(MeshKit::MK_BAD_INPUT,
        "The rule edge is not bound to any neighborhood edge.");
  }
  else
  {
    edgesInUse.erase(itr->second);
    edgeBindMap.erase(itr);
    release(ruleEdgePtr->getEnd(), false);
    release(ruleEdgePtr->getStart(), false);
  }
}

void AF2Binding::release(const AF2RuleExistVertex* ruleVertexPtr)
{
  release(ruleVertexPtr, true);
}

void AF2Binding::release(
    const AF2RuleExistVertex* ruleVertexPtr, bool isExplicit)
{
  typedef std::map<const AF2RuleExistVertex*,
      AF2Binding::VtxBindRec>::iterator VertexBindMapItr;
  VertexBindMapItr itr = vertexBindMap.find(ruleVertexPtr);
  if (itr == vertexBindMap.end())
  {
    // the rule vertex is not bound to anything
    throw MeshKit::Error(MeshKit::MK_BAD_INPUT,
        "The rule vertex is not bound to any neighborhood vertex.");
  }
  else
  {
    VtxBindRec* bindRecordPtr = &(itr->second);
    if (isExplicit)
    {
      if (!bindRecordPtr->xplctBound)
      {
        throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "The rule vertex is not explicitly bound to a neighborhood vertex.  The binding is implicit.  Release the edge(s) instead.");
      }
      bindRecordPtr->xplctBound = false;
    }
    else
    {
      --(bindRecordPtr->numBoundEdges);
    }

    if (!(bindRecordPtr->xplctBound) && bindRecordPtr->numBoundEdges == 0)
    {
      // after performing the requested release operation, the vertex
      // is no longer bound either explicitly or implicitly
      verticesInUse.erase(bindRecordPtr->pointPtr);
      vertexBindMap.erase(itr);
    }
  }
}
