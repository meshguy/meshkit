// C++
#include <stack>

// MeshKit
#include "meshkit/AF2Binding.hpp"
#include "meshkit/AF2Rule.hpp"
#include "meshkit/Error.hpp"

/**
 * A template method for allocating an array to hold the contents of a list
 * of pointers and copying the pointers from the list to the array.
 */
template<class T>
void aF2RuleCopyListToArray(std::list<T*> const & aListOfPtr,
    T** & anArrayOfPtr)
{
  anArrayOfPtr = new T*[aListOfPtr.size()];
  typedef typename std::list<T*>::const_iterator ItrType;
  int indx = 0;
  for (ItrType itr = aListOfPtr.begin(); itr != aListOfPtr.end(); ++itr)
  {
    if (*itr == NULL)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string(
          "AF2Rule constructor arguments may not contain any null pointers.");
      throw badArg;
    }
    anArrayOfPtr[indx] = *itr;
    ++indx;
  }
}

/**
 * A template method for deleting the objects that are pointed to by the
 * pointers in a dynamically allocated array of pointers and then deleting
 * the array itself.
 */
template<class T>
void aF2RuleDeepDeletePtrArray(T** & anArrayOfPtr, int arraySize)
{
  for (int indx = 0; indx < arraySize; ++indx)
  {
    delete anArrayOfPtr[indx];
  }
  delete[] anArrayOfPtr;
}

AF2Rule::AF2Rule(std::list<const AF2RuleExistVertex*> const & ruleVertices,
    std::list<const AF2RuleExistEdge*> const & ruleEdges,
    const AF2FreeZoneDef* freeZoneDef,
    std::list<const AF2RuleNewVertex*> const & ruleNewVertices,
    std::list<const AF2RuleNewEdge*> const & ruleNewEdges,
    std::list<const AF2RuleNewFace*> const & ruleNewFaces) :
    numExVertices(ruleVertices.size()),
    numExEdges(ruleEdges.size()),
    freeZoneDef(freeZoneDef),
    numNewVertices(ruleNewVertices.size()),
    numNewEdges(ruleNewEdges.size()),
    numNewFaces(ruleNewFaces.size())
{
  if (numExEdges < 1)
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("AF2Rule must define at least one existing edge.");
    throw badArg;
  }
  if (numExVertices < 2)
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("AF2Rule must define at least two existing vertices.");
    throw badArg;
  }
  aF2RuleCopyListToArray(ruleVertices, exVertices);
  aF2RuleCopyListToArray(ruleEdges, exEdges);
  aF2RuleCopyListToArray(ruleNewVertices, newVertices);
  aF2RuleCopyListToArray(ruleNewEdges, newEdges);
  aF2RuleCopyListToArray(ruleNewFaces, newFaces);
}

AF2Rule::~AF2Rule()
{
  aF2RuleDeepDeletePtrArray(newFaces, numNewFaces);
  aF2RuleDeepDeletePtrArray(newEdges, numNewEdges);
  aF2RuleDeepDeletePtrArray(newVertices, numNewVertices);
  aF2RuleDeepDeletePtrArray(exEdges, numExEdges);
  aF2RuleDeepDeletePtrArray(exVertices, numExVertices);
}

AF2Rule::AF2Rule(const AF2Rule & toCopy) :
    numExVertices(toCopy.numExVertices),
    numExEdges(toCopy.numExEdges),
    freeZoneDef(toCopy.freeZoneDef->clone()),
    numNewVertices(toCopy.numNewVertices),
    numNewEdges(toCopy.numNewEdges),
    numNewFaces(toCopy.numNewFaces)
{
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Rule copy construction is not supported.");
  throw notImpl;
}

AF2Rule& AF2Rule::operator=(const AF2Rule & rhs)
{
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Rule assignment operator is not supported.");
  throw notImpl;
}

std::map<const AF2RuleExistEdge*, std::list<AF2Edge2D*>*>*
    AF2Rule::findPotentialEdgeMatches(AF2Neighborhood const & ngbhd,
    int matchQuality) const
{
  // TODO: Match the first edge against only one base edge, presumably
  // the first AF2Edge2D in the neighborhood?
  const std::list<AF2Edge2D*>* ngbhdEdges = ngbhd.getEdges2D();
  std::map<const AF2RuleExistEdge*, std::list<AF2Edge2D*>*>* matchMap =
      new std::map<const AF2RuleExistEdge*, std::list<AF2Edge2D*>*>();
  for (unsigned indx = 0; indx < numExEdges; ++indx)
  {
    const AF2RuleExistEdge* ruleEdge = exEdges[indx];
    std::list<AF2Edge2D*>* possMatches = new std::list<AF2Edge2D*>();
    for (std::list<AF2Edge2D*>::const_iterator ngbEdgeItr =
        ngbhdEdges->begin(); ngbEdgeItr != ngbhdEdges->end(); ++ngbEdgeItr)
    {
      if (isMatchingEdge(**ngbEdgeItr, *ruleEdge, matchQuality))
      {
        possMatches->push_back(*ngbEdgeItr);
      }
    }
    (*matchMap)[ruleEdge] = possMatches;
  }

  return matchMap;
}

bool AF2Rule::isMatchingEdge(AF2Edge2D const & edge,
    AF2RuleExistEdge const & ruleEdge, int matchQuality) const
{
  if (!isMatchingVertex(*(edge.getStart()), *(ruleEdge.getStart()),
      matchQuality) || !isMatchingVertex(*(edge.getEnd()),
      *(ruleEdge.getEnd()), matchQuality))
  {
    return false;
  }

  double matchTol = 0.5 + 0.3 * matchQuality;
  return ruleEdge.isMatching(*(edge.getStart()), *(edge.getEnd()), matchTol);
}

bool AF2Rule::isMatchingVertex(AF2Point2D const & point,
    AF2RuleExistVertex const & ruleVertex, int matchQuality) const
{
  double matchTol = 0.5 + 0.3 * matchQuality;
  return ruleVertex.isMatching(point, matchTol);
}

void AF2Rule::applyRule(AF2Neighborhood const & ngbhd, int matchQuality) const
{
  std::map<const AF2RuleExistEdge*, std::list<AF2Edge2D*>*>* matchingEdgesMap =
    findPotentialEdgeMatches(ngbhd, matchQuality);

  // TODO: Check that there is at least one potential match for each edge?

  AF2Binding binding;

  std::stack<std::list<AF2Edge2D*>::const_iterator> edgeMatchItrStack;
  unsigned int edgeToMatchIndx = 0;
  --edgeToMatchIndx;
  bool  consistentMatch = true;
  const AF2RuleExistEdge* edgeToMatch = NULL;
  while (true)
  {
    if (consistentMatch)
    {
      ++edgeToMatchIndx;
      if (edgeToMatchIndx == numExEdges)
      {
        // all edges have been matched
        // TODO: Match free vertices, check free zone, and,
        // if successful, record and assess
      }
      else
      {
        // not all edges have been matched yet, so set up to attempt to
        // match the next edge given the current binding
        edgeToMatch = exEdges[edgeToMatchIndx];
        edgeMatchItrStack.push((*matchingEdgesMap)[edgeToMatch]->begin());
      }
    }
    else
    {
      --edgeToMatchIndx;
      if (edgeToMatchIndx < 0)
      {
        // all edges beyond the base edge have examined all possible matches
        break;
      }
      edgeToMatch = exEdges[edgeToMatchIndx];
      // release the edge binding from the last time edgeToMatch was bound
      binding.release(edgeToMatch);
    }

    consistentMatch = false;
    for (std::list<AF2Edge2D*>::const_iterator itr = edgeMatchItrStack.top();
        itr != (*matchingEdgesMap)[edgeToMatch]->end(); ++itr)
    {
      // pop off the top element of the stack (returns void)
      edgeMatchItrStack.pop();
      // check consistency
      consistentMatch = true;
      if (!binding.isConsistent(edgeToMatch, *itr))
      {
        consistentMatch = false;
      }
      else
      {
        // bind the edge and store the next position of the iterator
        // on the iterator stack for later examination
        binding.bind(edgeToMatch, *itr);
        edgeMatchItrStack.push(++itr);

        // break out of this loop and proceed to bind other edges or vertices
        break;
      }
    }
  }
}
