#include "meshkit/AF2Front.hpp"

// C++
#include <algorithm>
#include <cmath>

// temporary for debugging
#include <iostream>

// MeshKit
#include "meshkit/Error.hpp"
#include "meshkit/AF2LocalTransform.hpp"

bool EndPointLess::operator()(const AF2Edge3D* const & oneEdge,
    const AF2Edge3D* const & otherEdge) const
{
  if (oneEdge->getStart() < otherEdge->getStart())
  {
    return true;
  }
  if (oneEdge->getStart() > otherEdge->getStart())
  {
    return false;
  }
  // otherwise the start endpoints are equual
  if (oneEdge->getEnd() < otherEdge->getEnd())
  {
    return true;
  }
  return false;
}

AF2Front::AF2Front() : points(), edges(), qualityCount(4, 0u)
{
}

AF2Front::~AF2Front()
{
  typedef std::set<AF2Edge3D*>::iterator EdgeSetItr;
  for (EdgeSetItr edgeItr = edges.begin(); edgeItr != edges.end(); ++edgeItr)
  {
    delete *edgeItr;
  }
}

AF2Front::AF2Front(const AF2Front & toCopy)
{
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Front copy construction is not supported.");
  throw notImpl;
}

AF2Front& AF2Front::operator=(const AF2Front & rhs)
{
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Front assignment operator is not supported.");
  throw notImpl;
}

void AF2Front::qualityDecreased(const AF2Edge3D* const & anEdge)
{
  unsigned int curQual = anEdge->getQualityLevel();
  if (curQual > qualityCount.size())
  {
    qualityCount.resize(curQual, 0u);
  }
  ++qualityCount[curQual - 1];
  --qualityCount[curQual - 2];
}

void AF2Front::addPoint(AF2Point3D* pointPtr)
{
  if (points.find(pointPtr) != points.end())
  {
    MeshKit::Error duplicate(MeshKit::ErrorCode::MK_BAD_INPUT);
    duplicate.set_string("The point is already on the advancing front.");
    throw duplicate;
  }
  points[pointPtr] = 0;
}

void AF2Front::advanceFront(std::list<AF2Edge3D*> edgeList)
{
  typedef std::map<AF2Point3D*, int>::iterator PointCountItr;
  typedef std::set<AF2Edge3D*>::iterator EdgeSetItr;
  typedef std::list<AF2Edge3D*>::iterator EdgeListItr;

  // verify validity of endpoints, which could throw exception,
  // before doing any other processing
  for (EdgeListItr itr = edgeList.begin(); itr != edgeList.end(); ++itr)
  {
    PointCountItr startPntItr = points.find((*itr)->getStart());
    PointCountItr endPntItr = points.find((*itr)->getEnd());
    if (startPntItr == points.end() || endPntItr == points.end())
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string(
          "One or both endpoints of an edge are not on the advancing front.");
      throw badArg;
    }
    if ((*itr)->getQualityLevel() != 1u)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string(
          "Edges added to the advancing front must have quality level 1.");
      throw badArg;
    }
  }

  // decide which half edges need to be added
  // and which (reverse) half edges need to be removed
  std::list<AF2Edge3D*> edgesToAdd;
  std::list<EdgeSetItr> edgesToRemove;

  for (EdgeListItr itr = edgeList.begin(); itr != edgeList.end(); ++itr)
  {
    AF2Edge3D reverseEdge((*itr)->getEnd(), (*itr)->getStart());
    EdgeSetItr revEdgeItr = edges.find(&reverseEdge);
    if (revEdgeItr != edges.end())
    {
      // prepare for later processing of the reverse half edge that was
      // allocated on the heatp to remove it from the edge set and delete it
      // don't delete it now, since that would make it more difficult to
      //   determine whether its endpoints should be deleted
      edgesToRemove.push_back(revEdgeItr);
      // right now delete the forward half edge that was allocated on the heap
      delete *itr;
    }
    else
    {
      // prepare for later processing of the forward half edge
      // don't add it now, since that would prevent adding both half edges
      //   of an "isolated" edge
      edgesToAdd.push_back(*itr);
    }
  }

  // add the half edges that need to be added
  for (EdgeListItr itr = edgesToAdd.begin(); itr != edgesToAdd.end(); ++itr)
  {
    // add it to the set of edges
    edges.insert(*itr);
    // start tracking the edge's quality level
    ++qualityCount[0];
    (*itr)->setObserver(this);
    // increase the count of the number of edges incident to each endpoint
    ++points[(*itr)->getStart()];
    ++points[(*itr)->getEnd()];
    // update the distance to the boundary (locally)
    // * Note * The distance to the boundary may be globally inaccurate
    // after this update, but this agrees with how NetGen did the update.
    // It would be possible to cascade the distance to boundary updates
    // and ensure globally accuracy, but that might affect a lot of points.
    unsigned int maxDist = 1 + std::max(
        (*itr)->getStart()->getDistanceToBoundary(),
        (*itr)->getEnd()->getDistanceToBoundary());
    (*itr)->getStart()->limitDistanceToBoundary(maxDist);
    (*itr)->getEnd()->limitDistanceToBoundary(maxDist);
  }

  // remove the half edges that need to be removed and any endpoints of
  // those edges that have no incident edges after the edge is removed
  for (std::list<EdgeSetItr>::iterator itr = edgesToRemove.begin();
      itr != edgesToRemove.end(); ++itr)
  {
    EdgeSetItr edgeItr = *itr;
    // decrease the count of the number of edges incident to each endpoint
    --points[(*edgeItr)->getStart()];
    --points[(*edgeItr)->getEnd()];
    // remove endpoints from points if they now have zero incident edges
    if (points[(*edgeItr)->getStart()] == 0u)
    {
      points.erase((*edgeItr)->getStart());
    }
    if (points[(*edgeItr)->getEnd()] == 0u)
    {
      points.erase((*edgeItr)->getEnd());
    }
    // stop tracking the quality level of the edge
    --qualityCount[(*edgeItr)->getQualityLevel() - 1u];
    // delete the edge object that was allocated on the heap
    delete *edgeItr;
    // remove the (now dangling) pointer to the edge from the edge set
    edges.erase(edgeItr);
  }
}

bool AF2Front::isEmpty() const
{
  return edges.empty() && points.empty();
}

unsigned int AF2Front::getMaximumQuality() const
{
  for (unsigned int qualityLevel = 1;
      qualityLevel <= qualityCount.size(); ++qualityLevel)
  {
    if (qualityCount[qualityLevel - 1] > 0)
    {
      return qualityLevel;
    }
  }

  return 0u;
}

AF2Neighborhood* AF2Front::selectNeighborhood(
    const AF2LocalTransformMaker* const & transformMaker) const
{
  typedef std::set<AF2Edge3D*>::const_iterator EdgeSetItr;
  typedef std::map<AF2Point3D*, int>::const_iterator PointCountItr;

  // Select a baseline edge that minimizes a metric value for the metric
  // edge quality level + sum of endpoint distances to boundary.
  // The current implementation loops through all of the edges
  // to find an edge with minimal metric value.  This could be made faster
  // using a C++ priority queue, updating each time the metric value changes.
  AF2Edge3D* baselineEdge = NULL;
  unsigned int minMetricVal;
  for (EdgeSetItr edgeItr = edges.begin(); edgeItr != edges.end(); ++edgeItr)
  {
    if (baselineEdge == NULL)
    {
      minMetricVal = (*edgeItr)->getQualityLevel() +
          (*edgeItr)->getStart()->getDistanceToBoundary() +
          (*edgeItr)->getEnd()->getDistanceToBoundary();
      baselineEdge = *edgeItr;
    }
    else
    {
      unsigned int metricVal = (*edgeItr)->getQualityLevel() +
          (*edgeItr)->getStart()->getDistanceToBoundary() +
          (*edgeItr)->getEnd()->getDistanceToBoundary();
      if (metricVal < minMetricVal)
      {
        metricVal = minMetricVal;
        baselineEdge = *edgeItr;
      }
    }
  }

  // throw an exception if no baseline edge was found
  if (baselineEdge == NULL)
  {
    MeshKit::Error fail(MeshKit::ErrorCode::MK_FAILURE);
    fail.set_string("The advancing front has no edges.");
    throw fail;
  }

  // determine the neighborhood size from the length and quality
  // level of the baseline edge
  AF2Point3D* baseStart = baselineEdge->getStart();
  AF2Point3D* baseEnd = baselineEdge->getEnd();
  double xDiff = baseEnd->getX() - baseStart->getX();
  double yDiff = baseEnd->getY() - baseStart->getY();
  double zDiff = baseEnd->getZ() - baseStart->getZ();
  double sqLen = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff;
  double ngbhdSize = 4.0 * (3u + baselineEdge->getQualityLevel()) *
      (3u + baselineEdge->getQualityLevel()) * sqLen;

  // create lists to hold the points and (other) edges in the neighborhood
  std::list<AF2Point3D*> ngbhdPoints;
  std::list<const AF2Edge3D*> ngbhdEdges;
  // create a set to track points included in the neighorhood and avoid
  // duplicating them if they are endpoints of multiple neighborhood edges
  std::set<AF2Point3D*> ngbhdPointSet;

  // find edges that are close enough to the baseline edge
  // to be included in the neighborhood
  AF2Point3D* edgeStart = NULL;
  AF2Point3D* edgeEnd = NULL;
  double edgeXDiff = 0.0;
  double edgeYDiff = 0.0;
  double edgeZDiff = 0.0;
  double edgeSqLen = 0.0;
  double dotProd = 0.0;
  for (EdgeSetItr edgeItr = edges.begin(); edgeItr != edges.end(); ++edgeItr)
  {
    // skip the baseline edge to avoid duplicating it in the neighborhood
    if ((*edgeItr) == baselineEdge)
    {
      continue;
    }

    // calculate components of vector along the direction of the edge
    // and the squared length of the edge
    edgeStart = (*edgeItr)->getStart();
    edgeEnd = (*edgeItr)->getEnd();
    edgeXDiff = edgeEnd->getX() - edgeStart->getX();
    edgeYDiff = edgeEnd->getY() - edgeStart->getY();
    edgeZDiff = edgeEnd->getZ() - edgeStart->getZ();
    edgeSqLen = edgeXDiff * edgeXDiff +
        edgeYDiff * edgeYDiff + edgeZDiff * edgeZDiff;

    // calculate components of vector from edge start point to base start point
    xDiff = baseStart->getX() - edgeStart->getX();
    yDiff = baseStart->getY() - edgeStart->getY();
    zDiff = baseStart->getZ() - edgeStart->getZ();

    // calculate the dot product of the two vectors
    dotProd = edgeXDiff * xDiff + edgeYDiff * yDiff + edgeZDiff * zDiff;

    // calculate the components of the vector from the baseline edge
    // start point to the closest point on the edge, placing the results
    // in (and reusing the variables) xDiff, yDiff, and zDiff
    if (dotProd < 0)
    {
      // the orthogonal projection of the baseline edge start point onto
      // the line that contains the current edge lies outside the edge
      // and the edge start point is closest to the baseline edge start point

      // the vector from the edge start point to the baseline edge start
      // point is correct, and xDiff, yDiff, and zDiff already has that vector
      // --> do nothing
    }
    else if (dotProd > edgeSqLen)
    {
      // the orthogonal projection of the baseline edge start point onto
      // the line that contains the current edge lies outside the edge
      // and the edge end point is closest to the baseline edge start point

      // the vector from the edge end point to the baseline edge start
      // point is the one that is needed
      xDiff = baseStart->getX() - edgeEnd->getX();
      yDiff = baseStart->getY() - edgeEnd->getY();
      zDiff = baseStart->getZ() - edgeEnd->getZ();
    }
    else
    {
      // the orthogonal projection of the baseline edge start point onto
      // the line that contains the current edge lies on the edge
      // and is the closest point on the edge to the baseline edge start point

      // calculate the components of the projection of the vector from the
      // edge start point to the baseline edge start point onto the vector
      // in the direction of the edge, reusing edge difference variables
      edgeXDiff = dotProd * edgeXDiff / edgeSqLen;
      edgeYDiff = dotProd * edgeYDiff / edgeSqLen;
      edgeZDiff = dotProd * edgeZDiff / edgeSqLen;

      // calculate the components of the projection of the vector from the
      // edge start point to the baseline edge start point onto the plane
      // orthogonal to the direction of the edge
      xDiff = xDiff - edgeXDiff;
      yDiff = yDiff - edgeYDiff;
      zDiff = zDiff - edgeZDiff;
    }

    // calculate the square of the length of said vector, which is also the
    // square of the distance from the edge to the baseline edge start point
    sqLen = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff;

    // if the edge is close enough, add it and its endpoints
    // to the neighborhood
    if (sqLen < ngbhdSize)
    {
      ngbhdEdges.push_back(*edgeItr);
      if (ngbhdPointSet.find(edgeStart) == ngbhdPointSet.end())
      {
        ngbhdPointSet.insert(edgeStart);
        ngbhdPoints.push_back(edgeStart);
      }
      if (ngbhdPointSet.find(edgeEnd) == ngbhdPointSet.end())
      {
        ngbhdPointSet.insert(edgeEnd);
        ngbhdPoints.push_back(edgeEnd);
      }
    }
  }

  // find points that are close enough to the baseline edge
  // to be included in the neighborhood
  for (PointCountItr itr = points.begin(); itr != points.end(); ++itr)
  {
    AF2Point3D* frontPoint = itr->first;

    // calculate the square of the distance from the baseline edge start point
    xDiff = frontPoint->getX() - baseStart->getX();
    yDiff = frontPoint->getY() - baseStart->getY();
    zDiff = frontPoint->getZ() - baseStart->getZ();
    sqLen = xDiff * xDiff + yDiff * yDiff + zDiff * zDiff;

    // if the distance is small enough and the point has not been added
    // to the neighborhood yet, add it to the neighborhood
    if (sqLen < ngbhdSize &&
        ngbhdPointSet.find(frontPoint) == ngbhdPointSet.end())
    {
      ngbhdPointSet.insert(frontPoint);
      ngbhdPoints.push_back(frontPoint);
    }
  }

  // use the transform maker to make a local transform appropriate
  // to the neighborhood
  AF2LocalTransform* localTransform =
      transformMaker->makeLocalTransform(ngbhdPoints, baselineEdge, ngbhdEdges);

  // build and return the neighborhood object
  return new AF2Neighborhood(
      ngbhdPoints, baselineEdge, ngbhdEdges, localTransform);
}
