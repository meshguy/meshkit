#include "meshkit/AF2Algorithm.hpp"

// C++
#include <cstddef>
#include <vector>

// MeshKit
#include "meshkit/AF2DfltRuleAppVisitor.hpp"
#include "meshkit/AF2Edge3D.hpp"
#include "meshkit/AF2RuleApplication.hpp"
#include "meshkit/Error.hpp"

AF2Algorithm::AF2Algorithm(const std::list<const AF2Rule*> & ruleListArg) :
    ruleList(ruleListArg)
{
  // do nothing beyond copying the list of rules in initialization
}

AF2Algorithm::~AF2Algorithm()
{
  // do nothing, the standard deletion of the rule list is sufficient,
  // since this does not take ownership of the rules
}

AF2Algorithm::AF2Algorithm(const AF2Algorithm & toCopy)
{
  // Note: the default implementation would work at this point, but
  //   might not in the future, and there shouldn't be much need to
  //   copy this.
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Algorithm copy construction is not supported.");
  throw notImpl;
}

AF2Algorithm& AF2Algorithm::operator=(const AF2Algorithm & rhs)
{
  // Note: the default implementation would work at this point, but
  //   might not in the future, and there shouldn't be much need to
  //   assign this.
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Algorithm assignment operator is not supported.");
  throw notImpl;
}

AF2AlgorithmResult* AF2Algorithm::execute(
    const AF2LocalTransformMaker* const & transformMaker,
    const double* coords, unsigned int numPoints,
    const unsigned int* edges, unsigned int numEdges) const
{
  typedef std::list<const AF2Rule*>::const_iterator RuleConstItr;

  std::list<AF2Point3D*> allPoints;
  std::list<const AF2Polygon3D*> allFaces;
  AF2Front front;

  // initialize the front
  initFront(front, allPoints, coords, numPoints, edges, numEdges);

  // while the front is not empty and there is still progress
  while (!front.isEmpty() && front.getMaximumQuality() < 25u)
  {
    // select a neighborhood on the advancing front
    AF2Neighborhood* ngbhd = front.selectNeighborhood(transformMaker);
    AF2Edge3D* baselineEdge = ngbhd->getBaselineEdge3D();

    // attempt to apply each of the rules to the neighborhood
    AF2DfltRuleAppVisitor ruleAppVisitor;
    for (RuleConstItr itr = ruleList.begin(); itr != ruleList.end(); ++itr)
    {
      (*itr)->applyRule(*ngbhd,
          baselineEdge->getQualityLevel(), ruleAppVisitor);
    }

    // process the results of attempting to apply rules
    const AF2RuleApplication* bestRuleApp =
        ruleAppVisitor.getBestRuleApplication();
    if (bestRuleApp == NULL)
    {
      // there were no successful rule applications, so decrease the quality
      baselineEdge->decreaseQuality();
    }
    else
    {
      // there was a successful rule application, so process it and advance

      // add any new points added by the rule application
      std::map<const AF2Point2D*, AF2Point3D*> newPointsMap;
      for (unsigned int npi = 0; npi < bestRuleApp->getNumNewPoints(); ++npi)
      {
        processNewPoint(bestRuleApp->getNewPoint(npi),
            ngbhd, newPointsMap, allPoints, front);
      }

      // add the faces added by the rule application
      for (unsigned int nfi = 0; nfi < bestRuleApp->getNumNewFaces(); ++nfi)
      {
        processNewFace(bestRuleApp->getNewFace(nfi),
            ngbhd, newPointsMap, allFaces, front);
      }
    }
  }

  if (front.isEmpty())
  {
    // the advancing front algorithm successfully completed
    return new AF2AlgorithmResult(allPoints, allFaces);
  }

  // the advancing front algorithm failed
  release(allPoints, allFaces);
  return new AF2AlgorithmResult();
}

void AF2Algorithm::initFront(AF2Front & front, std::list<AF2Point3D*> & pntList,
    const double* coords, unsigned int numPoints,
    const unsigned int* edges, unsigned int numEdges) const
{
  // make the point objects and add them to the front
  std::vector<AF2Point3D*> pntVector;
  pntVector.reserve(numPoints);
  for (unsigned int pi = 0; pi < numPoints; ++pi)
  {
    AF2Point3D* point =
        new AF2Point3D(coords[3*pi], coords[3*pi + 1], coords[3*pi + 2]);
    pntVector.push_back(point);
    pntList.push_back(point);
    front.addPoint(point);
  }

  // make the edge objects and gather them in a list of edges
  std::list<AF2Edge3D*> edgeList;
  for (unsigned int ei = 0; ei < numEdges; ++ei)
  {
    if (edges[2*ei] >= numPoints || edges[2*ei + 1] >= numPoints)
    {
      MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
      badArg.set_string("An edge index exceeds the number of points.");
      throw badArg;
    }
    AF2Point3D* edgeStart = pntVector[edges[2*ei]];
    AF2Point3D* edgeEnd = pntVector[edges[2*ei + 1]];
    // mark the endpoints of the edges as on the initial boundary
    // Note: This is not done for isolated points, since the front
    //   does not start to advance from them.
    edgeStart->limitDistanceToBoundary(0);
    edgeEnd->limitDistanceToBoundary(0);
    AF2Edge3D* edge = new AF2Edge3D(edgeStart, edgeEnd);
    edgeList.push_back(edge);
  }

  // initialize the front with the edges
  front.advanceFront(edgeList);
}

void AF2Algorithm::processNewFace(const AF2Polygon2D* newFace2D,
    AF2Neighborhood* & ngbhd,
    std::map<const AF2Point2D*, AF2Point3D*> & newPointsMap,
    std::list<const AF2Polygon3D*> & allFaces, AF2Front & front) const
{
  std::list<const AF2Point3D*> facePoints3D;
  std::list<AF2Edge3D*> edgeList;

  AF2Point3D* firstVertex = NULL;
  AF2Point3D* prevVertex = NULL;
  AF2Point3D* curVertex = NULL;
  for (unsigned int fvi = 0; fvi < newFace2D->getNumVertices(); ++fvi)
  {
    // get the 3-D point corresponding to the 2-D point . . .
    const AF2Point2D* faceVertex2D = newFace2D->getVertex(fvi);
    // . . . either from the neighborhood . . .
    curVertex = ngbhd->getCorrespondingPoint(faceVertex2D);
    // . . . or from the map of new points that were added.
    if (curVertex == NULL)
    {
      curVertex = newPointsMap[faceVertex2D];
    }

    // add it to the list of points for future creation of the 3-D face
    facePoints3D.push_back(curVertex);

    // process the edge
    if (firstVertex == NULL)
    {
      // at the first vertex just store the first vertex and what will be
      // the previous vertex
      firstVertex = curVertex;
      prevVertex = curVertex;
    }
    else
    {
      // create and store an edge from the current to the previous vertex
      edgeList.push_back(new AF2Edge3D(curVertex, prevVertex));
      // update what will be the previous vertex
      prevVertex = curVertex;
    }
  }
  // create and store an edge from the first vertex to the last vertex
  edgeList.push_back(new AF2Edge3D(firstVertex, curVertex));

  // create the three-dimensional face and add it to the list of faces
  allFaces.push_back(new AF2Polygon3D(facePoints3D));

  // advance the front with the edges
  front.advanceFront(edgeList);
}

void AF2Algorithm::processNewPoint(const AF2Point2D* newPoint2D,
    AF2Neighborhood* & ngbhd,
    std::map<const AF2Point2D*, AF2Point3D*> & newPointsMap,
    std::list<AF2Point3D*> & allPoints, AF2Front & front) const
{
  AF2Point3D* newPoint3D = ngbhd->transformPoint(newPoint2D);
  newPointsMap[newPoint2D] = newPoint3D;
  allPoints.push_back(newPoint3D);
  front.addPoint(newPoint3D);
}

void AF2Algorithm::release(std::list<AF2Point3D*> & allPoints,
    std::list<const AF2Polygon3D*> & allFaces) const
{
  typedef std::list<const AF2Polygon3D*>::iterator FaceItr;
  typedef std::list<AF2Point3D*>::iterator PointItr;

  for (FaceItr itr = allFaces.begin(); itr != allFaces.end(); ++itr)
  {
    delete (*itr);
  }

  for (PointItr itr = allPoints.begin(); itr != allPoints.end(); ++itr)
  {
    delete (*itr);
  }
}
