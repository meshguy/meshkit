#include "meshkit/AF2DfltTriangleRules.hpp"

// C++
#include <list>

// MeshKit
#include "meshkit/AF2FreeZoneDef.hpp"
#include "meshkit/AF2FreeZoneDefLCQualLim.hpp"
#include "meshkit/AF2FreeZoneDefSimple.hpp"
#include "meshkit/AF2PntTrnsfrmLnrV.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2PointTransformNone.hpp"
#include "meshkit/AF2RuleExistEdge.hpp"
#include "meshkit/AF2RuleExistVertex.hpp"
#include "meshkit/AF2RuleNewEdge.hpp"
#include "meshkit/AF2RuleNewFace.hpp"
#include "meshkit/AF2RuleNewTriangle.hpp"
#include "meshkit/AF2RuleNewVertex.hpp"
#include "meshkit/Error.hpp"

AF2DfltTriangleRules::AF2DfltTriangleRules()
{
  ruleList.push_back(make180DegreeRuleQ1());
  ruleList.push_back(make180DegreeRuleQ5());
  ruleList.push_back(make180DegreeRuleQ10());
  ruleList.push_back(make180DegreeRuleQ20());
  ruleList.push_back(make60DegreeAngleRightRule());
  ruleList.push_back(make60DegreeAngleLeftRule());
  ruleList.push_back(make120DegreeAngleRightRule());
  ruleList.push_back(make120DegreeAngleLeftRule());
  ruleList.push_back(make120DegreeAngleBothRule());
  ruleList.push_back(makeFillTriangleRule());
  ruleList.push_back(makeConnectVertexRule());
  ruleList.push_back(makeConnectEdgeRule());
}

AF2DfltTriangleRules::~AF2DfltTriangleRules()
{
  typedef std::list<const AF2Rule*>::iterator RuleListItr;
  for (RuleListItr itr = ruleList.begin(); itr != ruleList.end(); ++itr)
  {
    delete *itr;
  }
}

AF2DfltTriangleRules::AF2DfltTriangleRules(const AF2DfltTriangleRules & toCopy)
{
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string(
      "AF2DfltTriangleRules copy constructor is not supported.");
  throw notImpl;
}

AF2DfltTriangleRules& AF2DfltTriangleRules::operator=(const AF2DfltTriangleRules & rhs)
{
  MeshKit::Error notImpl(MeshKit::ErrorCode::MK_NOT_IMPLEMENTED);
  notImpl.set_string(
      "AF2DfltTriangleRules assignment operator is not supported.");
  throw notImpl;
}

std::list<const AF2Rule*> AF2DfltTriangleRules::getRules() const
{
  return ruleList;
}

AF2PntTrnsfrmLnrV* AF2DfltTriangleRules::makeLinearTransformX(
    const AF2RuleExistVertex* ruleExVert, double xCoeff) const
{
  std::list<const AF2RuleExistVertex*> vertexList;
  std::list<double> xCoeffList;
  std::list<double> yCoeffList;

  vertexList.push_back(ruleExVert);
  xCoeffList.push_back(xCoeff);
  xCoeffList.push_back(0.0);
  yCoeffList.push_back(0.0);
  yCoeffList.push_back(0.0);

  return new AF2PntTrnsfrmLnrV(vertexList, xCoeffList, yCoeffList);
}

AF2PntTrnsfrmLnrV* AF2DfltTriangleRules::makeTranslation(
    const AF2RuleExistVertex* ruleExVert) const
{
  std::list<const AF2RuleExistVertex*> vertexList;
  std::list<double> xCoeffList;
  std::list<double> yCoeffList;

  vertexList.push_back(ruleExVert);
  xCoeffList.push_back(1.0);
  xCoeffList.push_back(0.0);
  yCoeffList.push_back(0.0);
  yCoeffList.push_back(1.0);

  return new AF2PntTrnsfrmLnrV(vertexList, xCoeffList, yCoeffList);
}

const AF2Rule* AF2DfltTriangleRules::make180DegreeRuleQ1() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist

  // assemble the non-baseline edges that must exist into a list
  std::list<const AF2RuleExistEdge*> existEdges;

  // define the reference preferred boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> prefBndryPnts;
  double prefXCoord[5] = {0.0, 1.0, 1.5, 0.5, -0.5};
  double prefYCoord[5] = {0.0, 0.0, 0.7, 1.5, 0.7};
  for (unsigned int pbpi = 0u; pbpi < 5u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    prefBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the preferred
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> prefBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* halfBaseXTransform =
      makeLinearTransformX(baseExist, 0.5);
  prefBndryTrnsfrms.push_back(noTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(halfBaseXTransform);
  prefBndryTrnsfrms.push_back(halfBaseXTransform);
  prefBndryTrnsfrms.push_back(halfBaseXTransform);

  // define the reference limiting boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> limBndryPnts;
  double limXCoord[5] = {0.0, 1.0, 0.5, 0.5, 0.5};
  double limYCoord[5] = {0.0, 0.0, 0.866, 0.866, 0.866};
  for (unsigned int lbpi = 0u; lbpi < 5u; ++lbpi)
  {
    AF2Point2D bndryPnt(limXCoord[lbpi], limYCoord[lbpi]);
    limBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the limiting
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> limBndryTrnsfrms;
  limBndryTrnsfrms.push_back(noTransform);
  limBndryTrnsfrms.push_back(baseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(prefBndryPnts,
      prefBndryTrnsfrms, limBndryPnts, limBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.866);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, halfBaseXTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 2));
  newEdges.push_back(new AF2RuleNewEdge(2, 1));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Triangle Off Line, Quality Level 1", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete halfBaseXTransform;
  delete baseXTransform;
  delete noTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make180DegreeRuleQ5() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist

  // assemble the non-baseline edges that must exist into a list
  std::list<const AF2RuleExistEdge*> existEdges;

  // define the reference preferred boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> prefBndryPnts;
  double prefXCoord[4] = {0.0, 1.0, 1.0, 0.0};
  double prefYCoord[4] = {0.0, 0.0, 0.7, 0.7};
  for (unsigned int pbpi = 0u; pbpi < 4u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    prefBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the preferred
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> prefBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* halfBaseXTransform =
      makeLinearTransformX(baseExist, 0.5);
  prefBndryTrnsfrms.push_back(noTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(noTransform);

  // define the reference limiting boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> limBndryPnts;
  double limXCoord[4] = {0.0, 1.0, 0.5, 0.5};
  double limYCoord[4] = {0.0, 0.0, 0.5, 0.5};
  for (unsigned int lbpi = 0u; lbpi < 4u; ++lbpi)
  {
    AF2Point2D bndryPnt(limXCoord[lbpi], limYCoord[lbpi]);
    limBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the limiting
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> limBndryTrnsfrms;
  limBndryTrnsfrms.push_back(noTransform);
  limBndryTrnsfrms.push_back(baseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(prefBndryPnts,
      prefBndryTrnsfrms, limBndryPnts, limBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.5);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, halfBaseXTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 2));
  newEdges.push_back(new AF2RuleNewEdge(2, 1));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Triangle Off Line, Quality Level 5", 5u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete halfBaseXTransform;
  delete baseXTransform;
  delete noTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make180DegreeRuleQ10() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist

  // assemble the non-baseline edges that must exist into a list
  std::list<const AF2RuleExistEdge*> existEdges;

  // define the reference preferred boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> prefBndryPnts;
  double prefXCoord[4] = {0.0, 1.0, 1.0, 0.0};
  double prefYCoord[4] = {0.0, 0.0, 0.5, 0.5};
  for (unsigned int pbpi = 0u; pbpi < 4u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    prefBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the preferred
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> prefBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* halfBaseXTransform =
      makeLinearTransformX(baseExist, 0.5);
  prefBndryTrnsfrms.push_back(noTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(noTransform);

  // define the reference limiting boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> limBndryPnts;
  double limXCoord[4] = {0.0, 1.0, 0.5, 0.5};
  double limYCoord[4] = {0.0, 0.0, 0.3, 0.3};
  for (unsigned int lbpi = 0u; lbpi < 4u; ++lbpi)
  {
    AF2Point2D bndryPnt(limXCoord[lbpi], limYCoord[lbpi]);
    limBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the limiting
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> limBndryTrnsfrms;
  limBndryTrnsfrms.push_back(noTransform);
  limBndryTrnsfrms.push_back(baseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(prefBndryPnts,
      prefBndryTrnsfrms, limBndryPnts, limBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.3);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, halfBaseXTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 2));
  newEdges.push_back(new AF2RuleNewEdge(2, 1));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Triangle Off Line, Quality Level 10", 10u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete halfBaseXTransform;
  delete baseXTransform;
  delete noTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make180DegreeRuleQ20() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist

  // assemble the non-baseline edges that must exist into a list
  std::list<const AF2RuleExistEdge*> existEdges;

  // define the reference preferred boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> prefBndryPnts;
  double prefXCoord[4] = {0.0, 1.0, 1.0, 0.0};
  double prefYCoord[4] = {0.0, 0.0, 0.2, 0.2};
  for (unsigned int pbpi = 0u; pbpi < 4u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    prefBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the preferred
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> prefBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* halfBaseXTransform =
      makeLinearTransformX(baseExist, 0.5);
  prefBndryTrnsfrms.push_back(noTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(noTransform);

  // define the reference limiting boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> limBndryPnts;
  double limXCoord[4] = {0.0, 1.0, 0.5, 0.5};
  double limYCoord[4] = {0.0, 0.0, 0.1, 0.1};
  for (unsigned int lbpi = 0u; lbpi < 4u; ++lbpi)
  {
    AF2Point2D bndryPnt(limXCoord[lbpi], limYCoord[lbpi]);
    limBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the limiting
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> limBndryTrnsfrms;
  limBndryTrnsfrms.push_back(noTransform);
  limBndryTrnsfrms.push_back(baseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);
  limBndryTrnsfrms.push_back(halfBaseXTransform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(prefBndryPnts,
      prefBndryTrnsfrms, limBndryPnts, limBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.1);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, halfBaseXTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 2));
  newEdges.push_back(new AF2RuleNewEdge(2, 1));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Triangle Off Line, Quality Level 20", 20u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete halfBaseXTransform;
  delete baseXTransform;
  delete noTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make60DegreeAngleRightRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);
  const AF2RuleExistVertex* peakExist = new AF2RuleExistVertex(0.5, 0.866);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(peakExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;
  existEdges.push_back(new AF2RuleExistEdge(baseExist, peakExist));

  // define a linear point transformation that balances the bound
  //   value of the base vertex and the peak vertex
  std::list<const AF2RuleExistVertex*> t0VertexList;
  std::list<double> t0XCoeffList;
  std::list<double> t0YCoeffList;

  t0VertexList.push_back(baseExist);
  t0VertexList.push_back(peakExist);
  t0XCoeffList.push_back(-0.5);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(0.75);
  t0XCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(-0.5);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.75);

  AF2PointTransform* prefV3Transform =
      new AF2PntTrnsfrmLnrV(t0VertexList, t0XCoeffList, t0YCoeffList);

  // define a linear point transformation that uses the bound
  //   value of the peak vertex, but has coefficients of 0.5 instead of 1.0
  std::list<const AF2RuleExistVertex*> t1VertexList;
  std::list<double> t1XCoeffList;
  std::list<double> t1YCoeffList;

  t1VertexList.push_back(peakExist);
  t1XCoeffList.push_back(0.5);
  t1XCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.5);

  AF2PointTransform* limV3Transform =
      new AF2PntTrnsfrmLnrV(t1VertexList, t1XCoeffList, t1YCoeffList);

  // define the reference preferred boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> prefBndryPnts;
  double prefXCoord[4] = {0.0, 1.0, 0.5, -0.125};
  double prefYCoord[4] = {0.0, 0.0, 0.866, 0.6495};
  for (unsigned int pbpi = 0u; pbpi < 4u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    prefBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the preferred
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> prefBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* translateToPeak = makeTranslation(peakExist);
  prefBndryTrnsfrms.push_back(noTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(translateToPeak);
  prefBndryTrnsfrms.push_back(prefV3Transform);

  // define the reference limiting boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> limBndryPnts;
  double limXCoord[4] = {0.0, 1.0, 0.5, 0.25};
  double limYCoord[4] = {0.0, 0.0, 0.866, 0.433};
  for (unsigned int lbpi = 0u; lbpi < 4u; ++lbpi)
  {
    AF2Point2D bndryPnt(limXCoord[lbpi], limYCoord[lbpi]);
    limBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the limiting
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> limBndryTrnsfrms;
  limBndryTrnsfrms.push_back(noTransform);
  limBndryTrnsfrms.push_back(baseXTransform);
  limBndryTrnsfrms.push_back(translateToPeak);
  limBndryTrnsfrms.push_back(limV3Transform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(prefBndryPnts,
      prefBndryTrnsfrms, limBndryPnts, limBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 2));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Angle at Right Vertex of 60 Degrees", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToPeak;
  delete baseXTransform;
  delete noTransform;
  delete limV3Transform;
  delete prefV3Transform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make60DegreeAngleLeftRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);
  const AF2RuleExistVertex* peakExist = new AF2RuleExistVertex(0.5, 0.866);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(peakExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;
  existEdges.push_back(new AF2RuleExistEdge(peakExist, originExist));

  // define a linear point transformation that balances the bound
  //   value of the base vertex and the peak vertex
  std::list<const AF2RuleExistVertex*> t0VertexList;
  std::list<double> t0XCoeffList;
  std::list<double> t0YCoeffList;

  t0VertexList.push_back(baseExist);
  t0VertexList.push_back(peakExist);
  t0XCoeffList.push_back(0.75);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(0.75);
  t0XCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.75);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.75);

  AF2PointTransform* prefV2Transform =
      new AF2PntTrnsfrmLnrV(t0VertexList, t0XCoeffList, t0YCoeffList);

  // define a linear point transformation that translates by the difference
  //   between the bound value and reference value of the midpoint
  //   between the base vertex and peak vertex
  std::list<const AF2RuleExistVertex*> t1VertexList;
  std::list<double> t1XCoeffList;
  std::list<double> t1YCoeffList;

  t1VertexList.push_back(baseExist);
  t1VertexList.push_back(peakExist);
  t1XCoeffList.push_back(0.5);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(0.5);
  t1XCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.5);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.5);

  AF2PointTransform* limV2Transform =
      new AF2PntTrnsfrmLnrV(t1VertexList, t1XCoeffList, t1YCoeffList);

  // define the reference preferred boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> prefBndryPnts;
  double prefXCoord[4] = {0.0, 1.0, 1.125, 0.5};
  double prefYCoord[4] = {0.0, 0.0, 0.6495, 0.866};
  for (unsigned int pbpi = 0u; pbpi < 4u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    prefBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the preferred
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> prefBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* translateToPeak = makeTranslation(peakExist);
  prefBndryTrnsfrms.push_back(noTransform);
  prefBndryTrnsfrms.push_back(baseXTransform);
  prefBndryTrnsfrms.push_back(prefV2Transform);
  prefBndryTrnsfrms.push_back(translateToPeak);

  // define the reference limiting boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> limBndryPnts;
  double limXCoord[4] = {0.0, 1.0, 0.75, 0.5};
  double limYCoord[4] = {0.0, 0.0, 0.433, 0.866};
  for (unsigned int lbpi = 0u; lbpi < 4u; ++lbpi)
  {
    AF2Point2D bndryPnt(limXCoord[lbpi], limYCoord[lbpi]);
    limBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the limiting
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> limBndryTrnsfrms;
  limBndryTrnsfrms.push_back(noTransform);
  limBndryTrnsfrms.push_back(baseXTransform);
  limBndryTrnsfrms.push_back(limV2Transform);
  limBndryTrnsfrms.push_back(translateToPeak);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(prefBndryPnts,
      prefBndryTrnsfrms, limBndryPnts, limBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(2, 1));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Angle at Left Vertex of 60 Degrees", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToPeak;
  delete baseXTransform;
  delete noTransform;
  delete limV2Transform;
  delete prefV2Transform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make120DegreeAngleRightRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);
  const AF2RuleExistVertex* highRightExist =
      new AF2RuleExistVertex(1.5, 0.866);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(highRightExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;
  existEdges.push_back(new AF2RuleExistEdge(baseExist, highRightExist));

  // define a linear point transformation that depends on
  //   the base vertex and highRight vertex
  std::list<const AF2RuleExistVertex*> t0VertexList;
  std::list<double> t0XCoeffList;
  std::list<double> t0YCoeffList;

  t0VertexList.push_back(baseExist);
  t0VertexList.push_back(highRightExist);
  t0XCoeffList.push_back(-1.0);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(1.0);
  t0XCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(-1.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(1.0);

  AF2PointTransform* newVertexTransform =
      new AF2PntTrnsfrmLnrV(t0VertexList, t0XCoeffList, t0YCoeffList);

  // define another linear point transformation that depends on
  //   the base vertex and highRight vertex
  std::list<const AF2RuleExistVertex*> t1VertexList;
  std::list<double> t1XCoeffList;
  std::list<double> t1YCoeffList;

  t1VertexList.push_back(baseExist);
  t1VertexList.push_back(highRightExist);
  t1XCoeffList.push_back(-2.0);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(2.0);
  t1XCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(-2.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(2.0);

  AF2PointTransform* fzV3Transform =
      new AF2PntTrnsfrmLnrV(t1VertexList, t1XCoeffList, t1YCoeffList);

  // define a third linear point transformation that depends on
  //   the base vertex and highRight vertex
  std::list<const AF2RuleExistVertex*> t2VertexList;
  std::list<double> t2XCoeffList;
  std::list<double> t2YCoeffList;

  t2VertexList.push_back(baseExist);
  t2VertexList.push_back(highRightExist);
  t2XCoeffList.push_back(-3.0);
  t2XCoeffList.push_back(0.0);
  t2XCoeffList.push_back(2.0);
  t2XCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(-3.0);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(2.0);

  AF2PointTransform* fzV4Transform =
      new AF2PntTrnsfrmLnrV(t2VertexList, t2XCoeffList, t2YCoeffList);

  // define a fourth linear point transformation that depends on
  //   the base vertex and highRight vertex
  std::list<const AF2RuleExistVertex*> t3VertexList;
  std::list<double> t3XCoeffList;
  std::list<double> t3YCoeffList;

  t3VertexList.push_back(baseExist);
  t3VertexList.push_back(highRightExist);
  t3XCoeffList.push_back(-2.0);
  t3XCoeffList.push_back(0.0);
  t3XCoeffList.push_back(1.0);
  t3XCoeffList.push_back(0.0);
  t3YCoeffList.push_back(0.0);
  t3YCoeffList.push_back(-2.0);
  t3YCoeffList.push_back(0.0);
  t3YCoeffList.push_back(1.0);

  AF2PointTransform* fzV5Transform =
      new AF2PntTrnsfrmLnrV(t3VertexList, t3XCoeffList, t3YCoeffList);

  // define the reference boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> fzBndryPnts;
  double prefXCoord[6] = {0.0, 1.0, 1.5, 1.0, 0.0, -0.5};
  double prefYCoord[6] = {0.0, 0.0, 0.866, 1.732, 1.732, 0.866};
  for (unsigned int pbpi = 0u; pbpi < 6u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    fzBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the boundary points
  //   of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> fzBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* translateToHighRight = makeTranslation(highRightExist);
  fzBndryTrnsfrms.push_back(noTransform);
  fzBndryTrnsfrms.push_back(baseXTransform);
  fzBndryTrnsfrms.push_back(translateToHighRight);
  fzBndryTrnsfrms.push_back(fzV3Transform);
  fzBndryTrnsfrms.push_back(fzV4Transform);
  fzBndryTrnsfrms.push_back(fzV5Transform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef =
      new AF2FreeZoneDefSimple(fzBndryPnts, fzBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.866);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, newVertexTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 3));
  newEdges.push_back(new AF2RuleNewEdge(3, 2));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 3));
  newTriangles.push_back(new AF2RuleNewTriangle(1, 2, 3));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Angle at Right Vertex of 120 Degrees", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToHighRight;
  delete baseXTransform;
  delete noTransform;
  delete fzV5Transform;
  delete fzV4Transform;
  delete fzV3Transform;
  delete newVertexTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make120DegreeAngleLeftRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);
  const AF2RuleExistVertex* highLeftExist =
      new AF2RuleExistVertex(-0.5, 0.866);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(highLeftExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;
  existEdges.push_back(new AF2RuleExistEdge(highLeftExist, originExist));

  // define a linear point transformation that depends on
  //   the base vertex and highLeft vertex
  std::list<const AF2RuleExistVertex*> t0VertexList;
  std::list<double> t0XCoeffList;
  std::list<double> t0YCoeffList;

  t0VertexList.push_back(baseExist);
  t0VertexList.push_back(highLeftExist);
  t0XCoeffList.push_back(1.0);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(1.0);
  t0XCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(1.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(1.0);

  AF2PointTransform* newVertexTransform =
      new AF2PntTrnsfrmLnrV(t0VertexList, t0XCoeffList, t0YCoeffList);

  // define another linear point transformation that depends on
  //   the base vertex and highLeft vertex
  std::list<const AF2RuleExistVertex*> t1VertexList;
  std::list<double> t1XCoeffList;
  std::list<double> t1YCoeffList;

  t1VertexList.push_back(baseExist);
  t1VertexList.push_back(highLeftExist);
  t1XCoeffList.push_back(2.0);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(1.0);
  t1XCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(2.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(1.0);

  AF2PointTransform* fzV2Transform =
      new AF2PntTrnsfrmLnrV(t1VertexList, t1XCoeffList, t1YCoeffList);

  // define a third linear point transformation that depends on
  //   the base vertex and highLeft vertex
  std::list<const AF2RuleExistVertex*> t2VertexList;
  std::list<double> t2XCoeffList;
  std::list<double> t2YCoeffList;

  t2VertexList.push_back(baseExist);
  t2VertexList.push_back(highLeftExist);
  t2XCoeffList.push_back(2.0);
  t2XCoeffList.push_back(0.0);
  t2XCoeffList.push_back(2.0);
  t2XCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(2.0);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(2.0);

  AF2PointTransform* fzV3Transform =
      new AF2PntTrnsfrmLnrV(t2VertexList, t2XCoeffList, t2YCoeffList);

  // define a fourth linear point transformation that depends on
  //   the base vertex and highLeft vertex
  std::list<const AF2RuleExistVertex*> t3VertexList;
  std::list<double> t3XCoeffList;
  std::list<double> t3YCoeffList;

  t3VertexList.push_back(baseExist);
  t3VertexList.push_back(highLeftExist);
  t3XCoeffList.push_back(1.0);
  t3XCoeffList.push_back(0.0);
  t3XCoeffList.push_back(2.0);
  t3XCoeffList.push_back(0.0);
  t3YCoeffList.push_back(0.0);
  t3YCoeffList.push_back(1.0);
  t3YCoeffList.push_back(0.0);
  t3YCoeffList.push_back(2.0);

  AF2PointTransform* fzV4Transform =
      new AF2PntTrnsfrmLnrV(t3VertexList, t3XCoeffList, t3YCoeffList);

  // define the reference boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> fzBndryPnts;
  double prefXCoord[6] = {0.0, 1.0, 1.5, 1.0, 0.0, -0.5};
  double prefYCoord[6] = {0.0, 0.0, 0.866, 1.732, 1.732, 0.866};
  for (unsigned int pbpi = 0u; pbpi < 6u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    fzBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the boundary points
  //   of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> fzBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* translateToHighLeft = makeTranslation(highLeftExist);
  fzBndryTrnsfrms.push_back(noTransform);
  fzBndryTrnsfrms.push_back(baseXTransform);
  fzBndryTrnsfrms.push_back(fzV2Transform);
  fzBndryTrnsfrms.push_back(fzV3Transform);
  fzBndryTrnsfrms.push_back(fzV4Transform);
  fzBndryTrnsfrms.push_back(translateToHighLeft);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef =
      new AF2FreeZoneDefSimple(fzBndryPnts, fzBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.866);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, newVertexTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(2, 3));
  newEdges.push_back(new AF2RuleNewEdge(3, 1));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 3));
  newTriangles.push_back(new AF2RuleNewTriangle(0, 3, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Angle at Left Vertex of 120 Degrees", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToHighLeft;
  delete baseXTransform;
  delete noTransform;
  delete fzV4Transform;
  delete fzV3Transform;
  delete fzV2Transform;
  delete newVertexTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::make120DegreeAngleBothRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);
  const AF2RuleExistVertex* highLeftExist =
      new AF2RuleExistVertex(-0.5, 0.866);
  const AF2RuleExistVertex* highRightExist =
      new AF2RuleExistVertex(1.5, 0.866);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(highLeftExist);
  existVertices.push_back(highRightExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;
  existEdges.push_back(new AF2RuleExistEdge(highLeftExist, originExist));
  existEdges.push_back(new AF2RuleExistEdge(baseExist, highRightExist));

  // define a linear point transformation that translates along the vector
  //   from the reference location to the bound value location of the
  //   midpoint between the highLeft and highRight vertices
  std::list<const AF2RuleExistVertex*> t0VertexList;
  std::list<double> t0XCoeffList;
  std::list<double> t0YCoeffList;

  t0VertexList.push_back(highLeftExist);
  t0VertexList.push_back(highRightExist);
  t0XCoeffList.push_back(0.5);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(0.5);
  t0XCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.5);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.5);

  AF2PointTransform* newVertexTransform =
      new AF2PntTrnsfrmLnrV(t0VertexList, t0XCoeffList, t0YCoeffList);

  // define a linear point transformation that depends on
  //   the base vertex, the highLeft vertex, and the highRight vertex
  std::list<const AF2RuleExistVertex*> t1VertexList;
  std::list<double> t1XCoeffList;
  std::list<double> t1YCoeffList;

  t1VertexList.push_back(baseExist);
  t1VertexList.push_back(highLeftExist);
  t1VertexList.push_back(highRightExist);
  t1XCoeffList.push_back(-0.5);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(0.375);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(1.125);
  t1XCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(-0.5);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.375);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(1.125);

  AF2PointTransform* fzV3Transform =
      new AF2PntTrnsfrmLnrV(t1VertexList, t1XCoeffList, t1YCoeffList);

  // define a linear point transformation that depends on only
  //   the highLeft vertex and the highRight vertex
  std::list<const AF2RuleExistVertex*> t2VertexList;
  std::list<double> t2XCoeffList;
  std::list<double> t2YCoeffList;

  t2VertexList.push_back(highLeftExist);
  t2VertexList.push_back(highRightExist);
  t2XCoeffList.push_back(1.125);
  t2XCoeffList.push_back(0.0);
  t2XCoeffList.push_back(0.375);
  t2XCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(1.125);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.375);

  AF2PointTransform* fzV4Transform =
      new AF2PntTrnsfrmLnrV(t2VertexList, t2XCoeffList, t2YCoeffList);

  // define the reference boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> fzBndryPnts;
  double prefXCoord[6] = {0.0, 1.0, 1.5, 1.0, 0.0, -0.5};
  double prefYCoord[6] = {0.0, 0.0, 0.866, 1.299, 1.299, 0.866};
  for (unsigned int pbpi = 0u; pbpi < 6u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    fzBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the boundary points
  //   of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> fzBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseXTransform =
      makeLinearTransformX(baseExist, 1.0);
  AF2PointTransform* translateToHighRight = makeTranslation(highRightExist);
  AF2PointTransform* translateToHighLeft = makeTranslation(highLeftExist);
  fzBndryTrnsfrms.push_back(noTransform);
  fzBndryTrnsfrms.push_back(baseXTransform);
  fzBndryTrnsfrms.push_back(translateToHighRight);
  fzBndryTrnsfrms.push_back(fzV3Transform);
  fzBndryTrnsfrms.push_back(fzV4Transform);
  fzBndryTrnsfrms.push_back(translateToHighLeft);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef =
      new AF2FreeZoneDefSimple(fzBndryPnts, fzBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.866);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, newVertexTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(2, 4));
  newEdges.push_back(new AF2RuleNewEdge(4, 3));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 4));
  newTriangles.push_back(new AF2RuleNewTriangle(0, 4, 2));
  newTriangles.push_back(new AF2RuleNewTriangle(1, 3, 4));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Angles at Both Vertices of 120 Degrees", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToHighLeft;
  delete translateToHighRight;
  delete baseXTransform;
  delete noTransform;
  delete fzV4Transform;
  delete fzV3Transform;
  delete newVertexTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::makeFillTriangleRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist =
      new AF2RuleExistVertex(1.0, 0.0, 0.5, 0.0, 0.5);
  const AF2RuleExistVertex* peakExist =
      new AF2RuleExistVertex(0.5, 0.866, 0.25, 0.0, 0.25);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(peakExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;
  existEdges.push_back(new AF2RuleExistEdge(baseExist, peakExist));
  existEdges.push_back(new AF2RuleExistEdge(peakExist, originExist));

  // define the reference boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> fzBndryPnts;
  double prefXCoord[3] = {0.0, 1.0, 0.5};
  double prefYCoord[3] = {0.0, 0.0, 0.866};
  for (unsigned int pbpi = 0u; pbpi < 3u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    fzBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the boundary points
  //   of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> fzBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseTransform = makeTranslation(baseExist);
  AF2PointTransform* translateToPeak = makeTranslation(peakExist);
  fzBndryTrnsfrms.push_back(noTransform);
  fzBndryTrnsfrms.push_back(baseTransform);
  fzBndryTrnsfrms.push_back(translateToPeak);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef =
      new AF2FreeZoneDefSimple(fzBndryPnts, fzBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Fill Triangle", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToPeak;
  delete baseTransform;
  delete noTransform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::makeConnectVertexRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);
  const AF2RuleExistVertex* peakExist = new AF2RuleExistVertex(0.5, 0.866);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(peakExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;

  // define a linear point transformation that balances the bound
  //   value of the base vertex and the peak vertex
  std::list<const AF2RuleExistVertex*> t0VertexList;
  std::list<double> t0XCoeffList;
  std::list<double> t0YCoeffList;

  t0VertexList.push_back(baseExist);
  t0VertexList.push_back(peakExist);
  t0XCoeffList.push_back(0.8);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(0.8);
  t0XCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.8);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.8);

  AF2PointTransform* prefV2Transform =
      new AF2PntTrnsfrmLnrV(t0VertexList, t0XCoeffList, t0YCoeffList);

  // define another linear point transformation that balances the bound
  //   value of the base vertex and the peak vertex
  std::list<const AF2RuleExistVertex*> t1VertexList;
  std::list<double> t1XCoeffList;
  std::list<double> t1YCoeffList;

  t1VertexList.push_back(baseExist);
  t1VertexList.push_back(peakExist);
  t1XCoeffList.push_back(-0.6);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(0.8);
  t1XCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(-0.6);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.8);

  AF2PointTransform* prefV4Transform =
      new AF2PntTrnsfrmLnrV(t1VertexList, t1XCoeffList, t1YCoeffList);

  // define a linear point transformation that translates along a
  //   vector from the reference location to the bound value location
  //   of the midpoint between the base vertex and peak vertex
  std::list<const AF2RuleExistVertex*> t2VertexList;
  std::list<double> t2XCoeffList;
  std::list<double> t2YCoeffList;

  t2VertexList.push_back(baseExist);
  t2VertexList.push_back(peakExist);
  t2XCoeffList.push_back(0.5);
  t2XCoeffList.push_back(0.0);
  t2XCoeffList.push_back(0.5);
  t2XCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.5);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.5);

  AF2PointTransform* limV2Transform =
      new AF2PntTrnsfrmLnrV(t2VertexList, t2XCoeffList, t2YCoeffList);

  // define a linear point transformation that translates along a
  //   vector from the reference location to the bound value location
  //   of the midpoint between the origin vertex and peak vertex
  //   (assuming the origin vertex is bound at (0, 0))
  std::list<const AF2RuleExistVertex*> t3VertexList;
  std::list<double> t3XCoeffList;
  std::list<double> t3YCoeffList;

  t3VertexList.push_back(peakExist);
  t3XCoeffList.push_back(0.5);
  t3XCoeffList.push_back(0.0);
  t3YCoeffList.push_back(0.0);
  t3YCoeffList.push_back(0.5);

  AF2PointTransform* limV4Transform =
      new AF2PntTrnsfrmLnrV(t3VertexList, t3XCoeffList, t3YCoeffList);

  // define the reference preferred boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> prefBndryPnts;
  double prefXCoord[5] = {0.0, 1.0, 1.2, 0.5, -0.2};
  double prefYCoord[5] = {0.0, 0.0, 0.693, 0.866, 0.693};
  for (unsigned int pbpi = 0u; pbpi < 5u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    prefBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the preferred
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> prefBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseTransform = makeTranslation(baseExist);
  AF2PointTransform* translateToPeak = makeTranslation(peakExist);
  prefBndryTrnsfrms.push_back(noTransform);
  prefBndryTrnsfrms.push_back(baseTransform);
  prefBndryTrnsfrms.push_back(prefV2Transform);
  prefBndryTrnsfrms.push_back(translateToPeak);
  prefBndryTrnsfrms.push_back(prefV4Transform);

  // define the reference limiting boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> limBndryPnts;
  double limXCoord[5] = {0.0, 1.0, 0.75, 0.5, 0.25};
  double limYCoord[5] = {0.0, 0.0, 0.433, 0.866, 0.433};
  for (unsigned int lbpi = 0u; lbpi < 5u; ++lbpi)
  {
    AF2Point2D bndryPnt(limXCoord[lbpi], limYCoord[lbpi]);
    limBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the limiting
  //   boundary points of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> limBndryTrnsfrms;
  limBndryTrnsfrms.push_back(noTransform);
  limBndryTrnsfrms.push_back(baseTransform);
  limBndryTrnsfrms.push_back(limV2Transform);
  limBndryTrnsfrms.push_back(translateToPeak);
  limBndryTrnsfrms.push_back(limV4Transform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(prefBndryPnts,
      prefBndryTrnsfrms, limBndryPnts, limBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 2));
  newEdges.push_back(new AF2RuleNewEdge(2, 1));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 2));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Connect to Opposite Vertex", 1u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToPeak;
  delete baseTransform;
  delete noTransform;
  delete limV4Transform;
  delete limV2Transform;
  delete prefV4Transform;
  delete prefV2Transform;

  return rule;
}

const AF2Rule* AF2DfltTriangleRules::makeConnectEdgeRule() const
{
  // define the vertices that must exist
  const AF2RuleExistVertex* originExist = new AF2RuleExistVertex(0.0, 0.0);
  const AF2RuleExistVertex* baseExist = new AF2RuleExistVertex(1.0, 0.0);
  const AF2RuleExistVertex* oppRightExist = new AF2RuleExistVertex(1.0, 1.732);
  const AF2RuleExistVertex* oppLeftExist = new AF2RuleExistVertex(0.0, 1.732);

  // assemble the vertices that must exist into a list
  std::list<const AF2RuleExistVertex*> existVertices;
  existVertices.push_back(originExist);
  existVertices.push_back(baseExist);
  existVertices.push_back(oppRightExist);
  existVertices.push_back(oppLeftExist);

  // define the baseline edge that must exist
  const AF2RuleExistEdge* baseEdge =
      new AF2RuleExistEdge(originExist, baseExist);

  // define the other edges that must exist and assemble them into a list
  std::list<const AF2RuleExistEdge*> existEdges;
  existEdges.push_back(new AF2RuleExistEdge(oppRightExist, oppLeftExist));

  // define a linear point transformation that balances the bound values of
  //   the base vertex, the opposite right vertex, and the opposite left vertex
  std::list<const AF2RuleExistVertex*> t0VertexList;
  std::list<double> t0XCoeffList;
  std::list<double> t0YCoeffList;

  t0VertexList.push_back(baseExist);
  t0VertexList.push_back(oppRightExist);
  t0VertexList.push_back(oppLeftExist);
  t0XCoeffList.push_back(0.25);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(0.25);
  t0XCoeffList.push_back(0.0);
  t0XCoeffList.push_back(0.25);
  t0XCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.25);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.25);
  t0YCoeffList.push_back(0.0);
  t0YCoeffList.push_back(0.25);

  AF2PointTransform* newVertexTransform =
      new AF2PntTrnsfrmLnrV(t0VertexList, t0XCoeffList, t0YCoeffList);

  // define another linear point transformation that balances bound values of
  //   the base vertex, the opposite right vertex, and the opposite left vertex
  std::list<const AF2RuleExistVertex*> t1VertexList;
  std::list<double> t1XCoeffList;
  std::list<double> t1YCoeffList;

  t1VertexList.push_back(baseExist);
  t1VertexList.push_back(oppRightExist);
  t1VertexList.push_back(oppLeftExist);
  t1XCoeffList.push_back(0.75);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(0.75);
  t1XCoeffList.push_back(0.0);
  t1XCoeffList.push_back(-0.25);
  t1XCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.75);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(0.75);
  t1YCoeffList.push_back(0.0);
  t1YCoeffList.push_back(-0.25);


  AF2PointTransform* fzV2Transform =
      new AF2PntTrnsfrmLnrV(t1VertexList, t1XCoeffList, t1YCoeffList);

  // define a third linear point transformation that balances bound values of
  //   the base vertex, the opposite right vertex, and the opposite left vertex
  std::list<const AF2RuleExistVertex*> t2VertexList;
  std::list<double> t2XCoeffList;
  std::list<double> t2YCoeffList;

  t2VertexList.push_back(baseExist);
  t2VertexList.push_back(oppRightExist);
  t2VertexList.push_back(oppLeftExist);
  t2XCoeffList.push_back(-0.25);
  t2XCoeffList.push_back(0.0);
  t2XCoeffList.push_back(-0.25);
  t2XCoeffList.push_back(0.0);
  t2XCoeffList.push_back(0.75);
  t2XCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(-0.25);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(-0.25);
  t2YCoeffList.push_back(0.0);
  t2YCoeffList.push_back(0.75);


  AF2PointTransform* fzV5Transform =
      new AF2PntTrnsfrmLnrV(t2VertexList, t2XCoeffList, t2YCoeffList);

  // define the reference boundary points of the free zone
  //   and assemble them into a list
  std::list<AF2Point2D> fzBndryPnts;
  double prefXCoord[6] = {0.0, 1.0, 1.5, 1.0, 0.0, -0.5};
  double prefYCoord[6] = {0.0, 0.0, 0.866, 1.732, 1.732, 0.866};
  for (unsigned int pbpi = 0u; pbpi < 6u; ++pbpi)
  {
    AF2Point2D bndryPnt(prefXCoord[pbpi], prefYCoord[pbpi]);
    fzBndryPnts.push_back(bndryPnt);
  }

  // define the point transforms that apply to the boundary points
  //   of the free zone and assemble them into a list
  std::list<const AF2PointTransform*> fzBndryTrnsfrms;
  AF2PointTransform* noTransform = new AF2PointTransformNone();
  AF2PointTransform* baseTransform = makeTranslation(baseExist);
  AF2PointTransform* translateToOppRight = makeTranslation(oppRightExist);
  AF2PointTransform* translateToOppLeft = makeTranslation(oppLeftExist);
  fzBndryTrnsfrms.push_back(noTransform);
  fzBndryTrnsfrms.push_back(baseTransform);
  fzBndryTrnsfrms.push_back(fzV2Transform);
  fzBndryTrnsfrms.push_back(translateToOppRight);
  fzBndryTrnsfrms.push_back(translateToOppLeft);
  fzBndryTrnsfrms.push_back(fzV5Transform);

  // define the free zone definition
  AF2FreeZoneDef* freeZoneDef =
      new AF2FreeZoneDefSimple(fzBndryPnts, fzBndryTrnsfrms);

  // define the new vertices that the rule would create and assemble them
  //   into a list
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D nvRefPoint(0.5, 0.866);
  newVertices.push_back(new AF2RuleNewVertex(nvRefPoint, newVertexTransform));

  // define the new edges and assemble them into a list
  std::list<const AF2RuleNewEdge*> newEdges;
  newEdges.push_back(new AF2RuleNewEdge(0, 4));
  newEdges.push_back(new AF2RuleNewEdge(2, 4));
  newEdges.push_back(new AF2RuleNewEdge(4, 1));
  newEdges.push_back(new AF2RuleNewEdge(4, 3));

  // define the new triangles and assemble them into a list
  std::list<const AF2RuleNewFace*> newTriangles;
  newTriangles.push_back(new AF2RuleNewTriangle(0, 1, 4));
  newTriangles.push_back(new AF2RuleNewTriangle(2, 3, 4));

  // construct the rule itself
  AF2Rule* rule =
      new AF2Rule("Connect to Opposite Edge", 3u, existVertices,
      baseEdge, existEdges, freeZoneDef, newVertices, newEdges, newTriangles);

  // delete transforms that were allocated with new but cloned when
  //   made part of the rule definition
  delete translateToOppLeft;
  delete translateToOppRight;
  delete baseTransform;
  delete noTransform;
  delete fzV5Transform;
  delete fzV2Transform;
  delete newVertexTransform;

  return rule;
}
