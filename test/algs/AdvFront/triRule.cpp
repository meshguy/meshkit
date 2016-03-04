/**
 * \file triRule.cpp \test
 *
 * Test advancing front rule for meshing with triangles.
 *
 */

// C++
#include <iostream>

// MeshKit
#include "meshkit/AF2FreeZoneDefLCQualLim.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2PointTransformNone.hpp"
#include "meshkit/AF2Rule.hpp"
#include "meshkit/AF2RuleNewTriangle.hpp"

// MeshKit testing
#include "TestUtil.hpp"

AF2Rule* makeFreeTriangleRule();
void testFreeTriangleRule();

int main(int argc, char **argv)
{

  // start up MK and load the geometry
  int num_fail = 0;

  num_fail += RUN_TEST(testFreeTriangleRule);

  return num_fail;
}

AF2Rule* makeFreeTriangleRule()
{
  // existing vertices
  std::list<const AF2RuleExistVertex*> exVertices;
  AF2RuleExistVertex* originVertexPtr = new AF2RuleExistVertex(0, 0);
  exVertices.push_back(originVertexPtr);
  AF2RuleExistVertex* baseVertexPtr = new AF2RuleExistVertex(1, 0);
  exVertices.push_back(baseVertexPtr);

  // existing edge
  std::list<const AF2RuleExistEdge*> exEdges;
  AF2RuleExistEdge* exEdgePtr =
      new AF2RuleExistEdge(originVertexPtr, baseVertexPtr);
  exEdges.push_back(exEdgePtr);

  // free zone definition
  // free zone definition lists
  std::list<AF2Point2D> prefPnts;
  std::list<const AF2PointTransform*> prefPntTransforms;
  std::list<AF2Point2D> limitPnts;
  std::list<const AF2PointTransform*> limitPntTransforms;

  // first element free zone definition lists
  AF2Point2D alpha(0, 0);
  AF2PointTransform* alphaTransform = new AF2PointTransformNone();
  prefPnts.push_back(alpha);
  prefPntTransforms.push_back(alphaTransform);
  limitPnts.push_back(alpha);
  limitPntTransforms.push_back(alphaTransform);

  // second element free zone definition lists
  AF2Point2D bravo(1, 0);
  // TODO: Replace this with the appropriate transform +1.0*x_{1}
  AF2PointTransform* bravoTransform = new AF2PointTransformNone();
  prefPnts.push_back(bravo);
  prefPntTransforms.push_back(bravoTransform);
  limitPnts.push_back(bravo);
  limitPntTransforms.push_back(bravoTransform);

  // third element free zone definition lists
  AF2Point2D charliePref(1.5, 0.7);
  // TODO: Replace this with the appropriate transform +0.5*x_{1}
  AF2PointTransform* charlieTransform = new AF2PointTransformNone();
  AF2Point2D charlieLimit(0.5, 0.866);
  prefPnts.push_back(charliePref);
  prefPntTransforms.push_back(charlieTransform);
  limitPnts.push_back(charlieLimit);
  limitPntTransforms.push_back(charlieTransform);

  // fourth element free zone definition lists
  AF2Point2D deltaPref(0.5, 1.5);
  // TODO: Replace this with the appropriate transform 0.5*x_{1}
  AF2PointTransform* deltaTransform = new AF2PointTransformNone();
  AF2Point2D deltaLimit(0.5, 0.866);
  prefPnts.push_back(deltaPref);
  prefPntTransforms.push_back(deltaTransform);
  limitPnts.push_back(deltaLimit);
  limitPntTransforms.push_back(deltaTransform);

  // fifth and final element free zone definition lists
  AF2Point2D echoPref(-0.5, 0.7);
  // TODO: Replace this with the appropriate transform 0.5*x_{1}
  AF2PointTransform* echoTransform = new AF2PointTransformNone();
  AF2Point2D echoLimit(0.5, 0.866);
  prefPnts.push_back(echoPref);
  prefPntTransforms.push_back(echoTransform);
  limitPnts.push_back(echoLimit);
  limitPntTransforms.push_back(echoTransform);

  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefLCQualLim(
      prefPnts, prefPntTransforms, limitPnts, limitPntTransforms);
  delete alphaTransform;
  delete bravoTransform;
  delete charlieTransform;
  delete deltaTransform;
  delete echoTransform;

  // new vertex
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D newVertexLoc(0.5, 0.866);
  // TODO: Replace this with the appropriate transform +0.5*x_{1}
  AF2PointTransform* newVertexTransform = new AF2PointTransformNone();
  AF2RuleNewVertex* newVertexPtr =
      new AF2RuleNewVertex(newVertexLoc, newVertexTransform);
  delete newVertexTransform;
  newVertices.push_back(newVertexPtr);

  // new edges
  std::list<const AF2RuleNewEdge*> newEdges;
  // TODO: Change constructor to use vertex pointers
  AF2RuleNewEdge* newEdgePtr = new AF2RuleNewEdge(0, 2);
  newEdges.push_back(newEdgePtr);
  // TODO: Change constructor to use vertex pointers
  newEdgePtr = new AF2RuleNewEdge(2, 1);
  newEdges.push_back(newEdgePtr);

  // new face
  std::list<const AF2RuleNewFace*> newFaces;
  // TODO: Change constructor to use vertex pointers
  AF2RuleNewFace* newFacePtr = new AF2RuleNewTriangle(0, 1, 2);
  newFaces.push_back(newFacePtr);

  AF2Rule* rulePtr = new AF2Rule(exVertices, exEdges, freeZoneDef,
      newVertices, newEdges, newFaces);

  return rulePtr;
}

void testFreeTriangleRule()
{
  AF2Rule* freeTriRule = makeFreeTriangleRule();
  CHECK(freeTriRule != NULL);
  delete freeTriRule;
}
