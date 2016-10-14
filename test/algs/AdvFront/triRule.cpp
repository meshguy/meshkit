/**
 * \file triRule.cpp \test
 *
 * Test advancing front rule for meshing with triangles.
 *
 */

// C++
#include <iostream>
#include <vector>

// MeshKit
#include "meshkit/AF2FreeZoneDefLCQualLim.hpp"
#include "meshkit/AF2FreeZoneDefSimple.hpp"
#include "meshkit/AF2Edge3D.hpp"
#include "meshkit/AF2Neighborhood.hpp"
#include "meshkit/AF2PlaneProjection.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/AF2PointTransformNone.hpp"
#include "meshkit/AF2PntTrnsfrmLnrV.hpp"
#include "meshkit/AF2Rule.hpp"
#include "meshkit/AF2RuleAppVisitor.hpp"
#include "meshkit/AF2RuleNewTriangle.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"

// MeshKit testing
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model
#if HAVE_OCC
#define FILE_EXT "stp"
#else
#define FILE_EXT "facet"
#endif

AF2Rule* makeFillTriangleRule();
AF2Rule* makeFreeTriangleRule();
AF2Rule* makeFreeTriangleRuleMaxQual5();
AF2Rule* makeMultiPointTriangleRule();
void testFillTriangleRule();
void testFreeTriangleRule();
void testMultiPointTriangleRule();

class TestVisitor : public AF2RuleAppVisitor
{
  private:
    unsigned int count;
    double firstX;
  public:
    TestVisitor()
    {
      count = 0u;
    }

    void visit(AF2RuleApplication const & ruleApp)
    {
      ++count;
      std::cout << "Visit #" << count << std::endl;
      std::cout << "  Number of new faces: " << ruleApp.getNumNewFaces()
          << "\n  Number of new points: " << ruleApp.getNumNewPoints() 
          << std::endl;
      if (ruleApp.getNumNewPoints() > 0)
      {
        firstX = ruleApp.getNewPoint(0u)->getX();
      }
    }

    double getFirstX() const
    {
      return firstX;
    }

    unsigned int getVisitCount() const
    {
      return count;
    }
};

// This variable is at global scope because (1) calling deleteAll on
// the MKCore geometry instance appears to cause memory inconsistencies
// with later use of the geometry instance and (2) it is more efficient
// to load the geometry model only once
MeshKit::ModelEnt* square = NULL;

int main(int argc, char **argv)
{
  // This variable is defined and used in main because a new MKCore
  // instance cannot be easily constructed after another MKCore
  // instance is deleted; there are problems with a tag left behind in
  // iGeom.
  MeshKit::MKCore* mk = new MeshKit::MKCore();

  // load a square in plane z = 0.5 with -1.0 <= x <= 0 and -0.5 <= y <= 0.5
  std::string file_name = TestDir + "/squaresurf." + FILE_EXT;
  mk->load_geometry_mesh(file_name.c_str(), file_name.c_str());
  MeshKit::MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);
  square = *surfs.begin();

  // start up MK and load the geometry
  int num_fail = 0;

  num_fail += RUN_TEST(testFreeTriangleRule);
  num_fail += RUN_TEST(testFillTriangleRule);
  num_fail += RUN_TEST(testMultiPointTriangleRule);

  delete mk;

  return num_fail;
}

AF2Rule* makeFillTriangleRule()
{
  // existing vertices
  std::list<const AF2RuleExistVertex*> exVertices;
  AF2RuleExistVertex* originVertexPtr = new AF2RuleExistVertex(0, 0);
  exVertices.push_back(originVertexPtr);
  AF2RuleExistVertex* baseVertexPtr = new AF2RuleExistVertex(1, 0);
  exVertices.push_back(baseVertexPtr);
  AF2RuleExistVertex* peakVertexPtr =
      new AF2RuleExistVertex(0.5, 0.8);
  exVertices.push_back(peakVertexPtr);

  // existing edges
  AF2RuleExistEdge* baseEdgePtr =
      new AF2RuleExistEdge(originVertexPtr, baseVertexPtr);
  std::list<const AF2RuleExistEdge*> exEdges;
  exEdges.push_back(new AF2RuleExistEdge(baseVertexPtr, peakVertexPtr));
  exEdges.push_back(new AF2RuleExistEdge(peakVertexPtr, originVertexPtr));

  // free zone definition
  // free zone definition lists
  std::list<AF2Point2D> bndryPnts;
  std::list<const AF2PointTransform*> bndryPntTransforms;

  // first element free zone definition lists
  AF2Point2D alpha(0, 0);
  AF2PointTransform* alphaTransform = new AF2PointTransformNone();
  bndryPnts.push_back(alpha);
  bndryPntTransforms.push_back(alphaTransform);

  // second element free zone definition lists
  AF2Point2D bravo(1, 0);
  AF2PointTransform* bravoTransform = new AF2PointTransformNone();
  bndryPnts.push_back(bravo);
  bndryPntTransforms.push_back(bravoTransform);

  // third element free zone definition lists
  AF2Point2D charlie(0.5, 0.8);
  AF2PointTransform* charlieTransform = new AF2PointTransformNone();
  bndryPnts.push_back(charlie);
  bndryPntTransforms.push_back(charlieTransform);

  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefSimple(
      bndryPnts, bndryPntTransforms);
  delete alphaTransform;
  delete bravoTransform;
  delete charlieTransform;

  // no new vertices
  std::list<const AF2RuleNewVertex*> newVertices;

  // no new edges
  std::list<const AF2RuleNewEdge*> newEdges;

  // new face
  std::list<const AF2RuleNewFace*> newFaces;
  AF2RuleNewFace* newFacePtr = new AF2RuleNewTriangle(0, 1, 2);
  newFaces.push_back(newFacePtr);

  AF2Rule* rulePtr = new AF2Rule("triRule Fill", 1u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
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
  AF2RuleExistEdge* baseEdgePtr =
      new AF2RuleExistEdge(originVertexPtr, baseVertexPtr);
  std::list<const AF2RuleExistEdge*> exEdges;

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
  std::list<const AF2RuleExistVertex*> trnsfrmVertices;
  trnsfrmVertices.push_back(baseVertexPtr);
  std::list<double> xCoeff;
  xCoeff.push_back(1.0);
  xCoeff.push_back(0.0);
  std::list<double> yCoeff;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* bravoTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  prefPnts.push_back(bravo);
  prefPntTransforms.push_back(bravoTransform);
  limitPnts.push_back(bravo);
  limitPntTransforms.push_back(bravoTransform);

  // third element free zone definition lists
  AF2Point2D charliePref(1.5, 0.7);
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* charlieTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  AF2Point2D charlieLimit(0.5, 0.866);
  prefPnts.push_back(charliePref);
  prefPntTransforms.push_back(charlieTransform);
  limitPnts.push_back(charlieLimit);
  limitPntTransforms.push_back(charlieTransform);

  // fourth element free zone definition lists
  AF2Point2D deltaPref(0.5, 1.5);
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* deltaTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  AF2Point2D deltaLimit(0.5, 0.866);
  prefPnts.push_back(deltaPref);
  prefPntTransforms.push_back(deltaTransform);
  limitPnts.push_back(deltaLimit);
  limitPntTransforms.push_back(deltaTransform);

  // fifth and final element free zone definition lists
  AF2Point2D echoPref(-0.5, 0.7);
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* echoTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
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
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* newVertexTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  AF2RuleNewVertex* newVertexPtr =
      new AF2RuleNewVertex(newVertexLoc, newVertexTransform);
  delete newVertexTransform;
  newVertices.push_back(newVertexPtr);

  // new edges
  std::list<const AF2RuleNewEdge*> newEdges;
  AF2RuleNewEdge* newEdgePtr = new AF2RuleNewEdge(0, 2);
  newEdges.push_back(newEdgePtr);
  newEdgePtr = new AF2RuleNewEdge(2, 1);
  newEdges.push_back(newEdgePtr);

  // new face
  std::list<const AF2RuleNewFace*> newFaces;
  AF2RuleNewFace* newFacePtr = new AF2RuleNewTriangle(0, 1, 2);
  newFaces.push_back(newFacePtr);

  AF2Rule* rulePtr = new AF2Rule("triRule Free", 1u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

AF2Rule* makeFreeTriangleRuleMaxQual5()
{
  // existing vertices
  std::list<const AF2RuleExistVertex*> exVertices;
  AF2RuleExistVertex* originVertexPtr = new AF2RuleExistVertex(0, 0);
  exVertices.push_back(originVertexPtr);
  AF2RuleExistVertex* baseVertexPtr = new AF2RuleExistVertex(1, 0);
  exVertices.push_back(baseVertexPtr);

  // existing edge
  AF2RuleExistEdge* baseEdgePtr =
      new AF2RuleExistEdge(originVertexPtr, baseVertexPtr);
  std::list<const AF2RuleExistEdge*> exEdges;

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
  std::list<const AF2RuleExistVertex*> trnsfrmVertices;
  trnsfrmVertices.push_back(baseVertexPtr);
  std::list<double> xCoeff;
  xCoeff.push_back(1.0);
  xCoeff.push_back(0.0);
  std::list<double> yCoeff;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* bravoTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  prefPnts.push_back(bravo);
  prefPntTransforms.push_back(bravoTransform);
  limitPnts.push_back(bravo);
  limitPntTransforms.push_back(bravoTransform);

  // third element free zone definition lists
  AF2Point2D charliePref(1.5, 0.7);
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* charlieTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  AF2Point2D charlieLimit(0.5, 0.866);
  prefPnts.push_back(charliePref);
  prefPntTransforms.push_back(charlieTransform);
  limitPnts.push_back(charlieLimit);
  limitPntTransforms.push_back(charlieTransform);

  // fourth element free zone definition lists
  AF2Point2D deltaPref(0.5, 1.5);
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* deltaTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  AF2Point2D deltaLimit(0.5, 0.866);
  prefPnts.push_back(deltaPref);
  prefPntTransforms.push_back(deltaTransform);
  limitPnts.push_back(deltaLimit);
  limitPntTransforms.push_back(deltaTransform);

  // fifth and final element free zone definition lists
  AF2Point2D echoPref(-0.5, 0.7);
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* echoTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
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
  trnsfrmVertices.clear();
  trnsfrmVertices.push_back(baseVertexPtr);
  xCoeff.clear();
  xCoeff.push_back(0.5);
  xCoeff.push_back(0.0);
  yCoeff.clear();;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PointTransform* newVertexTransform =
      new AF2PntTrnsfrmLnrV(trnsfrmVertices, xCoeff, yCoeff);
  AF2RuleNewVertex* newVertexPtr =
      new AF2RuleNewVertex(newVertexLoc, newVertexTransform);
  delete newVertexTransform;
  newVertices.push_back(newVertexPtr);

  // new edges
  std::list<const AF2RuleNewEdge*> newEdges;
  AF2RuleNewEdge* newEdgePtr = new AF2RuleNewEdge(0, 2);
  newEdges.push_back(newEdgePtr);
  newEdgePtr = new AF2RuleNewEdge(2, 1);
  newEdges.push_back(newEdgePtr);

  // new face
  std::list<const AF2RuleNewFace*> newFaces;
  AF2RuleNewFace* newFacePtr = new AF2RuleNewTriangle(0, 1, 2);
  newFaces.push_back(newFacePtr);

  AF2Rule* rulePtr = new AF2Rule("triRule Free (maxq5)", 5u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

AF2Rule* makeMultiPointTriangleRule()
{
  // existing vertices
  std::list<const AF2RuleExistVertex*> exVertices;
  AF2RuleExistVertex* originVertexPtr = new AF2RuleExistVertex(0, 0);
  exVertices.push_back(originVertexPtr);
  AF2RuleExistVertex* baseVertexPtr = new AF2RuleExistVertex(1, 0);
  exVertices.push_back(baseVertexPtr);
  AF2RuleExistVertex* rightPeakVertexPtr = new AF2RuleExistVertex(0.75, 0.8);
  exVertices.push_back(rightPeakVertexPtr);
  AF2RuleExistVertex* leftPeakVertexPtr = new AF2RuleExistVertex(-0.25, 0.8);
  exVertices.push_back(leftPeakVertexPtr);

  // existing edges
  AF2RuleExistEdge* baseEdgePtr =
      new AF2RuleExistEdge(originVertexPtr, baseVertexPtr);
  std::list<const AF2RuleExistEdge*> exEdges;

  // free zone definition
  // free zone definition lists
  std::list<AF2Point2D> bndryPnts;
  std::list<const AF2PointTransform*> bndryPntTransforms;

  // first element free zone definition lists
  AF2Point2D alpha(0, 0);
  AF2PointTransform* alphaTransform = new AF2PointTransformNone();
  bndryPnts.push_back(alpha);
  bndryPntTransforms.push_back(alphaTransform);

  // second element free zone definition lists
  AF2Point2D bravo(1.0, 0);
  AF2PointTransform* bravoTransform = new AF2PointTransformNone();
  bndryPnts.push_back(bravo);
  bndryPntTransforms.push_back(bravoTransform);

  // third element free zone definition lists
  AF2Point2D charlie(1.0, 0.79999);
  AF2PointTransform* charlieTransform = new AF2PointTransformNone();
  bndryPnts.push_back(charlie);
  bndryPntTransforms.push_back(charlieTransform);

  // fourth element free zone definition lists
  AF2Point2D delta(-0.5, 0.79999);
  AF2PointTransform* deltaTransform = new AF2PointTransformNone();
  bndryPnts.push_back(delta);
  bndryPntTransforms.push_back(deltaTransform);

  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefSimple(
      bndryPnts, bndryPntTransforms);
  delete alphaTransform;
  delete bravoTransform;
  delete charlieTransform;
  delete deltaTransform;

  // no new vertices
  std::list<const AF2RuleNewVertex*> newVertices;

  // new edges
  std::list<const AF2RuleNewEdge*> newEdges;
  AF2RuleNewEdge* newEdgePtr = new AF2RuleNewEdge(0, 3);
  newEdges.push_back(newEdgePtr);
  newEdgePtr = new AF2RuleNewEdge(3, 2);
  newEdges.push_back(newEdgePtr);
  newEdgePtr = new AF2RuleNewEdge(2, 1);
  newEdges.push_back(newEdgePtr);

  // new faces
  std::list<const AF2RuleNewFace*> newFaces;
  AF2RuleNewFace* newFacePtr = new AF2RuleNewTriangle(0, 1, 2);
  newFaces.push_back(newFacePtr);
  newFacePtr = new AF2RuleNewTriangle(0, 2, 3);
  newFaces.push_back(newFacePtr);

  AF2Rule* rulePtr = new AF2Rule("triRule MultiPoint", 1u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

void testFreeTriangleRule()
{
  AF2Rule* freeTriRule = makeFreeTriangleRule();

  CHECK(freeTriRule != NULL);

  MeshKit::Vector<3> origin;
  origin[0] = -0.75;
  origin[1] = -0.5;
  origin[2] = 0.0;
  MeshKit::Vector<3> normal;
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 1.0;
  MeshKit::Vector<3> xDir;
  xDir[0] = 0.25;
  xDir[1] = 0.0;
  xDir[2] = 0.0;
  AF2PlaneProjection* xyPlaneProj =
      new AF2PlaneProjection(square->igeom_instance(),
      square->geom_handle(), origin, normal, xDir, 0.21875);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -1.0, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.75, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.5, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.25, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(4, 0.0, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -1.0, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -1.0, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -1.0, 0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(8, 0.0, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(9, 0.0, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(10, 0.0, 0.25, 0.5));

  std::list<AF2Point3D*> points;
  AF2Edge3D* baselineEdge;
  std::list<const AF2Edge3D*> edges;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  edges.push_back(new AF2Edge3D(pointsVector[3], pointsVector[4]));
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[0]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[5]));
  edges.push_back(new AF2Edge3D(pointsVector[7], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[8]));
  edges.push_back(new AF2Edge3D(pointsVector[8], pointsVector[9]));
  edges.push_back(new AF2Edge3D(pointsVector[9], pointsVector[10]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, xyPlaneProj);
  TestVisitor visitor;
  freeTriRule->applyRule(ngbhd, 1u, visitor);
  // The exact answer should be 8/14 because the position
  // of the match for x_{1} should be 8/7
  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the free triangle rule."
      << std::endl;
  CHECK_REAL_EQUAL(8.0/14.0, visitor.getFirstX(), 1e-14);
  std::cout << "Good: The new point is located at about 8.0/14.0 in the\n"
      << "  scaled two-dimensional coordinate system." << std::endl;

  AF2Rule* freeTriRuleMxQ5 = makeFreeTriangleRuleMaxQual5();
  TestVisitor mxQ5Visitor;
  freeTriRuleMxQ5->applyRule(ngbhd, 1u, mxQ5Visitor);
  CHECK_EQUAL(0u, mxQ5Visitor.getVisitCount());
  std::cout << "Good: A rule is not applied if the match quality is too high."
      << std::endl;
  freeTriRuleMxQ5->applyRule(ngbhd, 5u, mxQ5Visitor);
  CHECK_EQUAL(1u, mxQ5Visitor.getVisitCount());
  std::cout << "Good: The low-quality free triangle rule can be applied when\n"
      << "  the match quality is low enough." << std::endl;
  CHECK_REAL_EQUAL(8.0/14.0, mxQ5Visitor.getFirstX(), 1e-14);
  std::cout << "Good: The new point is located at about 8.0/14.0 in the\n"
      << "  scaled two-dimensional coordinate system." << std::endl;

  for (std::list<const AF2Edge3D*>::iterator itr = edges.begin();
      itr != edges.end(); ++itr)
  {
    delete *itr;
  }
  delete baselineEdge;
  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    delete pointsVector[pi];
  }
  // projection deleted by neighborhood
  delete freeTriRule;
  delete freeTriRuleMxQ5;
}

void testFillTriangleRule()
{
  AF2Rule* fillTriRule = makeFillTriangleRule();

  CHECK(fillTriRule != NULL);

  MeshKit::Vector<3> origin;
  origin[0] = -0.75;
  origin[1] = -0.25;
  origin[2] = 0.5;
  MeshKit::Vector<3> normal;
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 1.0;
  MeshKit::Vector<3> xDir;
  xDir[0] = 1.0;
  xDir[1] = 0.0;
  xDir[2] = 0.0;
  AF2PlaneProjection* xyPlaneProj =
      new AF2PlaneProjection(square->igeom_instance(),
      square->geom_handle(), origin, normal, xDir, 0.25);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.75, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.625, -0.05, 0.5));

  std::list<AF2Point3D*> points;
  AF2Edge3D* baselineEdge;
  std::list<const AF2Edge3D*> edges;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  baselineEdge = new AF2Edge3D(pointsVector[0], pointsVector[1]);

  edges.push_back(new AF2Edge3D(pointsVector[1], pointsVector[2]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[0]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, xyPlaneProj);
  TestVisitor visitor;
  fillTriRule->applyRule(ngbhd, 1u, visitor);

  for (std::list<const AF2Edge3D*>::iterator itr = edges.begin();
      itr != edges.end(); ++itr)
  {
    delete *itr;
  }
  delete baselineEdge;
  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    delete pointsVector[pi];
  }
  // projection deleted by neighborhood
  delete fillTriRule;
}

void testMultiPointTriangleRule()
{
  AF2Rule* multiPointTriRule = makeMultiPointTriangleRule();

  CHECK(multiPointTriRule != NULL);

  MeshKit::Vector<3> origin;
  origin[0] = -0.75;
  origin[1] = -0.25;
  origin[2] = 0.0;
  MeshKit::Vector<3> normal;
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 1.0;
  MeshKit::Vector<3> xDir;
  xDir[0] = 0.25;
  xDir[1] = 0.0;
  xDir[2] = 0.0;
  AF2PlaneProjection* xyPlaneProj =
      new AF2PlaneProjection(square->igeom_instance(),
      square->geom_handle(), origin, normal, xDir, 0.25);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.75, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.5625, -0.05, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.8125, -0.05, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.875, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.375, -0.25, 0.5));

  std::list<AF2Point3D*> points;
  AF2Edge3D* baselineEdge;
  std::list<const AF2Edge3D*> edges;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  baselineEdge = new AF2Edge3D(pointsVector[0], pointsVector[1]);

  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[0]));
  edges.push_back(new AF2Edge3D(pointsVector[1], pointsVector[5]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, xyPlaneProj);
  TestVisitor visitor;
  multiPointTriRule->applyRule(ngbhd, 1u, visitor);

  for (std::list<const AF2Edge3D*>::iterator itr = edges.begin();
      itr != edges.end(); ++itr)
  {
    delete *itr;
  }
  delete baselineEdge;
  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    delete pointsVector[pi];
  }
  // projection deleted by neighborhood
  delete multiPointTriRule;
}
