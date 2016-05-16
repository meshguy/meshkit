/**
 * \file basicAlg2D.cpp \test
 *
 * \brief Test the AF2Algorithm with two basic rules
 *
 */

// C++
#include <cmath>
#include <iostream>
#include <vector>

// MeshKit
#include "meshkit/AF2Algorithm.hpp"
#include "meshkit/AF2DfltPlaneProjMaker.hpp"
#include "meshkit/AF2FreeZoneDefSimple.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2PointTransformNone.hpp"
#include "meshkit/AF2Rule.hpp"
#include "meshkit/AF2RuleNewTriangle.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"

// MeshKit testing
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model
#if HAVE_OCC
#define FILE_EXT "stp"
#else
#define FILE_EXT "sat"
#endif

// This variable is at global scope because (1) calling deleteAll on
// the MKCore geometry instance appears to cause memory inconsistencies
// with later use of the geometry instance and (2) it is more efficient
// to load the geometry model only once
MeshKit::ModelEnt* square = NULL;

AF2Rule* makeAddPeakVertexRule();
AF2Rule* makeConnectToVertexRule();
AF2Rule* makeCloseTriangleRule();
AF2Rule* makeFillTriangleRule();
void testAlgorithmFail();
void testAlgorithmSucceed();
void testAlgorithmSucceedIsoPoint();
void testAlgorithmSucceedAddPoint();

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

  int num_fail = 0;

  num_fail += RUN_TEST(testAlgorithmFail);
  num_fail += RUN_TEST(testAlgorithmSucceed);
  num_fail += RUN_TEST(testAlgorithmSucceedIsoPoint);
  num_fail += RUN_TEST(testAlgorithmSucceedAddPoint);

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
      new AF2RuleExistVertex(0.5, sqrt(3.0)/2.0);
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
  AF2Point2D charlie(0.5, sqrt(3.0)/2.0);
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

  AF2Rule* rulePtr = new AF2Rule("Fill Triangle", 1u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

AF2Rule* makeCloseTriangleRule()
{
  // existing vertices
  std::list<const AF2RuleExistVertex*> exVertices;
  AF2RuleExistVertex* originVertexPtr = new AF2RuleExistVertex(0, 0);
  exVertices.push_back(originVertexPtr);
  AF2RuleExistVertex* baseVertexPtr = new AF2RuleExistVertex(1, 0);
  exVertices.push_back(baseVertexPtr);
  AF2RuleExistVertex* peakVertexPtr =
      new AF2RuleExistVertex(0.5, sqrt(3.0)/2.0);
  exVertices.push_back(peakVertexPtr);

  // existing edges
  AF2RuleExistEdge* baseEdgePtr =
      new AF2RuleExistEdge(originVertexPtr, baseVertexPtr);
  std::list<const AF2RuleExistEdge*> exEdges;
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
  AF2Point2D charlie(0.5, sqrt(3.0)/2.0);
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

  // new edge
  std::list<const AF2RuleNewEdge*> newEdges;
  AF2RuleNewEdge* newEdgePtr = new AF2RuleNewEdge(2, 1);
  newEdges.push_back(newEdgePtr);

  // new face
  std::list<const AF2RuleNewFace*> newFaces;
  AF2RuleNewFace* newFacePtr = new AF2RuleNewTriangle(0, 1, 2);
  newFaces.push_back(newFacePtr);

  AF2Rule* rulePtr = new AF2Rule("Close Triangle", 1u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

AF2Rule* makeConnectToVertexRule()
{
  // existing vertices
  std::list<const AF2RuleExistVertex*> exVertices;
  AF2RuleExistVertex* originVertexPtr = new AF2RuleExistVertex(0, 0);
  exVertices.push_back(originVertexPtr);
  AF2RuleExistVertex* baseVertexPtr = new AF2RuleExistVertex(1, 0);
  exVertices.push_back(baseVertexPtr);
  AF2RuleExistVertex* peakVertexPtr =
      new AF2RuleExistVertex(0.5, sqrt(3.0)/2.0);
  exVertices.push_back(peakVertexPtr);

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
  AF2Point2D bravo(1, 0);
  AF2PointTransform* bravoTransform = new AF2PointTransformNone();
  bndryPnts.push_back(bravo);
  bndryPntTransforms.push_back(bravoTransform);

  // third element free zone definition lists
  AF2Point2D charlie(1.4, 0.75);
  AF2PointTransform* charlieTransform = new AF2PointTransformNone();
  bndryPnts.push_back(charlie);
  bndryPntTransforms.push_back(charlieTransform);

  // fourth element free zone definition lists
  AF2Point2D delta(0.5, sqrt(3.0)/2.0);
  AF2PointTransform* deltaTransform = new AF2PointTransformNone();
  bndryPnts.push_back(delta);
  bndryPntTransforms.push_back(deltaTransform);

  // fifth element free zone definition lists
  AF2Point2D echo(-0.4, 0.75);
  AF2PointTransform* echoTransform = new AF2PointTransformNone();
  bndryPnts.push_back(echo);
  bndryPntTransforms.push_back(echoTransform);

  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefSimple(
      bndryPnts, bndryPntTransforms);
  delete alphaTransform;
  delete bravoTransform;
  delete charlieTransform;
  delete deltaTransform;
  delete echoTransform;

  // no new vertices
  std::list<const AF2RuleNewVertex*> newVertices;

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

  AF2Rule* rulePtr = new AF2Rule("Connect To Vertex", 1u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

AF2Rule* makeAddPeakVertexRule()
{
  // existing vertices
  std::list<const AF2RuleExistVertex*> exVertices;
  AF2RuleExistVertex* originVertexPtr = new AF2RuleExistVertex(0, 0);
  exVertices.push_back(originVertexPtr);
  AF2RuleExistVertex* baseVertexPtr = new AF2RuleExistVertex(1, 0);
  exVertices.push_back(baseVertexPtr);

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
  AF2Point2D bravo(1, 0);
  AF2PointTransform* bravoTransform = new AF2PointTransformNone();
  bndryPnts.push_back(bravo);
  bndryPntTransforms.push_back(bravoTransform);

  // third element free zone definition lists
  AF2Point2D charlie(1.4, 0.75);
  AF2PointTransform* charlieTransform = new AF2PointTransformNone();
  bndryPnts.push_back(charlie);
  bndryPntTransforms.push_back(charlieTransform);

  // fourth element free zone definition lists
  AF2Point2D delta(0.5, 1.25);
  AF2PointTransform* deltaTransform = new AF2PointTransformNone();
  bndryPnts.push_back(delta);
  bndryPntTransforms.push_back(deltaTransform);

  // fifth element free zone definition lists
  AF2Point2D echo(-0.4, 0.75);
  AF2PointTransform* echoTransform = new AF2PointTransformNone();
  bndryPnts.push_back(echo);
  bndryPntTransforms.push_back(echoTransform);

  AF2FreeZoneDef* freeZoneDef = new AF2FreeZoneDefSimple(
      bndryPnts, bndryPntTransforms);
  delete alphaTransform;
  delete bravoTransform;
  delete charlieTransform;
  delete deltaTransform;
  delete echoTransform;

  // new vertex
  std::list<const AF2RuleNewVertex*> newVertices;
  AF2Point2D newVertexLoc(0.5, sqrt(3.0)/2.0);
  AF2PointTransform* newVertexTransform = new AF2PointTransformNone();
  AF2RuleNewVertex* newVertex =
      new AF2RuleNewVertex(newVertexLoc, newVertexTransform);
  delete newVertexTransform;
  newVertices.push_back(newVertex);

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

  AF2Rule* rulePtr = new AF2Rule("Add Peak Vertex", 1u, exVertices,
      baseEdgePtr, exEdges, freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

void testAlgorithmFail()
{
  AF2Rule* closeTriRule = makeCloseTriangleRule();
  AF2Rule* fillTriRule = makeFillTriangleRule();
  AF2LocalTransformMaker* transformMaker = new AF2DfltPlaneProjMaker(
      square->igeom_instance(), square->geom_handle());
  std::list<const AF2Rule*> ruleList;
  ruleList.push_back(closeTriRule);
  ruleList.push_back(fillTriRule);

  AF2Algorithm alg(ruleList);

  // coordinates of a regular hexagon with side length 0.1
  std::vector<double> coordinates;
  coordinates.push_back(-0.5);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.35);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.5);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.55);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);

  // edges of the regular hexagon
  std::vector<unsigned int> edges;
  edges.push_back(0u);
  edges.push_back(1u);
  edges.push_back(1u);
  edges.push_back(2u);
  edges.push_back(2u);
  edges.push_back(3u);
  edges.push_back(3u);
  edges.push_back(4u);
  edges.push_back(4u);
  edges.push_back(5u);
  edges.push_back(5u);
  edges.push_back(0u);

  unsigned int numPoints = coordinates.size() / 3;
  unsigned int numEdges = edges.size() / 2;
  AF2AlgorithmResult* result = alg.execute(
      transformMaker, &(coordinates[0]), numPoints, &(edges[0]), numEdges);

  CHECK(!result->isSuccessful());
  std::cout << "Good: Algorithm result is not successful." << std::endl;

  ruleList.clear();
  delete result;
  delete transformMaker;
  delete fillTriRule;
  delete closeTriRule;

  std::cout << "PASS: The algorithm appears to run properly on a hexagon\n"
      << "  with no rules that allow adding points and no free zones that\n"
      << "  have point transforms."  << std::endl;
}

void testAlgorithmSucceed()
{
  AF2Rule* closeTriRule = makeCloseTriangleRule();
  AF2Rule* fillTriRule = makeFillTriangleRule();
  AF2LocalTransformMaker* transformMaker = new AF2DfltPlaneProjMaker(
      square->igeom_instance(), square->geom_handle());
  std::list<const AF2Rule*> ruleList;
  ruleList.push_back(closeTriRule);
  ruleList.push_back(fillTriRule);

  AF2Algorithm alg(ruleList);

  // coordinates of a regular hexagon with side length 0.1
  // plus the coordinates of a central vertex
  std::vector<double> coordinates;
  coordinates.push_back(-0.5);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.35);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.5);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.55);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.45);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);

  // edges of the regular hexagon, plus a hanging edge to the central vertex
  std::vector<unsigned int> edges;
  edges.push_back(0u);
  edges.push_back(6u);
  edges.push_back(6u);
  edges.push_back(0u);
  edges.push_back(0u);
  edges.push_back(1u);
  edges.push_back(1u);
  edges.push_back(2u);
  edges.push_back(2u);
  edges.push_back(3u);
  edges.push_back(3u);
  edges.push_back(4u);
  edges.push_back(4u);
  edges.push_back(5u);
  edges.push_back(5u);
  edges.push_back(0u);

  unsigned int numPoints = coordinates.size() / 3;
  unsigned int numEdges = edges.size() / 2;
  AF2AlgorithmResult* result = alg.execute(
      transformMaker, &(coordinates[0]), numPoints, &(edges[0]), numEdges);

  CHECK(result->isSuccessful());
  std::cout << "Good: Algorithm result is successful." << std::endl;
  std::list<AF2Point3D*>::size_type sizeSeven(7u);
  CHECK_EQUAL(sizeSeven, result->getPoints()->size());
  std::cout << "Good: Algorithm result has seven points." << std::endl;
  std::list<const AF2Polygon3D*>::size_type sizeSix(6u);
  CHECK_EQUAL(sizeSix, result->getFaces()->size());
  std::cout << "Good: Algorithm result has six faces." << std::endl;

  ruleList.clear();
  delete result;
  delete transformMaker;
  delete fillTriRule;
  delete closeTriRule;

  std::cout << "PASS: The algorithm appears to run properly on a hexagon\n"
      << "  with a hanging edge." << std::endl;
}

void testAlgorithmSucceedIsoPoint()
{
  AF2Rule* closeTriRule = makeCloseTriangleRule();
  AF2Rule* fillTriRule = makeFillTriangleRule();
  AF2Rule* connectToVertexRule = makeConnectToVertexRule();
  AF2LocalTransformMaker* transformMaker = new AF2DfltPlaneProjMaker(
      square->igeom_instance(), square->geom_handle());
  std::list<const AF2Rule*> ruleList;
  ruleList.push_back(closeTriRule);
  ruleList.push_back(fillTriRule);
  ruleList.push_back(connectToVertexRule);

  AF2Algorithm alg(ruleList);

  // coordinates of a regular hexagon with side length 0.1
  // plus the coordinates of a central vertex
  std::vector<double> coordinates;
  coordinates.push_back(-0.5);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.35);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.5);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.55);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.45);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);

  // edges of the regular hexagon
  std::vector<unsigned int> edges;
  edges.push_back(0u);
  edges.push_back(1u);
  edges.push_back(1u);
  edges.push_back(2u);
  edges.push_back(2u);
  edges.push_back(3u);
  edges.push_back(3u);
  edges.push_back(4u);
  edges.push_back(4u);
  edges.push_back(5u);
  edges.push_back(5u);
  edges.push_back(0u);

  unsigned int numPoints = coordinates.size() / 3;
  unsigned int numEdges = edges.size() / 2;
  AF2AlgorithmResult* result = alg.execute(
      transformMaker, &(coordinates[0]), numPoints, &(edges[0]), numEdges);

  CHECK(result->isSuccessful());
  std::cout << "Good: Algorithm result is successful." << std::endl;
  std::list<AF2Point3D*>::size_type sizeSeven(7u);
  CHECK_EQUAL(sizeSeven, result->getPoints()->size());
  std::cout << "Good: Algorithm result has seven points." << std::endl;
  std::list<const AF2Polygon3D*>::size_type sizeSix(6u);
  CHECK_EQUAL(sizeSix, result->getFaces()->size());
  std::cout << "Good: Algorithm result has six faces." << std::endl;

  ruleList.clear();
  delete result;
  delete transformMaker;
  delete connectToVertexRule;
  delete fillTriRule;
  delete closeTriRule;

  std::cout << "PASS: The algorithm appears to run properly on a hexagon\n"
      << "  with an isolated central vertex." << std::endl;
}

void testAlgorithmSucceedAddPoint()
{
  AF2Rule* closeTriRule = makeCloseTriangleRule();
  AF2Rule* fillTriRule = makeFillTriangleRule();
  AF2Rule* connectToVertexRule = makeConnectToVertexRule();
  AF2Rule* addPeakVertexRule = makeAddPeakVertexRule();
  AF2LocalTransformMaker* transformMaker = new AF2DfltPlaneProjMaker(
      square->igeom_instance(), square->geom_handle());
  std::list<const AF2Rule*> ruleList;
  ruleList.push_back(closeTriRule);
  ruleList.push_back(fillTriRule);
  ruleList.push_back(connectToVertexRule);
  ruleList.push_back(addPeakVertexRule);

  AF2Algorithm alg(ruleList);

  // coordinates of a regular hexagon with side length 0.1
  std::vector<double> coordinates;
  coordinates.push_back(-0.5);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(0.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.35);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.4);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.5);
  coordinates.push_back(sqrt(3)/10.0);
  coordinates.push_back(-0.5);
  coordinates.push_back(-0.55);
  coordinates.push_back(sqrt(3)/20.0);
  coordinates.push_back(-0.5);

  // edges of the regular hexagon
  std::vector<unsigned int> edges;
  edges.push_back(0u);
  edges.push_back(1u);
  edges.push_back(1u);
  edges.push_back(2u);
  edges.push_back(2u);
  edges.push_back(3u);
  edges.push_back(3u);
  edges.push_back(4u);
  edges.push_back(4u);
  edges.push_back(5u);
  edges.push_back(5u);
  edges.push_back(0u);

  unsigned int numPoints = coordinates.size() / 3;
  unsigned int numEdges = edges.size() / 2;
  AF2AlgorithmResult* result = alg.execute(
      transformMaker, &(coordinates[0]), numPoints, &(edges[0]), numEdges);

  CHECK(result->isSuccessful());
  std::cout << "Good: Algorithm result is successful." << std::endl;
  std::list<AF2Point3D*>::size_type sizeSeven(7u);
  CHECK_EQUAL(sizeSeven, result->getPoints()->size());
  std::cout << "Good: Algorithm result has seven points." << std::endl;
  std::list<const AF2Polygon3D*>::size_type sizeSix(6u);
  CHECK_EQUAL(sizeSix, result->getFaces()->size());
  std::cout << "Good: Algorithm result has six faces." << std::endl;

  ruleList.clear();
  delete result;
  delete transformMaker;
  delete addPeakVertexRule;
  delete connectToVertexRule;
  delete fillTriRule;
  delete closeTriRule;

  std::cout << "PASS: The algorithm appears to run properly on a hexagon\n"
      << "  when it can add a central vertex." << std::endl;
}
