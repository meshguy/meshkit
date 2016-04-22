/**
 * \file ruleAppObj2D.cpp \test
 *
 * Test the two-dimensional polygon and two-dimensional rule application
 * object.
 *
 * Rule application objects are normally created when a rule is
 * successfully applied, but the the object can be tested without
 * the context of a rule.  This test creates points and polygons
 * like they might be created for a successful rule application
 * and uses them to test the rule application object.
 */

// C++
#include <iostream>
#include <list>

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Polygon2D.hpp"
#include "meshkit/AF2RuleApplication.hpp"
#include "meshkit/Error.hpp"

// MeshKit test utilities
#include "TestUtil.hpp"

void makeTriangle(AF2Point2D const & v0, AF2Point2D const & v1,
    AF2Point2D const & v2, AF2Polygon2D* & trianglePtr);
AF2Polygon2D* makeQuad(AF2Point2D const & v0, AF2Point2D const & v1,
    AF2Point2D const & v2, AF2Point2D const & v3);
AF2RuleApplication* makeZeroNewPointsRuleApplication(AF2Point2D const & v0,
    AF2Point2D const & v1, AF2Point2D const & v2);
AF2RuleApplication* makeTwoNewPointsRuleApplication(AF2Point2D const & v0,
    AF2Point2D const & v1, double v2x, double v2y, double v3x, double v3y);
void testPolygonCreateDestroy();
void testPolygonCopy();
void testPolygonAssign();
void testPolygonConstructorException();
void testPolygonRangeException();
void testRuleAppZeroPntsCreateDestroy();
void testRuleAppZeroPntsCopy();
void testRuleAppZeroPntsAssign();
void testRuleAppTwoPntsCreateDestroy();
void testRuleAppTwoPntsCopy();
void testRuleAppTwoPntsAssign();
void testRuleAppConstructorExceptions();
void testRuleAppRangeExceptions();

int main(int argc, char **argv)
{

  int num_fail = 0;

  num_fail += RUN_TEST(testPolygonCreateDestroy);
  num_fail += RUN_TEST(testPolygonCopy);
  num_fail += RUN_TEST(testPolygonAssign);
  num_fail += RUN_TEST(testPolygonConstructorException);
  num_fail += RUN_TEST(testPolygonRangeException);
  num_fail += RUN_TEST(testRuleAppZeroPntsCreateDestroy);
  num_fail += RUN_TEST(testRuleAppZeroPntsCopy);
  num_fail += RUN_TEST(testRuleAppZeroPntsAssign);
  num_fail += RUN_TEST(testRuleAppTwoPntsCreateDestroy);
  num_fail += RUN_TEST(testRuleAppTwoPntsCopy);
  num_fail += RUN_TEST(testRuleAppTwoPntsAssign);
  num_fail += RUN_TEST(testRuleAppConstructorExceptions);
  num_fail += RUN_TEST(testRuleAppRangeExceptions);

  return num_fail;
}

void makeTriangle(AF2Point2D const & v0, AF2Point2D const & v1,
    AF2Point2D const & v2, AF2Polygon2D* & trianglePtr)
{
  std::list<const AF2Point2D*> triangleVertices;
  triangleVertices.push_back(&v0);
  triangleVertices.push_back(&v1);
  triangleVertices.push_back(&v2);
  trianglePtr = new AF2Polygon2D(triangleVertices);
}

AF2Polygon2D* makeQuad(AF2Point2D const & v0, AF2Point2D const & v1,
    AF2Point2D const & v2, AF2Point2D const & v3)
{
  std::list<const AF2Point2D*> quadVertices;
  quadVertices.push_back(&v0);
  quadVertices.push_back(&v1);
  quadVertices.push_back(&v2);
  quadVertices.push_back(&v3);
  return new AF2Polygon2D(quadVertices);
}

/**
 * Make a rule application that has one triangle and depends
 * on existing vertices whose memory is managed outside the
 * rule application.
 */
AF2RuleApplication* makeZeroNewPointsRuleApplication(AF2Point2D const & v0,
    AF2Point2D const & v1, AF2Point2D const & v2)
{
  // make the triangle from the existing vertices
  AF2Polygon2D* ruleAppTriangle = NULL;
  makeTriangle(v0, v1, v2, ruleAppTriangle);

  // make an empty list of points -- there are no new points in the rule app
  std::list<const AF2Point2D*> ruleAppNewPoints;

  // make a list of polygons containing the triangle
  std::list<const AF2Polygon2D*> ruleAppFaces;
  ruleAppFaces.push_back(ruleAppTriangle);

  // construct the rule application
  AF2RuleApplication* ruleApplication =
      new AF2RuleApplication(ruleAppNewPoints, ruleAppFaces);

  // delete the triangle, which should have been copied in the rule
  // application constructor
  delete ruleAppTriangle;

  // return the rule application
  return ruleApplication;
}

/**
 * Make a rule application that has one quadrilateral.  The quadrilateral
 * depends on two existing vertices whose memory is managed outside the
 * rule application and two new vertices whose memory is managed
 * by the rule application.
 */
AF2RuleApplication* makeTwoNewPointsRuleApplication(AF2Point2D const & v0,
    AF2Point2D const & v1, double v2x, double v2y, double v3x, double v3y)
{
  // make the two new points that will be stored in the rule application
  AF2Point2D* v2 = new AF2Point2D(v2x, v2y);
  AF2Point2D* v3 = new AF2Point2D(v3x, v3y);

  // make the quadrilateral from the two existing and two new vertices
  AF2Polygon2D* ruleAppQuad = makeQuad(v0, v1, *v2, *v3);

  // make a list of the new points
  // the order should not matter
  std::list<const AF2Point2D*> ruleAppNewPoints;
  ruleAppNewPoints.push_back(v3);
  ruleAppNewPoints.push_back(v2);

  // make a list of polygons containing the quadrilateral
  std::list<const AF2Polygon2D*> ruleAppFaces;
  ruleAppFaces.push_back(ruleAppQuad);

  // construct the rule application
  AF2RuleApplication* ruleApplication =
      new AF2RuleApplication(ruleAppNewPoints, ruleAppFaces);

  // delete the new points and the quadrilateral, which all should
  // have been copied to new memory owned by the rule application
  delete ruleAppQuad;
  delete v3;
  delete v2;

  // return the rule application
  return ruleApplication;
}

void testPolygonCreateDestroy()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2Polygon2D* trianglePtr;
  makeTriangle(v0, v1, v2, trianglePtr);

  CHECK_EQUAL(trianglePtr->getNumVertices(), 3u);
  std::cout << "Good: The triangle has 3 vertices." << std::endl;
  CHECK(trianglePtr->getVertex(0) == &v0);
  std::cout << "Good: The vertex at index 0 matches expectation." << std::endl;
  CHECK(trianglePtr->getVertex(1) == &v1);
  std::cout << "Good: The vertex at index 1 matches expectation." << std::endl;
  CHECK(trianglePtr->getVertex(2) == &v2);
  std::cout << "Good: The vertex at index 2 matches expectation." << std::endl;

  delete trianglePtr;

  std::cout << "PASS: The triangle can be constructed and deleted."
      << std::endl;
}

void testPolygonCopy()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2Polygon2D* trianglePtr;
  makeTriangle(v0, v1, v2, trianglePtr);

  AF2Polygon2D copiedTriangle(*trianglePtr);
  delete trianglePtr;
  CHECK(copiedTriangle.getNumVertices() == 3 &&
      copiedTriangle.getVertex(0) == &v0 &&
      copiedTriangle.getVertex(1) == &v1 &&
      copiedTriangle.getVertex(2) == &v2);

  std::cout << "PASS: The copied triangle matches expectations." << std::endl;
}

void testPolygonAssign()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2Polygon2D* trianglePtr;
  makeTriangle(v0, v1, v2, trianglePtr);

  AF2Point2D v3(5, 5);
  AF2Point2D v4(6, 5);
  AF2Point2D v5(5.5, 5.8);
  AF2Polygon2D* otherTrianglePtr;
  makeTriangle(v3, v4, v5, otherTrianglePtr);

  AF2Polygon2D assignedTriangle(*otherTrianglePtr);
  delete otherTrianglePtr;
  assignedTriangle = *trianglePtr;
  delete trianglePtr;
  CHECK(assignedTriangle.getNumVertices() == 3 &&
      assignedTriangle.getVertex(0) == &v0 &&
      assignedTriangle.getVertex(1) == &v1 &&
      assignedTriangle.getVertex(2) == &v2);
  std::cout << "PASS: The triangle produced by assignment matches expectations."
      << std::endl;
}

void testPolygonConstructorException()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);

  std::list<const AF2Point2D*> badVertexList;
  badVertexList.push_back(&v0);
  badVertexList.push_back(NULL);
  badVertexList.push_back(&v1);
  bool exceptionThrown = false;
  try
  {
    AF2Polygon2D* shouldFail = new AF2Polygon2D(badVertexList);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(mkError.error_code(), MeshKit::MK_BAD_INPUT);
    std::cout << "Good: The error is a bad input error." << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "PASS: A null pointer in the vertex list passed to the\n"
      << "  constructor caused an exception." << std::endl;
}

void testPolygonRangeException()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2Polygon2D* trianglePtr;
  makeTriangle(v0, v1, v2, trianglePtr);

  bool exceptionThrown = false;
  try
  {
    trianglePtr->getVertex(3);
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(mkError.error_code(), MeshKit::MK_BAD_INPUT);
    std::cout << "Good: The error is a bad input error." << std::endl;
    exceptionThrown = true;
  }

  delete trianglePtr;
  CHECK(exceptionThrown);
  std::cout << "PASS: An index out of range passed to the getVertex() method\n"
      << "  caused an exception." << std::endl;
}

void testRuleAppZeroPntsCreateDestroy()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2RuleApplication* ruleApplication =
      makeZeroNewPointsRuleApplication(v0, v1, v2);

  CHECK_EQUAL(ruleApplication->getNumNewPoints(), 0u);
  std::cout << "Good: The rule application has 0 new points." << std::endl;
  CHECK_EQUAL(ruleApplication->getNumNewFaces(), 1u);
  std::cout << "Good: The rule application has 1 new face." << std::endl;
  const AF2Polygon2D* newFace = ruleApplication->getNewFace(0u);
  CHECK(newFace->getVertex(0u) == &v0);
  std::cout
      << "Good: The vertex at index 0 of the new face matches expectation."
      << std::endl;
  CHECK(newFace->getVertex(1u) == &v1);
  std::cout
      << "Good: The vertex at index 1 of the new face matches expectation."
      << std::endl;
  CHECK(newFace->getVertex(2u) == &v2);
  std::cout
      << "Good: The vertex at index 2 of the new face matches expectation."
      << std::endl;

  delete ruleApplication;

  std::cout << "PASS: The rule application with no new points\n"
      << "  can be constructed and deleted." << std::endl;
}

void testRuleAppZeroPntsCopy()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2RuleApplication* ruleApplication =
      makeZeroNewPointsRuleApplication(v0, v1, v2);

  AF2RuleApplication copiedRuleApplication(*ruleApplication);
  delete ruleApplication;
  const AF2Polygon2D* newFace = copiedRuleApplication.getNewFace(0u);
  CHECK(copiedRuleApplication.getNumNewPoints() == 0u &&
      copiedRuleApplication.getNumNewFaces() == 1u &&
      newFace->getVertex(0u) == &v0 &&
      newFace->getVertex(1u) == &v1 &&
      newFace->getVertex(2u) == &v2);

  std::cout << "PASS: The copied rule application matches expectations."
      << std::endl;
}

void testRuleAppZeroPntsAssign()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2RuleApplication* ruleApplication =
      makeZeroNewPointsRuleApplication(v0, v1, v2);

  AF2Point2D v3(5, 5);
  AF2Point2D v4(6, 5);
  AF2Point2D v5(5.5, 5.8);
  AF2RuleApplication* otherRuleApplication =
      makeZeroNewPointsRuleApplication(v3, v4, v5);

  AF2RuleApplication assignedRuleApplication(*otherRuleApplication);
  delete otherRuleApplication;
  assignedRuleApplication = *ruleApplication;
  delete ruleApplication;

  const AF2Polygon2D* newFace = assignedRuleApplication.getNewFace(0u);
  CHECK(assignedRuleApplication.getNumNewPoints() == 0u &&
      assignedRuleApplication.getNumNewFaces() == 1u &&
      newFace->getVertex(0u) == &v0 &&
      newFace->getVertex(1u) == &v1 &&
      newFace->getVertex(2u) == &v2);

  std::cout << "PASS: The rule application produced by assignment matches\n"
      << "  expectations." << std::endl;
}

void testRuleAppTwoPntsCreateDestroy()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2RuleApplication* ruleApplication =
      makeTwoNewPointsRuleApplication(v0, v1, 1.0625, 1.0, 0.125, 0.875);

  CHECK_EQUAL(ruleApplication->getNumNewPoints(), 2u);
  std::cout << "Good: The rule application has 2 new points." << std::endl;
  const AF2Point2D* v2 = ruleApplication->getNewPoint(1u);
  const AF2Point2D* v3 = ruleApplication->getNewPoint(0u);
  CHECK_REAL_EQUAL(1.0625, v2->getX(), 0.0);
  CHECK_REAL_EQUAL(1.0, v2->getY(), 0.0);
  CHECK_REAL_EQUAL(0.125, v3->getX(), 0.0);
  CHECK_REAL_EQUAL(0.875, v3->getY(), 0.0);
  std::cout << "Good: The new points have the expected coordinates."
      << std::endl;
  CHECK_EQUAL(ruleApplication->getNumNewFaces(), 1u);
  std::cout << "Good: The rule application has 1 new face." << std::endl;
  const AF2Polygon2D* newFace = ruleApplication->getNewFace(0u);
  CHECK(newFace->getVertex(0u) == &v0);
  std::cout
      << "Good: The vertex at index 0 of the new face matches expectation."
      << std::endl;
  CHECK(newFace->getVertex(1u) == &v1);
  std::cout
      << "Good: The vertex at index 1 of the new face matches expectation."
      << std::endl;
  CHECK(newFace->getVertex(2u) == v2);
  std::cout
      << "Good: The vertex at index 2 of the new face matches expectation."
      << std::endl;
  CHECK(newFace->getVertex(3u) == v3);
  std::cout
      << "Good: The vertex at index 3 of the new face matches expectation."
      << std::endl;

  delete ruleApplication;

  std::cout << "PASS: The rule application with two new points\n"
      << "  can be constructed and deleted." << std::endl;
}

void testRuleAppTwoPntsCopy()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2RuleApplication* ruleApplication =
      makeTwoNewPointsRuleApplication(v0, v1, 1.0625, 1.0, 0.125, 0.875);

  AF2RuleApplication copiedRuleApplication(*ruleApplication);
  CHECK(copiedRuleApplication.getNumNewPoints() == 2u &&
      copiedRuleApplication.getNumNewFaces() == 1u);
  std::cout << "Good: The copied rule application has the expected\n"
      << "  number of new points and new faces." << std::endl;
  const AF2Point2D* originalV2 = ruleApplication->getNewPoint(1u);
  const AF2Point2D* v2 = copiedRuleApplication.getNewPoint(1u);
  CHECK(originalV2 != v2);
  std::cout << "Good: The first new point in the copied rule is a different\n"
      << "  object than the first new point in the original rule."
      << std::endl;
  delete ruleApplication;
  v2 = copiedRuleApplication.getNewPoint(1u);
  const AF2Point2D* v3 = copiedRuleApplication.getNewPoint(0u);
  CHECK_REAL_EQUAL(1.0625, v2->getX(), 0.0);
  CHECK_REAL_EQUAL(1.0, v2->getY(), 0.0);
  CHECK_REAL_EQUAL(0.125, v3->getX(), 0.0);
  CHECK_REAL_EQUAL(0.875, v3->getY(), 0.0);
  const AF2Polygon2D* newFace = copiedRuleApplication.getNewFace(0u);
  CHECK(copiedRuleApplication.getNumNewPoints() == 2u &&
      copiedRuleApplication.getNumNewFaces() == 1u &&
      newFace->getVertex(0u) == &v0 &&
      newFace->getVertex(1u) == &v1 &&
      newFace->getVertex(2u) == v2 &&
      newFace->getVertex(3u) == v3);

  std::cout << "PASS: The copied rule application matches expectations."
      << std::endl;
}

void testRuleAppTwoPntsAssign()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2RuleApplication* ruleApplication =
      makeTwoNewPointsRuleApplication(v0, v1, 1.0625, 1.0, 0.125, 0.875);
  const AF2Point2D* originalV2 = ruleApplication->getNewPoint(1u);

  AF2Point2D v4(5, 5);
  AF2Point2D v5(6, 5);
  AF2RuleApplication* otherRuleApplication =
      makeTwoNewPointsRuleApplication(v4, v5, 6.0, 5.0625, 5.875, 5.125);

  AF2RuleApplication assignedRuleApplication(*otherRuleApplication);
  delete otherRuleApplication;
  assignedRuleApplication = *ruleApplication;
  // this should makes originalV2 a dangling pointer, but the code can
  // still check (later) that it is pointing to a different location than v2
  delete ruleApplication;

  const AF2Point2D* v2 = assignedRuleApplication.getNewPoint(1u);
  const AF2Point2D* v3 = assignedRuleApplication.getNewPoint(0u);
  CHECK_REAL_EQUAL(1.0625, v2->getX(), 0.0);
  CHECK_REAL_EQUAL(1.0, v2->getY(), 0.0);
  CHECK_REAL_EQUAL(0.125, v3->getX(), 0.0);
  CHECK_REAL_EQUAL(0.875, v3->getY(), 0.0);
  const AF2Polygon2D* newFace = assignedRuleApplication.getNewFace(0u);
  CHECK(assignedRuleApplication.getNumNewPoints() == 2u &&
      assignedRuleApplication.getNumNewFaces() == 1u &&
      originalV2 != v2 &&
      newFace->getVertex(0u) == &v0 &&
      newFace->getVertex(1u) == &v1 &&
      newFace->getVertex(2u) == v2 &&
      newFace->getVertex(3u) == v3);

  std::cout << "PASS: The rule application produced by assignment matches\n"
      << "  expectations." << std::endl;
}

void testRuleAppConstructorExceptions()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2Polygon2D* trianglePtr = NULL;
  makeTriangle(v0, v1, v2, trianglePtr);

  std::list<const AF2Point2D*> badNewPointsList;
  badNewPointsList.push_back(NULL);
  badNewPointsList.push_back(&v2);

  std::list<const AF2Polygon2D*> goodNewFacesList;
  goodNewFacesList.push_back(trianglePtr);

  bool exceptionThrown = false;
  try
  {
    AF2RuleApplication* shouldFail =
        new AF2RuleApplication(badNewPointsList, goodNewFacesList);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(mkError.error_code(), MeshKit::MK_BAD_INPUT);
    std::cout << "Good: The error is a bad input error." << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "Good: A null pointer in the list of new points passed to the\n"
      << "  constructor caused an exception." << std::endl;

  std::list<const AF2Point2D*> goodNewPointsList;
  goodNewPointsList.push_back(&v2);

  std::list<const AF2Polygon2D*> badNewFacesList;
  badNewFacesList.push_back(trianglePtr);
  badNewFacesList.push_back(NULL);

  exceptionThrown = false;
  try
  {
    AF2RuleApplication* shouldFail =
        new AF2RuleApplication(goodNewPointsList, badNewFacesList);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(mkError.error_code(), MeshKit::MK_BAD_INPUT);
    std::cout << "Good: The error is a bad input error." << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "Good: A null pointer in the list of new faces passed to the\n"
      << "  constructor caused an exception." << std::endl;
  std::cout << "PASS: Exceptions were thrown in both cases." << std::endl;
}

void testRuleAppRangeExceptions()
{
  AF2Point2D v0(0, 0);
  AF2Point2D v1(1, 0);
  AF2Point2D v2(0.5, 0.8);
  AF2RuleApplication* ruleApplication =
      makeZeroNewPointsRuleApplication(v0, v1, v2);

  bool exceptionThrown = false;
  try
  {
    ruleApplication->getNewPoint(0);
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(mkError.error_code(), MeshKit::MK_BAD_INPUT);
    std::cout << "Good: The error is a bad input error." << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "Good: An index out of range passed to the getNewPoint()\n"
      << "  method caused an exception." << std::endl;

  exceptionThrown = false;
  try
  {
    ruleApplication->getNewFace(1);
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(mkError.error_code(), MeshKit::MK_BAD_INPUT);
    std::cout << "Good: The error is a bad input error." << std::endl;
    exceptionThrown = true;
  }

  delete ruleApplication;
  CHECK(exceptionThrown);
  std::cout << "Good: An index out of range passed to the getNewFace()\n"
      << "  method caused an exception." << std::endl;
  std::cout << "PASS: Exceptions were thrown in both cases." << std::endl;
}
