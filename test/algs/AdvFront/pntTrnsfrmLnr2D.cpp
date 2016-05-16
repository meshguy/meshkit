/**
 * \file pntTrnsfrmLnr2D.cpp \test
 *
 * Test two-dimensional linear vertex-based point transform.
 */

// C++
#include <iostream>
#include <list>
#include <vector>

// MeshKit
#include "meshkit/AF2Binding.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2PntTrnsfrmLnrV.hpp"
#include "meshkit/AF2RuleExistVertex.hpp"
#include "meshkit/Error.hpp"

// MeshKit testing
#include "TestUtil.hpp"

void testListConstructErrorX();
void testListConstructErrorY();
void testVectorConstructErrorX();
void testVectorConstructErrorY();
void testUnboundError();
void testSingleVertexJustX();
void testSingleVertexJustY();
void testSingleVertexFull();
void testMultipleVertices();

int main(int argc, char **argv)
{
  int num_fail = 0;

  num_fail += RUN_TEST(testListConstructErrorX);
  num_fail += RUN_TEST(testListConstructErrorY);
  num_fail += RUN_TEST(testVectorConstructErrorX);
  num_fail += RUN_TEST(testVectorConstructErrorY);
  num_fail += RUN_TEST(testUnboundError);
  num_fail += RUN_TEST(testSingleVertexJustX);
  num_fail += RUN_TEST(testSingleVertexJustY);
  num_fail += RUN_TEST(testSingleVertexFull);
  num_fail += RUN_TEST(testMultipleVertices);

  return num_fail;
}

void testListConstructErrorX()
{
  std::list<const AF2RuleExistVertex*> emptyVertexList;
  std::list<double> xCoeff;
  std::list<double> yCoeff;

  xCoeff.push_back(1.0);
  xCoeff.push_back(0.0);

  bool exceptionThrown = false;
  try
  {
    AF2PntTrnsfrmLnrV* shouldFail =
        new AF2PntTrnsfrmLnrV(emptyVertexList, xCoeff, yCoeff);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(MeshKit::MK_BAD_INPUT, mkError.error_code());
    std::cout << "Good: The error is a bad input error." << std::endl;
    std::cout << "The error message is:\n" << mkError.what() << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "PASS: A mismatch in the size of the vertex list and the\n"
      << "  x-coefficient list in the constructor caused an exception."
      << std::endl;
}

void testListConstructErrorY()
{
  std::list<const AF2RuleExistVertex*> emptyVertexList;
  std::list<double> xCoeff;
  std::list<double> yCoeff;

  yCoeff.push_back(0.0);
  yCoeff.push_back(1.0);

  bool exceptionThrown = false;
  try
  {
    AF2PntTrnsfrmLnrV* shouldFail =
        new AF2PntTrnsfrmLnrV(emptyVertexList, xCoeff, yCoeff);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(MeshKit::MK_BAD_INPUT, mkError.error_code());
    std::cout << "Good: The error is a bad input error." << std::endl;
    std::cout << "The error message is:\n" << mkError.what() << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "PASS: A mismatch in the size of the vertex list and the\n"
      << "  y-coefficient list in the constructor caused an exception."
      << std::endl;
}

void testVectorConstructErrorX()
{
  std::vector<const AF2RuleExistVertex*> emptyVertexVector;
  std::vector<double> xCoeff;
  std::vector<double> yCoeff;

  xCoeff.push_back(1.0);
  xCoeff.push_back(0.0);

  bool exceptionThrown = false;
  try
  {
    AF2PntTrnsfrmLnrV* shouldFail =
        new AF2PntTrnsfrmLnrV(emptyVertexVector, xCoeff, yCoeff);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(MeshKit::MK_BAD_INPUT, mkError.error_code());
    std::cout << "Good: The error is a bad input error." << std::endl;
    std::cout << "The error message is:\n" << mkError.what() << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "PASS: A mismatch in the size of the vertex vector and the\n"
      << "  x-coefficient vector in the constructor caused an exception."
      << std::endl;
}

void testVectorConstructErrorY()
{
  std::vector<const AF2RuleExistVertex*> emptyVertexVector;
  std::vector<double> xCoeff;
  std::vector<double> yCoeff;

  yCoeff.push_back(0.0);
  yCoeff.push_back(1.0);

  bool exceptionThrown = false;
  try
  {
    AF2PntTrnsfrmLnrV* shouldFail =
        new AF2PntTrnsfrmLnrV(emptyVertexVector, xCoeff, yCoeff);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(MeshKit::MK_BAD_INPUT, mkError.error_code());
    std::cout << "Good: The error is a bad input error." << std::endl;
    std::cout << "The error message is:\n" << mkError.what() << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "PASS: A mismatch in the size of the vertex vector and the\n"
      << "  y-coefficient vector in the constructor caused an exception."
      << std::endl;
}

void testUnboundError()
{
  std::vector<const AF2RuleExistVertex*> vertexVector;
  AF2RuleExistVertex* existVertex = new AF2RuleExistVertex(0.0, 0.0);
  vertexVector.push_back(existVertex);
  std::vector<double> xCoeff;
  xCoeff.push_back(1.0);
  xCoeff.push_back(0.0);
  std::vector<double> yCoeff;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  
  AF2PntTrnsfrmLnrV pntTrnsfrm(vertexVector, xCoeff, yCoeff);
  AF2Binding emptyBinding;
  AF2Point2D aPoint(0.5, 0.5);

  bool exceptionThrown = false;
  try
  {
    pntTrnsfrm.transformPoint(aPoint, emptyBinding);
  }
  catch (MeshKit::Error& mkError)
  {
    CHECK_EQUAL(MeshKit::MK_BAD_INPUT, mkError.error_code());
    std::cout << "Good: The error is a bad input error." << std::endl;
    std::cout << "The error message is:\n" << mkError.what() << std::endl;
    exceptionThrown = true;
  }

  CHECK(exceptionThrown);
  std::cout << "PASS: A binding that does not have bound values for all of\n"
      << "  the vertices that the linear transformation uses causes an error."
      << std::endl;

  delete existVertex;
}

void testSingleVertexJustX()
{
  std::vector<const AF2RuleExistVertex*> vertexVector;
  AF2RuleExistVertex* existVertex = new AF2RuleExistVertex(0.0, 0.0);
  vertexVector.push_back(existVertex);
  std::vector<double> xCoeff;
  xCoeff.push_back(1.0);
  xCoeff.push_back(0.0);
  std::vector<double> yCoeff;
  yCoeff.push_back(0.0);
  yCoeff.push_back(0.0);
  AF2PntTrnsfrmLnrV pntTrnsfrm(vertexVector, xCoeff, yCoeff);

  AF2Binding binding;
  AF2Point2D boundPnt(0.75, 0.25);
  binding.bind(existVertex, &boundPnt);

  AF2Point2D pointArg(0.5, 0.5);
  AF2Point2D resultPoint = pntTrnsfrm.transformPoint(pointArg, binding);
  CHECK_REAL_EQUAL(1.25, resultPoint.getX(), 1e-14);
  std::cout << "Good: The resulting point has an x-coordinate of 1.25"
      << std::endl;
  CHECK_REAL_EQUAL(0.5, resultPoint.getY(), 1e-14);
  std::cout << "Good: The resulting point has a y-coordinate of 0.5"
      << std::endl;

  std::cout << "PASS: A translation in the x direction works."
      << std::endl;

  delete existVertex;
}

void testSingleVertexJustY()
{
  std::vector<const AF2RuleExistVertex*> vertexVector;
  AF2RuleExistVertex* existVertex = new AF2RuleExistVertex(1.0, 1.0);
  vertexVector.push_back(existVertex);
  std::vector<double> xCoeff;
  xCoeff.push_back(0.0);
  xCoeff.push_back(0.0);
  std::vector<double> yCoeff;
  yCoeff.push_back(0.0);
  yCoeff.push_back(1.0);
  AF2PntTrnsfrmLnrV pntTrnsfrm(vertexVector, xCoeff, yCoeff);

  AF2Binding binding;
  AF2Point2D boundPnt(1.75, 1.25);
  binding.bind(existVertex, &boundPnt);

  AF2Point2D pointArg(0.5, 0.5);
  AF2Point2D resultPoint = pntTrnsfrm.transformPoint(pointArg, binding);
  CHECK_REAL_EQUAL(0.5, resultPoint.getX(), 1e-14);
  std::cout << "Good: The resulting point has an x-coordinate of 0.5"
      << std::endl;
  CHECK_REAL_EQUAL(0.75, resultPoint.getY(), 1e-14);
  std::cout << "Good: The resulting point has a y-coordinate of 0.75"
      << std::endl;

  std::cout << "PASS: A translation in the y direction works."
      << std::endl;

  delete existVertex;
}

void testSingleVertexFull()
{
  std::vector<const AF2RuleExistVertex*> vertexVector;
  AF2RuleExistVertex* existVertex = new AF2RuleExistVertex(0.0, 0.0);
  vertexVector.push_back(existVertex);
  std::vector<double> xCoeff;
  xCoeff.push_back(1.0);
  xCoeff.push_back(-2.0);
  std::vector<double> yCoeff;
  yCoeff.push_back(-3.0);
  yCoeff.push_back(4.0);
  AF2PntTrnsfrmLnrV pntTrnsfrm(vertexVector, xCoeff, yCoeff);

  AF2Binding binding;
  AF2Point2D boundPnt(1.0, 10.0);
  binding.bind(existVertex, &boundPnt);

  AF2Point2D pointArg(0.0, 0.0);
  AF2Point2D resultPoint = pntTrnsfrm.transformPoint(pointArg, binding);
  CHECK_REAL_EQUAL(-19.0, resultPoint.getX(), 1e-14);
  std::cout << "Good: The resulting point has an x-coordinate of 8.0"
      << std::endl;
  CHECK_REAL_EQUAL(37.0, resultPoint.getY(), 1e-14);
  std::cout << "Good: The resulting point has a y-coordinate of 37.0"
      << std::endl;

  std::cout << "PASS: A translation involving a single vertex and all\n"
      << "  four coefficients works." << std::endl;

  delete existVertex;
}

void testMultipleVertices()
{
  std::vector<const AF2RuleExistVertex*> vertexVector;
  AF2RuleExistVertex firstExistVertex(0.0, 0.0);
  AF2RuleExistVertex secondExistVertex(1.0, 0.0);
  vertexVector.push_back(&firstExistVertex);
  vertexVector.push_back(&secondExistVertex);
  std::vector<double> xCoeff;
  xCoeff.push_back(1.0);
  xCoeff.push_back(0.0);
  xCoeff.push_back(0.0);
  xCoeff.push_back(3.0);
  std::vector<double> yCoeff;
  yCoeff.push_back(0.0);
  yCoeff.push_back(5.0);
  yCoeff.push_back(1.0);
  yCoeff.push_back(0.0);
  AF2PntTrnsfrmLnrV pntTrnsfrm(vertexVector, xCoeff, yCoeff);

  AF2Binding binding;
  AF2Point2D firstBoundPnt(0.125, 0.25);
  AF2Point2D secondBoundPnt(1.5, 0.0625);
  binding.bind(&firstExistVertex, &firstBoundPnt);
  binding.bind(&secondExistVertex, &secondBoundPnt);

  AF2Point2D pointArg(0.0, 0.0);
  AF2Point2D resultPoint = pntTrnsfrm.transformPoint(pointArg, binding);
  CHECK_REAL_EQUAL(0.3125, resultPoint.getX(), 1e-14);
  std::cout << "Good: The resulting point has an x-coordinate of 0.3125"
      << std::endl;
  CHECK_REAL_EQUAL(1.75, resultPoint.getY(), 1e-14);
  std::cout << "Good: The resulting point has a y-coordinate of 1.75"
      << std::endl;

  std::cout << "PASS: A translation involving multiple vertices works."
      << std::endl;
}
