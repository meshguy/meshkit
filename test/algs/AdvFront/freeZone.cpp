/** 
 * \file freeZone.cpp \test
 *
 * Test the two-dimensional free zone object.
 *
 * The free zone has four key capabilities --
 * (1) construction,
 * (2) verification of convexity,
 * (3) detection of point containment, and
 * (4) detection of intersection with a line segment.
 */

#include <iostream>
#include <list>

#include "meshkit/AF2FreeZone.hpp"
#include "meshkit/Matrix.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

AF2FreeZone makeSquareFreeZone();
AF2FreeZone makePentagonFreeZone();
AF2FreeZone makeDiamondFreeZone();
void testSquareConvex();
void testPentagonConvex();
void testClockwiseSquare();
void testLocalNonConvex();
void testGlobalNonConvex();
void testSquareContains();
void testSquareNotContains();
void testPentagonContains();
void testPentagonNotContains();
void testDiamondNearContains();
void testSquareIntersect();
void testSquareNotIntersect();
void testPentagonIntersect();
void testPentagonNearIntersect();
void testPentagonNotIntersectBB();
void testPentagonNotIntersectCW();
void testPentagonNotIntersectZoneCW();
void testPentagonNotIntersectZoneCCW();

int main(int argc, char **argv)
{

  int num_fail = 0;

  num_fail += RUN_TEST(testSquareConvex);
  num_fail += RUN_TEST(testPentagonConvex);
  num_fail += RUN_TEST(testClockwiseSquare);
  num_fail += RUN_TEST(testLocalNonConvex);
  num_fail += RUN_TEST(testGlobalNonConvex);
  num_fail += RUN_TEST(testSquareContains);
  num_fail += RUN_TEST(testSquareNotContains);
  num_fail += RUN_TEST(testPentagonContains);
  num_fail += RUN_TEST(testPentagonNotContains);
  num_fail += RUN_TEST(testDiamondNearContains);
  num_fail += RUN_TEST(testSquareIntersect);
  num_fail += RUN_TEST(testSquareNotIntersect);
  num_fail += RUN_TEST(testPentagonIntersect);
  num_fail += RUN_TEST(testPentagonNearIntersect);
  num_fail += RUN_TEST(testPentagonNotIntersectBB);
  num_fail += RUN_TEST(testPentagonNotIntersectCW);
  num_fail += RUN_TEST(testPentagonNotIntersectZoneCW);
  num_fail += RUN_TEST(testPentagonNotIntersectZoneCCW);

  return num_fail;
}

AF2FreeZone makeSquareFreeZone()
{
  std::list<Vector<2> > squareBndryPnts;
  Vector<2> myPoint;
  myPoint[0] = 0;
  myPoint[1] = 0;
  squareBndryPnts.push_back(myPoint);
  myPoint[0] = 1;
  squareBndryPnts.push_back(myPoint);
  myPoint[1] = 1;
  squareBndryPnts.push_back(myPoint);
  myPoint[0] = 0;
  squareBndryPnts.push_back(myPoint);

  AF2FreeZone squareFreeZone(squareBndryPnts);
  return squareFreeZone;
}

AF2FreeZone makePentagonFreeZone()
{
  // make a free zone in roughtly the shape of a regular pentagon
  // first two points at (0, 0) and (1e-14, 0)
  std::list<Vector<2> > pentBndryPnts;
  Vector<2> myPoint;
  myPoint[0] = 0;
  myPoint[1] = 0;
  pentBndryPnts.push_back(myPoint);
  myPoint[0] = 1e-14;
  pentBndryPnts.push_back(myPoint);
  myPoint[0] = 1.39735982487e-14;
  myPoint[1] = 1.1091842589e-14;
  pentBndryPnts.push_back(myPoint);
  myPoint[0] = 0.5e-14;
  myPoint[1] = 1.74370722192e-14;
  pentBndryPnts.push_back(myPoint);
  myPoint[0] = -0.39735982487e-14;
  myPoint[1] = 1.1091842589e-14;

  AF2FreeZone pentFreeZone(pentBndryPnts);
  return pentFreeZone;
}

AF2FreeZone makeDiamondFreeZone()
{
  std::list<Vector<2> > diamondBndryPnts;
  Vector<2> myPoint;
  myPoint[0] = 0;
  myPoint[1] = 0;
  diamondBndryPnts.push_back(myPoint);
  myPoint[0] = 0.5e40;
  myPoint[1] = 0.5e40;
  diamondBndryPnts.push_back(myPoint);
  myPoint[0] = 0;
  myPoint[1] = 1e40;
  diamondBndryPnts.push_back(myPoint);
  myPoint[0] = -0.5e40;
  myPoint[1] = 0.5e40;
  diamondBndryPnts.push_back(myPoint);

  AF2FreeZone diamondFreeZone(diamondBndryPnts);
  return diamondFreeZone;
}

void testSquareConvex()
{
  AF2FreeZone convexFreeZone = makeSquareFreeZone();
  CHECK(convexFreeZone.isConvex());
  std::cout << "PASS: The square free zone is convex." << std::endl;
}

void testPentagonConvex()
{
  AF2FreeZone convexFreeZone = makePentagonFreeZone();
  CHECK(convexFreeZone.isConvex());
  std::cout << "PASS: The pentagon free zone is convex." << std::endl;
}

void testClockwiseSquare()
{
  std::list<Vector<2> > squareBndryPnts;
  Vector<2> myPoint;
  myPoint[0] = 0;
  myPoint[1] = 0;
  squareBndryPnts.push_back(myPoint);
  myPoint[1] = 1;
  squareBndryPnts.push_back(myPoint);
  myPoint[0] = 1;
  squareBndryPnts.push_back(myPoint);
  myPoint[1] = 0;
  squareBndryPnts.push_back(myPoint);

  AF2FreeZone clockwiseSquareFreeZone(squareBndryPnts);
  CHECK(!clockwiseSquareFreeZone.isConvex());
  std::cout <<
     "PASS: The clockwise traversal of the square is considered nonconvex."
     << std::endl;
}

void testLocalNonConvex()
{
  std::list<Vector<2> > nonConvexBndryPnts;
  Vector<2> myPoint;
  myPoint[0] = 0;
  myPoint[1] = 0;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 1;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 0.9999;
  myPoint[1] = 0.5;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[1] = 1;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 0;
  nonConvexBndryPnts.push_back(myPoint);

  AF2FreeZone nonConvexFreeZone(nonConvexBndryPnts);
  CHECK(!nonConvexFreeZone.isConvex());
  std::cout << "PASS: The counterclockwise traversal of a nonconvex region\n"
     << "  is considered nonconvex."
     << std::endl;
}

/**
 * Test that even though every sequence of three points along the boundary
 * forms a counterclockwise triangle, the convexity test detects that
 * a region is not convex because there are multiple cycles before
 * the counterclockwise traversal is closed.
 */
void testGlobalNonConvex()
{
  std::list<Vector<2> > nonConvexBndryPnts;
  Vector<2> myPoint;
  myPoint[0] = 0;
  myPoint[1] = 0;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 1;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 0.5;
  myPoint[1] = 1;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = -0.125;
  myPoint[1] = 0.25;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 0.375;
  myPoint[1] = -0.25;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 0.875;
  myPoint[1] = 0.75;
  nonConvexBndryPnts.push_back(myPoint);
  myPoint[0] = 0;
  myPoint[1] = 1;
  nonConvexBndryPnts.push_back(myPoint);

  AF2FreeZone nonConvexFreeZone(nonConvexBndryPnts);
  CHECK(!nonConvexFreeZone.isConvex());
  std::cout << "PASS: The counterclockwise traversal of a nonconvex region\n"
     << "  is considered nonconvex, even if every consecutive triple is\n"
     << "  counterclockwise."
     << std::endl;
}

void testSquareContains()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  Vector<2> myPoint;
  myPoint[0] = 0.75;
  myPoint[1] = 0.875;
  CHECK(squareFreeZone.nearContains(myPoint));
  std::cout << "PASS: The square free zone contains (0.75, 0.875)."
      << std::endl;
}

void testSquareNotContains()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  Vector<2> myPoint;
  myPoint[0] = 0;
  myPoint[1] = 1.0625;
  CHECK(!squareFreeZone.nearContains(myPoint));
  std::cout << "PASS: The square free zone does not contain (0, 1.0625)."
      << std::endl;
}

void testDiamondNearContains()
{
  AF2FreeZone diamondFreeZone = makeDiamondFreeZone();
  Vector<2> myPoint;
  myPoint[0] = 0.25e40;
  myPoint[1] = 0.249999999999999e40;
  CHECK(diamondFreeZone.nearContains(myPoint));
  std::cout << "PASS: The diamond free zone nearly contains\n"
      << "(0.25e40, 0.249999999999999e40)." << std::endl;
}

void testPentagonContains()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> myPoint;
  myPoint[0] = 0.625e-14;
  myPoint[1] = 1.25e-14;
  CHECK(pentFreeZone.nearContains(myPoint));
  std::cout << "PASS: The pentagon free zone contains (0.625e-14, 1.25e-14)."
      << std::endl;
}

void testPentagonNotContains()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> myPoint;
  myPoint[0] = 1.25e-14;
  myPoint[1] = 0.125e-14;
  CHECK(!pentFreeZone.nearContains(myPoint));
  std::cout
      << "PASS: The pentagon free zone does not contain (1.25e-14, 0.125e-14)."
      << std::endl;
}

void testSquareIntersect()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  Vector<2> startEdge;
  startEdge[0] = 1.2499999;
  startEdge[1] = 0.75;
  Vector<2> endEdge;
  endEdge[0] = 0.75;
  endEdge[1] = 1.25;
  CHECK(squareFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The square free zone intersects an edge from\n"
      << "  (1.2499999, 0.75) to (0.75, 1.25)" << std::endl;
}

void testSquareNotIntersect()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  Vector<2> startEdge;
  startEdge[0] = 1.2500001;
  startEdge[1] = 0.75;
  Vector<2> endEdge;
  endEdge[0] = 0.75;
  endEdge[1] = 1.25;
  CHECK(!squareFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The square free zone does not intersect an edge from\n"
      << "  (1.2500001, 0.75) to (0.75, 1.25)" << std::endl;
}

void testPentagonIntersect()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> startEdge;
  startEdge[0] = 1.0e-14;
  startEdge[1] = 1.75e-14;
  Vector<2> endEdge;
  endEdge[0] = 0.5e-14;
  endEdge[1] = 1.7437e-14;
  CHECK(pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone intersects an edge from\n"
      << "  (1.0e-14, 1.75e-14) to (0.5e-14, 1.7437e-14)" << std::endl;
}

void testPentagonNearIntersect()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> startEdge;
  startEdge[0] = 1.0e-14;
  startEdge[1] = 1.75e-14;
  Vector<2> endEdge;
  // As of 02/17/2016, computations in the code suggested that the
  // actual change from intersection to non-intersection is between
  // 0.5000102134378802e-14 and 0.5000102134378803e-14, but the
  // rounding error involved means this might be just noise.
  // Values between 0.50001021343788022e-14 and 0.50001021343788029e-14
  // inclusive produced a result of 0, i.e., on the boundary of intersection.
  endEdge[0] = 0.500010213437881e-14;
  endEdge[1] = 1.7437e-14;
  CHECK(pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone nearly intersects an edge from\n"
      << "  (1.0e-14, 1.75e-14) to (0.500010213437881e-14, 1.7437e-14)"
      << std::endl;
}

void testPentagonNotIntersectBB()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> startEdge;
  startEdge[0] = -1.0e-30;
  startEdge[1] = 1.0e-14;
  Vector<2> endEdge;
  endEdge[0] = -1.0e-30;
  endEdge[1] = -1.0e-14;
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (-1.0e-30, 1.0e-14) to (-1.0e-30, -1.0e-14)" << std::endl;
}

void testPentagonNotIntersectCW()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> startEdge;
  startEdge[0] = 1.0e-14;
  startEdge[1] = 1.75e-14;
  Vector<2> endEdge;
  endEdge[0] = 0.50002e-14;
  endEdge[1] = 1.7437e-14;
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (1.0e-14, 1.75e-14) to (0.50002e-14, 1.7437e-14)" << std::endl;
}

void testPentagonNotIntersectZoneCW()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> startEdge;
  startEdge[0] = 1.3e-14;
  startEdge[1] = 1.3e-14;
  Vector<2> endEdge;
  endEdge[0] = 1.5e-14;
  endEdge[1] = 0.95e-14;
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (1.3e-14, 1.3e-14) to (1.5e-14, 0.95e-14)" << std::endl;
}

void testPentagonNotIntersectZoneCCW()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  Vector<2> startEdge;
  startEdge[0] = 1.5e-14;
  startEdge[1] = 0.95e-14;
  Vector<2> endEdge;
  endEdge[0] = 1.3e-14;
  endEdge[1] = 1.3e-14;
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (1.5e-14, 0.95e-14) to (1.3e-14, 1.3e-14)" << std::endl;
}
