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
 *
 * With the detection of point containment and detection of
 * intersection, the free zone has the capability to not detect
 * a containment or intersection in cases where the query point
 * or line segment is approximately equal to some portion of
 * the boundary.
 */

// C++
#include <iostream>
#include <list>

// MeshKit
#include "meshkit/AF2FreeZone.hpp"
#include "meshkit/AF2Point2D.hpp"

// MeshKit test utilities
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
void testSquareContainsVertex();
void testSquareNotContains();
void testPentagonContains();
void testPentagonNotContains();
void testPentagonNotContainsVertex();
void testDiamondNearContains();
void testDiamondNotContainsVertex();
void testSquareIntersect();
void testSquareNotIntersect();
void testSquareNotIntersectSegment();
void testSquareIntersectOutwardRay();
void testPentagonIntersect();
void testPentagonNearIntersect();
void testPentagonNotIntersectBB();
void testPentagonNotIntersectCW();
void testPentagonNotIntersectZoneCW();
void testPentagonNotIntersectZoneCCW();
void testPentagonIntersectSegment();
void testPentagonNotIntersectOutwardRay();
void testPentagonIntersectInwardRay();
void testDiamondNotIntersectSegment();
void testDiamondNotIntersectOutwardRay();
void testDiamondIntersectInwardRay();

int main(int argc, char **argv)
{

  int num_fail = 0;

  num_fail += RUN_TEST(testSquareConvex);
  num_fail += RUN_TEST(testPentagonConvex);
  num_fail += RUN_TEST(testClockwiseSquare);
  num_fail += RUN_TEST(testLocalNonConvex);
  num_fail += RUN_TEST(testGlobalNonConvex);
  num_fail += RUN_TEST(testSquareContains);
  num_fail += RUN_TEST(testSquareContainsVertex);
  num_fail += RUN_TEST(testSquareNotContains);
  num_fail += RUN_TEST(testPentagonContains);
  num_fail += RUN_TEST(testPentagonNotContains);
  num_fail += RUN_TEST(testPentagonNotContainsVertex);
  num_fail += RUN_TEST(testDiamondNearContains);
  num_fail += RUN_TEST(testDiamondNotContainsVertex);
  num_fail += RUN_TEST(testSquareIntersect);
  num_fail += RUN_TEST(testSquareNotIntersect);
  num_fail += RUN_TEST(testSquareNotIntersectSegment);
  num_fail += RUN_TEST(testSquareIntersectOutwardRay);
  num_fail += RUN_TEST(testPentagonIntersect);
  num_fail += RUN_TEST(testPentagonNearIntersect);
  num_fail += RUN_TEST(testPentagonNotIntersectBB);
  num_fail += RUN_TEST(testPentagonNotIntersectCW);
  num_fail += RUN_TEST(testPentagonNotIntersectZoneCW);
  num_fail += RUN_TEST(testPentagonNotIntersectZoneCCW);
  num_fail += RUN_TEST(testPentagonIntersectSegment);
  num_fail += RUN_TEST(testPentagonNotIntersectOutwardRay);
  num_fail += RUN_TEST(testPentagonIntersectInwardRay);
  num_fail += RUN_TEST(testDiamondNotIntersectSegment);
  num_fail += RUN_TEST(testDiamondNotIntersectOutwardRay);
  num_fail += RUN_TEST(testDiamondIntersectInwardRay);

  return num_fail;
}

AF2FreeZone makeSquareFreeZone()
{
  std::list<AF2Point2D> squareBndryPnts;
  AF2Point2D alpha(0, 0);
  AF2Point2D bravo(1, 0);
  AF2Point2D charlie(1, 1);
  AF2Point2D delta(0, 1);
  squareBndryPnts.push_back(alpha);
  squareBndryPnts.push_back(bravo);
  squareBndryPnts.push_back(charlie);
  squareBndryPnts.push_back(delta);

  AF2FreeZone squareFreeZone(squareBndryPnts);
  return squareFreeZone;
}

AF2FreeZone makePentagonFreeZone()
{
  // make a free zone in roughtly the shape of a regular pentagon
  // first two points at (0, 0) and (1e-14, 0)
  std::list<AF2Point2D> pentBndryPnts;
  AF2Point2D alpha(0, 0);
  AF2Point2D bravo(1e-14, 0);
  AF2Point2D charlie(1.39735982487e-14, 1.1091842589e-14);
  AF2Point2D delta(0.5e-14, 1.74370722192e-14);
  AF2Point2D echo(-0.39735982487e-14, 1.1091842589e-14);
  pentBndryPnts.push_back(alpha);
  pentBndryPnts.push_back(bravo);
  pentBndryPnts.push_back(charlie);
  pentBndryPnts.push_back(delta);
  pentBndryPnts.push_back(echo);

  AF2FreeZone pentFreeZone(pentBndryPnts);
  return pentFreeZone;
}

AF2FreeZone makeDiamondFreeZone()
{
  std::list<AF2Point2D> diamondBndryPnts;
  AF2Point2D alpha(0, 0);
  AF2Point2D bravo(0.5e40, 0.5e40);
  AF2Point2D charlie(0, 1e40);
  AF2Point2D delta(-0.5e-40, 0.5e40);
  diamondBndryPnts.push_back(alpha);
  diamondBndryPnts.push_back(bravo);
  diamondBndryPnts.push_back(charlie);
  diamondBndryPnts.push_back(delta);

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
  std::list<AF2Point2D> squareBndryPnts;
  AF2Point2D alpha(0, 0);
  AF2Point2D bravo(0, 1);
  AF2Point2D charlie(1, 1);
  AF2Point2D delta(1, 0);
  squareBndryPnts.push_back(alpha);
  squareBndryPnts.push_back(bravo);
  squareBndryPnts.push_back(charlie);
  squareBndryPnts.push_back(delta);

  AF2FreeZone clockwiseSquareFreeZone(squareBndryPnts);
  CHECK(!clockwiseSquareFreeZone.isConvex());
  std::cout <<
     "PASS: The clockwise traversal of the square is considered nonconvex."
     << std::endl;
}

void testLocalNonConvex()
{
  std::list<AF2Point2D> nonConvexBndryPnts;
  AF2Point2D alpha(0, 0);
  AF2Point2D bravo(1, 0);
  AF2Point2D charlie(0.9999, 0.5);
  AF2Point2D delta(0.9999, 1);
  AF2Point2D echo(0, 1);
  nonConvexBndryPnts.push_back(alpha);
  nonConvexBndryPnts.push_back(bravo);
  nonConvexBndryPnts.push_back(charlie);
  nonConvexBndryPnts.push_back(delta);
  nonConvexBndryPnts.push_back(echo);

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
  std::list<AF2Point2D> nonConvexBndryPnts;
  AF2Point2D alpha(0, 0);
  AF2Point2D bravo(1, 0);
  AF2Point2D charlie(0.5, 1);
  AF2Point2D delta(-0.125, 0.25);
  AF2Point2D echo(0.375, -0.25);
  AF2Point2D foxtrot(0.875, 0.75);
  AF2Point2D golf(0, 1);
  nonConvexBndryPnts.push_back(alpha);
  nonConvexBndryPnts.push_back(bravo);
  nonConvexBndryPnts.push_back(charlie);
  nonConvexBndryPnts.push_back(delta);
  nonConvexBndryPnts.push_back(echo);
  nonConvexBndryPnts.push_back(foxtrot);
  nonConvexBndryPnts.push_back(golf);

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
  AF2Point2D myPoint(0.75, 0.875);
  CHECK(squareFreeZone.nearContains(myPoint));
  std::cout << "PASS: The square free zone contains " << myPoint
      << "." << std::endl;
}

void testSquareContainsVertex()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  AF2Point2D myPoint(1, 0);
  CHECK(squareFreeZone.nearContains(myPoint, true));
  std::cout << "PASS: The square free zone contains " << myPoint
      << "\n  when passed an argument to contain boundary points."
      << std::endl;
}

void testSquareNotContains()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  AF2Point2D myPoint(0, 1.0625);
  CHECK(!squareFreeZone.nearContains(myPoint));
  std::cout << "PASS: The square free zone does not contain " << myPoint
      << "." << std::endl;
}

void testDiamondNearContains()
{
  AF2FreeZone diamondFreeZone = makeDiamondFreeZone();
  AF2Point2D myPoint(0.25e40, 0.249999999999999e40);
  CHECK(diamondFreeZone.nearContains(myPoint));
  std::cout << "PASS: The diamond free zone nearly contains\n"
      << "  (0.25e40, 0.249999999999999e40)." << std::endl;
}

void testDiamondNotContainsVertex()
{
  AF2FreeZone diamondFreeZone = makeDiamondFreeZone();
  AF2Point2D myPoint(0.0, 9.9999999999999e39);
  CHECK(!diamondFreeZone.nearContains(myPoint));
  std::cout << "PASS: The diamond free zone says it does not contain\n"
      << "  (0, 9.9999999999999e39) because it is near a boundary vertex"
      << std::endl;
}

void testPentagonContains()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D myPoint(0.625e-14, 1.25e-14);
  CHECK(pentFreeZone.nearContains(myPoint));
  std::cout << "PASS: The pentagon free zone contains " << myPoint
      << "." << std::endl;
}

void testPentagonNotContains()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D myPoint(1.25e-14, 0.125e-14);
  CHECK(!pentFreeZone.nearContains(myPoint));
  std::cout
      << "PASS: The pentagon free zone does not contain " << myPoint
      << "." << std::endl;
}

void testPentagonNotContainsVertex()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D myPoint(1.39735982487e-14, 1.1091842589e-14);
  CHECK(!pentFreeZone.nearContains(myPoint));
  std::cout
      << "PASS: The pentagon free zone says it does not contain\n"
      << "  (1.39735982487e-14, 1.1091842589e-14)\n"
      << "  because it is a boundary vertex" << std::endl;
}

void testSquareIntersect()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  AF2Point2D startEdge(1.2499999, 0.75);
  AF2Point2D endEdge(0.75, 1.25);
  CHECK(squareFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The square free zone intersects an edge from\n"
      << "  (1.2499999, 0.75) to (0.75, 1.25)" << std::endl;
}

void testSquareNotIntersect()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  AF2Point2D startEdge(1.2500001, 0.75);
  AF2Point2D endEdge(0.75, 1.25);
  CHECK(!squareFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The square free zone does not intersect an edge from\n"
      << "  (1.2500001, 0.75) to (0.75, 1.25)" << std::endl;
}

void testSquareNotIntersectSegment()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  AF2Point2D startEdge(1, 1);
  AF2Point2D endEdge(0, 1);
  CHECK(!squareFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The square free zone says it does not intersect "
      << "an edge from\n  (1, 1) to (0, 1) because it is a boundary segment"
      << std::endl;
}

void testSquareIntersectOutwardRay()
{
  AF2FreeZone squareFreeZone = makeSquareFreeZone();
  AF2Point2D startEdge(1, 1);
  AF2Point2D endEdge(1.5, 1.5);
  CHECK(squareFreeZone.nearIntersects(startEdge, endEdge, true));
  std::cout << "PASS: The square free zone intersects an edge from\n"
      << "  (1, 1) to (1.5, 1.5) when passed an argument to contain its "
      << "boundary" << std::endl;
}

void testPentagonIntersect()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(1.0e-14, 1.75e-14);
  AF2Point2D endEdge(0.5e-14, 1.7437e-14);
  CHECK(pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone intersects an edge from\n"
      << "  (1.0e-14, 1.75e-14) to (0.5e-14, 1.7437e-14)" << std::endl;
}

void testPentagonNearIntersect()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(1.0e-14, 1.75e-14);
  AF2Point2D endEdge(0.500010213437881e-14, 1.7437e-14);
  // As of 02/17/2016, computations in the code suggested that the
  // actual change from intersection to non-intersection is between
  // 0.5000102134378802e-14 and 0.5000102134378803e-14, but the
  // rounding error involved means this might be just noise.
  // Values between 0.50001021343788022e-14 and 0.50001021343788029e-14
  // inclusive produced a result of 0, i.e., on the boundary of intersection.
  CHECK(pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone nearly intersects an edge from\n"
      << "  (1.0e-14, 1.75e-14) to (0.500010213437881e-14, 1.7437e-14)"
      << std::endl;
}

void testPentagonNotIntersectBB()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(-0.39735982488e-14, 1.0-14);
  AF2Point2D endEdge(-0.39735982488e-14, -1.0e-14);
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (-0.39735982488e-14, 1.0e-14) to "
      << "(-0.39735982488e-14, -1.0e-14)" << std::endl;
}

void testPentagonNotIntersectCW()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(1.0e-14, 1.75e-14);
  AF2Point2D endEdge(0.50002e-14, 1.7437e-14);
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (1.0e-14, 1.75e-14) to (0.50002e-14, 1.7437e-14)" << std::endl;
}

void testPentagonNotIntersectZoneCW()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(1.3e-14, 1.3e-14);
  AF2Point2D endEdge(1.5e-14, 0.95e-14);
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (1.3e-14, 1.3e-14) to (1.5e-14, 0.95e-14)" << std::endl;
}

void testPentagonNotIntersectZoneCCW()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(1.5e-14, 0.95e-14);
  AF2Point2D endEdge(1.3e-14, 1.3e-14);
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone does not intersect an edge from\n"
      << "  (1.5e-14, 0.95e-14) to (1.3e-14, 1.3e-14)" << std::endl;
}

void testPentagonIntersectSegment()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(1.39735982487e-14, 1.1091842589e-14);
  AF2Point2D endEdge(0.5e-14, 1.74370722192e-14);
  CHECK(pentFreeZone.nearIntersects(startEdge, endEdge, true));
  std::cout << "PASS: The pentagon free zone intersects an edge from\n"
      << "  (1.39735982487e-14, 1.1091842589e-14) to"
      << " (0.5e-14, 1.74370722192e-14)\n"
      << "  when passed an argument to contain boundary segments."
      << std::endl;
}

void testPentagonNotIntersectOutwardRay()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(5.0e-15, 1.74370722192e-14);
  AF2Point2D endEdge(7.0e-15, 3.0e-14);
  CHECK(!pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone says it does not intersect "
      << "an edge from\n"
      << "  (5.0e-15, 1.74370722192e-14) to (7.0e-15, 3.0e-14)\n"
      << "  because the only intersection is near an endpoint that is\n"
      << "  near a boundary point."
      << std::endl;
}

void testPentagonIntersectInwardRay()
{
  AF2FreeZone pentFreeZone = makePentagonFreeZone();
  AF2Point2D startEdge(5.0e-15, 1.74370722192e-14);
  AF2Point2D endEdge(6.0e-15, 1.1e-14);
  CHECK(pentFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The pentagon free zone intersects an edge from\n"
      << "  (5.0e-15, 1.74370722192e-14) to (6.0e-15, 1.1e-14)."
      << std::endl;
}

void testDiamondNotIntersectSegment()
{
  AF2FreeZone diamondFreeZone = makeDiamondFreeZone();
  AF2Point2D startEdge(0.0, 9.9999999999999e39);
  AF2Point2D endEdge(5.0e39, 5.0e39);
  CHECK(!diamondFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The diamond free zone says it does not intersect\n"
      << "  an edge from (0.0, 9.9999999999999e39) to (5.0e39, 5.0e39)\n"
      << "  because it is approximately equal to a boundary edge." << std::endl;
}

void testDiamondNotIntersectOutwardRay()
{
  AF2FreeZone diamondFreeZone = makeDiamondFreeZone();
  AF2Point2D startEdge(-1.0e39, 1.4e40);
  AF2Point2D endEdge(0.0, 9.9999999999999e39);
  CHECK(!diamondFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The diamond free zone says it does not intersect "
      << "an edge from\n  " << startEdge
      << " to (0.0, 9.9999999999999e39)\n"
      << "  because the only intersection is near an endpoint that is\n"
      << "  near a boundary point." << std::endl;
}

void testDiamondIntersectInwardRay()
{
  AF2FreeZone diamondFreeZone = makeDiamondFreeZone();
  AF2Point2D startEdge(-1.0e39, 6.1e39);
  AF2Point2D endEdge(0.0, 9.9999999999999e39);
  CHECK(diamondFreeZone.nearIntersects(startEdge, endEdge));
  std::cout << "PASS: The diamond free zone intersects an edge from\n"
      << startEdge << " to (0.0, 9.9999999999999e39)\n" << std::endl;
}
