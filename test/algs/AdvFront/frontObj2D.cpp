/**
 * \file frontObj2D.cpp \test
 *
 * Test AF2Front, i.e., the object that manages the points and edges
 * that are on the advancing front of the two-dimensional advancing
 * front algorithm, and AF2DfltPlaneProjMaker, the object that implements
 * the default method of creating AF2PlaneProjection local transforms.
 */

// C++
#include <cstddef>
#include <iostream>
#include <list>

// MeshKit
#include "meshkit/AF2Front.hpp"
#include "meshkit/AF2DfltPlaneProjMaker.hpp"
#include "meshkit/Error.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"

// MeshKit test utilities
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model
#if HAVE_OCC
#define FILE_EXT "stp"
#else
#define FILE_EXT "sat"
#endif

// This variable is at global scope because a new MKCore instance cannot
// be easily constructed after another MKCore instance is deleted
// There are problems with a tag left behind in iGeom.
// Thus this instance is shared across tests.
MeshKit::MKCore* mk = NULL;
// These variables are at global scope because (1) calling deleteAll on
// the MKCore geometry instance appears to cause memory inconsistencies
// with later use of the geometry instance and (2) it is more efficient
// to load the geometry model only once
MeshKit::MEntVector surfs;
MeshKit::ModelEnt* square = NULL;

AF2LocalTransformMaker* makeTransformBuilder();
void initHex(AF2Front* front, AF2Point3D** & pointsAry, AF2Edge3D** & edgesAry);
void deletePoints(AF2Point3D** & pointsAry, int pointsAryLength);
void testSelectEmpty();
void testInitAndDestruct();
void testQualityDecrease();
void testSelectNeighborhood();
void testSelectIsolatedPoint();
void testSelectHourglass();
void testAdvanceMissingPoints();
void testAdvanceInvalidQuality();
void testAdvanceFront();
void testAdvanceHangingEdge();
void testCompletedAdvance();

int main(int argc, char **argv)
{
  mk = new MeshKit::MKCore();

  int num_fail = 0;

  num_fail += RUN_TEST(testSelectEmpty);
  num_fail += RUN_TEST(testInitAndDestruct);
  num_fail += RUN_TEST(testQualityDecrease);
  num_fail += RUN_TEST(testSelectNeighborhood);
  num_fail += RUN_TEST(testSelectIsolatedPoint);
  num_fail += RUN_TEST(testSelectHourglass);
  num_fail += RUN_TEST(testAdvanceMissingPoints);
  num_fail += RUN_TEST(testAdvanceInvalidQuality);
  num_fail += RUN_TEST(testAdvanceFront);
  num_fail += RUN_TEST(testAdvanceHangingEdge);
  num_fail += RUN_TEST(testCompletedAdvance);

  delete mk;

  return num_fail;
}

AF2LocalTransformMaker* makeTransformBuilder()
{
  if (square == NULL)
  {
    // load a square in plane z = 0.5 with -1.0 <= x <= 0 and -0.5 <= y <= 0.5
    std::string file_name = TestDir + "/squaresurf." + FILE_EXT;
    mk->load_geometry_mesh(file_name.c_str(), file_name.c_str());

    mk->get_entities_by_dimension(2, surfs);
    square = *surfs.begin();
  }

  return new AF2DfltPlaneProjMaker(
      square->igeom_instance(), square->geom_handle());
}

void initHex(AF2Front* front, AF2Point3D** & pointsAry, AF2Edge3D** & edgesAry)
{
  pointsAry = new AF2Point3D*[6];
  pointsAry[0] = new AF2Point3D(0, -0.625, -0.25, 0.5);
  pointsAry[1] = new AF2Point3D(1, -0.375, -0.25, 0.5);
  pointsAry[2] = new AF2Point3D(2, -0.25, 0.0, 0.5);
  pointsAry[3] = new AF2Point3D(3, -0.375, 0.25, 0.5);
  pointsAry[4] = new AF2Point3D(4, -0.625, 0.25, 0.5);
  pointsAry[5] = new AF2Point3D(5, -0.75, 0.0, 0.5);

  for (int pi = 0; pi < 6; ++pi)
  {
    pointsAry[pi]->limitDistanceToBoundary(0u);
    front->addPoint(pointsAry[pi]);
  }

  edgesAry = new AF2Edge3D*[6];
  edgesAry[0] = new AF2Edge3D(pointsAry[0], pointsAry[1]);
  edgesAry[1] = new AF2Edge3D(pointsAry[1], pointsAry[2]);
  edgesAry[2] = new AF2Edge3D(pointsAry[2], pointsAry[3]);
  edgesAry[3] = new AF2Edge3D(pointsAry[3], pointsAry[4]);
  edgesAry[4] = new AF2Edge3D(pointsAry[4], pointsAry[5]);
  edgesAry[5] = new AF2Edge3D(pointsAry[5], pointsAry[0]);

  std::list<AF2Edge3D*> edges;
  for (int ei = 0; ei < 6; ++ei)
  {
    edges.push_back(edgesAry[ei]);
  }
  front->advanceFront(edges);
}

void deletePoints(AF2Point3D** & pointsAry, int pointsAryLength)
{
  for (int i = 0; i < pointsAryLength; ++i)
  {
    delete pointsAry[i];
  }
  delete[] pointsAry;
}

// test that attempting to select a neighborhood when there
// are no edges on the advancing front throws an exception
void testSelectEmpty()
{
  AF2LocalTransformMaker* transformBuilder = makeTransformBuilder();

  AF2Front* emptyFront = new AF2Front();
  CHECK(emptyFront->isEmpty());
  std::cout << "Good: At initial construction, the front is empty."
      << std::endl;

  bool exceptionThrown = false;
  try
  {
    AF2Neighborhood* shouldFail =
        emptyFront->selectNeighborhood(transformBuilder);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    std::cout << "Good: Attempting to select a neighborhood caused an error."
        << std::endl;
    exceptionThrown = true;
    CHECK_EQUAL(MeshKit::MK_FAILURE, mkError.error_code());
    std::cout << "Good: The error is a failure error." << std::endl;
  }
  CHECK(exceptionThrown);

  AF2Point3D* aPoint = new AF2Point3D(0, -0.5, 0.0, 0.5);
  aPoint->limitDistanceToBoundary(0);
  emptyFront->addPoint(aPoint);
  std::cout << "Added a point to the advancing front, but no edges."
      << std::endl;
  CHECK(!emptyFront->isEmpty());
  std::cout << "Good: The front is no longer empty." << std::endl;

  exceptionThrown = false;
  try
  {
    AF2Neighborhood* shouldFail =
        emptyFront->selectNeighborhood(transformBuilder);
    delete shouldFail;
  }
  catch (MeshKit::Error& mkError)
  {
    std::cout << "Good: Attempting to select a neighborhood caused an error."
        << std::endl;
    exceptionThrown = true;
    CHECK_EQUAL(MeshKit::MK_FAILURE, mkError.error_code());
    std::cout << "Good: The error is a failure error." << std::endl;
  }
  CHECK(exceptionThrown);

  delete aPoint;
  delete emptyFront;

  delete transformBuilder;

  std::cout << "PASS: Attempts to select a neighborhood when the front has\n"
      << "  no edges failed appropriately." << std::endl;
}

void testInitAndDestruct()
{
  AF2Front* hexFront = new AF2Front();
  AF2Point3D** pointsAry = NULL;
  AF2Edge3D** edgesAry = NULL;
  initHex(hexFront, pointsAry, edgesAry);
  CHECK(!hexFront->isEmpty());
  std::cout << "Good: The front is not empty after it has been initialized."
      << std::endl;
  CHECK_EQUAL(1u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality is 1 after the front has\n"
      << "  been initialized." << std::endl;
  delete hexFront;
  delete[] edgesAry;
  deletePoints(pointsAry, 6);
  std::cout << "PASS: The front is initialized and destructed correctly."
      << std::endl;
}

void testQualityDecrease()
{
  AF2Front* hexFront = new AF2Front();
  AF2Point3D** pointsAry = NULL;
  AF2Edge3D** edgesAry = NULL;
  initHex(hexFront, pointsAry, edgesAry);

  for (int ei = 0; ei < 6; ++ei)
  {
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
  }

  CHECK_EQUAL(3u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality is 3 after each edge has\n"
      << "  had its quality decreased twice." << std::endl;

  for (int ei = 0; ei < 4; ++ei)
  {
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
  }

  CHECK_EQUAL(3u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality is 3 after some but not all edges\n"
      << "  have their quality decreased three more times." << std::endl;

  for (int ei = 4; ei < 6; ++ei)
  {
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
  }

  CHECK_EQUAL(6u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality is 6 after some edges\n"
      << "  have their quality decreased five times and others\n"
      << "  have their quality decreased seven times." << std::endl;

  for (int ei = 0; ei < 3; ++ei)
  {
    edgesAry[ei]->decreaseQuality();
    edgesAry[ei]->decreaseQuality();
  }

  CHECK_EQUAL(6u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality should be and is 6." << std::endl;
  edgesAry[3]->decreaseQuality();
  CHECK_EQUAL(7u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality should be and is 7." << std::endl;
  edgesAry[3]->decreaseQuality();
  CHECK_EQUAL(8u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality should be and is 8." << std::endl;
  edgesAry[3]->decreaseQuality();
  CHECK_EQUAL(8u, hexFront->getMaximumQuality());
  std::cout << "Good: The maximum quality should be and is 8." << std::endl;

  delete hexFront;
  delete[] edgesAry;
  deletePoints(pointsAry, 6);

  std::cout << "PASS: Quality decreases are being correctly reported to\n"
      << "  and tracked by the front." << std::endl;
}

void testSelectNeighborhood()
{
  AF2Front* frontObj = new AF2Front();
  AF2Point3D** pointsAry = NULL;
  AF2Edge3D** edgesAry = NULL;
  initHex(frontObj, pointsAry, edgesAry);

  AF2LocalTransformMaker* transformBuilder = makeTransformBuilder();

  AF2Neighborhood* ngbhd = frontObj->selectNeighborhood(transformBuilder);
  std::list<const AF2Point2D*>::size_type sizeSix(6u);
  CHECK_EQUAL(sizeSix, ngbhd->getPoints2D()->size());
  std::cout << "Good: The selected neighborhood has 6 points." << std::endl;
  CHECK_EQUAL(sizeSix, ngbhd->getEdges2D()->size());
  std::cout << "Good: The selected neighborhood has 6 edges." << std::endl;
  delete ngbhd;

  for (int ei = 0; ei < 6; ++ei)
  {
    if (ei == 2)
    {
      continue;
    }
    edgesAry[ei]->decreaseQuality();
  }
  std::cout << "Decreased quality on all but one of the edges." << std::endl;

  ngbhd = frontObj->selectNeighborhood(transformBuilder);
  const AF2Edge2D* baseEdge2D = ngbhd->getBaselineEdge2D();
  CHECK_EQUAL(edgesAry[2]->getStart(),
      ngbhd->getCorrespondingPoint(baseEdge2D->getStart()));
  CHECK_EQUAL(edgesAry[2]->getEnd(),
      ngbhd->getCorrespondingPoint(baseEdge2D->getEnd()));
  std::cout << "Good: The baseline edge is the maximum quality edge."
      << std::endl;
  CHECK_EQUAL(sizeSix, ngbhd->getPoints2D()->size());
  std::cout << "Good: The selected neighborhood has 6 points." << std::endl;
  CHECK_EQUAL(sizeSix, ngbhd->getEdges2D()->size());
  std::cout << "Good: The selected neighborhood has 6 edges." << std::endl;

  delete ngbhd;
  delete frontObj;
  delete[] edgesAry;
  deletePoints(pointsAry, 6);
  delete transformBuilder;

  std::cout << "PASS: The basic case of selecting a neighborhood with a\n"
      << "  contiguous front of edges appears to work properly." << std::endl;
}

void testSelectIsolatedPoint()
{
  AF2Front* frontObj = new AF2Front();
  AF2Point3D** pointsAry = NULL;
  AF2Edge3D** edgesAry = NULL;
  initHex(frontObj, pointsAry, edgesAry);

  // add an isolated point
  AF2Point3D* isoPoint = new AF2Point3D(6, -0.5, 0.0625, 0.5);
  frontObj->addPoint(isoPoint);

  // decrease quality on all but one of the edges
  for (int ei = 1; ei < 6; ++ei)
  {
    edgesAry[ei]->decreaseQuality();
  }

  AF2LocalTransformMaker* transformBuilder = makeTransformBuilder();

  AF2Neighborhood* ngbhd = frontObj->selectNeighborhood(transformBuilder);
  const AF2Edge2D* baseEdge2D = ngbhd->getBaselineEdge2D();
  CHECK_EQUAL(edgesAry[0]->getStart(),
      ngbhd->getCorrespondingPoint(baseEdge2D->getStart()));
  CHECK_EQUAL(edgesAry[0]->getEnd(),
      ngbhd->getCorrespondingPoint(baseEdge2D->getEnd()));
  std::cout << "Good: The baseline edge is the maximum quality edge."
      << std::endl;
  std::list<const AF2Point2D*>::size_type sizeSeven(7u);
  CHECK_EQUAL(sizeSeven, ngbhd->getPoints2D()->size());
  std::cout << "Good: The selected neighborhood has 7 points." << std::endl;
  std::list<const AF2Point2D*>::size_type sizeSix(6u);
  CHECK_EQUAL(sizeSix, ngbhd->getEdges2D()->size());
  std::cout << "Good: The selected neighborhood has 6 edges." << std::endl;

  delete ngbhd;
  delete frontObj;
  delete[] edgesAry;
  delete isoPoint;
  deletePoints(pointsAry, 6);
  delete transformBuilder;

  std::cout << "PASS: The case of selecting a neighborhood with an\n"
      << "  isolated point in addition to edges appears to work properly."
      << std::endl;
}

void testSelectHourglass()
{
  AF2Front* frontObj = new AF2Front();
  AF2Point3D** pointsAry = NULL;
  AF2Edge3D** edgesAry = NULL;

  unsigned int numPoints = 54u;
  pointsAry = new AF2Point3D*[numPoints];
  pointsAry[0] = new AF2Point3D(0, -0.501, -0.0005, 0.5);
  pointsAry[1] = new AF2Point3D(1, -0.499, -0.0005, 0.5);
  pointsAry[2] = new AF2Point3D(2, -0.497, -0.0006, 0.5);
  pointsAry[3] = new AF2Point3D(3, -0.495, -0.0008, 0.5);
  pointsAry[4] = new AF2Point3D(4, -0.493, -0.0011, 0.5);
  pointsAry[5] = new AF2Point3D(5, -0.491, -0.0015, 0.5);
  pointsAry[6] = new AF2Point3D(6, -0.489, -0.0020, 0.5);
  pointsAry[7] = new AF2Point3D(7, -0.487, -0.0024, 0.5);
  pointsAry[8] = new AF2Point3D(8, -0.485, -0.0027, 0.5);
  pointsAry[9] = new AF2Point3D(9, -0.483, -0.0029, 0.5);
  pointsAry[10] = new AF2Point3D(10, -0.481, -0.0030, 0.5);
  pointsAry[11] = new AF2Point3D(11, -0.479, -0.0030, 0.5);
  pointsAry[12] = new AF2Point3D(12, -0.479, -0.0020, 0.5);
  pointsAry[13] = new AF2Point3D(13, -0.479, -0.0010, 0.5);
  pointsAry[14] = new AF2Point3D(14, -0.479, 0.0, 0.5);
  pointsAry[15] = new AF2Point3D(15, -0.479, 0.0010, 0.5);
  pointsAry[16] = new AF2Point3D(16, -0.479, 0.0020, 0.5);
  pointsAry[17] = new AF2Point3D(17, -0.479, 0.0030, 0.5);
  pointsAry[18] = new AF2Point3D(18, -0.481, 0.0030, 0.5);
  pointsAry[19] = new AF2Point3D(19, -0.483, 0.0029, 0.5);
  pointsAry[20] = new AF2Point3D(20, -0.485, 0.0027, 0.5);
  pointsAry[21] = new AF2Point3D(21, -0.487, 0.0024, 0.5);
  pointsAry[22] = new AF2Point3D(22, -0.489, 0.0020, 0.5);
  pointsAry[23] = new AF2Point3D(23, -0.491, 0.0015, 0.5);
  pointsAry[24] = new AF2Point3D(24, -0.493, 0.0011, 0.5);
  pointsAry[25] = new AF2Point3D(25, -0.495, 0.0008, 0.5);
  pointsAry[26] = new AF2Point3D(26, -0.497, 0.0006, 0.5);
  pointsAry[27] = new AF2Point3D(27, -0.499, 0.0005, 0.5);
  pointsAry[28] = new AF2Point3D(28, -0.501, 0.0005, 0.5);
  pointsAry[29] = new AF2Point3D(29, -0.503, 0.0006, 0.5);
  pointsAry[30] = new AF2Point3D(30, -0.505, 0.0008, 0.5);
  pointsAry[31] = new AF2Point3D(31, -0.507, 0.0011, 0.5);
  pointsAry[32] = new AF2Point3D(32, -0.509, 0.0015, 0.5);
  pointsAry[33] = new AF2Point3D(33, -0.511, 0.0020, 0.5);
  pointsAry[34] = new AF2Point3D(34, -0.513, 0.0024, 0.5);
  pointsAry[35] = new AF2Point3D(35, -0.515, 0.0027, 0.5);
  pointsAry[36] = new AF2Point3D(36, -0.517, 0.0029, 0.5);
  pointsAry[37] = new AF2Point3D(37, -0.519, 0.0030, 0.5);
  pointsAry[38] = new AF2Point3D(38, -0.521, 0.0030, 0.5);
  pointsAry[39] = new AF2Point3D(39, -0.521, 0.0020, 0.5);
  pointsAry[40] = new AF2Point3D(40, -0.521, 0.0010, 0.5);
  pointsAry[41] = new AF2Point3D(41, -0.521, 0.0, 0.5);
  pointsAry[42] = new AF2Point3D(42, -0.521, -0.0010, 0.5);
  pointsAry[43] = new AF2Point3D(43, -0.521, -0.0020, 0.5);
  pointsAry[44] = new AF2Point3D(44, -0.521, -0.0030, 0.5);
  pointsAry[45] = new AF2Point3D(45, -0.519, -0.0030, 0.5);
  pointsAry[46] = new AF2Point3D(46, -0.517, -0.0029, 0.5);
  pointsAry[47] = new AF2Point3D(47, -0.515, -0.0027, 0.5);
  pointsAry[48] = new AF2Point3D(48, -0.513, -0.0024, 0.5);
  pointsAry[49] = new AF2Point3D(49, -0.511, -0.0020, 0.5);
  pointsAry[50] = new AF2Point3D(50, -0.509, -0.0015, 0.5);
  pointsAry[51] = new AF2Point3D(51, -0.507, -0.0011, 0.5);
  pointsAry[52] = new AF2Point3D(52, -0.505, -0.0008, 0.5);
  pointsAry[53] = new AF2Point3D(53, -0.503, -0.0006, 0.5);

  for (unsigned int pi = 0; pi < numPoints; ++pi)
  {
    pointsAry[pi]->limitDistanceToBoundary(0u);
    frontObj->addPoint(pointsAry[pi]);
  }

  edgesAry = new AF2Edge3D*[numPoints];
  std::list<AF2Edge3D*> edges;
  for (unsigned int ei = 0; ei < numPoints; ++ei)
  {
    edgesAry[ei] =
        new AF2Edge3D(pointsAry[ei], pointsAry[(ei + 1) % numPoints]);
    edges.push_back(edgesAry[ei]);
  }

  frontObj->advanceFront(edges);

  // decrease quality on all but the first edge
  for (unsigned int ei = 1u; ei < numPoints; ++ei)
  {
    edgesAry[ei]->decreaseQuality();
  }

  AF2LocalTransformMaker* transformBuilder = makeTransformBuilder();

  AF2Neighborhood* ngbhd = frontObj->selectNeighborhood(transformBuilder);
  const AF2Edge2D* baseEdge2D = ngbhd->getBaselineEdge2D();
  CHECK_EQUAL(edgesAry[0]->getStart(),
      ngbhd->getCorrespondingPoint(baseEdge2D->getStart()));
  CHECK_EQUAL(edgesAry[0]->getEnd(),
      ngbhd->getCorrespondingPoint(baseEdge2D->getEnd()));
  std::cout << "Good: The baseline edge is the maximum quality edge."
      << std::endl;
  // The neighborhood should have two disconnected paths
  CHECK(ngbhd->getPoints2D()->size() == ngbhd->getEdges2D()->size() + 2);
  std::cout << "Good: The selected neighborhood has 2 more points than edges."
      << std::endl;
  // The neighborhood should exclude more than just the two caps
  CHECK(ngbhd->getEdges2D()->size() < 40u);
  std::cout << "Good: The selected neighborhood has fewer than 40 edges.\n"
      << "  (It has " << ngbhd->getEdges2D()->size() << " edges.)"
      << std::endl;

  delete ngbhd;
  delete frontObj;
  delete[] edgesAry;
  deletePoints(pointsAry, numPoints);
  delete transformBuilder;

  std::cout << "PASS: Neighborhood selection can select two disjoint paths."
      << std::endl;
}

void testAdvanceMissingPoints()
{
  AF2Point3D** pointsAry = new AF2Point3D*[6];
  pointsAry[0] = new AF2Point3D(0, -0.625, -0.25, 0.5);
  pointsAry[1] = new AF2Point3D(1, -0.375, -0.25, 0.5);
  pointsAry[2] = new AF2Point3D(2, -0.25, 0.0, 0.5);
  pointsAry[3] = new AF2Point3D(3, -0.375, 0.25, 0.5);
  pointsAry[4] = new AF2Point3D(4, -0.625, 0.25, 0.5);
  pointsAry[5] = new AF2Point3D(5, -0.75, 0.0, 0.5);

  // do not add the points to the front before attempting to 
  // advance the front and add the edges

  AF2Edge3D** edgesAry = new AF2Edge3D*[6];
  edgesAry[0] = new AF2Edge3D(pointsAry[0], pointsAry[1]);
  edgesAry[1] = new AF2Edge3D(pointsAry[1], pointsAry[2]);
  edgesAry[2] = new AF2Edge3D(pointsAry[2], pointsAry[3]);
  edgesAry[3] = new AF2Edge3D(pointsAry[3], pointsAry[4]);
  edgesAry[4] = new AF2Edge3D(pointsAry[4], pointsAry[5]);
  edgesAry[5] = new AF2Edge3D(pointsAry[5], pointsAry[0]);

  std::list<AF2Edge3D*> edges;
  for (int ei = 0; ei < 6; ++ei)
  {
    edges.push_back(edgesAry[ei]);
  }

  AF2Front* frontObj = new AF2Front();
  bool exceptionThrown = false;
  try
  {
    frontObj->advanceFront(edges);
  }
  catch (MeshKit::Error& mkError)
  {
    std::cout << "Good: Attempting to initialize/advance the front using\n"
        << "  edges whose endpoints have not been added causes an error."
        << std::endl;
    exceptionThrown = true;
    CHECK_EQUAL(MeshKit::MK_BAD_INPUT, mkError.error_code());
    std::cout << "Good: The error is a bad input error." << std::endl;
  }
  CHECK(exceptionThrown);

  // delete the edges (due to error, none of them should be owned by the front)
  for (int i = 0; i < 6; ++i)
  {
    delete edgesAry[i];
  }
  delete[] edgesAry;
  deletePoints(pointsAry, 6);
  delete frontObj;

  std::cout << "PASS: Advancing the front fails if there are missing points."
      << std::endl;
}

void testAdvanceInvalidQuality()
{
  AF2Front* frontObj = new AF2Front();

  AF2Point3D** pointsAry = new AF2Point3D*[6];
  pointsAry[0] = new AF2Point3D(0, -0.625, -0.25, 0.5);
  pointsAry[1] = new AF2Point3D(1, -0.375, -0.25, 0.5);
  pointsAry[2] = new AF2Point3D(2, -0.25, 0.0, 0.5);
  pointsAry[3] = new AF2Point3D(3, -0.375, 0.25, 0.5);
  pointsAry[4] = new AF2Point3D(4, -0.625, 0.25, 0.5);
  pointsAry[5] = new AF2Point3D(5, -0.75, 0.0, 0.5);

  for (int pi = 0; pi < 6; ++pi)
  {
    pointsAry[pi]->limitDistanceToBoundary(0u);
    frontObj->addPoint(pointsAry[pi]);
  }

  AF2Edge3D** edgesAry = new AF2Edge3D*[6];
  edgesAry[0] = new AF2Edge3D(pointsAry[0], pointsAry[1]);
  edgesAry[1] = new AF2Edge3D(pointsAry[1], pointsAry[2]);
  edgesAry[2] = new AF2Edge3D(pointsAry[2], pointsAry[3]);
  edgesAry[3] = new AF2Edge3D(pointsAry[3], pointsAry[4]);
  edgesAry[4] = new AF2Edge3D(pointsAry[4], pointsAry[5]);
  edgesAry[5] = new AF2Edge3D(pointsAry[5], pointsAry[0]);

  // decrease the quality of one of the edges before advancing the front
  edgesAry[4]->decreaseQuality();

  std::list<AF2Edge3D*> edges;
  for (int ei = 0; ei < 6; ++ei)
  {
    edges.push_back(edgesAry[ei]);
  }

  bool exceptionThrown = false;
  try
  {
    frontObj->advanceFront(edges);
  }
  catch (MeshKit::Error& mkError)
  {
    std::cout << "Good: Attempting to initialize/advance the front using\n"
        << "  edges whose quality has been decreased causes an error."
        << std::endl;
    exceptionThrown = true;
    CHECK_EQUAL(MeshKit::MK_BAD_INPUT, mkError.error_code());
    std::cout << "Good: The error is a bad input error." << std::endl;
  }
  CHECK(exceptionThrown);

  // delete the edges (due to error, none of them should be owned by the front)
  for (int i = 0; i < 6; ++i)
  {
    delete edgesAry[i];
  }
  delete[] edgesAry;
  deletePoints(pointsAry, 6);
  delete frontObj;

  std::cout << "PASS: Advancing the front fails if edges have decreased in\n"
      << "  quality before they are added to the front."
      << std::endl;
}

void testAdvanceFront()
{
  AF2Front* frontObj = new AF2Front();
  AF2Point3D** pointsAry = NULL;
  AF2Edge3D** edgesAry = NULL;
  initHex(frontObj, pointsAry, edgesAry);

  AF2Point3D* addedPoint = new AF2Point3D(6, -0.5, 0.0625, 0.5);
  frontObj->addPoint(addedPoint);

  std::list<AF2Edge3D*> addedEdges;
  addedEdges.push_back(new AF2Edge3D(pointsAry[0], addedPoint));
  addedEdges.push_back(new AF2Edge3D(addedPoint, pointsAry[1]));
  addedEdges.push_back(new AF2Edge3D(pointsAry[1], pointsAry[0]));
  frontObj->advanceFront(addedEdges);
  std::cout << "Good: The front has advanced without error." << std::endl;

  AF2LocalTransformMaker* transformBuilder = makeTransformBuilder();
  AF2Neighborhood* ngbhd = frontObj->selectNeighborhood(transformBuilder);
  std::list<const AF2Point2D*>::size_type sizeSeven(7u);
  CHECK_EQUAL(sizeSeven, ngbhd->getPoints2D()->size());
  std::cout << "Good: After advance, the selected neighborhood has 7 points."
      << std::endl;
  CHECK_EQUAL(sizeSeven, ngbhd->getEdges2D()->size());
  std::cout << "Good: After advance, the selected neighborhood has 7 edges."
      << std::endl;
  CHECK_EQUAL(1u, addedPoint->getDistanceToBoundary());
  std::cout << "Good: After advance, the distance from the added point\n"
      << "  to the boundary is 1." << std::endl;

  delete ngbhd;
  delete frontObj;
  delete[] edgesAry;
  delete addedPoint;
  deletePoints(pointsAry, 6);
  delete transformBuilder;
  std::cout << "PASS: The front appears to have been advanced correctly."
      << std::endl;
}

void testAdvanceHangingEdge()
{
  AF2Front* frontObj = new AF2Front();

  AF2Point3D** pointsAry = new AF2Point3D*[7];
  pointsAry[0] = new AF2Point3D(0, -0.625, -0.25, 0.5);
  pointsAry[1] = new AF2Point3D(1, -0.375, -0.25, 0.5);
  pointsAry[2] = new AF2Point3D(2, -0.25, 0.0, 0.5);
  pointsAry[3] = new AF2Point3D(3, -0.375, 0.25, 0.5);
  pointsAry[4] = new AF2Point3D(4, -0.625, 0.25, 0.5);
  pointsAry[5] = new AF2Point3D(5, -0.75, 0.0, 0.5);
  pointsAry[6] = new AF2Point3D(6, -0.5, 0.0625, 0.5);

  for (int pi = 0; pi < 7; ++pi)
  {
    pointsAry[pi]->limitDistanceToBoundary(0u);
    frontObj->addPoint(pointsAry[pi]);
  }

  AF2Edge3D** edgesAry = new AF2Edge3D*[8];
  // first edge
  edgesAry[0] = new AF2Edge3D(pointsAry[0], pointsAry[1]);
  // both half edges of a hanging edge
  edgesAry[1] = new AF2Edge3D(pointsAry[1], pointsAry[6]);
  edgesAry[2] = new AF2Edge3D(pointsAry[6], pointsAry[1]);
  // other edges
  edgesAry[3] = new AF2Edge3D(pointsAry[1], pointsAry[2]);
  edgesAry[4] = new AF2Edge3D(pointsAry[2], pointsAry[3]);
  edgesAry[5] = new AF2Edge3D(pointsAry[3], pointsAry[4]);
  edgesAry[6] = new AF2Edge3D(pointsAry[4], pointsAry[5]);
  edgesAry[7] = new AF2Edge3D(pointsAry[5], pointsAry[0]);

  std::list<AF2Edge3D*> edges;
  for (int ei = 0; ei < 8; ++ei)
  {
    edges.push_back(edgesAry[ei]);
  }
  frontObj->advanceFront(edges);

  AF2LocalTransformMaker* transformBuilder = makeTransformBuilder();
  AF2Neighborhood* ngbhd = frontObj->selectNeighborhood(transformBuilder);
  std::list<const AF2Point2D*>::size_type sizeSeven(7u);
  std::list<const AF2Point2D*>::size_type sizeEight(8u);
  CHECK_EQUAL(sizeSeven, ngbhd->getPoints2D()->size());
  std::cout << "Good: After init, the selected neighborhood has 7 points."
      << std::endl;
  CHECK_EQUAL(sizeEight, ngbhd->getEdges2D()->size());
  std::cout << "Good: After init, the selected neighborhood has 8 edges."
      << std::endl;
  delete ngbhd;

  std::list<AF2Edge3D*> addedEdges;
  addedEdges.push_back(new AF2Edge3D(pointsAry[0], pointsAry[6]));
  addedEdges.push_back(new AF2Edge3D(pointsAry[6], pointsAry[1]));
  addedEdges.push_back(new AF2Edge3D(pointsAry[1], pointsAry[0]));
  frontObj->advanceFront(addedEdges);
  std::cout << "Good: The front has advanced without error." << std::endl;

  ngbhd = frontObj->selectNeighborhood(transformBuilder);
  CHECK_EQUAL(sizeSeven, ngbhd->getPoints2D()->size());
  std::cout << "Good: After advance, the selected neighborhood has 7 points."
      << std::endl;
  CHECK_EQUAL(sizeSeven, ngbhd->getEdges2D()->size());
  std::cout << "Good: After advance, the selected neighborhood has 7 edges."
      << std::endl;

  delete ngbhd;
  delete frontObj;
  delete[] edgesAry;
  deletePoints(pointsAry, 7);
  delete transformBuilder;
  std::cout << "PASS: The front appears to have initialized and advanced\n"
      << "  correctly in the presence of a hanging edge."
      << std::endl;
}

void testCompletedAdvance()
{
  AF2Front* frontObj = new AF2Front();
  AF2Point3D** pointsAry = NULL;
  AF2Edge3D** edgesAry = NULL;
  initHex(frontObj, pointsAry, edgesAry);

  // add a point near the center of the hexagon
  AF2Point3D* addedPoint = new AF2Point3D(6, -0.5, 0.0625, 0.5);
  frontObj->addPoint(addedPoint);

  // fill in triangles with the added point one at a time
  // until all should be filled in
  std::list<AF2Edge3D*> addedEdges;
  for (unsigned int i = 0; i < 6u; ++i)
  {
    addedEdges.clear();
    addedEdges.push_back(new AF2Edge3D(pointsAry[i], addedPoint));
    addedEdges.push_back(new AF2Edge3D(addedPoint, pointsAry[(i + 1) % 6]));
    addedEdges.push_back(new AF2Edge3D(pointsAry[(i + 1) % 6], pointsAry[i]));
    frontObj->advanceFront(addedEdges);
  }

  CHECK(frontObj->isEmpty());
  std::cout << "Good: The front is empty after appropriate advances."
      << std::endl;

  delete frontObj;
  delete[] edgesAry;
  delete addedPoint;
  deletePoints(pointsAry, 6);
  std::cout << "PASS: The front is empty after appropriate advances."
      << std::endl;
}
