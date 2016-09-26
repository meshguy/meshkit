/**
 * \file dfltTriRules.cpp \test
 *
 * Test the default set of triangle rules for two-dimensional advancing front.
 *
 * Attempts to apply each rule in the default set of triangle rules to its
 * reference case, testing that it can be applied in that case.
 *
 * This test is a bit non-standard in that it depends on knowledge of the
 * implementation details rather than merely the API of AF2DfltTriangleRules.
 * The purpose, though, is to make sure that the implementation details
 * match expectations.
 */

// C++
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

// MeshKit
#include "meshkit/AF2DfltTriangleRules.hpp"
#include "meshkit/AF2Edge3D.hpp"
#include "meshkit/AF2Neighborhood.hpp"
#include "meshkit/AF2PlaneProjection.hpp"
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/AF2Rule.hpp"
#include "meshkit/AF2RuleAppVisitor.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"

// MeshKit testing
#include "TestUtil.hpp"

// define the geometry file extension depending on the geometry model
#if HAVE_OCC
#define FILE_EXT "stp"
#else
#define FILE_EXT "sat"
#endif

const AF2Rule* findRule(
    const AF2DfltTriangleRules & ruleSet, const std::string & ruleName);
AF2LocalTransform* makePlaneProj();
void test180DegreeRuleQ1();
void test180DegreeRuleQ5();
void test180DegreeRuleQ10();
void test180DegreeRuleQ20();
void test60DegreeAngleRightRule();
void test60DegreeAngleLeftRule();
void test120DegreeAngleRightRule();
void test120DegreeAngleLeftRule();
void test120DegreeAngleBothRule();
void testFillTriangleRule();
void testConnectVertexRule();
void testConnectEdgeRule();

class SaveLastVisitor : public AF2RuleAppVisitor
{
  private:
    unsigned int count;
    AF2RuleApplication* lastRuleApp;
  public:
    SaveLastVisitor()
    {
      count = 0u;
      lastRuleApp = NULL;
    }

    ~SaveLastVisitor()
    {
      if (lastRuleApp != NULL)
      {
        delete lastRuleApp;
      }
    }

    void visit(AF2RuleApplication const & ruleApp)
    {
      ++count;
      if (lastRuleApp != NULL)
      {
        delete lastRuleApp;
      }
      lastRuleApp = new AF2RuleApplication(ruleApp);
    }

    AF2RuleApplication* getLast() const
    {
      return lastRuleApp;
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

  num_fail += RUN_TEST(test180DegreeRuleQ1);
  num_fail += RUN_TEST(test180DegreeRuleQ5);
  num_fail += RUN_TEST(test180DegreeRuleQ10);
  num_fail += RUN_TEST(test180DegreeRuleQ20);
  num_fail += RUN_TEST(test60DegreeAngleRightRule);
  num_fail += RUN_TEST(test60DegreeAngleLeftRule);
  num_fail += RUN_TEST(test120DegreeAngleRightRule);
  num_fail += RUN_TEST(test120DegreeAngleLeftRule);
  num_fail += RUN_TEST(test120DegreeAngleBothRule);
  num_fail += RUN_TEST(testFillTriangleRule);
  num_fail += RUN_TEST(testConnectVertexRule);
  num_fail += RUN_TEST(testConnectEdgeRule);

  delete mk;

  return num_fail;
}

const AF2Rule* findRule(
    const AF2DfltTriangleRules & ruleSet, const std::string & ruleName)
{
  typedef std::list<const AF2Rule*>::const_iterator RuleListItr;

  std::list<const AF2Rule*> ruleList = ruleSet.getRules();
  for (RuleListItr itr = ruleList.begin(); itr != ruleList.end(); ++itr)
  {
    if (ruleName == (*itr)->getName())
    {
      return *itr;
    }
  }

  return NULL;
}

AF2LocalTransform* makePlaneProj()
{
  MeshKit::Vector<3> origin;
  origin[0] = -0.5;
  origin[1] = 0.0;
  origin[2] = 0.5;
  MeshKit::Vector<3> normal;
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 1.0;
  MeshKit::Vector<3> xDir;
  xDir[0] = 1.0;
  xDir[1] = 0.0;
  xDir[2] = 0.0;

  return new AF2PlaneProjection(square->igeom_instance(),
      square->geom_handle(), origin, normal, xDir, 0.125);
}

void test180DegreeRuleQ1()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule180DegreeQ1 =
      findRule(triRules, "Triangle Off Line, Quality Level 1");

  CHECK(rule180DegreeQ1 != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.3123, 0.0876, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.4375, 0.1877, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.5627, 0.0876, 0.5));
  pointsVector.push_back(new AF2Point3D(8, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic line
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[5]));
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));
  edges.push_back(new AF2Edge3D(pointsVector[7], pointsVector[8]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule180DegreeQ1->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule180DegreeQ1->getName() << "\"" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.866, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.866)" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(9, -0.3752, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ1->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(10, -0.3127, 0.0874, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ1->applyRule(failNgbhdBravo, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3127, 0.0874, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(11, -0.4375, 0.1873, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ1->applyRule(failNgbhdCharlie, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4375, 0.1873, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(12, -0.5623, 0.0874, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ1->applyRule(failNgbhdDelta, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5623, 0.0874, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointEcho(13, -0.4998, 0.0001, 0.5);
  points.push_back(&failPointEcho);
  AF2Neighborhood failNgbhdEcho(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ1->applyRule(failNgbhdEcho, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.0001, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 180DegreeRuleQ1 passes the test." << std::endl;
}

void test180DegreeRuleQ5()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule180DegreeQ5 =
      findRule(triRules, "Triangle Off Line, Quality Level 5");

  CHECK(rule180DegreeQ5 != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.4248, 0.0676, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.4502, 0.0676, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic line
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[5]));
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule180DegreeQ5->applyRule(ngbhd, 5u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule180DegreeQ5->getName() << "\"" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.5)" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(8, -0.3752, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ5->applyRule(failNgbhdAlpha, 5u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(9, -0.4252, 0.0674, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ5->applyRule(failNgbhdBravo, 5u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4252, 0.0674, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(10, -0.4498, 0.0674, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ5->applyRule(failNgbhdCharlie, 5u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4498, 0.0674, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(11, -0.4998, 0.0001, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ5->applyRule(failNgbhdDelta, 5u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.0001, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 180DegreeRuleQ5 passes the test." << std::endl;
}

void test180DegreeRuleQ10()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule180DegreeQ10 =
      findRule(triRules, "Triangle Off Line, Quality Level 10");

  CHECK(rule180DegreeQ10 != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.43105, 0.0401, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.44395, 0.0401, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic line
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[5]));
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule180DegreeQ10->applyRule(ngbhd, 10u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule180DegreeQ10->getName() << "\"" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.3, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.3)" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(8, -0.3752, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ10->applyRule(failNgbhdAlpha, 10u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(9, -0.43145, 0.0399, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ10->applyRule(failNgbhdBravo, 10u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.43145, 0.0399, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(10, -0.44355, 0.0399, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ10->applyRule(failNgbhdCharlie, 10u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.44355, 0.0399, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(11, -0.4998, 0.0001, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ10->applyRule(failNgbhdDelta, 10u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.0001, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 180DegreeRuleQ10 passes the test." << std::endl;
}

void test180DegreeRuleQ20()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule180DegreeQ20 =
      findRule(triRules, "Triangle Off Line, Quality Level 20");

  CHECK(rule180DegreeQ20 != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.434175, 0.013225, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.440825, 0.013225, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic line
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[5]));
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule180DegreeQ20->applyRule(ngbhd, 20u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule180DegreeQ20->getName() << "\"" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.1, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.1)" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(8, -0.3756, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ20->applyRule(failNgbhdAlpha, 20u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3756, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(9, -0.434575, 0.013025, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule180DegreeQ20->applyRule(failNgbhdBravo, 20u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.434575, 0.013025, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(10, -0.440425, 0.013025, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ20->applyRule(failNgbhdCharlie, 20u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.440425, 0.013025, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(11, -0.4994, 0.0001, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  rule180DegreeQ20->applyRule(failNgbhdDelta, 20u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4994, 0.0001, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 180DegreeRuleQ20 passes the test." << std::endl;
}

void test60DegreeAngleRightRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule60DegreeRight =
      findRule(triRules, "Angle at Right Vertex of 60 Degrees");

  CHECK(rule60DegreeRight != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.4375, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.3125, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.4375, 0.10845, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.515725, 0.0812875, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  edges.push_back(new AF2Edge3D(pointsVector[3], pointsVector[4]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule60DegreeRight->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule60DegreeRight->getName() << "\"" << std::endl;
  CHECK_EQUAL(0u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There are no new points" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(8, -0.4375, 0.10805, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule60DegreeRight->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4375, 0.10805, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(9, -0.515525, 0.0810875, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule60DegreeRight->applyRule(failNgbhdBravo, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.515525, 0.0810875, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(10, -0.4998, 0.0001, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule60DegreeRight->applyRule(failNgbhdCharlie, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.0001, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 60DegreeAngleRightRule passes the test." << std::endl;
}

void test60DegreeAngleLeftRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule60DegreeLeft =
      findRule(triRules, "Angle at Left Vertex of 60 Degrees");

  CHECK(rule60DegreeLeft != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.5625, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.4375, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.359275, 0.0812875, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.4375, 0.10845, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[2], pointsVector[3]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[1], pointsVector[2]));
  edges.push_back(new AF2Edge3D(pointsVector[3], pointsVector[4]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule60DegreeLeft->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule60DegreeLeft->getName() << "\"" << std::endl;
  CHECK_EQUAL(0u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There are no new points" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(8, -0.3752, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule60DegreeLeft->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(9, -0.359475, 0.0810875, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule60DegreeLeft->applyRule(failNgbhdBravo, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.359475, 0.0810875, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(10, -0.4375, 0.10805, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule60DegreeLeft->applyRule(failNgbhdCharlie, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4375, 0.10805, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 60DegreeAngleLeftRule passes the test." << std::endl;
}

void test120DegreeAngleRightRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule120DegreeRight =
      findRule(triRules, "Angle at Right Vertex of 120 Degrees");

  CHECK(rule120DegreeRight != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.3125, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.1875, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.3124, 0.10835, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.375, 0.2167, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.5, 0.2167, 0.5));
  pointsVector.push_back(new AF2Point3D(8, -0.5626, 0.10835, 0.5));
  pointsVector.push_back(new AF2Point3D(9, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  edges.push_back(new AF2Edge3D(pointsVector[3], pointsVector[4]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));
  edges.push_back(new AF2Edge3D(pointsVector[7], pointsVector[8]));
  edges.push_back(new AF2Edge3D(pointsVector[8], pointsVector[9]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule120DegreeRight->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule120DegreeRight->getName() << "\"" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There is one new point" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.866, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.866)" << std::endl;
  CHECK_EQUAL(2u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There are two new faces" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(10, -0.3126, 0.10815, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule120DegreeRight->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3126, 0.10815, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(11, -0.375, 0.2163, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule120DegreeRight->applyRule(failNgbhdBravo, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.375, 0.2163, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(12, -0.5, 0.2163, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeRight->applyRule(failNgbhdCharlie, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5, 0.2163, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(13, -0.5624, 0.10815, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeRight->applyRule(failNgbhdDelta, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5624, 0.10815, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointEcho(14, -0.4998, 0.0001, 0.5);
  points.push_back(&failPointEcho);
  AF2Neighborhood failNgbhdEcho(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeRight->applyRule(failNgbhdEcho, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.0001, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 120DegreeAngleRightRule passes the test."
      << std::endl;
}

void test120DegreeAngleLeftRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule120DegreeLeft =
      findRule(triRules, "Angle at Left Vertex of 120 Degrees");

  CHECK(rule120DegreeLeft != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.6875, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5625, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.3124, 0.10835, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.375, 0.2167, 0.5));
  pointsVector.push_back(new AF2Point3D(8, -0.5, 0.2167, 0.5));
  pointsVector.push_back(new AF2Point3D(9, -0.5626, 0.10835, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[2], pointsVector[3]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[1], pointsVector[2]));
  edges.push_back(new AF2Edge3D(pointsVector[3], pointsVector[4]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));
  edges.push_back(new AF2Edge3D(pointsVector[7], pointsVector[8]));
  edges.push_back(new AF2Edge3D(pointsVector[8], pointsVector[9]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule120DegreeLeft->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule120DegreeLeft->getName() << "\"" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There is one new point" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.866, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.866)" << std::endl;
  CHECK_EQUAL(2u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There are two new faces" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(10, -0.3752, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule120DegreeLeft->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(11, -0.3126, 0.10815, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule120DegreeLeft->applyRule(failNgbhdBravo, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3126, 0.10815, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(12, -0.375, 0.2163, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeLeft->applyRule(failNgbhdCharlie, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.375, 0.2163, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(13, -0.5, 0.2163, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeLeft->applyRule(failNgbhdDelta, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5, 0.2163, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointEcho(14, -0.5624, 0.10815, 0.5);
  points.push_back(&failPointEcho);
  AF2Neighborhood failNgbhdEcho(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeLeft->applyRule(failNgbhdEcho, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5624, 0.10815, 0.5)" << std::endl;
  points.pop_back();


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

  std::cout << "PASS: The 120DegreeAngleLeftRule passes the test."
      << std::endl;
}

void test120DegreeAngleBothRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* rule120DegreeBoth =
      findRule(triRules, "Angles at Both Vertices of 120 Degrees");

  CHECK(rule120DegreeBoth != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.6875, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5625, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.3125, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.1875, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.3124, 0.10835, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.375, 0.162575, 0.5));
  pointsVector.push_back(new AF2Point3D(8, -0.5, 0.162575, 0.5));
  pointsVector.push_back(new AF2Point3D(9, -0.5626, 0.10835, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[2], pointsVector[3]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[1], pointsVector[2]));
  edges.push_back(new AF2Edge3D(pointsVector[3], pointsVector[4]));
  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[5]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));
  edges.push_back(new AF2Edge3D(pointsVector[7], pointsVector[8]));
  edges.push_back(new AF2Edge3D(pointsVector[8], pointsVector[9]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  rule120DegreeBoth->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << rule120DegreeBoth->getName() << "\"" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There is one new point" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.866, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.866)" << std::endl;
  CHECK_EQUAL(3u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There are three new faces" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(10, -0.3126, 0.10815, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  rule120DegreeBoth->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3126, 0.10815, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(11, -0.375, 0.162175, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  rule120DegreeBoth->applyRule(failNgbhdBravo, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.375, 0.162175, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(12, -0.5, 0.162175, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeBoth->applyRule(failNgbhdCharlie, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5, 0.162175, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(13, -0.5624, 0.10815, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  rule120DegreeBoth->applyRule(failNgbhdDelta, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5624, 0.10815, 0.5)" << std::endl;
  points.pop_back();

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

  std::cout << "PASS: The 120DegreeAngleBothRule passes the test."
      << std::endl;
}

void testFillTriangleRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* fillTriangleRule =
      findRule(triRules, "Fill Triangle");

  CHECK(fillTriangleRule != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.4375, 0.10825, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[0], pointsVector[1]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[1], pointsVector[2]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[0]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  fillTriangleRule->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << fillTriangleRule->getName() << "\"" << std::endl;
  CHECK_EQUAL(0u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There are no new points" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(3, -0.4375, 0.10815, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  fillTriangleRule->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4375, 0.10815, 0.5)" << std::endl;
  points.pop_back();

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

  std::cout << "PASS: The fillTriangleRule passes the test."
      << std::endl;
}

void testConnectVertexRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* connectVertexRule =
      findRule(triRules, "Connect to Opposite Vertex");

  CHECK(connectVertexRule != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.4375, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.3499, 0.086725, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.4375, 0.10845, 0.5));
  pointsVector.push_back(new AF2Point3D(8, -0.5251, 0.086725, 0.5));
  pointsVector.push_back(new AF2Point3D(9, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[5], pointsVector[6]));
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));
  edges.push_back(new AF2Edge3D(pointsVector[7], pointsVector[8]));
  edges.push_back(new AF2Edge3D(pointsVector[8], pointsVector[9]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  connectVertexRule->applyRule(ngbhd, 1u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << connectVertexRule->getName() << "\"" << std::endl;
  CHECK_EQUAL(0u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There are no new points" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There is one new face" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(10, -0.3752, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  connectVertexRule->applyRule(failNgbhdAlpha, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(11, -0.3501, 0.086525, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  connectVertexRule->applyRule(failNgbhdBravo, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3501, 0.086525, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(12, -0.5249, 0.086525, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  connectVertexRule->applyRule(failNgbhdCharlie, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5249, 0.086525, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(13, -0.4998, 0.0001, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  connectVertexRule->applyRule(failNgbhdDelta, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  // This vertex of the free zone cannot be tested by adding a point at
  //   (-0.4375, 0.10805) because the rule will just choose the vertex at
  //   (-0.4375, 0.10805) as the one it should connect to
  pointsVector.push_back(new AF2Point3D(14, -0.125, 0.10805, 0.5));
  pointsVector.push_back(new AF2Point3D(15, -0.875, 0.10805, 0.5));
  points.push_back(pointsVector[10]);
  points.push_back(pointsVector[11]);
  edges.push_back(new AF2Edge3D(pointsVector[10], pointsVector[11]));
  AF2Neighborhood failNgbhdEcho(points,
      baselineEdge, edges, makePlaneProj());
  connectVertexRule->applyRule(failNgbhdEcho, 1u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has an edge from (-0.125, 0.10805) to (-0.875, 0.10805, 0.5)"
      << std::endl;
  points.pop_back();

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

  std::cout << "PASS: The connectVertexRule passes the test." << std::endl;
}

void testConnectEdgeRule()
{
  AF2DfltTriangleRules triRules;
  const AF2Rule* connectEdgeRule =
      findRule(triRules, "Connect to Opposite Edge");

  CHECK(connectEdgeRule != NULL);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(0, -0.625, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(1, -0.5, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(2, -0.375, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(3, -0.25, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(4, -0.375, 0.2165, 0.5));
  pointsVector.push_back(new AF2Point3D(5, -0.5, 0.2165, 0.5));
  pointsVector.push_back(new AF2Point3D(6, -0.3748, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(7, -0.3123, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(8, -0.3748, 0.2165, 0.5));
  pointsVector.push_back(new AF2Point3D(9, -0.5002, 0.2165, 0.5));
  pointsVector.push_back(new AF2Point3D(10, -0.5627, 0.10825, 0.5));
  pointsVector.push_back(new AF2Point3D(11, -0.5002, 0.0, 0.5));

  std::list<AF2Point3D*> points;

  for (unsigned int pi = 0; pi < pointsVector.size(); ++pi)
  {
    points.push_back(pointsVector[pi]);
  }

  AF2Edge3D* baselineEdge = new AF2Edge3D(pointsVector[1], pointsVector[2]);

  std::list<const AF2Edge3D*> edges;
  // basic neighborhood edges
  edges.push_back(new AF2Edge3D(pointsVector[0], pointsVector[1]));
  edges.push_back(new AF2Edge3D(pointsVector[2], pointsVector[3]));
  edges.push_back(new AF2Edge3D(pointsVector[4], pointsVector[5]));
  // lines bounding the expected free zone location
  edges.push_back(new AF2Edge3D(pointsVector[6], pointsVector[7]));
  edges.push_back(new AF2Edge3D(pointsVector[7], pointsVector[8]));
  edges.push_back(new AF2Edge3D(pointsVector[9], pointsVector[10]));
  edges.push_back(new AF2Edge3D(pointsVector[10], pointsVector[11]));

  AF2Neighborhood ngbhd(points, baselineEdge, edges, makePlaneProj());

  SaveLastVisitor visitor;
  connectEdgeRule->applyRule(ngbhd, 3u, visitor);

  CHECK_EQUAL(1u, visitor.getVisitCount());
  std::cout << "Good: There is one way to apply the default triangle rule\n"
      << "  named \"" << connectEdgeRule->getName() << "\"" << std::endl;
  CHECK_EQUAL(1u, visitor.getLast()->getNumNewPoints());
  std::cout << "Good: There is one new point" << std::endl;
  CHECK_REAL_EQUAL(0.5, visitor.getLast()->getNewPoint(0u)->getX(), 1e-14);
  CHECK_REAL_EQUAL(0.866, visitor.getLast()->getNewPoint(0u)->getY(), 1e-14);
  std::cout << "Good: The new point is located near (0.5, 0.866)" << std::endl;
  CHECK_EQUAL(2u, visitor.getLast()->getNumNewFaces());
  std::cout << "Good: There are two new faces" << std::endl;

  SaveLastVisitor failVisitor;
  AF2Point3D failPointAlpha(12, -0.3752, 0.0001, 0.5);
  points.push_back(&failPointAlpha);
  AF2Neighborhood failNgbhdAlpha(points, baselineEdge, edges, makePlaneProj());
  connectEdgeRule->applyRule(failNgbhdAlpha, 3u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.0001, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointBravo(13, -0.3127, 0.10825, 0.5);
  points.push_back(&failPointBravo);
  AF2Neighborhood failNgbhdBravo(points, baselineEdge, edges, makePlaneProj());
  connectEdgeRule->applyRule(failNgbhdBravo, 3u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3127, 0.10825, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointCharlie(14, -0.3752, 0.2164, 0.5);
  points.push_back(&failPointCharlie);
  AF2Neighborhood failNgbhdCharlie(points,
      baselineEdge, edges, makePlaneProj());
  connectEdgeRule->applyRule(failNgbhdCharlie, 3u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.3752, 0.2164, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointDelta(15, -0.4998, 0.2164, 0.5);
  points.push_back(&failPointDelta);
  AF2Neighborhood failNgbhdDelta(points,
      baselineEdge, edges, makePlaneProj());
  connectEdgeRule->applyRule(failNgbhdDelta, 3u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.2164, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointEcho(16, -0.5623, 0.10825, 0.5);
  points.push_back(&failPointEcho);
  AF2Neighborhood failNgbhdEcho(points,
      baselineEdge, edges, makePlaneProj());
  connectEdgeRule->applyRule(failNgbhdEcho, 3u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.5623, 0.10825, 0.5)" << std::endl;
  points.pop_back();

  AF2Point3D failPointFoxtrot(17, -0.4998, 0.0001, 0.5);
  points.push_back(&failPointFoxtrot);
  AF2Neighborhood failNgbhdFoxtrot(points,
      baselineEdge, edges, makePlaneProj());
  connectEdgeRule->applyRule(failNgbhdFoxtrot, 3u, failVisitor);
  CHECK_EQUAL(0u, failVisitor.getVisitCount());
  std::cout << "Good: The rule fails to apply when the neighborhood\n"
      << "  has a point at (-0.4998, 0.0001, 0.5)" << std::endl;
  points.pop_back();

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

  std::cout << "PASS: The connectEdgeRule passes the test." << std::endl;
}
