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
#include "meshkit/AF2Edge3D.hpp"
#include "meshkit/AF2Neighborhood.hpp"
#include "meshkit/AF2PlaneProjection.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/AF2PointTransformNone.hpp"
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
#define FILE_EXT "sat"
#endif

AF2Rule* makeFreeTriangleRule();
void testFreeTriangleRule();

class TestVisitor : public AF2RuleAppVisitor
{
  private:
    int count;
  public:
    TestVisitor() {
      count = 0;
    }

    void visit(AF2RuleApplication const & ruleApp)
    {
      ++count;
      std::cout << "Visit #" << count << std::endl;
      std::cout << "  Number of new faces: " << ruleApp.getNumNewFaces()
          << "\n  Number of new points: " << ruleApp.getNumNewPoints() 
          << std::endl;
    }
};

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
  AF2RuleNewEdge* newEdgePtr = new AF2RuleNewEdge(0, 2);
  newEdges.push_back(newEdgePtr);
  newEdgePtr = new AF2RuleNewEdge(2, 1);
  newEdges.push_back(newEdgePtr);

  // new face
  std::list<const AF2RuleNewFace*> newFaces;
  AF2RuleNewFace* newFacePtr = new AF2RuleNewTriangle(0, 1, 2);
  newFaces.push_back(newFacePtr);

  AF2Rule* rulePtr = new AF2Rule(exVertices, baseEdgePtr, exEdges,
      freeZoneDef, newVertices, newEdges, newFaces);

  return rulePtr;
}

void testFreeTriangleRule()
{
  AF2Rule* freeTriRule = makeFreeTriangleRule();

  CHECK(freeTriRule != NULL);

  MeshKit::MKCore* mk = new MeshKit::MKCore();

  // load a square in plane z = 0.5 with -1.0 <= x <= 0 and -0.5 <= y <= 0.5
  std::string file_name = TestDir + "/squaresurf." + FILE_EXT;
  mk->load_geometry_mesh(file_name.c_str(), file_name.c_str());

  MeshKit::MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);
  MeshKit::ModelEnt* square = (*surfs.begin());

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
      square->geom_handle(), origin, normal, xDir, 0.25);

  std::vector<AF2Point3D*> pointsVector;
  pointsVector.push_back(new AF2Point3D(-1.0, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(-0.75, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(-0.5, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(-0.25, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(0.0, -0.5, 0.5));
  pointsVector.push_back(new AF2Point3D(-1.0, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(-1.0, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(-1.0, 0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(0.0, -0.25, 0.5));
  pointsVector.push_back(new AF2Point3D(0.0, 0.0, 0.5));
  pointsVector.push_back(new AF2Point3D(0.0, 0.25, 0.5));

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
  freeTriRule->applyRule(ngbhd, 1, visitor);

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
}
