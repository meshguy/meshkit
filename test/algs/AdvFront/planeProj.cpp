/**
 * \file planeProj.cpp \test
 *
 * Test the plane projection.
 *
 * Loads a geometry file representing a piece of a torus, and tests that
 * points are correctly projected from the torus onto various planes
 * and back from the plane onto the piece of the torus.
 */

#include <iostream>

#include "meshkit/AF2PlaneProjection.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

#if HAVE_OCC
#define FILE_EXT "stp"
#else
#define FILE_EXT "sat"
#endif

MKCore *mk = NULL;

void test_plane_projection();
bool testPlaneProj();

int main(int argc, char **argv)
{

  // start up MK and load the geometry
  int num_fail = 0;

  num_fail += RUN_TEST(test_plane_projection);

  return num_fail;

}

void test_plane_projection()
{
  CHECK(testPlaneProj());
}

bool testPlaneProj()
{
  int failCount = 0;

  mk = new MKCore();

  std::string file_name = TestDir + "/pieceOfTorus01." + FILE_EXT;
  mk->load_geometry_mesh(file_name.c_str(), file_name.c_str());

  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);
  ModelEnt *this_surf = (*surfs.begin());

  MeshKit::Vector<3> origin;
  // point on surface is at approx (2.82842712474619, 2.82842712474619, 0)
  // = (4*sin(pi/4), 4*sin(pi/4), 0) = (4*sqrt(2)/2, 4*sqrt(2)/2, 0)
  origin[0] = 2.8;
  origin[1] = 2.8;
  origin[2] = 0.0;
  MeshKit::Vector<3> normal;
  normal[0] = 0.7071067811865476; // sqrt(2)/2
  normal[1] = 0.7071067811865476;
  normal[2] = 0.0;
  MeshKit::Vector<3> xDir;
  xDir[0] = -0.7071067811865476;
  xDir[1] = 0.7071067811865476;
  xDir[2] = 0.0;
  AF2PlaneProjection planeProj(this_surf->igeom_instance(),
      this_surf->geom_handle(), origin, normal, xDir, 1.0);
  std::cout << "Created plane passing through (2.8, 2.8, 0) "
      << "normal to (1, 1, 0)." << std::endl;
  std::cout << "The y-coordinate in the plane is in the direction of "
      << "the z-axis." << std::endl;


  // 2.8059735404411 is (sqrt(2)/2) * (3 + cos(asin(0.25)))
  AF2Point3D testPnt(0, 2.8059735404411, 2.8059735404411, 0.25);
  AF2Point2D* planePnt = planeProj.transformFromSurface(testPnt);
  std::cout << "Transform (" << testPnt.getX() << ", " << testPnt.getY() << ", "
      << testPnt.getZ() << ") from surface to plane." << std::endl;
  std::cout << "Coordinates in plane: (" << planePnt->getX() << ", "
      << planePnt->getY() << ")" << std::endl;
  if ((fabs(planePnt->getX() - 0.0) + fabs(planePnt->getY() - 0.25)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates in the plane should be (0, 0.25)."
        << std::endl;
    ++failCount;
  }
  delete planePnt;

  // 0.1539268506556547 is sin(acos(2.82 * sqrt(2) - 3))
  AF2Point2D planeTestPntA(0, 0.1539268506556547);
  AF2Point3D* surfacePnt = planeProj.transformToSurface(planeTestPntA, 1);
  std::cout << "Transform (" << planeTestPntA.getX() << ", "
      << planeTestPntA.getY() << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt->getX() << ", "
      << surfacePnt->getY() << ", " << surfacePnt->getZ() << ")" << std::endl;
  if ((fabs(surfacePnt->getX() - 2.82) + fabs(surfacePnt->getY() - 2.82) +
      fabs(surfacePnt->getZ() - 0.15392685065547)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (2.82, 2.82, 0.15392685065547)."
        << std::endl;
    ++failCount;
  }
  delete surfacePnt;

  // 0.574682269122808 is sin(acos(2.7 * sqrt(2) - 3))
  AF2Point2D planeTestPntB(0, 0.574682269122808);
  surfacePnt = planeProj.transformToSurface(planeTestPntB, 2);
  std::cout << "Transform (" << planeTestPntB.getX() << ", "
      << planeTestPntB.getY() << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt->getX() << ", "
      << surfacePnt->getY() << ", " << surfacePnt->getZ() << ")" << std::endl;
  if ((fabs(surfacePnt->getX() - 2.7) + fabs(surfacePnt->getY() - 2.7) +
      fabs(surfacePnt->getZ() - 0.574682269122808)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (2.7, 2.7, 0.574682269122808)."
        << std::endl;
    ++failCount;
  }
  delete surfacePnt;


  // create a plane projection (with a normal that does not agree with the
  // normal of any plane tangent to the surface) such that there will
  // be points on the surface on both sides of the plane
  origin[0] = 3.5;
  origin[1] = 0.0;
  origin[2] = 0.0;
  normal[0] = -0.7071067811865476; // -sqrt(2)/2
  normal[1] = 0.0;
  normal[2] = 0.7071067811865476;
  xDir[0] = 0.7071067811865476;
  xDir[1] = 0.0;
  xDir[2] = 0.7071067811865476;
  AF2PlaneProjection wackyProj(this_surf->igeom_instance(),
      this_surf->geom_handle(), origin, normal, xDir, 1.0);
  std::cout << "Created plane passing through (3.5, 0, 0) "
      << "normal to (-1, 0, 1)." << std::endl;
  std::cout << "The y-coordinate in the plane is in the direction of "
      << "the y-axis." << std::endl;

  // test projecting a point from this plane onto the surface
  // given x = 3.95, y = 0.1,
  // z = sqrt(-8 + 6 * sqrt (x * x + y * y) - x * x - y * y)
  // ~ 0.3083726968400695 is the coordinate of the point on the surface
  // this projects to roughly (3.87918634842, 0.1, 0.37918634842) as a
  // 3-dimensional pooint on the plane, computing z-coordinate there as
  // z + (0.45 - z)/2.  Thus 2-dimensional planar coordinates are as
  // follows, with the first coordinat as 3-d z-coord divided by sqrt(2)/2
  AF2Point2D planeTestPntC(0.53625047660234298437463798744018, 0.1);
  surfacePnt = wackyProj.transformToSurface(planeTestPntC, 3);
  std::cout << "Transform (" << planeTestPntC.getX() << ", "
      << planeTestPntC.getY() << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt->getX() << ", "
      << surfacePnt->getY() << ", " << surfacePnt->getZ() << ")" << std::endl;
  // expect ~ (3.95, 0.1, 0.3083726968400695)
  if ((fabs(surfacePnt->getX() - 3.95) + fabs(surfacePnt->getY() - 0.1) +
      fabs(surfacePnt->getZ() - 0.3083726968400695)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (3.95, 0.1, 0.3083726968400695)."
        << std::endl;
    ++failCount;
  }
  // The line also intersects the surface near (3.30723, 0.1, 0.951148),
  // but that point is farther away from the plane than the expected result.
  delete surfacePnt;


  // Perform a similar test with the plane origin at x = 2.6
  origin[0] = 2.6;
  AF2PlaneProjection wackyProjShift01(this_surf->igeom_instance(),
      this_surf->geom_handle(), origin, normal, xDir, 1.0);
  std::cout << "Created plane passing through (2.6, 0, 0) "
      << "normal to (-1, 0, 1)." << std::endl;
  std::cout << "The y-coordinate in the plane is in the direction of "
      << "the y-axis." << std::endl;

  AF2Point2D planeTestPntD(1.17264657967023277, 0.1);
  surfacePnt = wackyProjShift01.transformToSurface(planeTestPntD, 4);
  std::cout << "Transform (" << planeTestPntD.getX() << ", "
      << planeTestPntD.getY() << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt->getX() << ", "
      << surfacePnt->getY() << ", " << surfacePnt->getZ() << ")" << std::endl;
  if ((fabs(surfacePnt->getX() - 3.3072251285719) +
      fabs(surfacePnt->getY() - 0.1) +
      fabs(surfacePnt->getZ() - 0.9511475682682)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (3.30722512857195, 0.1, 0.9511475682682)."
        << std::endl;
    ++failCount;
  }
  // The line also intersects the surface near (3.95, 0.1, 0.3083726968400695)
  // but that point is farther away from the plane than the expected result.
  delete surfacePnt;


  // Perform a similar test with the plane origin at x = 2.1
  origin[0] = 2.1;
  AF2PlaneProjection wackyProjShift02(this_surf->igeom_instance(),
      this_surf->geom_handle(), origin, normal, xDir, 1.0);
  std::cout << "Created plane passing through (2.1, 0, 0) "
      << "normal to (-1, 0, 1)." << std::endl;
  std::cout << "The y-coordinate in the plane is in the direction of "
      << "the y-axis." << std::endl;

  AF2Point2D planeTestPntE(1.5261999702635065, 0.1);
  surfacePnt = wackyProjShift02.transformToSurface(planeTestPntE, 5);
  std::cout << "Transform (" << planeTestPntE.getX() << ", "
      << planeTestPntE.getY() << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt->getX() << ", "
      << surfacePnt->getY() << ", " << surfacePnt->getZ() << ")" << std::endl;
  if ((fabs(surfacePnt->getX() - 3.3072251285719) +
      fabs(surfacePnt->getY() - 0.1) +
      fabs(surfacePnt->getZ() - 0.9511475682682)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (3.30722512857195, 0.1, 0.9511475682682)."
        << std::endl;
    ++failCount;
  }
  // The line also intersects the surface near (3.95, 0.1, 0.3083726968400695)
  // but that point is farther away from the plane than the expected result.
  // Both points are in the negative normal direction from the plane.
  delete surfacePnt;

  delete mk;

  return ( failCount == 0 );
}
