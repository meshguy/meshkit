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
      this_surf->geom_handle(), origin, normal, xDir);
  std::cout << "Created plane passing through (2.8, 2.8, 0) "
      << "normal to (1, 1, 0)." << std::endl;
  std::cout << "The y-coordinate in the plane is in the direction of "
      << "the z-axis." << std::endl;


  MeshKit::Vector<3> testPnt;
  testPnt[0] = 2.8059735404411; // (sqrt(2)/2) * (3 + cos(asin(0.25)))
  testPnt[1] = 2.8059735404411;
  testPnt[2] = 0.25;
  MeshKit::Vector<2> planePnt;
  planeProj.transformFromSurface(testPnt, planePnt);
  std::cout << "Transform (" << testPnt[0] << ", " << testPnt[1] << ", "
      << testPnt[2] << ") from surface to plane." << std::endl;
  std::cout << "Coordinates in plane: (" << planePnt[0] << ", "
      << planePnt[1] << ")" << std::endl;
  if ((fabs(planePnt[0] - 0.0) + fabs(planePnt[1] - 0.25)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates in the plane should be (0, 0.25)."
        << std::endl;
    ++failCount;
  }

  MeshKit::Vector<2> planeTestPnt;
  planeTestPnt[0] = 0;
  planeTestPnt[1] = 0.1539268506556547; // sin(acos(2.82 * sqrt(2) - 3))
  MeshKit::Vector<3> surfacePnt;
  planeProj.transformToSurface(planeTestPnt, surfacePnt);
  std::cout << "Transform (" << planeTestPnt[0] << ", " << planeTestPnt[1]
      << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt[0] << ", "
      << surfacePnt[1] << ", " << surfacePnt[2] << ")" << std::endl;
  if ((fabs(surfacePnt[0] - 2.82) + fabs(surfacePnt[1] - 2.82) +
      fabs(surfacePnt[2] - 0.15392685065547)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (2.82, 2.82, 0.15392685065547)."
        << std::endl;
    ++failCount;
  }

  planeTestPnt[0] = 0;
  planeTestPnt[1] = 0.574682269122808; // sin(acos(2.7 * sqrt(2) - 3))
  planeProj.transformToSurface(planeTestPnt, surfacePnt);
  std::cout << "Transform (" << planeTestPnt[0] << ", " << planeTestPnt[1]
      << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt[0] << ", "
      << surfacePnt[1] << ", " << surfacePnt[2] << ")" << std::endl;
  if ((fabs(surfacePnt[0] - 2.7) + fabs(surfacePnt[1] - 2.7) +
      fabs(surfacePnt[2] - 0.574682269122808)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (2.7, 2.7, 0.574682269122808)."
        << std::endl;
    ++failCount;
  }


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
      this_surf->geom_handle(), origin, normal, xDir);
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
  planeTestPnt[0] = 0.53625047660234298437463798744018;
  planeTestPnt[1] = 0.1;
  wackyProj.transformToSurface(planeTestPnt, surfacePnt);
  std::cout << "Transform (" << planeTestPnt[0] << ", " << planeTestPnt[1]
      << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt[0] << ", "
      << surfacePnt[1] << ", " << surfacePnt[2] << ")" << std::endl;
  // expect ~ (3.95, 0.1, 0.3083726968400695)
  if ((fabs(surfacePnt[0] - 3.95) + fabs(surfacePnt[1] - 0.1) +
      fabs(surfacePnt[2] - 0.3083726968400695)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (3.95, 0.1, 0.3083726968400695)."
        << std::endl;
    ++failCount;
  }
  // The line also intersects the surface near (3.30723, 0.1, 0.951148),
  // but that point is farther away from the plane than the expected result.


  // Perform a similar test with the plane origin at x = 2.6
  origin[0] = 2.6;
  AF2PlaneProjection wackyProjShift01(this_surf->igeom_instance(),
      this_surf->geom_handle(), origin, normal, xDir);
  std::cout << "Created plane passing through (2.6, 0, 0) "
      << "normal to (-1, 0, 1)." << std::endl;
  std::cout << "The y-coordinate in the plane is in the direction of "
      << "the y-axis." << std::endl;

  planeTestPnt[0] = 1.17264657967023277;
  wackyProjShift01.transformToSurface(planeTestPnt, surfacePnt);
  std::cout << "Transform (" << planeTestPnt[0] << ", " << planeTestPnt[1]
      << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt[0] << ", "
      << surfacePnt[1] << ", " << surfacePnt[2] << ")" << std::endl;
  if ((fabs(surfacePnt[0] - 3.3072251285719) + fabs(surfacePnt[1] - 0.1) +
      fabs(surfacePnt[2] - 0.9511475682682)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (3.30722512857195, 0.1, 0.9511475682682)."
        << std::endl;
    ++failCount;
  }
  // The line also intersects the surface near (3.95, 0.1, 0.3083726968400695)
  // but that point is farther away from the plane than the expected result.


  // Perform a similar test with the plane origin at x = 2.1
  origin[0] = 2.1;
  AF2PlaneProjection wackyProjShift02(this_surf->igeom_instance(),
      this_surf->geom_handle(), origin, normal, xDir);
  std::cout << "Created plane passing through (2.1, 0, 0) "
      << "normal to (-1, 0, 1)." << std::endl;
  std::cout << "The y-coordinate in the plane is in the direction of "
      << "the y-axis." << std::endl;

  planeTestPnt[0] = 1.5261999702635065;
  wackyProjShift02.transformToSurface(planeTestPnt, surfacePnt);
  std::cout << "Transform (" << planeTestPnt[0] << ", " << planeTestPnt[1]
      << ") from plane to surface." << std::endl;
  std::cout << "Coordinates on surface: (" << surfacePnt[0] << ", "
      << surfacePnt[1] << ", " << surfacePnt[2] << ")" << std::endl;
  if ((fabs(surfacePnt[0] - 3.3072251285719) + fabs(surfacePnt[1] - 0.1) +
      fabs(surfacePnt[2] - 0.9511475682682)) > 1e-12)
  {
    std::cout << "FAIL: The coordinates on the surface should be\n"
        << "  (3.30722512857195, 0.1, 0.9511475682682)."
        << std::endl;
    ++failCount;
  }
  // The line also intersects the surface near (3.95, 0.1, 0.3083726968400695)
  // but that point is farther away from the plane than the expected result.
  // Both points are in the negative normal direction from the plane.

  delete mk;

  return ( failCount == 0 );
}
