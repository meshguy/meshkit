#include "meshkit/AF2FreeZone.hpp"

#include <algorithm>

// temporary
#include <iostream>

AF2FreeZone::AF2FreeZone(
    std::list<MeshKit::Vector<2> > const& bndryPoints) :
    vertices(bndryPoints.begin(), bndryPoints.end())
{
  typedef std::list<MeshKit::Vector<2> >::const_iterator ItrType;

  // determine the bounding box
  minX = 0;
  maxX = 0;
  minY = 0;
  maxY = 0;
  for (ItrType itr = vertices.begin(); itr != vertices.end(); ++itr)
  {
    minX = std::min(minX, (*itr)[0]);
    maxX = std::max(maxX, (*itr)[0]);
    minY = std::min(minY, (*itr)[1]);
    maxY = std::max(maxY, (*itr)[1]);
  }

  // determine the scale
  scale = std::max(maxX - minX, maxY - minY);
  scale = std::max(scale, fabs(maxX));
  scale = std::max(scale, fabs(minX));
  scale = std::max(scale, fabs(maxY));
  scale = std::max(scale, fabs(minY));
}

bool AF2FreeZone::nearContains(MeshKit::Vector<2> const & testPnt) const
{
  // Remark: This method could execute faster if the class precomputed
  // and stored the coefficients needed for the comparisons.  That is
  // how netgen version 5.2 did it.
  // The test point x-coefficient would be -rayYDiff.
  // The test point y-coefficient would be rayXDiff.
  // The constant would be ( -rayXDiff * (*itr)[1] ) + ( rayYDiff * (*itr)[0] )

  typedef std::list<MeshKit::Vector<2> >::const_iterator ItrType;

  // check whether the test point lies outside of the free zone bounding box
  if (testPnt[0] < minX || testPnt[0] > maxX ||
      testPnt[1] < minY || testPnt[1] > maxY)
  {
    return false;
  }

  // check whether the test point is clockwise (i.e., outside) relative
  // to one of the bounding line segments of the free zone
  for (ItrType itr = vertices.begin(); itr != vertices.end(); ++itr)
  {
    // get an iterator pointing to the point after the current point,
    // wrapping around if necessary
    ItrType next = itr;
    ++next;
    if (next == vertices.end())
    {
      next = vertices.begin();
    }

    // calculate the ray direction
    double rayXDiff = (*next)[0] - (*itr)[0];
    double rayYDiff = (*next)[1] - (*itr)[1];

    // test whether the test point is clockwise from the ray from the
    // current point to the next point
    double testXDiff = (testPnt[0] - (*itr)[0]);
    double testYDiff = (testPnt[1] - (*itr)[1]);
    if (testYDiff*rayXDiff - rayYDiff*testXDiff < -1.0e-14*scale*scale)
    {
      // if this condition is true (test would be < 0 in perfect arithmetic,
      // but want to guarantee the point is outside here), then the
      // test point is clockwise from the ray
      // (== 0 would mean collinear in perfect arithmetic)
      return false;
    }
  }

  return true;
}

bool AF2FreeZone::nearIntersects(MeshKit::Vector<2> const & startPoint,
    MeshKit::Vector<2> const & endPoint) const
{
  // Remark: This method could execute faster if the class precomputed
  // and stored the coefficients needed for the comparisons.  That is
  // how netgen version 5.2 did it.
  // The edge endpoints x-coefficient would be -rayYDiff.
  // The edge endpoints y-coefficient would be rayXDiff.
  // The constant would be ( -rayXDiff * (*itr)[1] ) + ( rayYDiff * (*itr)[0] )

  typedef std::list<MeshKit::Vector<2> >::const_iterator ItrType;

  // check whether both endpoints lie outside of the free zone bounding box
  if ((startPoint[0] < minX && endPoint[0] < minX) ||
      (startPoint[0] > maxX && endPoint[0] > maxX) ||
      (startPoint[1] < minY && endPoint[1] < minY) ||
      (startPoint[1] > maxY && endPoint[1] > maxY))
  {
    return false;
  }

  // check whether both endpoints of the edge are clockwise (i.e., outside)
  // relative to one of the bounding line segments of the free zone
  for (ItrType itr = vertices.begin(); itr != vertices.end(); ++itr)
  {
    // get an iterator pointing to the point after the current point,
    // wrapping around if necessary
    ItrType next = itr;
    ++next;
    if (next == vertices.end())
    {
      next = vertices.begin();
    }

    // calculate the ray direction
    double rayXDiff = (*next)[0] - (*itr)[0];
    double rayYDiff = (*next)[1] - (*itr)[1];

    // test whether the start point and end point are both clockwise from
    // the ray from the current point to the next point
    double startXDiff = (startPoint[0] - (*itr)[0]);
    double startYDiff = (startPoint[1] - (*itr)[1]);
    double endXDiff = (endPoint[0] - (*itr)[0]);
    double endYDiff = (endPoint[1] - (*itr)[1]);
    if ((startYDiff*rayXDiff - rayYDiff*startXDiff < -1.0e-14*scale*scale) &&
        (endYDiff*rayXDiff - rayYDiff*endXDiff < -1.0e-14*scale*scale))
    {
      // if these condition are true (test would be < 0 in perfect arithmetic,
      // but want to guarantee the point is outside here), then the
      // whole edge is clockwise from the ray
      // (== 0 would mean collinear in perfect arithmetic)
      return false;
    }
  }

  // Check whether the whole free zone is either clockwise (right) or
  // counterclockwise (left) of the test line segment
  bool allCW = true;
  bool allCCW = true;
  // First calculate the x and y difference of the test line segment
  double testSegXDiff = (endPoint[0] - startPoint[0]);
  double testSegYDiff = (endPoint[1] - startPoint[1]);
  for (ItrType itr = vertices.begin();
      itr != vertices.end() && (allCW || allCCW); ++itr)
  {

    // calculate the difference from the start point to the
    // free zone vertex
    double fzvXDiff = (*itr)[0] - startPoint[0];
    double fzvYDiff = (*itr)[1] - startPoint[1];

    // test whether the free zone vertex is clockwise from (right of)
    // the test line segment
    if (fzvYDiff*testSegXDiff - testSegYDiff*fzvXDiff > -1.0e-14*scale*scale)
    {
      // if this condition is true (test would be >= 0 in perfect arithmetic,
      // but want to guarantee the point is clockwise if the condition
      // is false), then the free zone vertex is (probably) not clockwise
      // from the test line segment.
      allCW = false;
    }

    // test whether the free zone vertex is counterclockwise from (left of)
    // the test line segment
    if (fzvYDiff*testSegXDiff - testSegYDiff*fzvXDiff < 1.0e-14*scale*scale)
    {
      // if this condition is true (test would be <= 0 in perfect arithmetic,
      // but want to guarantee the point is counterclockwise if the condition
      // is false), then the free zone vertex is (probably) not
      // counterclockwise from the test line segment.
      allCCW = false;
    }
  }

  return ( !allCW && !allCCW );
}

bool AF2FreeZone::isConvex() const
{
  typedef std::list<MeshKit::Vector<2> >::const_iterator ItrType;

  for (ItrType itr = vertices.begin(); itr != vertices.end(); ++itr)
  {
    // get an iterator pointing to the point after the current point,
    // wrapping around if necessary
    ItrType next = itr;
    ++next;
    if (next == vertices.end())
    {
      next = vertices.begin();
    }

    // calculate the ray direction
    double rayXDiff = (*next)[0] - (*itr)[0];
    double rayYDiff = (*next)[1] - (*itr)[1];

    // test whether all other points are counterclockwise
    for (ItrType testPnt = vertices.begin();
        testPnt != vertices.end(); ++testPnt)
    {
      // skip the test if the test point is the current point or the next point
      if ((*testPnt) == (*itr) || (*testPnt) == (*next))
      {
        continue;
      }

      // test that the test point is counterclockwise from the ray from the
      // current point to the next point
      double testXDiff = ((*testPnt)[0] - (*itr)[0]);
      double testYDiff = ((*testPnt)[1] - (*itr)[1]);
      if (testYDiff*rayXDiff - rayYDiff*testXDiff < -1.0e-14*scale*scale)
      {
        // if this condition is true (test would be < 0 in perfect arithmetic),
        // then the test point is clockwise from the ray
        // (== 0 would mean collinear in perfect arithmetic)
        return false;
      }
    }
  }

  return true;
}
