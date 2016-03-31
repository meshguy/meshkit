#include "meshkit/AF2FreeZone.hpp"

#include <algorithm>

AF2FreeZone::AF2FreeZone(
    std::list<AF2Point2D> const & bndryPoints) :
    vertices(bndryPoints.begin(), bndryPoints.end())
{
  typedef std::list<AF2Point2D>::const_iterator ItrType;

  // determine the bounding box
  minX = 0;
  maxX = 0;
  minY = 0;
  maxY = 0;
  bool isFirst = true;
  for (ItrType itr = vertices.begin(); itr != vertices.end(); ++itr)
  {
    if (isFirst)
    {
      minX = itr->getX();
      maxX = itr->getX();
      minY = itr->getY();
      maxY = itr->getY();
      isFirst = false;
    }
    else
    {
      minX = std::min(minX, itr->getX());
      maxX = std::max(maxX, itr->getX());
      minY = std::min(minY, itr->getY());
      maxY = std::max(maxY, itr->getY());
    }
  }

  // determine the scale
  scale = std::max(maxX - minX, maxY - minY);
  scale = std::max(scale, fabs(maxX));
  scale = std::max(scale, fabs(minX));
  scale = std::max(scale, fabs(maxY));
  scale = std::max(scale, fabs(minY));
}

bool AF2FreeZone::nearContains(AF2Point2D const & testPnt,
    bool const & containsBndry) const
{
  // Remark: This method could execute faster if the class precomputed
  // and stored the coefficients needed for the comparisons.  That is
  // how netgen version 5.2 did it.
  // The test point x-coefficient would be -rayYDiff.
  // The test point y-coefficient would be rayXDiff.
  // The constant would be ( -rayXDiff * itr->getY() ) + . . .
  //     ( rayYDiff * itr->getX() )

  typedef std::list<AF2Point2D>::const_iterator ItrType;

  // check whether the test point lies outside of the free zone bounding box
  if (testPnt.getX() < minX || testPnt.getX() > maxX ||
      testPnt.getY() < minY || testPnt.getY() > maxY)
  {
    return false;
  }

  // check whether the test point is clockwise (i.e., outside) relative
  // to one of the bounding line segments of the free zone
  for (ItrType itr = vertices.begin(); itr != vertices.end(); ++itr)
  {
    // if the parameter is defined such that the free zone should not
    // be considered to contain points approximately equal to its boundary
    // points, explicitly test whether the test point is approximately
    // equal to the boundary point, and, if so, count it as clockwise
    if (!containsBndry && nearEqual(testPnt, *itr))
    {
      return false;
    }
 
    // get an iterator pointing to the point after the current point,
    // wrapping around if necessary
    ItrType next = itr;
    ++next;
    if (next == vertices.end())
    {
      next = vertices.begin();
    }

    // calculate the ray direction
    double rayXDiff = next->getX() - itr->getX();
    double rayYDiff = next->getY() - itr->getY();

    // test whether the test point is clockwise from the ray from the
    // current point to the next point
    double testXDiff = (testPnt.getX() - itr->getX());
    double testYDiff = (testPnt.getY() - itr->getY());
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

bool AF2FreeZone::nearEqual(AF2Point2D const & pointAlpha,
    AF2Point2D const & pointBravo) const
{
  double xDiff = pointAlpha.getX() - pointBravo.getX();
  double yDiff = pointAlpha.getY() - pointBravo.getY();
  return xDiff*xDiff + yDiff*yDiff < 2.0e-28*scale*scale;
}

bool AF2FreeZone::nearIntersects(AF2Point2D const & startPoint,
    AF2Point2D const & endPoint, bool const & containsBndry) const
{
  // Remark: This method could execute faster if the class precomputed
  // and stored the coefficients needed for the comparisons.  That is
  // how netgen version 5.2 did it.
  // The edge endpoints x-coefficient would be -rayYDiff.
  // The edge endpoints y-coefficient would be rayXDiff.
  // The constant would be ( -rayXDiff * itr->getY() ) + . . .
  //     ( rayYDiff * itr->getX() )

  typedef std::list<AF2Point2D>::const_iterator ItrType;

  // check whether both endpoints lie outside of the free zone bounding box
  if ((startPoint.getX() < minX && endPoint.getX() < minX) ||
      (startPoint.getX() > maxX && endPoint.getX() > maxX) ||
      (startPoint.getY() < minY && endPoint.getY() < minY) ||
      (startPoint.getY() > maxY && endPoint.getY() > maxY))
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

    // establish offsets that can be used to avoid detecting intersections
    // with the boundary as intersections
    double startXOffset = 0.0;
    double startYOffset = 0.0;
    double endXOffset = 0.0;
    double endYOffset = 0.0;
    bool offsetStart = false;
    bool offsetEnd = false;

    if (!containsBndry)
    {
      if (nearEqual(startPoint, *itr))
      {
        if (nearEqual(endPoint, *next))
        {
          // the query line segment matches a boundary line segment
          return false;
        }
        // the start point is approximately equal to a boundary point
        // so it should be offset in the direction of the end point
        offsetStart = true;
      }
      else if (nearEqual(startPoint, *next))
      {
        if (nearEqual(endPoint, *itr))
        {
          // the query line segment matches a boundary line segment
          return false;
        }
        // the start point is approximately equal to a boundary point
        // so it should be offset in the direction of the end point
        offsetStart = true;
      }
      else if (nearEqual(endPoint, *itr) || nearEqual(endPoint, *next))
      {
        offsetEnd = true;
      }

      if (offsetStart || offsetEnd)
      {
        double testSegXDiff = endPoint.getX() - startPoint.getX();
        double testSegYDiff = endPoint.getY() - startPoint.getY();
        double testSegLength = sqrt(testSegXDiff * testSegXDiff  +
            testSegYDiff * testSegYDiff);
        if (testSegLength == 0)
        {
          testSegLength = 1.0;
        }
        if (offsetStart)
        {
          startXOffset = 8.0e-14 * scale * testSegXDiff / testSegLength;
          startYOffset = 8.0e-14 * scale * testSegYDiff / testSegLength;
        }
        if (offsetEnd)
        {
          endXOffset = -8.0e-14 * scale * testSegXDiff / testSegLength;
          endYOffset = -8.0e-14 * scale * testSegYDiff / testSegLength;
        }
      }
    }

    // calculate the ray direction
    double rayXDiff = next->getX() - itr->getX();
    double rayYDiff = next->getY() - itr->getY();

    // test whether the start point and end point are both clockwise from
    // the ray from the current point to the next point
    double startXDiff = (startPoint.getX() + startXOffset - itr->getX());
    double startYDiff = (startPoint.getY() + startYOffset - itr->getY());
    double endXDiff = (endPoint.getX() + endXOffset - itr->getX());
    double endYDiff = (endPoint.getY() + endYOffset - itr->getY());
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
  double testSegXDiff = (endPoint.getX() - startPoint.getX());
  double testSegYDiff = (endPoint.getY() - startPoint.getY());
  for (ItrType itr = vertices.begin();
      itr != vertices.end() && (allCW || allCCW); ++itr)
  {

    // calculate the difference from the start point to the
    // free zone vertex
    double fzvXDiff = itr->getX() - startPoint.getX();
    double fzvYDiff = itr->getY() - startPoint.getY();

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
  typedef std::list<AF2Point2D>::const_iterator ItrType;

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
    double rayXDiff = next->getX() - itr->getX();
    double rayYDiff = next->getY() - itr->getY();

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
      double testXDiff = (testPnt->getX() - itr->getX());
      double testYDiff = (testPnt->getY() - itr->getY());
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
