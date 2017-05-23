#include "meshkit/AF2Neighborhood.hpp"

// MeshKit
#include "meshkit/Error.hpp"

// C++
#include <set>

AF2Neighborhood::AF2Neighborhood(const std::list<AF2Point3D*> & points,
    AF2Edge3D* baselineEdge,
    const std::list<const AF2Edge3D*> & otherEdges,
    const AF2LocalTransform* localTransformArg)
{
  typedef std::list<AF2Point3D*>::const_iterator ConstPoint3DItr;
  typedef std::list<const AF2Edge3D*>::const_iterator ConstEdge3DItr;
  typedef std::map<AF2Point3D*, const AF2Point2D*>::const_iterator MapItr;

  std::set<AF2Point3D*> illegalPoints;

  baseEdge3D = baselineEdge;
  localTransform = localTransformArg;

  for (ConstPoint3DItr itr = points.begin(); itr != points.end(); ++itr)
  {
    bool legal = true;
    AF2Point2D* point2D = localTransform->transformFromSurface(**itr, legal);
    if (legal)
    {
      points2D.push_back(point2D);
      map2DTo3D[point2D] = *itr;
    }
    else
    {
      illegalPoints.insert(*itr);
    }
    map3DTo2D[*itr] = point2D;
  }

  MapItr baseStartItr = map3DTo2D.find(baselineEdge->getStart());
  MapItr baseEndItr = map3DTo2D.find(baselineEdge->getEnd());
  if (baseStartItr == map3DTo2D.end() || baseEndItr == map3DTo2D.end())
  {
    MeshKit::Error badArg(MeshKit::MK_BAD_INPUT);
    badArg.set_string(
        "A baseline edge endpoint is not listed in the neighborhood points.");
    throw badArg;
  }
  baseEdge2D = new AF2Edge2D(baseStartItr->second, baseEndItr->second);
  edges2D.push_back(baseEdge2D);

  for (ConstEdge3DItr itr = otherEdges.begin();
      itr != otherEdges.end(); ++itr)
  {
    if (*itr == baselineEdge)
    {
      // the baseline edge should be listed only once (and listed first)
      // in the list of edges
      continue;
    }
    MapItr startItr = map3DTo2D.find((*itr)->getStart());
    MapItr endItr = map3DTo2D.find((*itr)->getEnd());
    if (startItr == map3DTo2D.end() || endItr == map3DTo2D.end())
    {
      MeshKit::Error badArg(MeshKit::MK_BAD_INPUT);
      badArg.set_string(
          "An edge endpoint is not listed in the neighborhood points.");
      throw badArg;
    }
    if ((illegalPoints.find((*itr)->getStart()) != illegalPoints.end()) ||
        (illegalPoints.find((*itr)->getEnd()) != illegalPoints.end()))
    {
      // Don't create an edge if its endpoints are illegal
      // TODO: Decide whether there is something better to do if one endpoint
      // is legal and the other is not legal
      continue;
    }
    const AF2Edge2D* edge2D =
        new AF2Edge2D(startItr->second, endItr->second);
    edges2D.push_back(edge2D);
  }
}

AF2Neighborhood::~AF2Neighborhood()
{
  typedef std::list<const AF2Point2D*>::const_iterator ConstPoint2DItr;
  typedef std::list<const AF2Edge2D*>::const_iterator ConstEdge2DItr;

  for (ConstEdge2DItr itr = edges2D.begin(); itr != edges2D.end(); ++itr)
  {
    delete *itr;
  }
  for (ConstPoint2DItr itr = points2D.begin(); itr != points2D.end(); ++itr)
  {
    delete *itr;
  }
  delete localTransform;
}

AF2Neighborhood::AF2Neighborhood(const AF2Neighborhood & toCopy)
{
  MeshKit::Error notImpl(MeshKit::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Neighborhood copy construction is not supported.");
  throw notImpl;
}

AF2Neighborhood& AF2Neighborhood::operator=(const AF2Neighborhood & rhs)
{
  MeshKit::Error notImpl(MeshKit::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2Neighborhood assignment operator is not supported.");
  throw notImpl;
}

const AF2Edge2D* AF2Neighborhood::getBaselineEdge2D() const
{
  return baseEdge2D;
}

AF2Edge3D* AF2Neighborhood::getBaselineEdge3D() const
{
  return baseEdge3D;
}

AF2Point3D* AF2Neighborhood::getCorrespondingPoint(
    const AF2Point2D* const & ngbhdPoint2D) const
{
  typedef std::map<const AF2Point2D*, AF2Point3D*>::const_iterator MapItr;
  MapItr ngbhdPntItr = map2DTo3D.find(ngbhdPoint2D);
  if (ngbhdPntItr == map2DTo3D.end())
  {
    return NULL;
  }
  return ngbhdPntItr->second;
}

const std::list<const AF2Edge2D*>* AF2Neighborhood::getEdges2D() const
{
  return &edges2D;
}

const std::list<const AF2Point2D*>* AF2Neighborhood::getPoints2D() const
{
  return &points2D;
}

AF2Point3D* AF2Neighborhood::transformPoint(
    const AF2Point2D* const & point2D, unsigned long const & pntId) const
{
  return localTransform->transformToSurface(*point2D, pntId);
}
