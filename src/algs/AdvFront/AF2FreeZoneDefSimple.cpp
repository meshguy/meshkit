#include "meshkit/AF2FreeZoneDefSimple.hpp"
#include "meshkit/Error.hpp"

AF2FreeZoneDefSimple::AF2FreeZoneDefSimple(
        std::list<AF2Point2D> const & rfrncBndryPnts,
        std::list<const AF2PointTransform*> const & bndryPntTrnsfrms)
{
  // confirm that the lists are the same size
  if (rfrncBndryPnts.size() != bndryPntTrnsfrms.size())
  {
    MeshKit::Error badArg(MeshKit::MK_BAD_INPUT);
    badArg.set_string("The lists of reference points and point transforms are not the same size");
    throw badArg;
  }

  numPoints = rfrncBndryPnts.size();
  bndryPoints = new AF2Point2D[numPoints];
  pointTransforms = new const AF2PointTransform*[numPoints];
  int pIndx = 0;
  for (std::list<AF2Point2D>::const_iterator itr =
      rfrncBndryPnts.begin(); itr != rfrncBndryPnts.end(); ++itr)
  {
    bndryPoints[pIndx] = *itr;
    ++pIndx;
  }

  pIndx = 0;
  for (std::list<const AF2PointTransform*>::const_iterator itr =
      bndryPntTrnsfrms.begin(); itr != bndryPntTrnsfrms.end(); ++itr)
  {
    pointTransforms[pIndx] = (*itr)->clone();
    ++pIndx;
  }
}

AF2FreeZoneDefSimple::~AF2FreeZoneDefSimple()
{
  delete[] bndryPoints;
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    delete pointTransforms[pIndx];
  }
  delete[] pointTransforms;
}

AF2FreeZoneDefSimple::AF2FreeZoneDefSimple(const AF2FreeZoneDefSimple & toCopy)
{
  numPoints = toCopy.numPoints;
  bndryPoints = new AF2Point2D[numPoints];
  pointTransforms = new const AF2PointTransform*[numPoints];
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    bndryPoints[pIndx] = toCopy.bndryPoints[pIndx];
    pointTransforms[pIndx] = toCopy.pointTransforms[pIndx]->clone();
  }
}

AF2FreeZoneDefSimple& AF2FreeZoneDefSimple::operator=(
    const AF2FreeZoneDefSimple & rhs)
{
  // copy constructor functionality,
  // but to other parts of memory, not yet to this
  AF2Point2D* otherBndryPoints = new AF2Point2D[numPoints];
  const AF2PointTransform** otherPointTransforms =
      new const AF2PointTransform*[numPoints];
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    otherBndryPoints[pIndx] = rhs.bndryPoints[pIndx];
    otherPointTransforms[pIndx] = rhs.pointTransforms[pIndx]->clone();
  }

  // destructor functionality
  delete[] bndryPoints;
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    delete pointTransforms[pIndx];
  }
  delete[] pointTransforms;

  // transfer ownership from other parts of memory to this object
  numPoints = rhs.numPoints;
  bndryPoints = otherBndryPoints;
  otherBndryPoints = NULL; // not necessary, but to be explicit
  pointTransforms = otherPointTransforms;
  otherPointTransforms = NULL; // not necessary, but to be explicit

  // return this
  return *this;
}

AF2FreeZoneDefSimple* AF2FreeZoneDefSimple::clone() const
{
  return new AF2FreeZoneDefSimple(*this);
}

AF2FreeZone* AF2FreeZoneDefSimple::makeFreeZone(
    AF2Binding const & vertexBinding, unsigned int qualityClass) const
{
  std::list<AF2Point2D> freeZoneBndry;
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    freeZoneBndry.push_back(pointTransforms[pIndx]->transformPoint(
        bndryPoints[pIndx], vertexBinding));
  }
  AF2FreeZone* freeZone = new AF2FreeZone(freeZoneBndry);
  return freeZone;
}
