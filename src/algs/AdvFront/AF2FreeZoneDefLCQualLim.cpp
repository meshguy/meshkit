// MeshhKit
#include "meshkit/AF2FreeZoneDefLCQualLim.hpp"
#include "meshkit/Error.hpp"

AF2FreeZoneDefLCQualLim::AF2FreeZoneDefLCQualLim(
        std::list<AF2Point2D> const & preferBndryPnts,
        std::list<const AF2PointTransform*> const & preferPntTrnsfrms,
        std::list<AF2Point2D> const & limitBndryPnts,
        std::list<const AF2PointTransform*> const & limitPntTrnsfrms)
{
  // confirm that the lists are the same size
  if ((preferBndryPnts.size() != preferPntTrnsfrms.size()) ||
      (preferBndryPnts.size() != limitBndryPnts.size()) ||
      (preferBndryPnts.size() != limitPntTrnsfrms.size()))
  {
    MeshKit::Error badArg(MeshKit::MK_BAD_INPUT);
    badArg.set_string(
        "The lists passed to the constructor are not all the same size");
    throw badArg;
  }

  numPoints = preferBndryPnts.size();
  prefBndryPoints = new AF2Point2D[numPoints];
  prefPointTransforms = new const AF2PointTransform*[numPoints];
  limBndryPoints = new AF2Point2D[numPoints];
  limPointTransforms = new const AF2PointTransform*[numPoints];
  int pIndx = 0;
  for (std::list<AF2Point2D>::const_iterator itr =
      preferBndryPnts.begin(); itr != preferBndryPnts.end(); ++itr)
  {
    prefBndryPoints[pIndx] = *itr;
    ++pIndx;
  }

  pIndx = 0;
  for (std::list<const AF2PointTransform*>::const_iterator itr =
      preferPntTrnsfrms.begin(); itr != preferPntTrnsfrms.end(); ++itr)
  {
    prefPointTransforms[pIndx] = (*itr)->clone();
    ++pIndx;
  }

  pIndx = 0;
  for (std::list<AF2Point2D>::const_iterator itr =
      limitBndryPnts.begin(); itr != limitBndryPnts.end(); ++itr)
  {
    limBndryPoints[pIndx] = *itr;
    ++pIndx;
  }

  pIndx = 0;
  for (std::list<const AF2PointTransform*>::const_iterator itr =
      limitPntTrnsfrms.begin(); itr != limitPntTrnsfrms.end(); ++itr)
  {
    limPointTransforms[pIndx] = (*itr)->clone();
    ++pIndx;
  }
}

AF2FreeZoneDefLCQualLim::~AF2FreeZoneDefLCQualLim()
{
  delete[] prefBndryPoints;
  delete[] limBndryPoints;
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    delete prefPointTransforms[pIndx];
    delete limPointTransforms[pIndx];
  }
  delete[] prefPointTransforms;
  delete[] limPointTransforms;
}

AF2FreeZoneDefLCQualLim::AF2FreeZoneDefLCQualLim(
    const AF2FreeZoneDefLCQualLim & toCopy)
{
  numPoints = toCopy.numPoints;
  prefBndryPoints = new AF2Point2D[numPoints];
  prefPointTransforms = new const AF2PointTransform*[numPoints];
  limBndryPoints = new AF2Point2D[numPoints];
  limPointTransforms = new const AF2PointTransform*[numPoints];
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    prefBndryPoints[pIndx] = toCopy.prefBndryPoints[pIndx];
    prefPointTransforms[pIndx] = toCopy.prefPointTransforms[pIndx]->clone();
    limBndryPoints[pIndx] = toCopy.limBndryPoints[pIndx];
    limPointTransforms[pIndx] = toCopy.limPointTransforms[pIndx]->clone();
  }
}

AF2FreeZoneDefLCQualLim& AF2FreeZoneDefLCQualLim::operator=(
    const AF2FreeZoneDefLCQualLim & rhs)
{
  // copy constructor functionality,
  // but to other parts of memory, not yet to this
  AF2Point2D* otherPrefBndryPoints = new AF2Point2D[numPoints];
  const AF2PointTransform** otherPrefPointTransforms =
      new const AF2PointTransform*[numPoints];
  AF2Point2D* otherLimBndryPoints = new AF2Point2D[numPoints];
  const AF2PointTransform** otherLimPointTransforms =
      new const AF2PointTransform*[numPoints];
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    otherPrefBndryPoints[pIndx] = rhs.prefBndryPoints[pIndx];
    otherPrefPointTransforms[pIndx] = rhs.prefPointTransforms[pIndx]->clone();
    otherLimBndryPoints[pIndx] = rhs.limBndryPoints[pIndx];
    otherLimPointTransforms[pIndx] = rhs.limPointTransforms[pIndx]->clone();
  }

  // destructor functionality
  delete[] prefBndryPoints;
  delete[] limBndryPoints;
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    delete prefPointTransforms[pIndx];
    delete limPointTransforms[pIndx];
  }
  delete[] prefPointTransforms;
  delete[] limPointTransforms;

  // transfer ownership from other parts of memory to this object
  numPoints = rhs.numPoints;
  prefBndryPoints = otherPrefBndryPoints;
  otherPrefBndryPoints = NULL; // not necessary, but to be explicit
  prefPointTransforms = otherPrefPointTransforms;
  otherPrefPointTransforms = NULL; // not necessary, but to be explicit
  limBndryPoints = otherLimBndryPoints;
  otherLimBndryPoints = NULL; // not necessary, but to be explicit
  limPointTransforms = otherLimPointTransforms;
  otherLimPointTransforms = NULL; // not necessary, but to be explicit

  // return this
  return *this;
}

AF2FreeZoneDefLCQualLim* AF2FreeZoneDefLCQualLim::clone() const
{
  return new AF2FreeZoneDefLCQualLim(*this);
}

AF2FreeZone* AF2FreeZoneDefLCQualLim::makeFreeZone(
    AF2Binding const & vertexBinding, unsigned int qualityClass) const
{
  // Check that the quality class is greater than zero
  if (qualityClass == 0)
  {
    MeshKit::Error badArg(MeshKit::MK_BAD_INPUT);
    badArg.set_string("The quality class is not greater than zero.");
    throw badArg;
  }

  // compute the coefficients to use in the linear combination
  double prefCoeff = 1.0/qualityClass;
  double limCoeff = 1.0 - prefCoeff;

  // compute the list of actual free zone boundary points
  std::list<AF2Point2D> freeZoneBndry;
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    AF2Point2D prefPnt = prefPointTransforms[pIndx]->transformPoint(
        prefBndryPoints[pIndx], vertexBinding);
    AF2Point2D limPnt = limPointTransforms[pIndx]->transformPoint(
        limBndryPoints[pIndx], vertexBinding);
    AF2Point2D fzBndryPnt(
        prefCoeff * prefPnt.getX()  +  limCoeff * limPnt.getX(),
        prefCoeff * prefPnt.getY()  +  limCoeff * limPnt.getY());
    freeZoneBndry.push_back(fzBndryPnt);
  }

  // construct and return the free zone
  AF2FreeZone* freeZone = new AF2FreeZone(freeZoneBndry);
  return freeZone;
}
