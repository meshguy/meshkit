#include "meshkit/AF2FreeZoneDefLCQualLim.hpp"
#include "meshkit/Error.hpp"

AF2FreeZoneDefLCQualLim::AF2FreeZoneDefLCQualLim(
        std::list<MeshKit::Vector<2> > const & preferBndryPnts,
        std::list<const AF2PointTransform*> const & preferPntTrnsfrms,
        std::list<MeshKit::Vector<2> > const & limitBndryPnts,
        std::list<const AF2PointTransform*> const & limitPntTrnsfrms)
{
  // confirm that the lists are the same size
  if ((preferBndryPnts.size() != preferPntTrnsfrms.size()) ||
      (preferBndryPnts.size() != limitBndryPnts.size()) ||
      (preferBndryPnts.size() != limitPntTrnsfrms.size()))
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string(
        "The lists passed to the constructor are not all the same size");
    throw badArg;
  }

  numPoints = preferBndryPnts.size();
  prefBndryPoints = new MeshKit::Vector<2>[numPoints];
  prefPointTransforms = new const AF2PointTransform*[numPoints];
  limBndryPoints = new MeshKit::Vector<2>[numPoints];
  limPointTransforms = new const AF2PointTransform*[numPoints];
  int pIndx = 0;
  for (std::list<MeshKit::Vector<2> >::const_iterator itr =
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
  for (std::list<MeshKit::Vector<2> >::const_iterator itr =
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
  prefBndryPoints = new MeshKit::Vector<2>[numPoints];
  prefPointTransforms = new const AF2PointTransform*[numPoints];
  limBndryPoints = new MeshKit::Vector<2>[numPoints];
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
  MeshKit::Vector<2>* otherPrefBndryPoints = new MeshKit::Vector<2>[numPoints];
  const AF2PointTransform** otherPrefPointTransforms =
      new const AF2PointTransform*[numPoints];
  MeshKit::Vector<2>* otherLimBndryPoints = new MeshKit::Vector<2>[numPoints];
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
    AF2VertexBinding & vertexBinding, int qualityClass) const
{
  // Check that the quality class is greater than or equal to zero
  if (qualityClass < 0)
  {
    MeshKit::Error badArg(MeshKit::ErrorCode::MK_BAD_INPUT);
    badArg.set_string("The quality class is less than zero.");
    throw badArg;
  }

  // compute the coefficients to use in the linear combination
  double prefCoeff = 1.0;
  if (qualityClass > 0)
  {
    prefCoeff = 1.0/qualityClass;
  }
  double limCoeff = 1.0 - prefCoeff;

  // compute the list of actual free zone boundary points
  std::list<MeshKit::Vector<2> > freeZoneBndry;
  for (int pIndx = 0; pIndx < numPoints; ++pIndx)
  {
    MeshKit::Vector<2> prefPnt = prefPointTransforms[pIndx]->transformPoint(
        prefBndryPoints[pIndx], vertexBinding);
    MeshKit::Vector<2> limPnt = limPointTransforms[pIndx]->transformPoint(
        limBndryPoints[pIndx], vertexBinding);
    prefPnt *= prefCoeff;
    limPnt *= limCoeff;
    freeZoneBndry.push_back(prefPnt + limPnt);
  }

  // construct and return the free zone
  AF2FreeZone* freeZone = new AF2FreeZone(freeZoneBndry);
  return freeZone;
}
