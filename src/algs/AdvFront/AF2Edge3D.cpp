#include "meshkit/AF2Edge3D.hpp"

// C++
#include <cstddef>

// MeshKit
#include "meshkit/AF2Front.hpp"

AF2Edge3D::AF2Edge3D(AF2Point3D* start, AF2Point3D* end) :
    startPnt(start), endPnt(end), qualityLevel(1u), observer(NULL)
{
  // no work to do beyond the member initializers
}

void AF2Edge3D::decreaseQuality()
{
  ++qualityLevel;
  if (observer != NULL)
  {
    observer->qualityDecreased(this);
  }
}

AF2Point3D* AF2Edge3D::getStart() const
{
  return startPnt;
}

AF2Point3D* AF2Edge3D::getEnd() const
{
  return endPnt;
}

unsigned int AF2Edge3D::getQualityLevel() const
{
  return qualityLevel;
}

void AF2Edge3D::setObserver(QualityDecreaseObserver* observerArg)
{
  observer = observerArg;
}
