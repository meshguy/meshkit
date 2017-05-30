#include "meshkit/AF2AlgorithmResult.hpp"

// C++
#include <cstddef>

// MeshKit
#include "meshkit/Error.hpp"


AF2AlgorithmResult::AF2AlgorithmResult() :
     succeeded(false), points(), faces()
{
  // nothing more to do
}

AF2AlgorithmResult::AF2AlgorithmResult(
    const std::list<AF2Point3D*> & pointsArg,
    const std::list<const AF2Polygon3D*> & facesArg) :
    succeeded(true), points(pointsArg), faces(facesArg)
{
  // nothing more to do
}

AF2AlgorithmResult::~AF2AlgorithmResult()
{
  typedef std::list<const AF2Polygon3D*>::const_iterator FaceConstItr;
  typedef std::list<AF2Point3D*>::const_iterator PointConstItr;

  for (FaceConstItr itr = faces.begin(); itr != faces.end(); ++itr)
  {
    delete (*itr);
  }

  for (PointConstItr itr = points.begin(); itr != points.end(); ++itr)
  {
    delete (*itr);
  }
}

AF2AlgorithmResult::AF2AlgorithmResult(const AF2AlgorithmResult & toCopy)
{
  MeshKit::Error notImpl(MeshKit::MK_NOT_IMPLEMENTED);
  notImpl.set_string("AF2AlgorithmResult copy construction is not supported.");
  throw notImpl;
}

AF2AlgorithmResult& AF2AlgorithmResult::operator=(
    const AF2AlgorithmResult & rhs)
{
  MeshKit::Error notImpl(MeshKit::MK_NOT_IMPLEMENTED);
  notImpl.set_string(
      "AF2AlgorithmResult assignment operator is not supported.");
  throw notImpl;
}

const std::list<const AF2Polygon3D*>* AF2AlgorithmResult::getFaces() const
{
  if (succeeded)
  {
    return &faces;
  }
  return NULL;
}

const std::list<AF2Point3D*>* AF2AlgorithmResult::getPoints() const
{
  if (succeeded)
  {
    return &points;
  }
  return NULL;
}

bool AF2AlgorithmResult::isSuccessful() const
{
  return succeeded;
}
