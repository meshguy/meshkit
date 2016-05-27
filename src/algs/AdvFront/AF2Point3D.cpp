#include "meshkit/AF2Point3D.hpp"

// C++
#include <climits>

// MeshKit
#include "meshkit/Error.hpp"

AF2Point3D::AF2Point3D(double xVal, double yVal, double zVal)
{
  x = xVal;
  y = yVal;
  z = zVal;
  distToBndry = UINT_MAX;
  committed = false;
}

unsigned int AF2Point3D::getDistanceToBoundary() const
{
  return distToBndry;
}

moab::EntityHandle AF2Point3D::getVertexHandle() const
{
  if (!committed)
  {
    MeshKit::Error badState(MeshKit::ErrorCode::MK_FAILURE);
    badState.set_string("The point has not been committed to the mesh yet.");
    throw badState;
  }
  return vertexHandle;
}

double AF2Point3D::getX() const
{
  return x;
}

double AF2Point3D::getY() const
{
  return y;
}

double AF2Point3D::getZ() const
{
  return z;
}

bool AF2Point3D::isCommitted() const
{
  return committed;
}

void AF2Point3D::limitDistanceToBoundary(unsigned int upperBound)
{
  if (upperBound < distToBndry)
  {
    distToBndry = upperBound;
  }
}

void AF2Point3D::setCommittedHandle(const moab::EntityHandle & vertexHandleArg)
{
  vertexHandle = vertexHandleArg;
  committed = true;
}
