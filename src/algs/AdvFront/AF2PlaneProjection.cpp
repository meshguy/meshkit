#include "meshkit/AF2PlaneProjection.hpp"

// C++
#include <cmath>

// CGM
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "RefFace.hpp"
#include "Surface.hpp"

// MeshKit
#include "meshkit/Error.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/VecUtil.hpp"

AF2PlaneProjection::AF2PlaneProjection(iGeom* iGeomPtrArg,
    iGeom::EntityHandle srfcHandle,
    MeshKit::Vector<3> const & planeOrigin,
    MeshKit::Vector<3> const & planeNormal,
    MeshKit::Vector<3> const & planeXDir) :
    iGeomPtr(iGeomPtrArg),
    surface(srfcHandle),
    pOrigin(planeOrigin),
    pNormal(planeNormal),
    pXDir(planeXDir)
{
  // exercise the iGeom pointer and check that the entity is a face
  iBase_EntityType eTypeResult;
  MeshKit::IBERRCHK(iGeomPtr->getEntType(surface, eTypeResult), *iGeomPtr);
  if (eTypeResult != iBase_FACE)
  {
    throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "Geometry entity handle does not refer to a surface in AF2PlaneProjection constructor.");
  }

  // compute the squared length of the normal vector
  double nvNormSqrd = pNormal[0] * pNormal[0] +
      pNormal[1] * pNormal[1] + pNormal[2] * pNormal[2];

  // check that the length is sufficiently non-zero
  if (std::fabs(nvNormSqrd) < 1.0e-16)
  {
    throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "Normal vector length is near zero in AF2PlaneProjection constructor.");
  }

  // normalize the normal vector to have length one, if necessary
  if (std::fabs(nvNormSqrd - 1.0) > 1.0e-8)
  {
    pNormal *= 1.0 / std::sqrt(nvNormSqrd);
  }

  // compute the squared length of the x-direction vector
  double xdvNormSqrd = pXDir[0] * pXDir[0] +
      pXDir[1] * pXDir[1] + pXDir[2] * pXDir[2];

  // check that the length is sufficiently non-zero
  if (std::fabs(xdvNormSqrd) < 1.0e-16)
  {
    throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "X-direction vector length is near zero in AF2PlaneProjection constructor.");
  }

  // normalize the x-direction vector to have length one, if necessary
  if (std::fabs(xdvNormSqrd - 1.0) > 1.0e-8)
  {
    pXDir *= 1.0 / std::sqrt(xdvNormSqrd);
  }

  double nxDot = MeshKit::VecUtil::dot((double*)pNormal.data(),
      (double*)pXDir.data());
  if (std::fabs(nxDot) > 1.0e-12)
  {
    throw MeshKit::Error(MeshKit::MK_BAD_INPUT, "X-direction vector is not normal to the normal vector in AF2PlaneProjection constructor.");
  }

  // compute the y-direction vector
  MeshKit::VecUtil::cross((double*)pNormal.data(), (double*)pXDir.data(),
      (double*)pYDir.data());
}

void AF2PlaneProjection::transformFromSurface(
    MeshKit::Vector<3> const & srfcPnt, MeshKit::Vector<2> & planePnt) const
{
  MeshKit::Vector<3> diffVec = srfcPnt - pOrigin;
  MeshKit::Vector<3> normalCmpnnt(pNormal);
  normalCmpnnt *= MeshKit::VecUtil::dot((double*)diffVec.data(),
      (double*)pNormal.data());
  MeshKit::Vector<3> planeVec(diffVec);
  planeVec -= normalCmpnnt;
  planePnt[0] = MeshKit::VecUtil::dot((double*)planeVec.data(),
      (double*)pXDir.data());
  planePnt[1] = MeshKit::VecUtil::dot((double*)planeVec.data(),
      (double*)pYDir.data());
}

void AF2PlaneProjection::transformToSurface(
    MeshKit::Vector<2> const & planePnt, MeshKit::Vector<3> & srfcPnt) const
{
  MeshKit::Vector<3> rayOrigin(pOrigin);
  rayOrigin += planePnt[0] * pXDir + planePnt[1] * pYDir;
  CubitVector cvRayOrigin(rayOrigin[0], rayOrigin[1], rayOrigin[2]);
  CubitVector cvRayDir(pNormal[0], pNormal[1], pNormal[2]);
  CubitVector cvPointOnSrfc;
  CubitVector closestPointOnSrfc;

  // closest point in positive normal direction
  RefFace* cgmFacePtr =
      dynamic_cast<RefFace*>(reinterpret_cast<RefEntity*>(surface));
  Surface* cgmSrfcPtr = cgmFacePtr->get_surface_ptr();
  CubitStatus posDirResult = cgmSrfcPtr->closest_point_along_vector(
      cvRayOrigin, cvRayDir, cvPointOnSrfc);
  if (posDirResult == CUBIT_SUCCESS)
  {
    closestPointOnSrfc = cvPointOnSrfc;
  }

  cvRayDir *= -1;
  CubitStatus negDirResult = cgmSrfcPtr->closest_point_along_vector(
      cvRayOrigin, cvRayDir, cvPointOnSrfc);
  if (negDirResult == CUBIT_SUCCESS)
  {
    if (posDirResult == CUBIT_SUCCESS)
    {
      // REMARK: It might be better to measure the distance in
      // parametric coordinates from some source point . . .

      // compute square distance to point in positive direction
      double posXDiff = closestPointOnSrfc.x() - rayOrigin[0];
      double posYDiff = closestPointOnSrfc.y() - rayOrigin[1];
      double posZDiff = closestPointOnSrfc.z() - rayOrigin[2];
      double posSqDist = posXDiff * posXDiff +
          posYDiff * posYDiff + posZDiff * posZDiff;

      // compute square distance to point in negative direction
      double negXDiff = cvPointOnSrfc.x() - rayOrigin[0];
      double negYDiff = cvPointOnSrfc.y() - rayOrigin[1];
      double negZDiff = cvPointOnSrfc.z() - rayOrigin[2];
      double negSqDist = negXDiff * negXDiff +
          negYDiff * negYDiff + negZDiff * negZDiff;

      if (negSqDist < posSqDist)
      {
        closestPointOnSrfc = cvPointOnSrfc;
      }
    }
    else
    {
      closestPointOnSrfc = cvPointOnSrfc;
    }
  }

  srfcPnt[0] = closestPointOnSrfc.x();
  srfcPnt[1] = closestPointOnSrfc.y();
  srfcPnt[2] = closestPointOnSrfc.z();
}
