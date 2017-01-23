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
    MeshKit::Vector<3> const & planeXDir,
    double scaleFactor) :
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

  if (scaleFactor <= 0.0)
  {
    throw MeshKit::Error(MeshKit::MK_BAD_INPUT,
        "Scale factor is not positive.");
  }
  scale = 1.0 / scaleFactor;
}

AF2Point2D* AF2PlaneProjection::transformFromSurface(
    AF2Point3D const & srfcPnt, bool & legal) const
{
  double diffX = srfcPnt.getX() - pOrigin[0];
  double diffY = srfcPnt.getY() - pOrigin[1];
  double diffZ = srfcPnt.getZ() - pOrigin[2];

  double diffDotNormal = diffX * pNormal[0] +
      diffY * pNormal[1] + diffZ * pNormal[2];
  double normalCmpX = diffDotNormal * pNormal[0];
  double normalCmpY = diffDotNormal * pNormal[1];
  double normalCmpZ = diffDotNormal * pNormal[2];

  double planarCmpX = diffX - normalCmpX;
  double planarCmpY = diffY - normalCmpY;
  double planarCmpZ = diffZ - normalCmpZ;

  double planePntX = scale * (planarCmpX * pXDir[0] +
      planarCmpY * pXDir[1] + planarCmpZ * pXDir[2]);
  double planePntY = scale * (planarCmpX * pYDir[0] +
      planarCmpY * pYDir[1] + planarCmpZ * pYDir[2]);

  // compute normal at the surface point; if too different orientation, we should not
  // consider the point legal
  // maybe we should keep the normal to the points on the front available always, we repeat this too much
  double normalX, normalY, normalZ;
  iGeomPtr->getEntNrmlXYZ(surface, srfcPnt.getX(), srfcPnt.getY(),
      srfcPnt.getZ(), normalX, normalY, normalZ);
  // normalize?
  double sqrN= sqrt(normalX*normalX + normalY*normalY + normalZ*normalZ);
  if (sqrN<=1.e-12)
    legal= false;
  else
  {
    normalX/=sqrN; normalY/=sqrN; normalZ/=sqrN;
    double nrmlDot = normalX * pNormal[0] +
        normalY * pNormal[1] + normalZ * pNormal[2];

    if (nrmlDot < 0.0009765625) // 1.0/1024 == 0.0009765625)
      legal= false;
  }
  /*if (legal && (fabs(diffDotNormal) > maxDistanceFromPlane))
  {
    maxDistanceFromPlane = fabs(diffDotNormal);
  }*/
  return new AF2Point2D(planePntX, planePntY);
}

AF2Point3D* AF2PlaneProjection::transformToSurface(
    AF2Point2D const & planePnt, unsigned long const & pntId) const
{
  MeshKit::Vector<3> rayOrigin(pOrigin);
  rayOrigin += (planePnt.getX() / scale) * pXDir +
      (planePnt.getY() / scale) * pYDir;
 /* // we are looking for the point on the surface, within some distance away
  // from our plane;
  const double MYTOL = 1.e-7;
  // maxDistanceFromPlane is a good measure of how far is the surface from the projection plane
  double alfa = 1000*MYTOL + 3.0*maxDistanceFromPlane;
  if (alfa >2.0)
   alfa =2.0; // do not go overboard
  if (alfa < 0.5)
     alfa = 0.5;*/

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

  return new AF2Point3D(pntId, closestPointOnSrfc.x(),
      closestPointOnSrfc.y(), closestPointOnSrfc.z());
}
