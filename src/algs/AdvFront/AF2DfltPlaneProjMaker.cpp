#include "meshkit/AF2DfltPlaneProjMaker.hpp"

// C++
#include <cmath>

// MeshKit
#include "meshkit/Matrix.hpp"
#include "meshkit/AF2PlaneProjection.hpp"

AF2DfltPlaneProjMaker::AF2DfltPlaneProjMaker(
    iGeom* iGeomPtrArg, iGeom::EntityHandle surf) :
    iGeomPtr(iGeomPtrArg), surface(surf)
{
  // do nothing beyond constructing the members as above
}

AF2LocalTransform* AF2DfltPlaneProjMaker::makeLocalTransform(
    const std::list<AF2Point3D*> & ngbhdPoints,
    const AF2Edge3D* const & baselineEdge,
    const std::list<const AF2Edge3D*> & otherNgbhdEdges) const
{
  // In NetGen PLANESPACE projection type, the x-direction is the direction
  // of the vector along the baseline edge, the normal is the mean of the
  // normals at the endpoints of the baseline edge with the component
  // in the x-direction removed, and the y-direction is the cross product.
  // This implementation matches that, I think.
  AF2Point3D* baseStart = baselineEdge->getStart();
  AF2Point3D* baseEnd = baselineEdge->getEnd();

  // define the origin of the transformation to be the start of the
  // baseline edge
  MeshKit::Vector<3> planeOrigin;
  planeOrigin[0] = baseStart->getX(); 
  planeOrigin[1] = baseStart->getY(); 
  planeOrigin[2] = baseStart->getZ(); 

  // define the x-direction of the transformation as a unit vector in
  // the direction of the baseline edge
  MeshKit::Vector<3> planeXDir;
  planeXDir[0] = baseEnd->getX() - baseStart->getX(); 
  planeXDir[1] = baseEnd->getY() - baseStart->getY(); 
  planeXDir[2] = baseEnd->getZ() - baseStart->getZ(); 
  double scale = std::sqrt(planeXDir[0] * planeXDir[0] +
      planeXDir[1] * planeXDir[1] + planeXDir[2] * planeXDir[2]);
  planeXDir /= scale;

  // Get the components of the vectors that are normal to the surface
  // at the endpoints of the baseline edge
  double startNrmlX;
  double startNrmlY;
  double startNrmlZ;
  iGeomPtr->getEntNrmlXYZ(surface, baseStart->getX(), baseStart->getY(),
      baseStart->getZ(), startNrmlX, startNrmlY, startNrmlZ);
  double endNrmlX;
  double endNrmlY;
  double endNrmlZ;
  iGeomPtr->getEntNrmlXYZ(surface, baseEnd->getX(), baseEnd->getY(),
      baseEnd->getZ(), endNrmlX, endNrmlY, endNrmlZ);

  // Start defining the normal vector to the plane as the mean of the
  // normal vectors to the surface at the endpoints
  MeshKit::Vector<3> planeNormal;
  planeNormal[0] = 0.5 * (startNrmlX + endNrmlX);
  planeNormal[1] = 0.5 * (startNrmlY + endNrmlY);
  planeNormal[2] = 0.5 * (startNrmlZ + endNrmlZ);
  // subtract from the normal vector any component in the x-direction
  planeNormal -= MeshKit::inner_product(planeNormal, planeXDir) * planeXDir;
  // normalize the normal vector to unit length
  planeNormal /= std::sqrt(planeNormal[0] * planeNormal[0] +
      planeNormal[1] * planeNormal[1] + planeNormal[2] * planeNormal[2]);

  // TODO: Support using a sizing function here instead of
  // passing sizing that matches the scale of the baseline edge
  // In NetGen, the scale is the smaller of the global h set by the user
  // and the value of the sizing function at the 3-dimensional midpoint
  // of the baseline edge.
  return new AF2PlaneProjection(iGeomPtr, surface, planeOrigin,
      planeNormal, planeXDir, scale);
}
