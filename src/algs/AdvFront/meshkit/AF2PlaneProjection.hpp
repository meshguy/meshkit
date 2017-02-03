/*
 * AF2PlaneProjection.hpp
 *
 * A transformation between 2-dimensional space and 3-dimensional space
 * that transforms from three dimensions to two dimensions by subtracting
 * some origin point from the 3-dimensional point and removing the
 * component that is in the direction normal to the plane.  In other words,
 * the transformation projects the 3-dimensional point onto a plane.
 * After the projection, the point is scaled in the plane.
 * 
 * The inverse works by first reversing the scaling, and then firing a
 * ray from the plane in the direction normal to the plane and finding
 * the intersection of the ray with the surface.
 */

#ifndef AF2PLANEPROJECTION_HPP
#define AF2PLANEPROJECTION_HPP

// MeshKit
#include "meshkit/iGeom.hpp"
#include "meshkit/AF2LocalTransform.hpp"
#include "meshkit/Matrix.hpp"

class AF2PlaneProjection : public AF2LocalTransform
{
  private:

    iGeom* iGeomPtr;
    iGeom::EntityHandle surface;
    MeshKit::Vector<3> pOrigin;
    MeshKit::Vector<3> pNormal;
    MeshKit::Vector<3> pXDir;
    MeshKit::Vector<3> pYDir;
    double scale;

  public:

    /**
     * \brief Constructor
     *
     * It is the responsibility of the context that constructs/uses an
     * AF2PlaneProjection to ensure that the pointer to the iGeom instance
     * and the handle to the geometric surface remain valid as long as
     * the AF2PlaneProjection is in use.
     *
     * The normal vector must have length bounded away from zero.  Ideally
     * it would have length one, but the constructor will normalize it
     * to have length one if it does not.
     *
     * The x-direction vector must also have length bounded away from zero.
     * Ideally it would have length one, but the constructor will normalize it
     * to have length one if it does not.  The x-direction vector must be
     * orthogonal to the normal vector.
     *
     * The y-direction of the transformation will be defined such that
     * the right hand rule is followed.  In other words, the normal vector
     * to the plane will be the standard right-handed cross product of the
     * (normalized) x-direction vector and the y-direction.
     *
     * After projecting three-dimensional points into the plane, the
     * coordinates will be multiplied by the inverse of the scaling factor.
     * When used for the advancing front algorithm implementation, the
     * scale should be given as the approximate distance between two
     * adjacent points on the front near the baseline edge, since this
     * will produce points that are nearly unit distance apart after
     * the transformation is applied and the advancing front rules are
     * defined relative to a unit length baseline edge.
     *
     * \param iGeomPtrArg a pointer to an iGeom instance
     * \param srfcHandle a handle to a 2-dimensional surface embedded in
     *   3-dimensional space that is a valid handle to access the surface
     *   through the iGeom instance referenced by iGeomPtrArg
     * \param planeOrigin a point within the plane onto which this
     *   transformation will project points
     * \param planeNormal a direction vector that is normal to the plane
     *   onto which this transformation will project points 
     * \param planeXDir a direction vector (parallel to the plane onto which
     *   this transformation will project points) that defines the direction
     *   that will be the first coordinate of the 2-dimensional space
     * \param scaleFactor a positive scaling factor that approximates
     *   the distance between adjacent points on the advancing front
     */
    AF2PlaneProjection(iGeom* iGeomPtrArg,
        iGeom::EntityHandle srfcHandle,
        MeshKit::Vector<3> const & planeOrigin,
        MeshKit::Vector<3> const & planeNormal,
        MeshKit::Vector<3> const & planeXDir,
        double scaleFactor);

    AF2Point2D* transformFromSurface(AF2Point3D const & srfcPnt, bool & legal) const;

    AF2Point3D* transformToSurface(AF2Point2D const & planePnt,
        unsigned long const & pntId) const;
};

#endif
