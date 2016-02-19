/*
 * AF2PlaneProjection.hpp
 *
 * A transformation between 2-dimensional space and 3-dimensional space
 * that transforms from three dimensions to two dimensions by subtracting
 * some origin point from the 3-dimensional point and removing the
 * component that is in the direction normal to the plane.  In other words,
 * the transformation projects the 3-dimensional point onto a plane.
 * 
 * The inverse works by firing a ray from the plane in the direction
 * normal to the plane and finding the intersection of the ray with
 * the surface.
 *
 * NOTE: In the future this could include scaling in the x-direction,
 * but at the moment it does not support that.
 */

#ifndef AF2PLANEPROJECTION_HPP
#define AF2PLANEPROJECTION_HPP

#include "meshkit/iGeom.hpp"
#include "meshkit/AF2LocalTransform.hpp"

class AF2PlaneProjection : public AF2LocalTransform
{
  private:

    iGeom* iGeomPtr;
    iGeom::EntityHandle surface;
    MeshKit::Vector<3> pOrigin;
    MeshKit::Vector<3> pNormal;
    MeshKit::Vector<3> pXDir;
    MeshKit::Vector<3> pYDir;

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
     * (normaliezed) x-direction vector and the y-direction.
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
     */
    AF2PlaneProjection(iGeom* iGeomPtrArg,
        iGeom::EntityHandle srfcHandle,
        MeshKit::Vector<3> const & planeOrigin,
        MeshKit::Vector<3> const & planeNormal,
        MeshKit::Vector<3> const & planeXDir);

    void transformFromSurface(MeshKit::Vector<3> const & srfcPnt,
        MeshKit::Vector<2> & planePnt) const;

    void transformToSurface(MeshKit::Vector<2> const & planePnt,
        MeshKit::Vector<3> & srfcPnt) const;
};

#endif
