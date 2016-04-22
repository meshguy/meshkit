/*
 * AF2LocalTransform.hpp
 *
 * A local transformation between some 2-dimensional space and some
 * 2-dimensional subspace of a 3-dimensional space.  This transformation
 * is used by the 2-dimensional advancing front algorithm to transform
 * points between a surface embedded in 3 dimensions and a 2-dimensional
 * space.  The transformation should be a bijection between the
 * 2-dimensional space and some local patch on the surface, but
 * it does not need to be a global parametrization.
 */
#ifndef AF2LOCALTRANSFORM_HPP
#define AF2LOCALTRANSFORM_HPP

#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Point3D.hpp"

class AF2LocalTransform
{

  public:

    /**
     * \brief Transform from a 3-dimensional point on the surface to a
     * point in a 2-dimensional space.
     *
     * The returned point is returned by pointer.  It is allocated
     * on the heap by this method using new, and it is the responsibility
     * of the calling context to deallocate it with a call to delete.
     *
     * \param srfcPnt the input 3-dimensional point on the surface
     * \return a 2-dimensional point
     */
    virtual AF2Point2D* transformFromSurface(
        AF2Point3D const & srfcPnt) const = 0;

    /**
     * \brief Transform from a point in the 2-dimensional space of this
     * transformation to a 3-dimensional point on the surface.
     *
     * The returned point is returned by pointer.  It is allocated
     * on the heap by this method using new, and it is the responsibility
     * of the calling context to deallocate it with a call to delete.
     *
     * \param planePnt the input 2-dimensional point
     * \return a 3-dimensional point on the surface
     */
    virtual AF2Point3D* transformToSurface(
        AF2Point2D const & planePnt) const = 0;
};

#endif
