/*
 * AF2LocalTransform.hpp
 *
 * A local transformation between some 2-dimensional space and 3-dimensional
 * space.  This transformation is used by the 2-dimensional advancing
 * front algorithm to transform points between a surface embedded in 3
 * dimensions and a 2-dimensional space.  The transformation should
 * be a bijection between the 2-dimensional space and some local patch
 * on the surface, but it does not need to be a global parametrization.
 */
#ifndef AF2LOCALTRANSFORM_HPP
#define AF2LOCALTRANSFORM_HPP

#include "meshkit/Matrix.hpp"

class AF2LocalTransform
{

  public:

    /**
     * \brief Transform from a 3-dimensional point on the surface to a
     * point in a 2-dimensional space.
     *
     * \param srfcPnt the input 3-dimensional point on the surface
     * \param planePnt a 2-dimensional point to hold the output coordinates
     */
    virtual void transformFromSurface(MeshKit::Vector<3> const & srfcPnt,
        MeshKit::Vector<2> & planePnt) const = 0;

    /**
     * \brief Transform from a point in the 2-dimensional space of this
     * transformation to a 3-dimensional point on the surface.
     *
     * \param planePnt the input 2-dimensional point
     * \param srfcPnt a 3-dimensional point to hold the output coordinates
     */
    virtual void transformToSurface(MeshKit::Vector<2> const & planePnt,
        MeshKit::Vector<3> & srfcPnt) const = 0;
};

#endif
