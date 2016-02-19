/*
 * AF2PointTransformNone.hpp
 *
 * \brief A point transformation that does no actual transformation
 * to the point.
 *
 * A copy of the point passed in is returned without modification;  the
 * transformation does not make use of or depend on the vertex binding.
 */

#ifndef AF2POINTTRANSFORMNONE_HPP
#define AF2POINTTRANSFORMNONE_HPP

// MeshKit
#include "meshkit/AF2PointTransform.hpp"

class AF2PointTransformNone : AF2PointTransform
{
  public:

    /**
     * \brief Implements AF2PointTransform::transformPoint
     *
     * See the documentation at AF2PointTransform::transformPoint
     * for the general contract for this method.  In this implementation
     * there is no actual transformation.  The point returned by the
     * method is at the same location as the point passed into the
     * method regardless of the vertex binding.
     *
     * \param point the coordinates of some 2-dimensional point
     * \param vBinding a binding of reference vertices to points with actual
     *   coordinates
     */
    virtual MeshKit::Vector<2> transformPoint(MeshKit::Vector<2> const & point,
        AF2VertexBinding & vBinding) const;

    /**
     * \brief Makes and returns an independent copy of this
     * AF2PointTransformNone.
     *
     * Implements AF2PointTransform::clone.  See additional documentation
     * there.
     */
    virtual AF2PointTransformNone* clone() const;
};

#endif
