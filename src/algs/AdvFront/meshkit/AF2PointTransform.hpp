/*
 * AF2PointTransform.hpp
 *
 * An AF2PointTransform is an object that is able to modify the coordinates
 * of a reference point based on an AF2Binding.
 */

#ifndef AF2POINTTRANSFORM_HPP
#define AF2POINTTRANSFORM_HPP

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Binding.hpp"

class AF2PointTransform
{
  public:

    virtual ~AF2PointTransform() {}

    /**
     * \brief Transform a point based on an AF2Binding
     *
     * Given (the reference coordinates of) a 2-dimensional point and
     * a binding that assigns actual coordinates (AF2Point2D points) to
     * the reference AF2RuleExistingVertex vertices of a rule, transform
     * the point from its reference location to a location that is
     * appropriate to the binding.
     *
     * This method assumes that the binding is complete, i.e.,
     * that every reference AF2RuleExistingVertex on the AF2Binding
     * has been assigned an AF2Point2D.
     *
     * \param point the coordinates of some 2-dimensional point
     * \param vBinding a binding of reference vertices to points with actual
     *   coordinates
     * \return a 2-dimensional point transformed from the value of the first
     *   argument based on the binding
     */
    virtual AF2Point2D transformPoint(AF2Point2D const & point,
        AF2Binding const & vBinding) const = 0;

    /**
     * \brief Make and return an independent copy of yourself.
     *
     * This method supports making copies of the concrete derived class that
     * an AF2PointTransform pointer references.
     *
     * Deletion of the object pointed to by the pointer returned from this
     * method must be managed by the context that calls the method.
     */
    virtual AF2PointTransform* clone() const = 0;
};

#endif
