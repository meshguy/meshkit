/*
 * AF2PointTransform.hpp
 *
 */

#ifndef AF2POINTTRANSFORM_HPP
#define AF2POINTTRANSFORM_HPP

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2VertexBinding.hpp"

class AF2PointTransform
{
  public:

    /**
     * \brief Transform a point based on an AF2VertexBinding
     *
     * Given (the reference coordinates of) a 2-dimensional point and
     * a vertex binding that assigns actual coordinates to the reference
     * AF2RuleExistingVertex vertices of a rule, transform the point
     * from its reference location to a location that is appropriate
     * to the vertex binding.
     *
     * This method assumed that the vertex binding is complete, i.e.,
     * that every reference AF2RuleExistingVertex on the AF2VertexBinding
     * has been assigned actual coordinates.  The method is allowed to
     * modify the vertex binding to improve computational efficiency by
     * updating internal caches in the AF2VertexBinding, but is guaranteed
     * not to change the coordinates that reference vertices are bound to.
     *
     * \param point the coordinates of some 2-dimensional point
     * \param vBinding a binding of reference vertices to points with actual
     *   coordinates
     * \return a 2-dimensional point transformed from point based on the
     *   vertex binding
     */
    virtual AF2Point2D transformPoint(AF2Point2D const & point,
        AF2VertexBinding & vBinding) const = 0;

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
