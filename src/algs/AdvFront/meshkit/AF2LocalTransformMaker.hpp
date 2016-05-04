/*
 * AF2LocalTransformMaker.hpp
 *
 * An AF2LocalTransformMaker is an object that can create an
 * AF2LocalTransform, given some information about the neighborhood
 * in which the local transformation will be applied.
 */
#ifndef AF2LOCALTRANSFORMMAKER_HPP
#define AF2LOCALTRANSFORMMAKER_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2LocalTransform.hpp"
#include "meshkit/AF2Edge3D.hpp"
#include "meshkit/AF2Point3D.hpp"

class AF2LocalTransformMaker
{

  public:

    /**
     * \brief Make a local transformation
     *
     * Make a local transformation to be used in a neighborhood, given
     * access to the information about the three-dimensional points and edges
     * that will be included in the neighborhood.
     *
     * The AF2LocalTransform returned from this method will be allocated
     * on the heap with the new operator.  The transform belongs to the
     * context that calls the method, and it is the responsibility
     * of the calling context to call delete.
     *
     * \param ngbhdPoints the points that will be included in the neighborhood
     * \param baselineEdge the baseline edge of the neighborhood
     * \param otherNgbhdEdges the other edges that will be included in
     *   the neighborhood
     * \return a local transformation appropriate for use with a neighborhood
     *   consisting of the specified points and edges
     */
    virtual AF2LocalTransform* makeLocalTransform(
        const std::list<AF2Point3D*> & ngbhdPoints,
        const AF2Edge3D* const & baselineEdge,
        const std::list<const AF2Edge3D*> & otherNgbhdEdges) const = 0;
};

#endif
