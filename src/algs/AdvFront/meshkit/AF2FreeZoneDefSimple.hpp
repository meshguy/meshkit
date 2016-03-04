/*
 * AF2FreeZoneDefSimple.hpp
 *
 * A free zone definition that does not depend on the quality level
 * of the rule application.  This free zone definition consists of
 * the reference locations of the free zone boundary vertices and point
 * transforms that will transform the free zone based on the vertex
 * binding that maps reference vertices to actual existing vertices.
 */

#ifndef AF2FREEZONEDEFSIMPLE_HPP
#define AF2FREEZONEDEFSIMPLE_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2FreeZoneDef.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2PointTransform.hpp"

class AF2FreeZoneDefSimple : public AF2FreeZoneDef
{
  private:

    int numPoints;
    AF2Point2D* bndryPoints;
    const AF2PointTransform** pointTransforms;

  public:

    /**
     * \brief Constructor
     *
     * The constructor requires a list of boundary points and a list
     * of point transforms.  The lists must be the same length.  Each
     * point transform will be applied to the boundary point at the
     * corresponding position in the list of boundary points.
     *
     * Points should be listed in order such that after the transforms
     * are applied they will define a counterclockwise traversal of the
     * vertices of a convex polygon.
     *
     * This class will clone the point transforms that are passed into
     * the method, so the calling context is responsible for managing
     * memory for the AF2PointTransform objects in the list that is
     * passed into this method.
     *
     * \param rfrncBndryPnts the reference locations of the boundary points
     *   of the free zone, i.e., the locations that the boundary points
     *   will have if the existing points of the rule are all
     *   placed in their ideal positions
     * \param bndryPntTrnsfrms the list of point transformations
     */
    AF2FreeZoneDefSimple(
        std::list<AF2Point2D> const & rfrncBndryPnts,
        std::list<const AF2PointTransform*> const & bndryPntTrnsfrms);

    ~AF2FreeZoneDefSimple();

    /**
     * \brief Copy constructor
     */
    AF2FreeZoneDefSimple(const AF2FreeZoneDefSimple & toCopy);

    /**
     * \brief Assignment operator
     */
    AF2FreeZoneDefSimple& operator=(const AF2FreeZoneDefSimple & rhs);

    AF2FreeZoneDefSimple* clone() const;

    /**
     * \brief Implement the makeFreeZone method defined in the parent
     * pure virtual class AF2FreeZoneDef.
     *
     * See the documentation in the parent class for the general contract
     * of the makeFreeZone method.
     *
     * This implementation works by taking each reference boundary point
     * in the list of points passed to the constructor and transforming
     * it to a different location using the corresponding point transform
     * passed to the constructor together with the vertex binding passed
     * into this method.  The quality class is ignored in this implementation.
     */
    AF2FreeZone* makeFreeZone(
        AF2VertexBinding & vertexBinding, int qualityClass) const;
};

#endif
