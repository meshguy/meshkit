/*
 * AF2FreeZone.hpp
 *
 * A region of 2-dimensional space, supposed to be convex, defined by
 * the points of which the region is supposed to be the convex hull,
 * listed in a counter-clockwise order around the region.  This region
 * can be checked (1) to verify that the points defining it are, indeed,
 * all on the convex hull of the region, (2) to see whether it contains 
 * other points, and (3) to see whether it intersects an edge defined
 * by its two endpoints.
 *
 * Within an advancing front algorithm, a free zone may be used to
 * determine whether insertion of vertices and edges can be accomplished
 * without forcing inverted or poor quality elements to appear in the mesh.
 */
#ifndef AF2FREEZONE_HPP
#define AF2FREEZONE_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2Point2D.hpp"

class AF2FreeZone
{
  private:

    // the vertices on the convex hull of the free zone,
    // in counterclockwise order
    std::list<AF2Point2D> vertices;

    // the bounding box of the free zone
    double minX;
    double maxX;
    double minY;
    double maxY;

    // approximate scale of the coordinates that define the free zone
    double scale;

    /**
     * \brief Check whether two points are equal or nearly equal
     *
     * Near equality is measured relative to the scale of the free zone.
     *
     * \param pointAlpha one of the two points to test
     * \param pointBravo the other of the two points to test
     * \return true if the points are equal or nearly equal, false otherwise
     */
    bool nearEqual(AF2Point2D const & pointAlpha,
        AF2Point2D const & pointBravo) const;

  public:

    /**
     * \brief Constructor
     */
    AF2FreeZone(std::list<AF2Point2D> const& bndryPoints);

    /**
     * \brief Check whether the free zone contains (or nearly contains)
     * a specified point
     *
     * \param point the point to check
     * \param containsBndry if false, the method treats query points that
     *   are nearly equal to a free zone boundary point as not contained,
     *   even if they are actually mathematically on or in the free zone
     * \return true if the free zone contains the point or nearly
     *   contains the point, false otherwise
     */
    bool nearContains(AF2Point2D const & testPnt,
        bool const & containsBndry = false) const;

    /**
     * \brief Check whether the free zone intersects (or nearly intersects)
     * a line segment given the two endpoints of the line segment
     *
     * \param startPoint one endpoint of the line segment
     * \param endPoint the other endpoint of the line segment
     * \param containsBndry if false, the method treats query line segments
     *   that are nearly equal to free zone boundary line segments as not
     *   intersecting, even if they actually mathematically do intersect
     *   the free zone or lie on its boundary.  This applies only if the
     *   query line segment is nearly equal to the full boundary line
     *   segment, i.e., if each endpoint of the query line segment
     *   is nearly equal to one of a pair of consective free zone
     *   boundary points.  Also, if false, a query line segment
     *   that has one endpoint that is nearly equal to a free
     *   zone boundary point will be treated as nonintersecting if
     *   the only intersection or near intersection is in the
     *   neighborhood of that free zone boundary point.
     * \return true if the free zone intersects or nearly intersects the
     *   line segment, false otherwise
     */
    bool nearIntersects(AF2Point2D const & startPoint,
        AF2Point2D const & endPoint, bool const & containsBndry = false) const;

    /**
     * \brief Verify that the free zone is convex, and that the
     * boundary points are given in a counter-clockwise ordering.
     *
     * If the free zone is not actually convex or the boundary points
     * were not given in a counter-clockwise ordering, the behavior
     * of other methods in this class is not defined.
     *
     * \return true if the free zone is convex, false otherwise
     */
    bool isConvex() const;
};

#endif
