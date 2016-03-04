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
     * \return true if the free zone contains the point or nearly
     *   contains the point, false otherwise
     */
    bool nearContains(AF2Point2D const & testPnt) const;

    /**
     * \brief Check whether the free zone intersects (or nearly intersects)
     * a line segment given the two endpoints of the line segment
     *
     * \param startPoint one endpoint of the line segment
     * \param endPoint the other endpoint of the line segment
     * \return true if the free zone intersects or nearly intersects the
     *   line segment, false otherwise
     */
    bool nearIntersects(AF2Point2D const & startPoint,
        AF2Point2D const & endPoint) const;

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
