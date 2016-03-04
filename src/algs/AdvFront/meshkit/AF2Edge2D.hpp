/*
 * AF2Edge2D.hpp
 *
 * An immutable edge that connects two AF2Point2D points, one labeled as
 * the start of the edge, and the other labeled as the end of the edge.
 *
 * This class holds on to its endpoints with pointers, but it does not
 * manage the memory for its endpoints.  It is the reponsibility of the
 * context that uses AF2Edge2D to ensure that pointers to endpoints of edges
 * remain valid as long as they may be referenced by an AF2Edge2D instance.
 */
#ifndef AF2EDGE2D_HPP
#define AF2EDGE2D_HPP

#include "meshkit/AF2Point2D.hpp"

class AF2Edge2D
{
  private:

    const AF2Point2D* startPnt;
    const AF2Point2D* endPnt;

  public:

    /**
     * \brief Constructor
     *
     * Construct an edge between the two specified points.
     * The context is responsible for maintaining the validity
     * of the pointers to the points as long as this Edge2D
     * (or somy copy of it) is in use.
     *
     * \param start the starting endpoint of the edge
     * \param end the ending endpoint of the edge.
     */
    AF2Edge2D(const AF2Point2D* start, const AF2Point2D* end);

    /**
     * \brief Get a pointer to the starting endpoint of the edge.
     */
    const AF2Point2D* getStart() const;

    /**
     * \brief Get a pointer to the ending endpoint of the edge.
     */
    const AF2Point2D* getEnd() const;
};

#endif
