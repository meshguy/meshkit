/*
 * AF2Edge3D.hpp
 *
 * An edge (or half-edge) that connects two AF2Point3D points, one labeled
 * as the start of the edge, and the other labeled as the end of the edge.
 *
 * This class holds on to its endpoints with pointers, but it does not
 * manage the memory for its endpoints.  It is the reponsibility of the
 * context that uses AF2Edge3D to ensure that pointers to endpoints of edges
 * remain valid as long as they may be referenced by an AF2Edge3D instance.
 */
#ifndef AF2EDGE3D_HPP
#define AF2EDGE3D_HPP

#include "meshkit/AF2Point3D.hpp"

class AF2Edge3D
{
  private:

    const AF2Point3D* startPnt;
    const AF2Point3D* endPnt;

  public:

    /**
     * \brief Constructor
     *
     * Construct an edge between the two specified points.
     * The context is responsible for maintaining the validity
     * of the pointers to the points as long as this AF2Edge3D
     * (or somy copy of it) is in use.
     *
     * \param start the starting endpoint of the edge
     * \param end the ending endpoint of the edge.
     */
    AF2Edge3D(const AF2Point3D* start, const AF2Point3D* end);

    /**
     * \brief Get a pointer to the starting endpoint of the edge.
     */
    const AF2Point3D* getStart() const;

    /**
     * \brief Get a pointer to the ending endpoint of the edge.
     */
    const AF2Point3D* getEnd() const;
};

#endif
