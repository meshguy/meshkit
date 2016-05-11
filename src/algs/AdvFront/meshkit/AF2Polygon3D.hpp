/*
 * AF2Polygon3D.hpp
 *
 * \brief An immutable not-necessarily flat polygon whose vertices are
 *   AF2Point3D points.
 *
 * Perhaps a better name would be a closed wire rather than a polygon,
 * for the way it is represented here, but the object is defined in
 * terms of its vertices and is intended to correspond to a
 * two-dimensional element of a mesh.  This is also a close analog
 * of the AF2Polygon2D, except that it has three-dimensional vertices.
 *
 * The vertices are numbered 0 through n - 1, where n is the number of
 * polygon vertices.  As the numbers increase, the polygon is traversed
 * counterclockwise.  Thus the interior of the polygon is to the left
 * of the edge from vertex i to vertex i + 1 modulo n.
 *
 * This class holds on to its vertices with pointers, but it does not
 * manage the memory for its vertices.  It is the reponsibility of the
 * context that uses AF2Polygon3D to ensure that pointers to vertices
 * remain valid as long as they may be referenced by an AF2Polygon3D instance.
 */
#ifndef AF2POLYGON3D_HPP
#define AF2POLYGON3D_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2Point3D.hpp"

class AF2Polygon3D
{
  private:

    unsigned int numVertices;
    const AF2Point3D** vertices;

  public:

    /**
     * \brief Constructor
     *
     * The vertices passed into the constructor are pointers to
     * polygon vertices.  The order of the polygon vertices within
     * the list must be counterclockwise.  In other words, the
     * interior of the polygon must be to the left of the edge
     * from each vertex to the next vertex in the list.
     *
     * This class does not take ownership of the vertices that are
     * passed into it (by pointer), and it is the reponsibility of the
     * context that uses AF2Polygon3D to ensure that pointers to vertices
     * remain valid as long as they are referenced by this object or
     * a copy of it.
     *
     * \param polygonVertices a list of AF2Point3D* pointing to the vertices
     *     of the polygon traversed in counterclockwise order
     */
    AF2Polygon3D(std::list<const AF2Point3D*> const & polygonVertices);

    /**
     * \brief Destructor
     */
    ~AF2Polygon3D();

    /**
     * \brief Copy constructor
     *
     * This is the standard copy constructor.
     *
     * \param toCopy an AF2Polygon3D that should be copied to construct
     *   a new AF2Polygon3D
     */
    AF2Polygon3D(const AF2Polygon3D & toCopy);

    /**
     * \brief Copy assignment
     *
     * This is the standard copy assignment operator.
     *
     * \param rhs an AF2Polygon3D that should be assigned to this AF2Polygon3D,
     *   overwriting (and destructing) whatever may be in this AF2Polygon3D
     */
    AF2Polygon3D& operator=(const AF2Polygon3D & rhs);

    /**
     * \brief Get the number of vertices that this polygon has
     *
     * \return the number of vertices that this polygon has
     */
    unsigned int getNumVertices() const;

    /**
     * \brief Get a vertex of this polygon
     *
     * The numbering of the vertices of the polygon begins with 0, so
     * the valid arguments to this method are 0 through n - 1, where
     * n is the value returned from getNumVertices().
     *
     * \param vtxNum the number of the vertex to return
     * \return a pointer to a vertex of the polygon
     */
    const AF2Point3D* getVertex(unsigned int vtxNum) const;
};

#endif
