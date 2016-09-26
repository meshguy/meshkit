/*
 * AF2Polygon2D.hpp
 *
 * An immutable 2-dimensional polygon whose vertices are AF2Point2D points.
 *
 * The vertices are numbered 0 through n - 1, where n is the number of
 * polygon vertices.  As the numbers increase, the polygon is traversed
 * counterclockwise.  Thus the interior of the polygon is to the left
 * of the edge from vertex i to vertex i + 1 modulo n.
 *
 * This class holds on to its vertices with pointers, but it does not
 * manage the memory for its vertices.  It is the reponsibility of the
 * context that uses AF2Polygon2D to ensure that pointers to vertices
 * remain valid as long as they may be referenced by an AF2Polygon2D instance.
 */
#ifndef AF2POLYGON2D_HPP
#define AF2POLYGON2D_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2Point2D.hpp"

class AF2Polygon2D
{
  private:

    unsigned int numVertices;
    const AF2Point2D** vertices;

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
     * \param polygonVertices a list of AF2Point2D* pointing to the vertices
     *     of the polygon traversed in counterclockwise order
     */
    AF2Polygon2D(std::list<const AF2Point2D*> const & polygonVertices);

    /**
     * \brief Destructor
     */
    ~AF2Polygon2D();

    /**
     * \brief Copy constructor
     *
     * This is the standard copy constructor.
     *
     * \param toCopy an AF2Polygon2D that should be copied to construct
     *   a new AF2Polygon2D
     */
    AF2Polygon2D(const AF2Polygon2D & toCopy);

    /**
     * \brief Copy assignment
     *
     * This is the standard copy assignment operator.
     *
     * \param rhs an AF2Polygon2D that should be assigned to this AF2Polygon2D,
     *   overwriting (and destructing) whatever may be in this AF2Polygon2D
     */
    AF2Polygon2D& operator=(const AF2Polygon2D & rhs);

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
     * When pointers returned from this method point to vertices that
     * belong to an AF2RuleApplication, i.e., new vertices that are
     * created as part of applying a rule, the user should be aware that
     * copying the AF2RuleApplication creates new instances of the
     * vertices, and the pointers will not reference the instances in
     * the copy.
     *
     * \param vtxNum the number of the vertex to return
     * \return a pointer to a vertex of the polygon
     */
    const AF2Point2D* getVertex(unsigned int vtxNum) const;
};

#endif
