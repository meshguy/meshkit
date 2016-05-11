/*
 * AF2Algorithm.hpp
 *
 * \brief An AF2Algorithm provides the capability to execute
 * a two-dimensional advancing front algorithm.
 *
 * // TODO: Complete the documentation
 */
#ifndef AF2ALGORITHM_HPP
#define AF2ALGORITHM_HPP

// C++
#include <list>
#include <map>

// MeshKit
#include "meshkit/AF2AlgorithmResult.hpp"
#include "meshkit/AF2LocalTransformMaker.hpp"
#include "meshkit/AF2Front.hpp"
#include "meshkit/AF2Neighborhood.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/AF2Polygon2D.hpp"
#include "meshkit/AF2Polygon3D.hpp"
#include "meshkit/AF2Rule.hpp"

class AF2Algorithm
{
  private:

    const std::list<const AF2Rule*> ruleList;

    /**
     * \brief Initialize the advancing front
     */
    void initFront(AF2Front & front, std::list<AF2Point3D*> & pntList,
        const double* coords, unsigned int numPoints,
        const unsigned int* edges, unsigned int numEdges) const;

    /**
     * \brief Process a new face that is being added to the advancing front
     */
    void processNewFace(const AF2Polygon2D* newFace2D,
        AF2Neighborhood* & ngbhd,
        std::map<const AF2Point2D*, AF2Point3D*> & newPointsMap,
        std::list<const AF2Polygon3D*> & allFaces, AF2Front & front) const;

    /**
     * \brief Process a new point that is being added to the advancing front
     */
    void processNewPoint(const AF2Point2D* newPoint2D,
        AF2Neighborhood* & ngbhd,
        std::map<const AF2Point2D*, AF2Point3D*> & newPointsMap,
        std::list<AF2Point3D*> & allPoints, AF2Front & front) const;

    /**
     * \brief Release the memory consumed by points and faces in the lists
     */
    void release(std::list<AF2Point3D*> & allPoints,
        std::list<const AF2Polygon3D*> & allFaces) const;

  public:

    /**
     * \brief Constructor
     *
     * // TODO:
     *
     * \param // TODO:
     */
    AF2Algorithm(const std::list<const AF2Rule*> & ruleListArg);

    /**
     * \brief Destructor
     */
    ~AF2Algorithm();

    /**
     * \brief Copy constructor (throws exception)
     *
     * This object does not currently support copying, so this method
     * is implemented to throw an exception.
     */
    AF2Algorithm(const AF2Algorithm & toCopy);

    /**
     * \brief Assignment operator (throws exception)
     *
     * This object does not currently support assignment, so this method
     * is implemented to throw an exception.
     */
    AF2Algorithm& operator=(const AF2Algorithm & rhs);

    /**
     * \brief Execute the advancing front algorithm
     *
     * TODO: Complete the documentation
     *
     * \param transformMaker An AF2LocalTransforMaker that encapsulates
     *   surface geometry information and is responsible for making local
     *   transformations between the three-dimensionsal surface geometry
     *   and two-dimensional coordinate spaces that locally approximate
     *   the surface
     * \param coords A pointer to an array of double values that contains
     *   the coordinates of mesh points that must be included in the mesh.
     *   The coordinates array should be ordered such that the x, y, and z
     *   coordinates of each point are listed together in that order.
     *   Coordinates of points that are endpoints of edges along the
     *   boundary must be included.  Points in the interior of the surface
     *   that must be in the mesh may also be included.
     * \param numPoints The number of points that are defined in the
     *   array of coordinates.  The coordinate array must have 3*numPoints
     *   double values (or more).
     * \param edges A pointer to an array of int values that define half-
     *   edges that must be included in the mesh.  Each edge is defined by
     *   two unsigned integer values.  The integer values are
     *   indices into the list of points defined in the coordinates
     *   array.  The index is 0-based, so the first three values of
     *   coords define point 0, the next three values of coords define
     *   point 1, etc.  The starting point of the edge is listed before the
     *   ending point of the edge.  When the edge is traversed from the
     *   starting point to the ending point, the surface that is to be
     *   meshed must be to the left.  Isolated edges in the interior of
     *   the surface, i.e., edges that have the surface both to the left
     *   and the right may be defined.  To require such an edge to exist
     *   in the mesh, both half-edges of the edge must be included in the
     *   array of edges.
     * \param numEdges The number of half-edges that are defined in
     *   the array of edges.  The array of edges must have 2*numEdges
     *   unsigned integer values (or more).
     */
    AF2AlgorithmResult* execute(
        const AF2LocalTransformMaker* const & transformMaker,
        const double* coords, unsigned int numPoints,
        const unsigned int* edges, unsigned int numEdges) const;
};

#endif
