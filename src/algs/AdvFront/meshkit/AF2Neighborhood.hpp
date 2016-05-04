/*
 * AF2Neighborhood.hpp
 *
 * An AF2Neighborhood represents a collection of points and edges
 * on the advancing front that are in the neighborhood of some specific
 * edge on the advancing front.  The specific edge is called the
 * baseline edge.  The advancing front algorithm uses an AF2Neighborhood
 * to assess whether it should add one or more additional faces and
 * perhaps also points and edges built off of the neighborhood's
 * baseline edge.
 *
 * The AF2Neighborhood manages conversion between an assumed
 * two-dimensional surface embedded in a three-dimensional
 * coordinate system and a local two-dimensional coordinate system
 * that is appropriate to the neighborhood.
 */
#ifndef AF2NEIGHBORHOOD_HPP
#define AF2NEIGHBORHOOD_HPP

// C++
#include <list>
#include <map>

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Point3D.hpp"
#include "meshkit/AF2Edge2D.hpp"
#include "meshkit/AF2Edge3D.hpp"
#include "meshkit/AF2LocalTransform.hpp"

class AF2Neighborhood
{
  private:

    const AF2LocalTransform* localTransform;
    std::list<const AF2Point2D*> points2D;
    const AF2Edge2D* baseEdge2D;
    std::list<const AF2Edge2D*> edges2D;
    std::map<AF2Point3D*, const AF2Point2D*> map3DTo2D;
    std::map<const AF2Point2D*, AF2Point3D*> map2DTo3D;

  public:

    /**
     * \brief Constructor
     *
     * Constructs a neighborhood, given the 3-dimensional points, the
     * 3-dimensional baseline edge, and the other 3-dimensional edges
     * that make up the neighborhood, along with a local transform
     * to manage transforming between the embedding of the surface
     * in three-dimensions and a local two-dimensional coordinate system
     * that is appropriate to the neighborhood.
     *
     * The baseline edge need not be listed in the list of other edges
     * and will be ignored if it is in that list.
     *
     * None of the pointers passed to the constructor may be null or invalid.
     *
     * This object takes ownership of the local transform that is passed
     * into the constructor.  The local transform will be deleted when
     * this object is deleted.  It is assumed that the local transform
     * was allocated with a call to new.
     *
     * This object does not take ownership of any of the three-dimensional
     * points or edges that are passed into the method.  The context that
     * calls the method is responsible for the memory of those objects.
     *
     * \param points3D a list of three-dimensional points that are
     *   in the neighborhood
     * \param baseLineEdge the three-dimensional baseline edge that
     *   is the basis of the neighborhood
     * \param edges3D a list of other three-dimensional edges that
     *   are in the neighborhod
     * \param localTransformArg a local transformation between points
     *   in the three-dimensional coordinate space that contains the
     *   neighborhood points and edges and an appropriate
     *   two-dimensional coordinate space
     */
    AF2Neighborhood(const std::list<AF2Point3D*> & points,
        const AF2Edge3D* baselineEdge,
        const std::list<const AF2Edge3D*> & otherEdges,
        const AF2LocalTransform* localTransformArg);

    /**
     * \brief Destructor
     */
    ~AF2Neighborhood();

    /**
     * \brief Copy constructor is not supported, so this throws an exception
     */
    AF2Neighborhood(const AF2Neighborhood & toCopy);

    /**
     * \brief Assignment operator is not supported, so this throws an exception
     */
    AF2Neighborhood & operator=(const AF2Neighborhood & rhs);

    /**
     * \brief Get the neighborhood's baseline edge
     *
     * Get the baseline edge of the neighborhood as it exists
     * in the two-dimensional coordinate system transformation
     * of the neighborhood.
     */
    const AF2Edge2D* getBaselineEdge2D() const;

    /**
     * \brief Get the three-dimensional point (if any) in the neighborhood
     *   that corresponds to the specified two-dimensional point
     *
     * If the two-dimensional point is a point in the neighborhood, then
     * this method will return a pointer to the three-dimensional point
     * in the neighborhood that generated the two-dimensional point.
     * If the two-dimensional point is not a point in the neighborhood, then
     * there will not be a corresponding three-dimensional point in the
     * neighborhood, and this method will return a NULL pointer.
     *
     * \param ngbhdPoint2D a two-dimensional point that may or may not be
     *   from the neighborhood
     * \return a three-dimensional point corresponding to the two-dimensional
     *   point, if there is one, or NULL
     */
    AF2Point3D* getCorrespondingPoint(
        const AF2Point2D* const & ngbhdPoint2D) const;

    /**
     * \brief Get all of the two-dimesional edges in the neighborhood
     *
     * Get the edges of the neighborhood as they have been transformed
     * into the two-dimensional coordinate system of the neighborhood.
     * The baseline edge will be the first edge in the list.
     */
    const std::list<const AF2Edge2D*>* getEdges2D() const;

    /**
     * \brief Get the two-dimesional points in the neighborhood
     *
     * Get the points of the neighborhood as they have been transformed
     * into the two-dimensional coordinate system of the neighborhood.
     */
    const std::list<const AF2Point2D*>* getPoints2D() const;

    /**
     * \brief Transform the specified two-dimensional point into a
     *   three-dimensional point using this neighborhood's local transform
     *
     * Without regard for whether the two-dimensional point is in this
     * neighborhood, this method applies the local transformation to the
     * two-dimensional point to construct a three-dimensional point.  The
     * three-dimensional point returned from this method is allocated
     * using the new operator.  It belongs to the calling context and
     * it is the responsibility of that context to call delete.
     *
     * A new point is allocated each time this method is called regardless
     * of whether the two-dimensional point has been passed to this
     * method previously.
     *
     * \param point2D a two-dimensional point
     * \return the local transformation of the two-dimensional point into
     *   a three-dimensional point
     */
    AF2Point3D* transformPoint(const AF2Point2D* const & point2D) const;
};

#endif
