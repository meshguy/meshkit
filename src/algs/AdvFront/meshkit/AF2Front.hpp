/*
 * AF2Front.hpp
 *
 * An AF2Front is an object that manages the front for a two-dimensional
 * advancing front algorithm.  It tracks which points and which half edges
 * are on the front.
 *
 * The main functions of an AF2Front are to initialize the front, to update
 * the front when a face is ready to be added, and to select a portion
 * of the front -- a neighborhood -- in which the algorithm should attempt
 * to advance the front.  The first two functions are accomplished by
 * the addPoint and advanceFront methods.  The third function is
 * implemented in the selectNeighborhood method.
 *
 * The AF2Front provides a method to check whether the front is emtpy,
 * which would indicate that the algorithm has completed if the front has
 * been in use.  It also tracks the quality of the edges that are on
 * the front and provides a convenient method to determine the quality
 * level of the highest quality edge that is on the front.
 */

#ifndef AF2FRONT_HPP
#define AF2FRONT_HPP

// C++
#include <list>
#include <map>
#include <set>
#include <vector>

// MeshKit
#include "meshkit/AF2Edge3D.hpp"
#include "meshkit/AF2LocalTransformMaker.hpp"
#include "meshkit/AF2Neighborhood.hpp"
#include "meshkit/AF2Point3D.hpp"

/**
 * An EndPointLess is a comparison object that compares AF2Edge3D pointers
 * as half edges according to the values of their endpoints.  This allows
 * different instances of AF2Edge3D that reference the same endpoints
 * (in the same order) to compare as equal.
 */
struct EndPointLess {
  /**
   * The comparison method.
   */
  bool operator()(const AF2Edge3D* const & oneEdge,
      const AF2Edge3D* const & otherEdge) const;
};

class AF2Front : QualityDecreaseObserver
{
  friend void AF2Edge3D::decreaseQuality();

  private:

    // a map from the points on the front to a count of the
    // number of front half edges that are incident to the point
    std::map<AF2Point3D*, int> points;
    // a set of half edges that are on the front
    std::set<AF2Edge3D*, EndPointLess> edges;
    // a vector tracking how many edges of each quality level are
    // currently on the front
    std::vector<unsigned int> qualityCount;

    /**
     * \brief Notify this advancing front that the quality has been decreased
     *   for the specified edge
     *
     * \param anEdge the edge for which the quality level number has increased
     */
    void qualityDecreased(const AF2Edge3D* const & anEdge);

  public:

    /**
     * \brief Constructor
     */
    AF2Front();

    /**
     * \brief Destructor
     *
     * Deletes any remaining edges that this front owns due to calls to
     * advanceFront.
     */
    virtual ~AF2Front();

    /**
     * \brief Copy constructor (throws exception)
     *
     * This object does not currently support copying, so this method
     * is implemented to throw an exception.
     */
    AF2Front(const AF2Front & toCopy);

    /**
     * \brief Assignment operator (throws exception)
     *
     * This object does not currently support assignment, so this method
     * is implemented to throw an exception.
     */
    AF2Front& operator=(const AF2Front & rhs);

    /**
     * \brief Add a point to the front.
     *
     * Adds a point to the advancing front.  The added point may be
     * an isolated point or an endpoint of an edge that the user is planning
     * to add to the front.
     *
     * This object does not take ownership of the point.  It is the
     * responsiblity of the context to manage the memory for the point
     * and keep the point available in memory as long as this AF2Front
     * is in memory.  This object requires permission to update the
     * point's distance to the boundary.
     *
     * \param pointPtr a pointer to the point that should be added to the
     *   front
     */
    void addPoint(AF2Point3D* pointPtr);

    /**
     * \brief Initialize or advance the front with a list of half edges.
     *
     * Each endpoint of each edge must have been added to the
     * to the advancing front already.  The half edges must be at quality
     * level 1, i.e., their quality must not have been decreased
     * after they were constructed.  If either of these conditions
     * is violated, an exception will be thrown.
     *
     * The half edges must be oriented such that the area that still
     * has to be meshed or the half edges that are already on the
     * advancing front are to the left when the edge is traversed from
     * its start toward its end.  Other than initialization, the targeted
     * use case is to advance the front with a list of edges that bound
     * a face.  In this use case the edges should be oriented such that
     * the bounded face is to the right.
     *
     * If the area that has to be meshed is on both sides of an edge,
     * both half edges of the edge should be passed to this method in the
     * same method call.
     *
     * If the opposite half edge of an edge in the list already exists
     * on the front, having been added by a previous call to advanceFront,
     * then that half edge will be removed from the front.  Also, if this
     * reduces the number of edges incident to a point on the front to zero,
     * that point will be removed from the front.  If the opposite half
     * edge of an edge in the list does not yet exist on the front, then
     * that half edge will be added to the advancing front.
     *
     * This object takes ownership of the edges that are passed to this
     * method and deletes the edges when they are no longer needed.  Edges
     * may be deleted during execution of this method.  It is assumed
     * that the edges were allocated with a call to new.
     *
     * This object manages the quality level of edges on the advancing
     * front, so it requires permission to modify the AF2Edge3D that
     * are passed to the method.
     *
     * \param edgeList a list of edges that will be used to initialize or
     *   advance the front
     */
    void advanceFront(std::list<AF2Edge3D*> edgeList);

    /**
     * \brief Check whether the front is empty
     *
     * The front is considered empty if it has no remaining points
     * and no remaining half edges
     *
     * \return true if the advancing front is empty, false if it is nonempty
     */
    bool isEmpty() const;

    /**
     * \brief Get the maximum quality of edges on the front.
     *
     * In other words, get the smallest integer value among the quality
     * level values for edges.
     *
     * If there are no edges on the front, this method will return
     * zero.
     */
    unsigned int getMaximumQuality() const;

    /**
     * \brief Select a baseline edge.  Then build and return the neighborhood
     *   of that baseline edge.
     *
     * The neighborhood returned from this method is allocated with a
     * call to new.  It is the responsibility of the user to delete the
     * neighborhood.  Since advancing the front may invalidate some of
     * the pointers to edges and points in a neighborhood, a neighborhood
     * returned by this method should be deleted before or in conjunction
     * with advancing the front.
     *
     * If there are no edges on the front, this method throws an exception.
     *
     * \param transformMaker an object responsible for making the local
     *   transform that will be used by the selected neighborhood
     * \return a neighborhood of an edge on the advancing front
     */
    AF2Neighborhood* selectNeighborhood(
        const AF2LocalTransformMaker* const & transformMaker) const;
};

#endif
