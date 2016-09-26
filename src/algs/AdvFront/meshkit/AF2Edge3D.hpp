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

// C++
#include <list>

// MeshKit
#include "meshkit/AF2Point3D.hpp"

class QualityDecreaseObserver;

class AF2Edge3D
{
  private:

    AF2Point3D* startPnt;
    AF2Point3D* endPnt;
    unsigned int qualityLevel;
    QualityDecreaseObserver* observer;

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
    AF2Edge3D(AF2Point3D* start, AF2Point3D* end);

    /**
     * \brief Decrease the required quality of future attempts
     *   to advance the front using this edge as the baseline edge.
     *
     * This method increments the quality level, since a higher value
     * for the quality level indicates lower quality requirements.
     */
    void decreaseQuality();

    /**
     * \brief Get a pointer to the starting endpoint of the edge.
     */
    AF2Point3D* getStart() const;

    /**
     * \brief Get a pointer to the ending endpoint of the edge.
     */
    AF2Point3D* getEnd() const;

    /**
     * \brief Get the current quality level of this advancing front edge.
     *
     * Higher numbers for the quality level correspond to lower
     * requirements on the quality of advancing the front when
     * this edge is the baseline edge.
     *
     * \return the current quality level
     */
    unsigned int getQualityLevel() const;

    /**
     * \brief Set this edge's quality decrease observer.
     *
     * The intended use of this method is for the advancing front
     * to observe quality decreases to the edges that are on the
     * advancing front.  The method is public as an implementation
     * artifact, but users should not call this method.
     *
     * This object does not take ownership of the observer.
     *
     * \param observerArg an observer that desires to be notified
     *   when this edge's quality decreases
     */
    void setObserver(QualityDecreaseObserver* observerArg);
};

class QualityDecreaseObserver
{
  friend void AF2Edge3D::decreaseQuality();

  protected:
    /**
     * \brief Notify the observer that the quality has been decreased
     *   for the specified edge
     *
     * \param anEdge the edge for which the quality level number has increased
     */
    virtual void qualityDecreased(const AF2Edge3D* const & anEdge) = 0;
};

#endif
