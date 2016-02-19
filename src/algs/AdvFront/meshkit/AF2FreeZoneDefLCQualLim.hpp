/*
 * AF2FreeZoneDefLCQualLim.hpp
 *
 * A free zone definition that uses the quality level of the rule
 * application to scale the location of the points from a preferred
 * location toward some limiting acceptable location.
 *
 * When applying this free zone definition at any quality level,
 * the preferred point locations and limiting point locations
 * (for the vertex binding that maps the rule's reference vertices
 * to actual coordinates) are computed by applying point transforms
 * to reference preferred and limiting locations.
 * After that, the actual location for each boundary point of the free
 * zone at the given quality level by takin a linear combination
 * of the preferred location for the point and the limiting acceptable
 * location for the point.  The coefficients of the linear combination
 * are defined such that as the qualityClass number q increases (i.e.,
 * the quality decreases), the preferred location's coefficient is 1/q.
 */

#ifndef AF2FREEZONEDEFLCQUALLIM_HPP
#define AF2FREEZONEDEFLCQUALLIM_HPP

// C++
#include <list>

// MeshKit
#include "meshkit/AF2FreeZoneDef.hpp"
#include "meshkit/AF2PointTransform.hpp"
#include "meshkit/Matrix.hpp"

class AF2FreeZoneDefLCQualLim : public AF2FreeZoneDef
{
  private:

    int numPoints;
    MeshKit::Vector<2>* prefBndryPoints;
    MeshKit::Vector<2>* limBndryPoints;
    const AF2PointTransform** prefPointTransforms;
    const AF2PointTransform** limPointTransforms;

  public:

    /**
     * \brief Constructor
     *
     * The constructor requires a list of preferred boundary point
     * locations, a list of point transforms for the preferred point
     * locations, a list of limiting boundary point locations, and a list
     * of point transforms for the limiting point locations.  The lists
     * must all be the same length.
     * The index within the four lists is used to associate the
     * preferred location of a point with its point transform, 
     * the limiting location of the point, and the limiting location's
     * point transform.
     *
     * Points should be listed in order such that after the transforms
     * and linear combination
     * are applied they will define a counterclockwise traversal of the
     * vertices of a convex polygon.
     *
     * This class will clone the point transforms that are passed into
     * the method, so the calling context is responsible for managing
     * memory for the AF2PointTransform objects in the lists that are
     * passed into this method.
     *
     * \param preferBndryPnts the preferred locations of the boundary points
     *   of the free zone relative to vertices at their reference positions
     * \param preferPntTrnsfrms the point transformations for the preferred
     *   point locations that will transform the points to locations
     *   appropriate to the actual vertex binding
     * \param limitBndryPnts the limiting locations of the boundary points
     *   of the free zone relative to vertices at their reference positions
     * \param limitPntTrnsfrms the point transformations for the limiting
     *   point locations that will transform the points to locations
     *   appropriate to the actual vertex binding
     */
    AF2FreeZoneDefLCQualLim(
        std::list<MeshKit::Vector<2> > const & preferBndryPnts,
        std::list<const AF2PointTransform*> const & preferPntTrnsfrms,
        std::list<MeshKit::Vector<2> > const & limitBndryPnts,
        std::list<const AF2PointTransform*> const & limitPntTrnsfrms);

    ~AF2FreeZoneDefLCQualLim();

    /**
     * \brief Copy constructor
     */
    AF2FreeZoneDefLCQualLim(const AF2FreeZoneDefLCQualLim & toCopy);

    /**
     * \brief Assignment operator
     */
    AF2FreeZoneDefLCQualLim& operator=(const AF2FreeZoneDefLCQualLim & rhs);

    AF2FreeZoneDefLCQualLim* clone() const;

    /**
     * \brief Implement the makeFreeZone method defined in the parent
     * pure virtual class AF2FreeZoneDef.
     *
     * See the documentation in the parent class for the general contract
     * of the makeFreeZone method.
     *
     * This implementation starts by using the vertex binding to compute
     * the preferred location of each free zone boundary point and the
     * limiting location of each free zone boundary point. Then, given
     * a quality class q, this implementation takes a linear combination
     * of (1/q) times each point's preferred location with (1 - (1/q))
     * times the point's limiting location to compute the actual
     * coordinates that will be used for the free zone boundary point.
     */
    AF2FreeZone* makeFreeZone(
        AF2VertexBinding & vertexBinding, int qualityClass) const;
};

#endif
