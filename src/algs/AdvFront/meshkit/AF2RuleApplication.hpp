/*
 * AF2RuleApplication.hpp
 *
 * An AF2RuleApplication is a representation of what is involved in applying
 * a particular AF2Rule is applied with a particular AF2Binding to points
 * and edges in some AF2Neighborhood.  The representation is independent
 * of the particular AF2Rule and AF2Binding so that it can be easily
 * compared with other AF2RuleApplications that involve other rules
 * or bind to other points and edges in the AF2Neighborhood.
 * 
 * The AF2RuleApplication tells which new points and new faces are to
 * be added by applying the particular rule with the particular binding.
 * The AF2Point2D objects referenced by the AF2Polygon2D faces are
 * either from the AF2Neighborhood or from the list of new points.  Since
 * the points are copied when the AF2RuleApplication is copied, a reference
 * to a point from a polygon from a AF2RuleApplication will not match
 * any AF2Point2D reference returned by the getNewPoint method of
 * a copy of the AF2RuleApplication.
 */
#ifndef AF2RULEAPPLICATION_HPP
#define AF2RULEAPPLICATION_HPP

// C++
#include <list>
#include <map>

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2Polygon2D.hpp"

class AF2RuleApplication
{
  private:

    unsigned int numNewPoints;
    const AF2Point2D** newPoints;
    unsigned int numNewFaces;
    const AF2Polygon2D** newFaces;

    /**
     * \brief construct a polygon with vertices mapped from a source
     *   polygon according to the defined map
     *
     * For use by constructors and the copy assignment operator.
     *
     * \param sourcePolygon the source polygon that is to be copied
     * \param targetPolygon a (currently null) pointer to the location
     *   in memory where the copy should be constructed using new
     * \param a map of the new vertices from the rule application instance
     *   that provided the context of the source polygon to the
     *   corresponding new vertices from this rule application instance
     */
    void mapCopyPolygon(const AF2Polygon2D & sourcePolygon,
        const AF2Polygon2D* & targetPolygon,
        const std::map<const AF2Point2D*, const AF2Point2D*> & vertexMap);

  public:

    /**
     * \brief Constructor
     *
     * Construct an AF2RuleApplication, defining the (two-dimensional)
     * points, if any, that would be added and the polygonal faces
     * that would be added by applying some particular AF2Rule with
     * some AF2Binding of the rule's required existing vertices and
     * edges to the two-dimensional points and edges in some AF2Neighborhood.
     * The vertices of the polygonal faces are assumed to be points that
     * come from either the AF2Neighborhood or the list of new points
     * that is passed to this constructor.
     *
     * Objects that are passed into this constructor by pointer are
     * copied, so the original instances are owned by the context
     * that calls the constructor and their memory should be managed
     * outside of this instance.
     *
     * \param newPointsList A list of pointers to AF2Point2D points
     *   that should be projected onto the surface and added to the
     *   mesh as part of accepting this AF2RuleApplication
     * \param newFacesList A list of polygonal faces (referencing
     *   two-dimensional points from the newPointsList or from
     *   an AF2Neighborhood) that should be added to the mesh
     *   as part of accepting this AF2RuleApplication
     */
    AF2RuleApplication(std::list<const AF2Point2D*> const & newPointsList,
        std::list<const AF2Polygon2D*> const & newFacesList);

    /**
     * \brief Destructor
     */
    ~AF2RuleApplication();

    /**
     * \brief Copy constructor
     *
     * This is the standard copy constructor.
     *
     * \param toCopy an AF2RuleApplication that should be copied to construct
     *   a new AF2RuleApplication
     */
    AF2RuleApplication(const AF2RuleApplication & toCopy);

    /**
     * \brief Copy assignment
     *
     * This is the standard copy assignment operator.
     *
     * \param rhs an AF2RuleApplication that should be assigned
     *   to this AF2RuleApplication, overwriting (and destructing)
     *   whatever may be in this AF2RuleApplication
     */
    AF2RuleApplication& operator=(const AF2RuleApplication & rhs);

    /**
     * \brief Get the number of new faces that would be produced by
     *   this rule application.
     *
     * \return the number of new faces
     */
    unsigned int getNumNewFaces() const;

    /**
     * \brief Get a pointer to one of the new faces that would be
     *   added by this rule application
     *
     * The numbering of the new faces that would be produced begins with 0,
     * so the valid arguments to this method are 0 through n - 1, where
     * n is the value returned from getNumNewFaces().
     *
     * \param newFaceIndex the number of the new face to return
     * \return a pointer to the requested new face that would be added
     *   by this rule application
     */
    const AF2Polygon2D* getNewFace(unsigned int newFaceIndex) const;

    /**
     * \brief Get the number of new points that would be added by
     *   this rule application.
     *
     * \return the number of new vertices
     */
    unsigned int getNumNewPoints() const;

    /**
     * \brief Get a pointer to one of the new (two-dimensional) points
     *   that would be added by this rule application
     *
     * The numbering of the new points that would be added begins with 0,
     * so the valid arguments to this method are 0 through n - 1, where
     * n is the value returned from getNumNewPoints().
     *
     * \param newPointIndex the number of the new point to return
     * \return a pointer to the requested new point that would be added
     *   by this rule application
     */
    const AF2Point2D* getNewPoint(unsigned int newPointIndex) const;
};

#endif
