/*
 * AF2Rule.hpp
 *
 * TODO: Document
 */
#ifndef AF2RULE_HPP
#define AF2RULE_HPP

// C++
#include <list>
#include <map>

// MeshKit
#include "meshkit/AF2Edge2D.hpp"
#include "meshkit/AF2FreeZoneDef.hpp"
#include "meshkit/AF2Neighborhood.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2RuleExistEdge.hpp"
#include "meshkit/AF2RuleExistVertex.hpp"
#include "meshkit/AF2RuleNewEdge.hpp"
#include "meshkit/AF2RuleNewFace.hpp"
#include "meshkit/AF2RuleNewVertex.hpp"

class AF2Rule
{
  private:

    const unsigned int numExVertices, numExEdges;
    const AF2FreeZoneDef* const freeZoneDef;
    const unsigned int numNewVertices, numNewEdges, numNewFaces;
    const AF2RuleExistVertex** exVertices;
    const AF2RuleExistEdge** exEdges;
    const AF2RuleNewVertex** newVertices;
    const AF2RuleNewEdge** newEdges;
    const AF2RuleNewFace** newFaces;

    /**
     * \brief Find which of the edges in a given neighborhood are potential
     *   matches for this rule's existing edges at a specified match quality.
     *
     * TODO: Complete the documentation
     */
    std::map<const AF2RuleExistEdge*, std::list<AF2Edge2D*>*>*
        findPotentialEdgeMatches(AF2Neighborhood const & ngbhd,
        int matchQuality) const;

    /**
     * \brief Check whether a specified AF2Edge2D is a match for
     *     a specified AF2RuleExistEdge at a given match quality level.
     * 
     * All quality levels are specified as strictly positive integers.
     * Lower match quality numbers correspond to higher quality matches.
     *
     * \param point an AF2Point2D point at defined coordinates
     * \param ruleVertex an AF2RuleExistVertex that is part of the definition
     *   of this rule
     * \param matchQuality an integer measuring how closely the coordinates
     *   of the point must match the coordinates of the rule existing vertex
     *   in order to be considered a match
     * \return true if the point matches the rule vertex to the desired
     *   of quality; false otherwise
     */
    bool isMatchingEdge(AF2Edge2D const & edge,
        AF2RuleExistEdge const & ruleEdge, int matchQuality) const;

    /**
     * \brief Check whether a specified AF2Point2D is a match for
     *     a specified AF2RuleExistVertex at a given match quality level.
     * 
     * All quality levels are specified as strictly positive integers.
     * Lower match quality numbers correspond to higher quality matches.
     *
     * \param point an AF2Point2D point at defined coordinates
     * \param ruleVertex an AF2RuleExistVertex that is part of the definition
     *   of this rule
     * \param matchQuality an integer measuring how closely the coordinates
     *   of the point must match the coordinates of the rule existing vertex
     *   in order to be considered a match
     * \return true if the point matches the rule vertex to the desired
     *   of quality; false otherwise
     */
    bool isMatchingVertex(AF2Point2D const & point,
        AF2RuleExistVertex const & ruleVertex, int matchQuality) const;

  public:

    /**
     * \brief Constructor
     *
     * This object takes ownership of all of the objects that are
     * passed into the constructor.  The objects will be destroyed
     * when this object is destroyed, and the context that constructs
     * this object must not delete the objects that it passes to
     * this constructor.
     *
     * The constructor requires a list of the rule's existing vertices,
     * existing edges, free zone definition, new vertices, new edges,
     * and new faces.  These objects may reference each other.  For example,
     * the existing edges reference existing vertices as their endpoints.
     * (This is one of the reasons that pointers are passed into this
     * method and this object takes ownership of the objects that are
     * passed in without supporting copy construction or assignment.)
     *
     * \param paramName paramDesc
     */
    AF2Rule(std::list<const AF2RuleExistVertex*> const & ruleVertices,
        std::list<const AF2RuleExistEdge*> const & ruleEdges,
        const AF2FreeZoneDef* freeZoneDef,
        std::list<const AF2RuleNewVertex*> const & ruleNewVertices,
        std::list<const AF2RuleNewEdge*> const & ruleNewEdges,
        std::list<const AF2RuleNewFace*> const & ruleNewFaces);

    /**
     * \brief Destructor
     */
    ~AF2Rule();

    /**
     * \brief Copy constructor (throws exception)
     *
     * This object does not currently support copying, so this method
     * is implemented to throw an exception.
     */
    AF2Rule(const AF2Rule & toCopy);

    /**
     * \brief Assignment operator (throws exception)
     *
     * This object does not currently support assignment, so this method
     * is implemented to throw an exception.
     */
    AF2Rule& operator=(const AF2Rule & rhs);

    /**
     * TODO: Document this method
     */
    void applyRule(AF2Neighborhood const & ngbhd, int matchQuality) const;
};

#endif
