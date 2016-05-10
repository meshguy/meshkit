/*
 * AF2Rule.hpp
 *
 * \brief An AF2Rule defines a two-dimensional pattern for advancing the front
 * of a two-dimensional advancing front algorithm.
 *
 * The pattern consists of a number of existing edges and existing vertices,
 * a free zone definition, new vertices, new edges, and new faces.
 * In order for the advancing algorithm to apply the pattern, there must
 * be edges and points on the advancing front in the neighborhood
 * of some baseline edge that match the defined existing edges and existing
 * vertices of the pattern.  Depending on the particular binding of
 * existing edges and vertices to edges and points in the neighborhood,
 * as well as the quality level, the free zone definition will define
 * a free zone.  This free zone must not contain any of the neighborhood's
 * points or intersect any of its edges, or the algorithm will not apply
 * the pattern.
 *
 * If the pattern can be successfully applied, then applying the pattern
 * involves adding the new vertices, new edges, and new faces of the pattern
 * to the mesh.  Since there may be several different ways to apply
 * the pattern or several different patterns that could be applied,
 * each potential successful application of the pattern is handed off
 * to a visitor for processing.  The rule makes no guarantee that
 * the new faces are not inverted, and the visitor should verify any
 * desired properties before committing the results to the mesh.
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
#include "meshkit/AF2RuleAppVisitor.hpp"
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
    unsigned int numExIsoVertices;
    const AF2RuleExistVertex** exIsoVertices;

    /**
     * \brief The second stage of applying a rule
     *
     * This stage matches any of the rule's existing vertices that
     * are not endpoints of edges with points in the neighborhood.
     *
     * \param ngbhd the neighborhood passed to applyRule
     * \param matchQuality the match quality passed to applyRule
     * \param visitor the visitor passsed to applyRule
     * \param matchingVerticesMap a map of the rule's existing
     *   vertices that are not endpoints of edges to a list
     *   of points in the neighborhood that are potential matches
     *   for those vertices
     * \param binding a binding of this rule's existing edges and
     *   vertices to edges and points in the neighborhood, understood
     *   to be in progress, containing consistent bindings for all
     *   of the rule's existing edges and (implicitly) the endpoints
     *   of those edges
     */
    void applyRuleStageTwo(AF2Neighborhood const & ngbhd,
        unsigned int matchQuality, AF2RuleAppVisitor & visitor,
        std::map<const AF2RuleExistVertex*,
        std::list<const AF2Point2D*>*>* const & matchingVerticesMap,
        AF2Binding & binding) const;

    /**
     * \brief The third stage of applying a rule
     *
     * This stage checks that the rule's free zone does not
     * contain any of the neighborhood's points or intersect
     * any of the neighborhood's edges except for points and
     * edges that coincide with the boundary of the free zone.
     *
     * If the free zone is appropriately empty, then this stage
     * also builds the rule application object and passes it to the visitor.
     *
     * \param ngbhd the neighborhood passed to applyRule
     * \param matchQuality the match quality passed to applyRule
     * \param visitor the visitor passsed to applyRule
     * \param binding a binding of this rule's existing edges and
     *   vertices to edges and points in the neighborhood, understood
     *   to be consistent and completed, i.e., mapping all of this
     *   rule's existing edges and vertices to edges and points in
     *   the neighborhood
     */
    void applyRuleStageThree(AF2Neighborhood const & ngbhd,
        unsigned int matchQuality, AF2RuleAppVisitor & visitor,
        AF2Binding const & binding) const;

    /**
     * \brief Find which of the edges in a given neighborhood are potential
     *   matches for this rule's existing edges at a specified match quality.
     *
     * \param ngbhd the neighborhood passed to applyRule
     * \param matchQuality the match quality passed to applyRule
     */
    std::map<const AF2RuleExistEdge*, std::list<const AF2Edge2D*>*>*
        findPotentialEdgeMatches(AF2Neighborhood const & ngbhd,
        unsigned int matchQuality) const;

    /**
     * \brief Find which of the vertices in a given neighborhood are potential
     *   matches for this rule's existing vertices at a specified match
     *   quality.
     *
     * \param ngbhd the neighborhood passed to applyRule
     * \param matchQuality the match quality passed to applyRule
     */
    std::map<const AF2RuleExistVertex*, std::list<const AF2Point2D*>*>*
        findPotentialVertexMatches(AF2Neighborhood const & ngbhd,
        unsigned int matchQuality) const;

    /**
     * \brief Check that the endpoints of the rule's existing edges
     *   are listed among the rule's existing vertices and find
     *   the any of the rule's existing vertices that are not
     *   endpoints of the rule's existing edges
     *
     * This method is intended to be used by the constructor.  It throws
     * an exception if there is an existing edge that has an endpoint
     * that is not listed in the existing vertices, and it fills the
     * data members that store the isolated vertices, i.e., the vertices
     * that are not endpoints of edges.
     */
    void checkExEndpointsAndFindIsolatedVertices();

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
        AF2RuleExistEdge const & ruleEdge, unsigned int matchQuality) const;

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
        AF2RuleExistVertex const & ruleVertex, unsigned int matchQuality) const;

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
     * The constructor will throw an exception if it detects an
     * inconsistency in the arguments.  For instance, each rule must have
     * a baseline existing edge and two existing vertices.  All of the
     * endpoints of the rule's existing edges must be listed in the rule's
     * existing vertices.  The baseline edge may not be listed among
     * the rule's other existing edges.
     *
     * \param ruleVertices the rule's existing vertices, i.e., the vertices
     *   that must match some vertex in the local neighborhood in order for
     *   the rule to be applied
     * \param baselineEdge the baseline edge that must match the baseline
     *   edge of the local neighborhood in order for the rule to be applied
     * \param otherRuleEdges the other existing edges of the rule, i.e.,
     *   the other edges that must match some edge in the local
     *   neighborhood in order for the rule to be applied
     * \param freeZoneDef a free zone definition that can be used to construct
     *   a free zone based on which local neighborhood vertices the rule's
     *   existing vertices are bound to.  The rule may be applied only if
     *   no local neighborhood vertices or edges may intersect the
     *   interior of the free zone
     * \param ruleNewVertices new vertices that will be created if the
     *   rule is applied
     * \param ruleNewEdges new edges that will be created if the rule
     *   is applied
     * \param ruleNewFaces new face elements, such as triangles or
     *   quadrilaterals, that will be created if the rule is applied
     */
    AF2Rule(std::list<const AF2RuleExistVertex*> const & ruleVertices,
        const AF2RuleExistEdge* baselineEdge,
        std::list<const AF2RuleExistEdge*> const & otherRuleEdges,
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
     * \brief Attempt to apply this rule to a neighborhood at the
     *   specified quality level, passing each successful rule
     *   application to the visitor for processing
     *
     * The quality level is a positive integer.  Lower values require
     * higher-quality matching, i.e., the rule's edges and vertices
     * must match items in the neighborhood more precisely when the
     * quality level number is a lower value.
     *
     * \param ngbhd an advancing front neighborhood
     * \param matchQuality the quality level for matching this rule's
     *   existing edges and vertices
     * \param visitor a visitor that will handle the processing
     *   of any possible applications of this rule within the
     *   neighborhood that meet the specified level of quality
     */
    void applyRule(AF2Neighborhood const & ngbhd, unsigned int matchQuality,
        AF2RuleAppVisitor & visitor) const;
};

#endif
