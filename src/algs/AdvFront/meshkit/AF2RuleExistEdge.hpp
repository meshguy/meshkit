/*
 * AF2RuleExistEdge.hpp
 *
 * A specification of an edge that must exist in order for a rule to
 * be applied.
 *
 * Note that the vector difference referred to in this class is the
 * second endpoint minus the first endpoint.  The order of the points
 * also matters when testing whether an actual edge matches this
 * specification.  The points should be ordered so that the region that
 * still needs to be meshed is to the left, i.e., ordered in a
 * counterclockwise traversal of the boundary of the unmeshed region.
 */
#ifndef AF2RULEEXISTEDGE_HPP
#define AF2RULEEXISTEDGE_HPP

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2RuleExistVertex.hpp"

class AF2RuleExistEdge
{
  private:

    const AF2RuleExistVertex* firstVertexPtr;
    const AF2RuleExistVertex* secondVertexPtr;
    double vecDiffX, vecDiffY;
    double a, b, c;

  public:

    /**
     * \brief Constructor
     *
     * The constructor requires pointers to the two endpoints of the
     * edge that must exist in order for the local neighborhood to
     * match the rule.  The context that calls this constructor is
     * responsible for ensuring that the pointers are valid at the
     * time the constructor is called and remain valid until this
     * object is destructed.
     *
     * The constructor also requires quadratic coefficients to use in
     * checking whether the vector difference between the two points that
     * match the endpoints of the edge (assuming such points exist) is
     * an acceptable match for the vector difference between points
     * at the ideal coordinates.
     *
     * The default values for the quadratic coefficients have the effect
     * that any vector difference is acceptable as long as the endpoints
     * of the edge match.
     *
     * \param firstVtx a pointer to the starting endpoint of the edge
     * \param secondVtx a pointer to the ending endpoint of the edge
     * \param coeffAlpha the coefficient of dx*dx to use when measuring
     *   deviation from the vector difference between the ideal coordinates
     * \param coeffBravo the coefficient of dx*dy to use when measuring
     *   deviation from the vector difference between the ideal coordinates
     * \param coeffCharlie the coefficient of dy*dy to use when measuring
     *   deviation from the vector difference between the ideal coordinates
     */
    AF2RuleExistEdge(const AF2RuleExistVertex* firstVtx,
        const AF2RuleExistVertex* secondVtx, double coeffAlpha = 0.0,
        double coeffBravo = 0.0, double coeffCharlie = 0.0);

    /**
     * \brief Get a constant pointer to the AF2RuleExistVertex that the
     * starting endpoint of this edge must match.
     */
    const AF2RuleExistVertex* getStart() const;

    /**
     * \brief Get a constant pointer to the AF2RuleExistVertex that the
     * ending endpoint of this edge must match.
     */
    const AF2RuleExistVertex* getEnd() const;

    /**
     * \brief Check whether a pair of points is a match for this
     * edge, assuming that the first point of the pair matches the
     * start endpoint of the edge and the second point of the pair
     * matches the end endpoint of the edge.
     * 
     * \param startPnt a point that is a potential match for the start
     *   endpoint of the edge and the first of a pair of points that is
     *   to be checked as a potential match for the edge
     * \param endPnt a point that is a potential match for the end
     *   endpoint of the edge and the second of a pair of points that is
     *   to be checked as a potential match for the edge
     * \param maxDeviation the maximum deviation that the vector difference
     *   between the points passed in may have from the vector difference
     *   between the reference coordinates (when measured with the quadratic
     *   coefficients) and still be counted as a match
     * \return true if the specified points have a vector difference
     *   within the specified deviation from the vector difference
     *   between the reference vertices, false otherwise
     */
    bool isMatching(AF2Point2D const & startPnt,
        AF2Point2D const & endPnt, double maxDeviation) const;
};

#endif
