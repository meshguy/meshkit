/*
 * AF2RuleExistEdge.hpp
 *
 * A specification of an edge that must exist in order for a rule to
 * be applied.
 *
 * Note that the vector difference referred to in this class is the
 * second endpoint minus the first endpoint.
 */
#ifndef AF2RULEEXISTEDGE_HPP
#define AF2RULEEXISTEDGE_HPP

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
     * \param firstVtx a pointer to one of the endpoints of the edge
     * \param secondVtx a pointer to the other endpoint of the edge
     * \param coeffAlpha the coefficient of dx*dx to use when measuring
     *   deviation from the vector difference between the ideal coordinates
     * \param coeffBravo the coefficient of dx*dy to use when measuring
     *   deviation from the vector difference between the ideal coordinates
     * \param coeffCharlie the coefficient of dy*dy to use when measuring
     *   deviation from the vector difference between the ideal coordinates
     */
    AF2RuleExistEdge(const AF2RuleExistVertex* firstVtx,
        const AF2RuleExistVertex* secondVtx,
        double coeffAlpha, double coeffBravo, double coeffCharlie);

    /**
     * \brief Get a constant pointer to the AF2RuleExistVertex that the first
     * endpoint of this edge must match.
     */
    const AF2RuleExistVertex* getFirstVertex() const;

    /**
     * \brief Get a constant pointer to the AF2RuleExistVertex that the second
     * endpoint of this edge must match.
     */
    const AF2RuleExistVertex* getSecondVertex() const;

    /**
     * \brief Check whether a pair of vertices is a match for this
     * edge, assuming that the first vertex of the pair matches the
     * first endpoint and the second vertex of the pair matches the
     * second endpoint.
     * 
     * \param matchFirstX the x-coordinate of the first vertex of the pair
     *   of vertices that is to be checked as a potential match for the edge
     * \param matchFirstY the y-coordinate of the first vertex of the pair
     *   of vertices that is to be checked as a potential match for the edge
     * \param matchSecondX the x-coordinate of the second vertex of the pair
     *   of vertices that is to be checked as a potential match for the edge
     * \param matchSecondY the y-coordinate of the second vertex of the pair
     *   of vertices that is to be checked as a potential match for the edge
     * \param maxDeviation the maximum deviation that the vector difference
     *   between the vertices passed in may have from the vector difference
     *   between the reference coordinates (when measured with the quadratic
     *   coefficients) and still be counted as a match
     * \return true if vertices at the specified coordinates have a vector
     *   difference within the specified deviation from the vector difference
     *   between the ideal vertices, false otherwise
     */
    bool isMatching(double matchFirstX, double matchFirstY,
        double matchSecondX, double matchSecondY, double maxDeviation) const;
};

#endif
