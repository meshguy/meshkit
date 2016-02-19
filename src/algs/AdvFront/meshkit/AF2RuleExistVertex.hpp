/*
 * AF2RuleExistVertex.hpp
 *
 * A specification of a vertex that must exist in order for a rule to
 * be applied.  A vertex with exact coordinates will not exist in general,
 * but a vertex within some neighborhood of given coordinates may be found.
 */
#ifndef AF2RULEEXISTVERTEX_HPP
#define AF2RULEEXISTVERTEX_HPP

class AF2RuleExistVertex
{
  private:

    double x, y;
    double a, b, c;

  public:

    /**
     * \brief Constructor
     *
     * A constructor that defines the ideal coordinates along with
     * quadratic coefficients to use in checking whether a vertex
     * is an acceptable match for this vertex.
     *
     * \param refXCoord the x-coordinate of the ideal matching vertex
     * \param refYCoord the y-coordinate of the ideal matching vertex
     * \param coeffAlpha the coefficient of dx*dx to use when measuring
     *   deviation from the ideal coordinates
     * \param coeffBravo the coefficient of dx*dy to use when measuring
     *   deviation from the ideal coordinates
     * \param coeffCharlie the coefficient of dy*dy to use when measuring
     *   deviation from the ideal coordinates
     */
    AF2RuleExistVertex(double refXCoord, double refYCoord,
        double coeffAlpha, double coeffBravo, double coeffCharlie);

    /**
     * \brief Get the ideal x-coordinate
     */
    double getX() const;

    /**
     * \brief Get the ideal y-coordinate
     */
    double getY() const;

    /**
     * \brief Check whether a vertex is a match for this ideal vertex
     * 
     * \param matchX the x-coordinate of the vertex that is to be checked
     *   as a potential match for the ideal vertex
     * \param matchY the y-coordinate of the vertex that is to be checked
     *   as a potential match for the ideal vertex
     * \param maxDeviation the maximum deviation from the reference
     *   coordinates (when measured with the quadratic coefficients)
     *   that is allowed to match this vertex
     * \return true if a vertex at the specified coordinates is within
     *   the specified deviation of the ideal vertex, false otherwise
     */
    bool isMatching(double matchX, double matchY, double maxDeviation) const;
};

#endif
