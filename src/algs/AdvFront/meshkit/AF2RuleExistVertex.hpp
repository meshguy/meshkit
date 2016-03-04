/*
 * AF2RuleExistVertex.hpp
 *
 * A specification of a vertex that must exist in order for a
 * rule to be applied.
 *
 * In general, a point having coordinates that exactly match the
 * reference coordinates will not exist, but a point within some
 * neighborhood of the ideal coordinates is still considered a
 * matching * vertex.
 */
#ifndef AF2RULEEXISTVERTEX_HPP
#define AF2RULEEXISTVERTEX_HPP

// MeshKit
#include "meshkit/AF2Point2D.hpp"

class AF2RuleExistVertex
{
  private:

    double x, y;
    double a, b, c;

  public:

    /**
     * \brief Constructor
     *
     * A constructor that defines the reference coordinates along with
     * quadratic coefficients to use in checking whether a vertex
     * is an acceptable match for this vertex.
     *
     * \param refXCoord the x-coordinate of the ideal matching vertex
     * \param refYCoord the y-coordinate of the ideal matching vertex
     * \param coeffAlpha the coefficient of dx*dx to use when measuring
     *   deviation from the reference coordinates
     * \param coeffBravo the coefficient of dx*dy to use when measuring
     *   deviation from the reference coordinates
     * \param coeffCharlie the coefficient of dy*dy to use when measuring
     *   deviation from the reference coordinates
     */
    AF2RuleExistVertex(double refXCoord, double refYCoord,
        double coeffAlpha = 1.0, double coeffBravo = 0.0,
        double coeffCharlie = 1.0);

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
     * \param matchPoint the point that is to be checked as a potential
     *   match for the ideal vertex
     * \param maxDeviation the maximum deviation from the reference
     *   coordinates (when measured with the quadratic coefficients)
     *   that is allowed to match this vertex
     * \return true if the coordinates of the specified point are within
     *   the specified deviation of the ideal vertex, false otherwise
     */
    bool isMatching(AF2Point2D const & matchPoint, double maxDeviation) const;
};

#endif
