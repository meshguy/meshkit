/*
 * AF2DfltRuleAppVisitor.hpp
 *
 * An AF2DfltRuleAppVisitor is a default implementation of an
 * AF2RuleAppVisitor appropriate to rules that produce triangular faces.
 * It stores at most one rule application, which is whichever rule
 * application it judges to be best among the rule applications it
 * has visited.
 *
 * Since AF2DfltRuleAppVisitor stores the best result and does not
 * provide a method to reset or clear what it has judged to be best,
 * it should be used to assess rule applications for at most one neighborhood.
 *
 * The evaluation of best is based on metrics measuring the new faces
 * that the rule application would introduce.  The evaluation favors
 * elements that are round, having an area to perimeter ratio similar
 * to that of an equilateral triangle and having edges that are not
 * too much longer or shorter than unit length.  When multiple faces
 * are introduced by a rule, the metric for the first face is taken
 * as the metric for the full rule application.
 */

#ifndef AF2DFLTRULEAPPVISITOR_HPP
#define AF2DFLTRULEAPPVISITOR_HPP

// C++
#include <cmath>

// MeshKit
#include "meshkit/AF2RuleAppVisitor.hpp"

class AF2DfltRuleAppVisitor : public AF2RuleAppVisitor
{
  private:

    static const double eqTriAreaPerimSqRatio;
    double bestMetricVal;
    AF2RuleApplication* bestRuleApp;

  public:

    /**
     * \brief Constructor
     */
    AF2DfltRuleAppVisitor();

    /**
     * \brief Destructor
     */
    ~AF2DfltRuleAppVisitor();

    /**
     * \brief Copy constructor
     */
    AF2DfltRuleAppVisitor(const AF2DfltRuleAppVisitor & toCopy);

    /**
     * \brief Assignment operator
     */
    AF2DfltRuleAppVisitor& operator=(const AF2DfltRuleAppVisitor & rhs);

    /**
     * \brief Get the best rule application
     *
     * Get the rule application, if any, that this AF2RuleAppVisitor
     * has visited, as assessed by its internal metrics.  If no
     * successful rule applications have been visited, this method
     * returns a NULL pointer.
     *
     * \return the best rule application, or null if no rule applications
     *   have been visited
     */
    const AF2RuleApplication* getBestRuleApplication() const;

    /**
     * \brief Evaluate an AF2RuleApplication to see if it is
     *   the best so far.
     *
     * This implements the pure virtual method of the superclass.
     * See additional documentation there.
     *
     * \param ruleApp A potential rule application of acceptable
     *   quality
     */
    void visit(AF2RuleApplication const & ruleApp);
};

#endif
