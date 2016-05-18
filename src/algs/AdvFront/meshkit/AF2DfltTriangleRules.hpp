/*
 * AF2DfltTriangleRules.hpp
 *
 * \brief An AF2DfltTriangleRules defines a default set of rules to
 * use for a two-dimensional advancing front triangle-meshing algorithm.
 *
 * Constructing an instance of AF2DfltTriangleRules constructs all of
 * associated default rules.  The AF2DfltTriangleRules retains ownership
 * of the rule objects until it is deleted, so if any subset of the rules
 * from an instance of AF2DfltTriangleRules are used in an AF2Algorithm,
 * the context should ensure that the instance of AF2DfltTriangleRules
 * remains in memory.
 */

#ifndef AF2DFLTTRIANGLERULES
#define AF2DFLTTRIANGLERULES

// C++
#include <list>

// MeshKit
#include "meshkit/AF2PntTrnsfrmLnrV.hpp"
#include "meshkit/AF2Rule.hpp"

class AF2DfltTriangleRules
{
  private:

    std::list<const AF2Rule*> ruleList;

    const AF2Rule* make180DegreeRuleQ1() const;
    const AF2Rule* make180DegreeRuleQ5() const;
    const AF2Rule* make180DegreeRuleQ10() const;
    const AF2Rule* make180DegreeRuleQ20() const;
    const AF2Rule* make60DegreeAngleRightRule() const;
    const AF2Rule* make60DegreeAngleLeftRule() const;
    const AF2Rule* make120DegreeAngleRightRule() const;
    const AF2Rule* make120DegreeAngleLeftRule() const;
    const AF2Rule* make120DegreeAngleBothRule() const;
    const AF2Rule* makeFillTriangleRule() const;
    const AF2Rule* makeConnectVertexRule() const;
    const AF2Rule* makeConnectEdgeRule() const;

    /**
     * \brief Make a linear point transformation based on a single
     *   vertex that translates in the x-direction
     *
     * The translation will depend on only the x-coordinate of the
     * bound value of the rule existing vertex.  It will affect only
     * the x-coordinate of the result of the linear transformation.
     * The coefficient specifies how far to translate relative to
     * the x-coordinate of the bound value.
     *
     * The returned value will be allocated on the heap with the new operator.
     * It belongs to the calling context and must be deleted there.
     *
     * \param ruleExVert the rule existing vertex whose bound value
     *   will be consulted
     * \param xCoeff the translation multiple of the difference between
     *   the x-coordinate of the bound value and the x-coordinate of the
     *   reference value
     */
    AF2PntTrnsfrmLnrV* makeLinearTransformX(
        const AF2RuleExistVertex* ruleExVert, double xCoeff) const;

    /**
     * \brief Make a linear point transformation based on a single vertex
     *   that will translate a point by the vector from the
     *   reference location to the bound value of the specified vertex.
     *
     * The translation affects both the x-coordinate and the y-coordinate.
     * If applied to a point at the reference location, as is frequently
     * the case, the effect will be to translate to the bound value location.
     *
     * The returned value will be allocated on the heap with the new operator.
     * It belongs to the calling context and must be deleted there.
     *
     * \param ruleExVert the rule existing vertex whose bound value
     *   will be used
     */
    AF2PntTrnsfrmLnrV* makeTranslation(
        const AF2RuleExistVertex* ruleExVert) const;

  public:

    /**
     * \brief Constructor
     */
    AF2DfltTriangleRules();

    /**
     * \brief Destructor
     */
    ~AF2DfltTriangleRules();

    /**
     * \brief Copy constructor (throws exception)
     *
     * This object does not currently support copying, so this method
     * is implemented to throw an exception.
     */
    AF2DfltTriangleRules(const AF2DfltTriangleRules & toCopy);

    /**
     * \brief Assignment operator (throws exception)
     *
     * This object does not currently support assignment, so this method
     * is implemented to throw an exception.
     */
    AF2DfltTriangleRules& operator=(const AF2DfltTriangleRules & rhs);

    /**
     * \brief Get a list of the default triangle rules for two-dimensional
     *   advancing front.
     *
     * The list returned from this method may be freely modified by the
     * calling context, but the AF2Rule objects that are returned via
     * pointer are owned by this object.
     *
     * \return a list of const AF2Rule* suitable for running a two-dimensional
     *   advancing front triangle-meshing algorithm
     */
    std::list<const AF2Rule*> getRules() const;
};

#endif
