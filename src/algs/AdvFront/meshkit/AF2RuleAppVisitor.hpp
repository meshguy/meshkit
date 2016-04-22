/*
 * AF2RuleAppVisitor.hpp
 *
 * An AF2RuleAppVisitor is an object that visits AF2RuleApplication
 * objects produced by successful potential rule applications.  It
 * is a pure virtual class that requires implementing the appropriate
 * visit method.
 */

#ifndef AF2RULEAPPVISITOR_HPP
#define AF2RULEAPPVISITOR_HPP

// MeshKit
#include "meshkit/AF2RuleApplication.hpp"

class AF2RuleAppVisitor
{
  public:

    /**
     * \brief Visit an AF2RuleApplication
     *
     * This pure virtual method provides an interface for an
     * AF2RuleAppVisitor to examine and process AF2RuleApplication
     * objects that represent ways that a rule could successfully
     * be applied.  When an AF2RuleAppVisitor is passed to the
     * applyRule method of an AF2Rule, each acceptable quality
     * potential application of the rule will be passed to the
     * AF2RuleAppVisitor as an AF2RuleApplication for processing.
     *
     * This method may change the state of the visitor.
     *
     * \param ruleApp A potential rule application of acceptable
     *   quality
     */
    virtual void visit(AF2RuleApplication const & ruleApp) = 0;
};

#endif
