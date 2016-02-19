/*
 * AF2RuleNewEdge.hpp
 *
 * A specification of a new edge that would be created by applying some
 * 2-dimensional advancing front rule.
 */
#ifndef AF2RULENEWEDGE_HPP
#define AF2RULENEWEDGE_HPP

class AF2RuleNewEdge
{
  private:

    int a, b;

  public:

    /**
     * \brief Constructor
     *
     * The constructor requires arguments for the indices of the two
     * endpoints of the edge that the rule would create.  The indices
     * must index into the list of vertices associated with whichever
     * rules this AF2RuleNewEdge is a member of.  The list of vertices
     * associated with the rule is the concatenation of the list of
     * existing vertices in the rule and the list of new vertices
     * that the rule would create.
     *
     * \param firstIndex the index of one endpoint of the edge
     * \param secondIndex the index of the other endpoint of the edge
     */
    AF2RuleNewEdge(int firstIndex, int secondIndex);

    /**
     * \brief Get the index of one of the two endpoints of the edge that
     * the rule would create.
     * 
     * The index will be an index into the list of vertices associated
     * with the rule, with the existing vertices listed before new
     * vertices that the rule would create.
     */
    int getFirstIndex() const;

    /**
     * \brief Get the index of the other of the two endpoints of the edge
     * that the rule would create.
     *
     * The index will be an index into the list of vertices associated
     * with the rule, with the existing vertices listed before new
     * vertices that the rule would create.
     */
    int getSecondIndex() const;
};

#endif
