/*
 * AF2RuleNewFace.hpp
 *
 * A specification of a new element that would be created by applying some
 * 2-dimensional advancing front rule.
 */
#ifndef AF2RULENEWFACE_HPP
#define AF2RULENEWFACE_HPP

class AF2RuleNewFace
{
  public:

    /**
     * \brief Get the number of vertices that the element this rule would
     * create would have.
     *
     * \return the number of vertices that the element would have
     */
    virtual unsigned int getNumVertices() const = 0;

    /**
     * \brief Get the index of one of the vertices of the element
     * that the rule would create.
     *
     * Indices returned by this method index into the list of vertices
     * associated with the rule, i.e., a list consisting of a list
     * of the rule's existing vertices followed by a list of the new
     * vertices that the rule would create.  The numbering of the vertices
     * of the face begins with 0, so the valid arguments to this method are
     * 0 through n - 1, where n is the value returned from getNumVertices().
     *
     * \param vtxNum the number of the vertex for which to get the index
     * \return the index in the list of the rule vertices of
     *   vertex number vtxNum of the element
     */
    virtual unsigned int getVertexIndex(unsigned int vtxNum) const = 0;
};

#endif
