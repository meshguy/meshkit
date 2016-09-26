/*
 * AF2RuleNewElement.hpp
 *
 * A specification of a new triangle that would be created by applying some
 * 2-dimensional advancing front rule.
 */
#ifndef AF2RULENEWTRIANGLE_HPP
#define AF2RULENEWTRIANGLE_HPP

#include "meshkit/AF2RuleNewFace.hpp"

class AF2RuleNewTriangle : public AF2RuleNewFace
{
  private:

    unsigned int triVtxIndices[3];

  public:

    /**
     * \brief Constructor
     *
     * The constructor requires arguments for the indices of the three
     * vertices of the triangle that the rule would create.  The indices
     * must index into the list of vertices associated with whichever
     * rules this AF2RuleNewTriangle is a member of.  The list of vertices
     * associated with the rule is the concatenation of the list of
     * existing vertices in the rule and the list of new vertices
     * that the rule would create.
     *
     * \param firstIndex the index of a vertex of the triangle
     * \param secondIndex the index of the next vertex of the triangle after
     *   the firstIndex traversing the vertices counterclockwise
     * \param thirdIndex the index of the remaining vertex of the element
     */
    AF2RuleNewTriangle(unsigned int firstIndex,
        unsigned int secondIndex, unsigned int thirdIndex);

    unsigned int getNumVertices() const;

    unsigned int getVertexIndex(unsigned int vtxNum) const;
};

#endif
