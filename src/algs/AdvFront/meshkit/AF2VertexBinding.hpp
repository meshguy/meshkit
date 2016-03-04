/*
 * AF2VertexBinding.hpp
 *
 * An AF2VertexBinding is an object that maintains a binding of
 * AF2RuleExistingVertex to a point located at coordinates.  It may
 * also cache somecomputations based on the difference between the
 * reference coordinates specified in the AF2RuleExistingVertex and
 * the actual coordinates to which those reference coordinates are bound.
 */

#ifndef AF2VERTEXBINDING_HPP
#define AF2VERTEXBINDING_HPP

// C++
#include <list>
#include <map>

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2RuleExistVertex.hpp"

class AF2VertexBinding
{
  private:

    std::map<const AF2RuleExistVertex*, const AF2Point2D*> bindingMap;
    mutable bool validCache;

  public:

    /**
     * \brief A sentinel value that an AF2RuleExistVertex will be bound
     *   to if it is not bound to anything more substantive.
     */
    static const AF2Point2D NOT_BOUND;

    /**
     * \brief Constructor
     *
     * This constructor takes a list of pointers to AF2RuleExistVertex
     * objects.  The list defines which vertices may be
     * bound to actual coordinates by this vertex binding.
     *
     * \param verticesToBind a list of (pointers to) the vertices that may
     *   be bound by this AF2VertexBinding
     */
    AF2VertexBinding(
        std::list<const AF2RuleExistVertex*> const & verticesToBind);

    /**
     * \brief Bind a reference vertex to actual coordinates.
     *
     * Stores a pointer to an AF2Point2D, binding the reference
     * vertex to the coordinates of that AF2Point2D.  The reference
     * vertex must have been in the list of vertices that was passed to
     * the constructor.  If it is not, the method will throw an exception.
     *
     * This method does not take ownership of the pointers passed into
     * it.  The calling context must ensure that the pointers remain valid
     * as long as this AF2VertexBinding (or any copy of it) is in use.
     *
     * \param vertexPtr a pointer to the vertex that should be bound
     *   to coordinates by this method call
     * \param point a pointer to the point to bind the vertex to
     */
    void bindVertex(const AF2RuleExistVertex* vertexPtr,
        const AF2Point2D* pointPtr);

    /**
     * \brief Retrieve the point to which a reference vertex is bound.
     *
     * The value returned will be a pointer to AF2VertexBinding::NOT_BOUND
     * if the reference vertex is not bound to any other point.  The reference
     * vertex must have been in the list of vertices that was passed to
     * the constructor.  If it is not, the method will throw an exception.
     *
     * \param vertexPtr a pointer to the rule existing vertex for which
     *   to retrieve the bound value
     * \return the point to which the rule existing vertex is bound
     */
    const AF2Point2D* getBoundValue(const AF2RuleExistVertex* vertexPtr) const;
};

#endif
