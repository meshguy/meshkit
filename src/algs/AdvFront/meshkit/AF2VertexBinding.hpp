/*
 * AF2VertexBinding.hpp
 *
 * An AF2VertexBinding is an object that maintains a binding of
 * AF2RuleExistingVertex to coordinates.  It also caches some
 * computations based on the difference between the reference coordinates
 * specified in the AF2RuleExistingVertex and the actual coordinates
 * to which those reference coordinates are bound.
 */

#ifndef AF2VERTEXBINDING_HPP
#define AF2VERTEXBINDING_HPP

// C++
#include <list>
#include <map>

// MeshKit
#include "meshkit/AF2RuleExistVertex.hpp"
#include "meshkit/Matrix.hpp"

class AF2VertexBinding
{
  private:

    std::map<const AF2RuleExistVertex*, MeshKit::Vector<2> > bindingMap;
    bool validCache;

  public:

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
     * Stores the coordinates of an actual vertex, binding the reference
     * vertex to those coordinates.  The reference vertex must be in
     * the list of vertices that was passed to the constructor.  If it
     * is not, the method will throw an exception.
     *
     * \param vertexPtr a pointer to the vertex that should be bound
     *   to coordinates by this method call
     * \param coordinates the coordinates to bind the vertex to
     */
    void bindVertex(const AF2RuleExistVertex* vertexPtr,
        MeshKit::Vector<2> const & coordinates);
};

#endif
