/*
 * AF2Binding.hpp
 *
 * An AF2Binding is an object that manages binding the edges and vertices
 * that a rule has specified must exist (i.e., AF2RuleExistingEdge and
 * AF2RuleExistingVertex) to edges and vertices (i.e., AF2Edge2D
 * and AF2Point2D) that actually do exist in some AF2Neighborhood.
 *
 * The AF2Binding provides methods to bind an edge or a vertex and to check
 * whether binding an edge or vertex is consistent with what has already
 * been bound.  Binding an edge implies binding the endpoint vertices
 * of that edge.  Attempting to bind an edge or vertex that is not
 * consistent with what is already bound will result in an exception.
 * Attempting to bind an edge or vertex that is already bound will succeed
 * if the binding is consistent.
 *
 * AF2Binding also provides methods to release edges or vertices that are
 * bound.  Attempting to release an edge or vertex that is not bound
 * will result in an exception.  The model for binding a vertex tracks
 * both whether the vertex has been explicitly bound with a call to the
 * method that binds vertices and whether it is implicitly bound because it
 * is an endpoint of an edge that is bound.  Attempting to release a vertex
 * that has been bound in both ways will release the explicit bound, but not
 * the implicit bound.  Attempting to release a vertex that is implicitly
 * bound but not explicitly bound will result in an exception.
 */

#ifndef AF2BINDING_HPP
#define AF2BINDING_HPP

// C++
#include <map>
#include <set>

// MeshKit
#include "meshkit/AF2Edge2D.hpp"
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2RuleExistEdge.hpp"
#include "meshkit/AF2RuleExistVertex.hpp"

class AF2Binding
{
  private:

    struct VtxBindRec
    {
      const AF2Point2D* pointPtr;
      unsigned int numBoundEdges;
      bool xplctBound;

      /**
       * \brief Constructor
       */
      VtxBindRec();
    };

    std::map<const AF2RuleExistVertex*, VtxBindRec> vertexBindMap;
    std::set<const AF2Point2D*> verticesInUse;
    std::map<const AF2RuleExistEdge*, const AF2Edge2D*> edgeBindMap;
    std::set<const AF2Edge2D*> edgesInUse;

    /**
     * \brief Bind a vertex.
     *
     * \param isExplicit differentiate between binding a vertex due
     *   to an explicit method call from a user and binding a vertex
     *   implicitly because it is an endpoint of an edge
     */
    void bind(const AF2RuleExistVertex* ruleVertexPtr,
        const AF2Point2D* ngbhdVertexPtr, bool isExplicit);

    /**
     * \brief Release a vertex binding.
     *
     * \param isExplicit differentiate between releasing a vertex due
     *   to an explicit method call from a user and releasing an implicit
     *   binding of the vertex due to the vertex being an endpoint of
     *   an edge that is being released
     */
    void release(const AF2RuleExistVertex* ruleVertexPtr, bool isExplicit);

  public:

    /**
     * \brief Bind a rule edge to a neighborhood edge.
     *
     * Binds a rule edge to a neighborhood edge and implicitly binds the
     * rule vertices that are endpoints of the rule edge to the neighborhood
     * vertices that are endpoints of the neighborhood edge.  If binding
     * the edge is not consistent, this method will throw an exception.
     * The consistency of the binding can be checked before attempting
     * the binding by calling the isConsistent method.
     *
     * This method does not take ownership of the pointers passed into
     * it.  The calling context must ensure that the pointers remain valid
     * as long as this AF2Binding (or any copy of it) is in use.
     *
     * \param ruleEdgePtr a pointer to the rule edge that should be bound
     * \param ngbhdEdgePtr a pointer to the neighborhood edge to which the
     *   rule edge should be bound
     */
    void bind(const AF2RuleExistEdge* ruleEdgePtr,
        const AF2Edge2D* ngbhdEdgePtr);

    /**
     * \brief Explicitly bind a rule vertex to a neighborhood vertex.
     *
     * Explicitly binds a rule vertex to a neighborhood vertex, even if
     * the rule vertex is already implicitly bound to the neighborhood
     * vertex.  If binding the vertex is not consistent, this method will
     * throw an exception.  The consistency of the binding can be checked
     * before attempting the binding by calling the isConsistent method.
     *
     * This method does not take ownership of the pointers passed into
     * it.  The calling context must ensure that the pointers remain valid
     * as long as this AF2Binding (or any copy of it) is in use.
     *
     * \param ruleVertexPtr a pointer to the rule vertex that should be bound
     * \param ngbhdVertexPtr a pointer to the neighborhood vertex to which the
     *   rule vertex should be bound
     */
    void bind(const AF2RuleExistVertex* ruleVertexPtr,
        const AF2Point2D* ngbhdVertexPtr);


    /**
     * \brief Retrieve the neighborhood vertex to which a rule vertex is bound.
     *
     * If the rule vertex is bound to a neighborhood vertex implicitly,
     * explicitly, or both implicitly and explicitly, this method will return
     * a pointer to the neighborhood vertex that the rule vertex is bound to.
     * If the rule vertex is not bound to any neighborhood vertex,
     * this method will return a null pointer.
     *
     * \param ruleVertexPtr a pointer to the rule existing vertex for which
     *   to retrieve the bound value
     * \return a pointer to the neighborhood vertex to which the rule
     *   existing vertex is bound, or a null pointer if rule veretx
     *   is not bound to any neighborhood vertex
     */
    const AF2Point2D* getBoundValue(
        const AF2RuleExistVertex* ruleVertexPtr) const;

    /**
     * \brief Check whether binding a rule edge to a neighborhood edge
     *   is consistent with this AF2Binding's current bindings.
     *
     * Binding an edge may be inconsistent for any of the following reasons.
     * (1) The AF2Binding may have the rule edge bound to some other
     *   neighborhood edge.
     * (2) The AF2Binding may have some other rule edge bound to the
     *   neighborhood edge.
     * (3) There may be an inconsistency in binding one of the endpoints of
     *   the edge.  See the isConsistent method for checking consistency of
     *   binding a rule vertex to a neighborhood vertex.
     *
     * \param ruleEdgePtr a pointer to a rule edge
     * \param ngbhdEdgePtr a pointer to a neighborhood edge
     * \return true if binding the rule edge to the neighborhood edge
     *   is consistent with this AF2Binding's current bindings; false otherwise
     */
    bool isConsistent(const AF2RuleExistEdge* ruleEdgePtr,
        const AF2Edge2D* ngbhdEdgePtr) const;

    /**
     * \brief Check whether binding a rule vertex to a neighborhood vertex
     *   is consistent with this AF2Binding's current bindings.
     *
     * Binding the specified rule vertex to the specified neighborhood
     * vertex may be inconsistent for any of the following reasons.
     * (1) The AF2Binding may have the rule vertex bound to some other
     *   neighborhood vertex.
     * (2) The AF2Binding may have some other rule vertex bound to the
     *   neighborhood vertex.
     *
     * \param ruleVertexPtr a pointer to a rule vertex
     * \param ngbhdVertexPtr a pointer to a neighborhood vertex
     * \return true if binding the rule vertex to the neighborhood vertex
     *   is consistent with this AF2Binding's current bindings; false otherwise
     */
    bool isConsistent(const AF2RuleExistVertex* ruleVertexPtr,
        const AF2Point2D* ngbhdVertexPtr) const;

    /**
     * \brief Release the binding of the specified rule edge to whatever
     *   neighborhood edge it is bound to.
     *
     * Release the binding of the specified rule edge and release the
     * implicit binding of the rule vertices that are endpoints of the
     * rule edge.  Releasing the endpoints will completely release an
     * endpoint's vertex binding only if the vertex is not explicitly
     * bound and there are no other bound edges that share the endpoint.
     *
     * This method will throw an exception if the rule edge is not bound
     * to a neighborhood edge.
     *
     * \param ruleEdgePtr a pointer to the rule edge for which to release
     *   the binding
     */
    void release(const AF2RuleExistEdge* ruleEdgePtr);

    /**
     * \brief Release the explicit binding of the specified rule vertex to
     *   whatever neighborhood vertex the rule vertex is bound to.
     *
     * This method will throw an exception if the specified rule vertex
     * is not bound to any neighborhood vertex or if the vertex is
     * only implicitly bound to a neighborhood vertex because it is
     * the endpoint of one or more rule edges that are bound to
     * neighborhood edges.
     *
     * \param ruleVertexPtr a pointer to the rule vertex for which to release
     *   the explicit binding
     */
    void release(const AF2RuleExistVertex* ruleVertexPtr);
};

#endif
