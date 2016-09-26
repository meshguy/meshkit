/*
 * AF2RuleNewVertex.hpp
 *
 * A specification of a vertex that will be created during application
 * of a rule.  The vertex should rarely be created with the coordinates
 * at their reference locations, so a point transform is provided that
 * will transform the location to a more appropriate location based
 * on the actual locations that the rules existing vertices are bound to.
 */

#ifndef AF2RULENEWVERTEX_HPP
#define AF2RULENEWVERTEX_HPP

// MeshKit
#include "meshkit/AF2Point2D.hpp"
#include "meshkit/AF2PointTransform.hpp"
#include "meshkit/AF2Binding.hpp"

class AF2RuleNewVertex
{
  private:

    AF2Point2D referencePoint;
    const AF2PointTransform* pointTransform;

  public:

    /**
     * \brief Constructor
     *
     * Construct an AF2RuleNewVertex from a reference point location
     * and a pointer to a point transform.  The point transform is
     * cloned, so the calling context retains ownership of the point
     * transform and is responsible for deleting it.
     *
     * \param rfrncPoint the reference location of the new vertex,
     *   i.e., the location that the new vertex will have if the
     *   existing points of the rule are all placed in their ideal positions
     * \param pntTrnsfrm the point transformation that will transform
     *   the point from its reference location to an actual location
     *   based on the vertex binding
     */
    AF2RuleNewVertex(AF2Point2D const & rfrncPoint,
        const AF2PointTransform* const & pntTrnsfrm);

    /**
     * \brief Destructor
     */
    ~AF2RuleNewVertex();

    /**
     * \brief Copy constructor
     */
    AF2RuleNewVertex(const AF2RuleNewVertex & toCopy);

    /**
     * \brief Assignment operator
     */
    AF2RuleNewVertex& operator=(const AF2RuleNewVertex & rhs);

    /**
     * \brief Return the coordinates at which the new vertex should
     * be placed based on the vertex binding.
     */
    AF2Point2D getLocation(AF2Binding const & vertexBinding) const;
};

#endif
