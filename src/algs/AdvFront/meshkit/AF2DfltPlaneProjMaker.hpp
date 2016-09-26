/*
 * AF2DfltPlaneProjMaker.hpp
 *
 * An AF2DfltPlaneProjMaker is an object that implements the default
 * method for making an AF2PlaneProjection.
 *
 * The method depends on having access to the geometric surface.  It uses the
 * start point of the neighborhood baseline edge as the origin of the plane.
 * It uses the direction of the neighborhood baseline edge as the x-direction
 * of the plane.  The mean of the normal vectors to the surface at the
 * endpoints of the baseline edge, after any component in the x-direction
 * is removed, is used as the normal to the plane.
 *
 * If a sizing function is provided, the size at the midpoint of the
 * (3-dimensional) baseline edge, which is a point that might not lie on
 * the surface, will be used as the scale for the plane projection,
 * provided that the size is larger than zero and is not too much different
 * than the actual length of the baseline edge.  In other words, a distance
 * in the plane of that length will be considered a distance of length
 * one when a rule is applied.  If no sizing function is provided or
 * the sizing function returns a value that is very near to zero relative
 * to the actual length of the baseline edge, then the actual length of
 * the baseline edge will be used as the scale.  If the sizing function
 * appears to return a reasonably defined value, but the value is much
 * smaller or larger than the actual length of the baseline edge, then
 * a scale closer to the actual length that attempts to change the size
 * towards the desired size will be used.
 */
#ifndef AF2DFLTPLANEPROJMAKER_HPP
#define AF2DFLTPLANEPROJMAKER_HPP

// C++
#include <cstddef>

// MeshKit
#include "meshkit/AF2LocalTransformMaker.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/iGeom.hpp"

class AF2DfltPlaneProjMaker: public AF2LocalTransformMaker
{
  private:

    iGeom* iGeomPtr;
    iGeom::EntityHandle surface;
    MeshKit::SizingFunction* sizing;

  public:

    /**
     * \brief Constructor
     *
     * This object does not own any of the objects passed into it.
     * The memory for these objects must be managed outside of this object.
     *
     * \param iGeomPtrArg an iGeom instance
     * \param surf a handle to a surface geometry representation within the
     *   iGeom instance
     * \param sizingArg a sizing function that, if provided, will define the
     *   preferred size
     */
     AF2DfltPlaneProjMaker(iGeom* iGeomPtrArg, iGeom::EntityHandle surf,
         MeshKit::SizingFunction* sizingArg = NULL);

    /**
     * \brief Make an AF2PlaneProjection
     *
     * Implements the method from the superclass AF2LocalTransformMaker
     * to make an AF2PlaneProjection.  See additional documentation
     * in the superclass.
     */
    AF2LocalTransform* makeLocalTransform(
        const std::list<AF2Point3D*> & ngbhdPoints,
        const AF2Edge3D* const & baselineEdge,
        const std::list<const AF2Edge3D*> & otherNgbhdEdges) const;
};

#endif
