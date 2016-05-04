/*
 * AF2PlaneProjMakerDef.hpp
 *
 * An AF2PlaneProjMakerDef is an object that implements the default
 * method for making an AF2PlaneProjection.
 *
 * The method depends on having access to the geometric surface.  It uses the
 * start point of the neighborhood baseline edge as the origin of the plane.
 * It uses the direction of the neighborhood baseline edge as the x-direction
 * of the plane.  The mean of the normal vectors to the surface at the
 * endpoints of the baseline edge, after any component in the x-direction
 * is removed, is used as the normal to the plane.
 */
#ifndef AF2PLANEPROJMAKERDEF_HPP
#define AF2PLANEPROJMAKERDEF_HPP

// MeshKit
#include "meshkit/iGeom.hpp"
#include "meshkit/AF2LocalTransformMaker.hpp"

class AF2PlaneProjMakerDef : public AF2LocalTransformMaker
{
  private:

    iGeom* iGeomPtr;
    iGeom::EntityHandle surface;
    // TODO: Sizing function

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
     */
     AF2PlaneProjMakerDef(iGeom* iGeomPtrArg, iGeom::EntityHandle surf);

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
