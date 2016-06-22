/*
 * AF2Point3D.hpp
 *
 * An AF2Point3D is a three-dimensional point used in the two-dimensional
 * advancing front algorithm.
 *
 * The coordinates are named x, y, and z.
 */
#ifndef AF2POINT_3D_HPP
#define AF2POINT_3D_HPP

// MOAB
#include "moab/Types.hpp"

class AF2Point3D
{
  private:

    unsigned long localId;
    double x, y, z;
    unsigned int distToBndry;
    bool committed;
    moab::EntityHandle vertexHandle;
    moab::EntityHandle dbgHandle; // this should be used only to output intermediate mesh

  public:

    /**
     * \brief Standard constructor
     *
     * Construct a point at the specified coordinates.
     *
     * \param pntId a number that uniquely identifies this point within
     *   the (local) context of a single execution of the AF2Algorithm
     * \param xVal the x coordinate of the point
     * \param yVal the y coordinate of the point
     * \param zVal the z coordinate of the point
     */
    AF2Point3D(unsigned long pntId, double xVal, double yVal, double zVal);

    /**
     * \brief Get this points current distance to the boundary
     *   based on the advancement of the advancing front.
     */
    unsigned int getDistanceToBoundary() const;

    /**
     * \brief Get the value of the handle to the vertex in the mesh.
     *
     * This will throw an exception if the point has not been committed.
     * Users can check whether the point has been committed using the
     * isCommitted() method.
     */
    moab::EntityHandle getVertexHandle() const;

    /**
     * \brief Get the number that uniquely identifies this point within
     *   the context of the current execution of the AF2Algorithm
     */
    unsigned long getLocalId() const;

    /** for intermediate output, use another handle
     */
    moab::EntityHandle getTmpVHandle() const { return dbgHandle;};
    void setTmpVHandle (moab::EntityHandle vhandle) { dbgHandle = vhandle;};
    /**
     * \brief Get the value of the x coordinate.
     */
    double getX() const;

    /**
     * \brief Get the value of the y coordinate.
     */
    double getY() const;

    /**
     * \brief Get the value of the z coordinate.
     */
    double getZ() const;

    /**
     * \brief Determine whether this point has been committed to the mesh
     */
    bool isCommitted() const;

    /**
     * \brief If the current distance to the boundary is greater than
     *   the specified upper bound, reduce the distance to the boundary
     *   to the specified upper bound.
     */
    void limitDistanceToBoundary(unsigned int upperBound);

    /**
     * \brief Set the value of the vertex handle in the mesh that
     *   corresponds to this point.
     *
     * This method should be called at most once.  If this method has
     * already been called, then the point is already noted as committed
     * to the mesh, and calling the method again will throw an exception.
     */
    void setCommittedHandle(const moab::EntityHandle & vertexHandleArg);
};

#endif
