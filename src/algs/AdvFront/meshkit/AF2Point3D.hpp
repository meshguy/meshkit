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

#include <iGeom.h>

class AF2Point3D
{
  private:

    double x, y, z;
    unsigned int distToBndry;
    bool committed;
    iBase_EntityHandle vertexHandle;

  public:

    /**
     * \brief Standard constructor
     *
     * Construct a point at the specified coordinates.
     *
     * \param xVal the x coordinate of the point
     * \param yVal the y coordinate of the point
     * \param zVal the z coordinate of the point
     */
    AF2Point3D(double xVal, double yVal, double zVal);

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
    iBase_EntityHandle getVertexHandle() const;

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
    void setCommittedHandle(iBase_EntityHandle & vertexHandleArg);
};

#endif
