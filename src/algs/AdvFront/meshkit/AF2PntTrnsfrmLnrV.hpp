/*
 * AF2PntTrnsfmLnrV.hpp
 *
 * \brief A point transformation that moves the point by some linear
 * transformation of the difference between the points that specified
 * reference vertices are bound to and the locations of the reference
 * vertices themselves.
 *
 * This implementation does not make use of a linear transformation matrix.
 * Instead it performs the low-level calculations of multiplying coefficients
 * by the vertex binding differences and summing the result.
 */

#ifndef AF2PNTTRNSFRMLNRV_HPP
#define AF2PNTTRNSFRMLNRV_HPP

// C++
#include <list>
#include <vector>

// MeshKit
#include "meshkit/AF2PointTransform.hpp"
#include "meshkit/AF2RuleExistVertex.hpp"

class AF2PntTrnsfrmLnrV : public AF2PointTransform
{
  private:

    std::vector<const AF2RuleExistVertex*> refVertices;
    std::vector<double> xCoeff;
    std::vector<double> yCoeff;

  public:

    /**
     * \brief Constructor with lists
     *
     * Construct an AF2PntTrnsfrmLnrV that is based on the differences
     * between bound locations of the specified vertices and transforms
     * these differences using a linear transformation with the specified
     * coefficients.
     *
     * The lists of coefficients should be twice as long as the list
     * of vertices.  For the vertex v at position i in the list of vertices
     * (where the first position is i = 0), the xDiffCoeff at position 2*i
     * defines the coefficient on the contribution of x-coordinate difference
     * for vertex v to the x translation, and the xDiffCoeff at
     * position 1 + 2*i defines the coefficient on the contribution of the
     * y-coordinate difference for vertex v to the x translation.  Similarly,
     * the yDiffCoeff at positions 2*i and 1 + 2*i define the coefficients
     * of the contributions of the x and y-coordinate differences for
     * vertex v to the y translation.
     *
     * For example, if a single vertex v is given the xDiffCoeff are 0.5, 0,
     * and the yDiffCoeff are 0, 0, then the transformation will translate
     * a point by half of the difference between the x-coordinate of the
     * vertex to which v is bound and the x-coordinate of the reference
     * location of vertex v.
     *
     * If the coefficient lists are not exactly twice as long as the list
     * of vertices, this method will throw an exception.
     *
     * This object does not take ownership of the AF2RuleExistVertex that
     * are passed into the constructor.  It is the responsibility of the
     * context to ensure that the pointers remain valid as long as this
     * object or a copy, assignment, or clone of it may be used to transform
     * a point.
     *
     * \param vertices a list of AF2RuleExistVertex
     * \param xDiffCoeff a list of double values twice as long as
     *   the list of vertices that defines coefficients affecting
     *   translation in the x direction
     * \param yDiffCoeff a list of double values twice as long as
     *   the list of vertices that defines coefficients affecting
     *   translation in the y direction
     */
    AF2PntTrnsfrmLnrV(std::list<const AF2RuleExistVertex*> vertices,
        std::list<double> xDiffCoeff, std::list<double> yDiffCoeff);

    /**
     * \brief Constructor with vectors
     *
     * This constructor is equivalent to the constructor with lists, but
     * allows passing in vectors instead.  See the constructor with
     * lists for additional documentation.
     *
     * \param vertices a vector of AF2RuleExistVertex
     * \param xDiffCoeff a vector of double values twice as long as
     *   the vector of vertices that defines coefficients affecting
     *   translation in the x direction
     * \param yDiffCoeff a vector of double values twice as long as
     *   the vector of vertices that defines coefficients affecting
     *   translation in the y direction
     */
    AF2PntTrnsfrmLnrV(std::vector<const AF2RuleExistVertex*> vertices,
        std::vector<double> xDiffCoeff, std::vector<double> yDiffCoeff);

    /**
     * \brief Makes and returns an independent copy of this
     * AF2PntTrnsfrmLnrV.
     *
     * Implements AF2PointTransform::clone.  See additional documentation
     * there.
     */
    virtual AF2PntTrnsfrmLnrV* clone() const;

    /**
     * \brief Implements AF2PointTransform::transformPoint
     *
     * See the documentation at AF2PointTransform::transformPoint
     * for the general contract for this method.
     *
     * This implementation consults the binding and computes the
     * differences between the actual bound location and reference
     * location of certain vertices.  The point is then translated
     * from its current position by a vector that is the result of a
     * linear transformation applied to the differences.
     *
     * If one of the vertices that this transformation is supposed to
     * use is not bound to a value in the binding, this method will throw
     * an exception.
     *
     * \param point the coordinates of some 2-dimensional point
     * \param vBinding a binding of reference vertices and edges to points
     *   with actual coordinates and edges between those points
     * \return a point offset from the original point determined by a linear
     *   transformation and the bound locations of certain vertices
     */
    virtual AF2Point2D transformPoint(AF2Point2D const & point,
        AF2Binding const & vBinding) const;
};

#endif
