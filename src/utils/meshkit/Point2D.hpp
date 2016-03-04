/*
 * Point2D.hpp
 *
 * An immutable point with two double coordinates.
 *
 * The coordinates are named x and y.
 */
#ifndef MESHKIT_UTILS_POINT2D_HPP
#define MESHKIT_UTILS_POINT2D_HPP

// C++
#include <ostream>

class Point2D
{
  private:

    double x, y;

  public:

    /**
     * \brief No-argument constructor
     *
     * Construct a point at (0, 0).
     */
    Point2D();

    /**
     * \brief Standard constructor
     *
     * Construct a point at the specified coordinates.
     *
     * \param xVal the x coordinate of the point
     * \param yVal the y coordinate of the point
     */
    Point2D(double xVal, double yVal);

    /**
     * \brief Get the value of the x coordinate.
     */
    double getX() const;

    /**
     * \brief Get the value of the y coordinate.
     */
    double getY() const;
};

bool operator==(const Point2D & onePnt, const Point2D & otherPnt);
bool operator!=(const Point2D & onePnt, const Point2D & otherPnt);
std::ostream& operator<<(std::ostream& outStream, const Point2D & point);

#endif
