#ifndef MESHKIT_VECUTIL_HPP
#define MESHKIT_VECUTIL_HPP

#include <math.h>

namespace MeshKit {

class ModelEnt;
class MKCore;

/** \class VecUtil VecUtil.hpp "meshkit/VecUtil.hpp"
 * \brief A class of vector utilities
 */

class VecUtil 
{
public:
    /** \brief Dot product
     *
     * \param a First vector
     * \param b Second vector
     * \return Dot project
     */
  static double dot(double *a, double *b);

    /** \brief Length squared of a vector
     * \param a Vector
     * \return Length squared
     */
  static double length_sq(double *a);

    /** \brief Distance between two vectors
     * \param a Vector 1
     * \param b Vector 2
     * \return Distance
     */
  static double dist2(double *a, double *b);

    /** \brief Normalize a vector
     * \param a Vector
     */
  static void normalize(double *a);

    /** \brief Cross product
     * \param a Vector 1
     * \param b Vector 2
     * \param c Cross product
     */
  static void cross(double *a, double *b, double *c);

    //! Constant PI
  static double PI;

    //! Constant 2*pi
  static double TWO_PI;
};
    
inline double VecUtil::dot(double *a, double *b) {return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];}

inline double VecUtil::length_sq(double *a) {return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
    
inline double VecUtil::dist2(double *a, double *b) 
{
  return sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1])
              + (b[2]-a[2])*(b[2]-a[2]));
}

inline void VecUtil::normalize(double *a) 
{
  double lsq = length_sq(a);
  if (lsq == 0.0) return;
  lsq = 1.0 / sqrt(lsq);
  for (int i = 0; i < 3; i++) a[i] *= lsq;
}

inline void VecUtil::cross(double *a, double *b, double *c) 
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

} // namespace MeshKit

#endif
