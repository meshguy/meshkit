#ifndef BASIC_MATH_H
#define BASIC_MATH_H

#include <math.h>
#include <values.h>
#include <iostream>
using namespace std;

#define ANGLE_IN_DEGREES  0
#define ANGLE_IN_RADIANS  1

#ifdef USE_BOOST_LIBS
#include <boost/utility.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include <boost/array.hpp>
typedef boost::array<double, 3> Vec3D;
typedef boost::array<double, 3> Point3D;
#endif

template<class DataType, int n>
class Array
{
public:
  double &operator[](int i)
  {
    return data[i];
  }
  double operator[](int i) const
  {
    return data[i];
  }
private:
  DataType  data[n];
};

typedef Array<double,3> Point3D;
typedef Array<double,4> Array4D;

namespace Math
{
inline Point3D create_vector( const Point3D &head, const Point3D &tail)
{
  Point3D xyz;
  xyz[0] = head[0] - tail[0];
  xyz[1] = head[1] - tail[1];
  xyz[2] = head[2] - tail[2];
  return xyz;
}


inline double length( const Point3D &A, const Point3D &B)
{
   double dx = A[0] - B[0];
   double dy = A[1] - B[1];
   double dz = A[2] - B[2];
   return sqrt( dx*dx + dy*dy + dz*dz );
}

inline double length2( const Point3D &A, const Point3D &B)
{
   double dx = A[0] - B[0];
   double dy = A[1] - B[1];
   double dz = A[2] - B[2];
   return dx*dx + dy*dy + dz*dz;
}

inline double magnitude( const Point3D &A )
{
  return sqrt( A[0]*A[0] + A[1]*A[1] + A[2]*A[2] );
}

inline double dot_product( const Point3D &A, const Point3D &B)
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

inline Point3D cross_product( const Point3D &A, const Point3D &B)
{
  Point3D C;
  C[0] = A[1]*B[2] - A[2]*B[1];
  C[1] = A[2]*B[0] - A[0]*B[2];
  C[2] = A[0]*B[1] - A[1]*B[0];
  return C;
}

///////////////////////////////////////////////////////////////////////////////

inline double getVectorAngle( const Point3D &A, const Point3D &B)
{
  double AB = dot_product(A,B);
  double Am = magnitude(A);
  double Bm = magnitude(B);

  if( Am < 1.0E-15 || Bm < 1.0E-15) return 0.0;

  double x = AB/(Am*Bm);

  if( x > 1.0) x = 1.0;
  if( x < -1.0) x = -1.0;

  return 180*acos(x)/M_PI;
}

//////////////////////////////////////////////////////////////////////////////

template <class T, size_t n>
inline double getAngle(const Array<T, n> &VecA, const Array<T, n> &VecB, 
                int unit_measure)
  {
      double Abar, Bbar, theta;
      Abar = magnitude(VecA);
      Bbar = magnitude(VecB);

      if (Abar < 1.0E-10 || Bbar < 1.0E-10) {
          cout << " Warning: Error in Angle calculation " << endl;
          cout << " Magnitude of Vector A is " << Abar << endl;
          cout << " Magnitude of Vector B is " << Bbar << endl;
          return 0.0;
      }

      double value = dot_product(VecA, VecB) / (Abar * Bbar);

      if (value > 1.0) value = 1.0;
      if (value < -1.0) value = -1.0;

      theta = acos(value);

      if (unit_measure == ANGLE_IN_DEGREES) theta *= (180.0 / M_PI);

      return theta;
  }

  //////////////////////////////////////////////////////////////////////////////

template <class T, size_t n>
inline T getAngle(const Array<T, n> &pa, const Array<T, n> &pb,
const Array<T, n> &pc, int unit_measure = 0)
{
      Array<T, n> VecA = create_vector(pb, pa);
      Array<T, n> VecB = create_vector(pc, pa);
      return getAngle(VecA, VecB, unit_measure);
}


}


#endif

