// IAIntParabolaNlp.hpp
// Interval Assignment for Meshkit
//
// use a cosine-"wave" function as a constraint, so that the only feasible solutions are integer
// use a cosine-"wave" with twice the period to enforce the evenality constraints

#ifndef MESHKIT_IA_IAINTPARABOLANLP_HP
#define MESHKIT_IA_IAINTPARABOLANLP_HP

#include "IAIntWaveNlp.hpp"

#include <math.h>
#include <limits.h>

namespace MeshKit 
{

class IAIntParabolaNlp : public IAIntWaveNlp
{
  // first set of functions required by TNLP
public:
  /** default constructor */
  IAIntParabolaNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr,
               const bool set_silent = true) :
    IAIntWaveNlp(data_ptr, ip_data_ptr, solution_ptr, set_silent)
  {}

  /** default destructor */
  virtual ~IAIntParabolaNlp() {}

protected:

  // nearest integer
  const double nearest_int( const double x )
  {
    const double xm = floor(x);
    const double xp = ceil(x);
    return (fabs(xm - x) < fabs(xp - x)) ? xm : xp;
  }
  const double nearest_even( const double s )
  {
    return 2. * nearest_int( s / 2. );
  }
  const double delta_x( const double x )
  {
    return x - nearest_int(x);
  }
  const double delta_s( const double s )
  {
    return s - nearest_even(s);
  }

  virtual double eval_g_int_x( const double x )
  {
    const double d = delta_x(x);
    return 1. - d * d;
  }
  virtual double eval_g_int_s( const double s )
  {
    const double d = delta_s(s);
    return 1. - d * d;
  }
  virtual double eval_jac_int_x( const double x )
  {
    const double d = delta_x( x );
    return -2. * d; // d' is 1
  }
  virtual double eval_jac_int_s( const double s )
  {
    double d = nearest_even(s);
    return -2. * d; // d' is 1
  }
  virtual double eval_hess_int_x( const double x )
  {
    return -2.; 
  }
  virtual double eval_hess_int_s( const double s )
  {
    return -2.;
  }
  
};


} // namespace MeshKit

#endif
