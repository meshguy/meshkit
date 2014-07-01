// IAIntCosNlp.hpp
// Interval Assignment for Meshkit
//
// use a cosine-"wave" function as a constraint, so that the only feasible solutions are integer
// use a cosine-"wave" with twice the period to enforce the evenality constraints

#ifndef MESHKIT_IA_IAINTCOSNLP_HP
#define MESHKIT_IA_IAINTCOSNLP_HP

#include "meshkit/IAIntWaveNlp.hpp"

#include <math.h>
#include <limits.h>

namespace MeshKit 
{

class IAIntCosNlp : public IAIntWaveNlp
{
  // first set of functions required by TNLP
public:
  /** default constructor */
  IAIntCosNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr,
               const bool set_silent = true) :
    IAIntWaveNlp(data_ptr, ip_data_ptr, solution_ptr, set_silent),
    PI( 3.1415926535897932384626433832795 )
  {}

  /** default destructor */
  virtual ~IAIntCosNlp() {}

protected:
  
  const double PI;
  // or PI( 2. * acos(0.0) )
  
  virtual double eval_g_int_x( const double x )
  {
    return cos( 2. * PI * x );
  }
  virtual double eval_g_int_s( const double s )
  {
    return cos( PI * s );
  }
  virtual double eval_jac_int_x( const double x )
  {
    return -2. * PI * sin( 2. * PI * x );    
  }
  virtual double eval_jac_int_s( const double s )
  {
    return -PI * sin( PI * s );
  }
  virtual double eval_hess_int_x( const double x )
  {
    return -4. * PI * PI * cos( 2. * PI * x );
  }
  virtual double eval_hess_int_s( const double s )
  {
    return -PI * PI * cos( PI * s );    
  }
  
};


} // namespace MeshKit

#endif
