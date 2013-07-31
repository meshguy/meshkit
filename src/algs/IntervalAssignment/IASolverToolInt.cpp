// IASolverToolInt.cpp

#include "IASolverToolInt.hpp"
#include "IAData.hpp"
#include "IASolution.hpp"
#include "IPData.hpp"
#include "IAWeights.hpp"

#include <vector>
#include <math.h>
#include <cstdio>
#include <assert.h>

namespace MeshKit {
 
IASolverToolInt::IASolverToolInt() 
  : IASolverTool(), 
  debuggingToolInt(false) 
  {ipData = new IPData;}
  
IASolverToolInt::IASolverToolInt( const IAData *ia_data, IASolution *ia_solution, bool debugging_set) 
  : IASolverTool(ia_data, ia_solution, debugging_set), 
  debuggingToolInt( debugging_set ) 
  {ipData = new IPData;}
  
IASolverToolInt::IASolverToolInt(const bool debugging_set) 
  : IASolverTool(debugging_set), 
  debuggingToolInt(debugging_set) 
  {ipData = new IPData;}

  
IASolverToolInt::~IASolverToolInt()
{
  ipData = NULL;
}
  
  
void IASolverToolInt::round_solution()
{
  assert( iaSolution );
  for (unsigned int i = 0; i < iaSolution->x_solution.size(); ++i )
  {
    assert(iaSolution->x_solution[i] >= 0.);
    iaSolution->x_solution[i] = floor( iaSolution->x_solution[i] + 0.5 );
  }
}
  
  
bool IASolverToolInt::solution_is_integer(const bool print_non_integer)
{
  assert( iaSolution );
  
  int num_non_int = 0;
  class NonInt
  {
  public:
    double x;
    int i;
    double dist () const
    {
      double xm = fabs( x - floor(x));
      double xp = fabs( floor(x) + 1. - x);
      // printf("dist: x %f floor %f xm %f xp %f\n", x, floor(x), xm, xp); // debug
      if (xm > xp)
        return xp;
      else
        return xm;
    }
    NonInt() : x(0.), i(-1) {}
    NonInt(const double x_val, const unsigned int i_val) : x(x_val), i( (int) i_val ) {}
  };
  
  NonInt worst_non_int;
  
  bool all_int = true;  
  for (unsigned int i = 0; i<iaSolution->x_solution.size(); ++i)
  {
    const double x = iaSolution->x_solution[i];
    if (!is_integer(x))
    {
      all_int = false;
      if (print_non_integer)
      {
        NonInt xnon(x,i);
        if (xnon.dist() > worst_non_int.dist())
          worst_non_int = xnon;
        if (num_non_int++ == 0)
        {
          printf("\nChecking solution integrality\n");
        }
        const double goal = iaData->I[i];
        printf(" x %d (g %f) is %f NON-INTEGER\n", i, goal, x);
      }
      // shortcut return if we're not printing anything
      else
        return false; 
    }
  }
  if (print_non_integer)
  {
    if (num_non_int == 0)
    {
      printf("all integer\n");
    }
    else
    {
      printf("%d of %lu NON-INTEGERS\n", num_non_int, iaSolution->x_solution.size());
      printf("worst was x_%d = %g\n", worst_non_int.i, worst_non_int.x );
    }
  }
  return all_int;
}


void IASolverToolInt::get_frac( double x, int &integer, double &frac) const
{
  assert(x>=0.);
  integer = floor(x); 
  frac = x - integer;
  assert( frac < 1. );
  assert( frac >= 0. );
}

double IASolverToolInt::get_xl( int i ) const
{
  return floor( ipData->relaxedSolution[i] );
}

double IASolverToolInt::get_xh( int i ) const
{
  return floor( ipData->relaxedSolution[i] )+1;
}

double IASolverToolInt::get_km( int i, double x ) const
{
  const double y =  get_xl(i) - x;
  return y>0. ? y : 0.;
}

void IASolverToolInt::get_km( int i, double x, int &km_integer, double &km_frac ) const
{
  get_frac( get_km(i,x), km_integer, km_frac);
}

double IASolverToolInt::get_kp( int i, double x ) const
{
  const double y = x - get_xh(i);
  return y > 0 ? y : 0.;
}

void IASolverToolInt::get_kp( int i, double x, int &kp_integer, double &kp_frac ) const
{
  get_frac( get_kp(i,x), kp_integer, kp_frac);
}
  
  
  
bool IASolverToolInt::randomize_weights_of_non_int(IAWeights* weights, const double rand_factor)
{
    
  for (size_t i=0; i<iaSolution->x_solution.size(); ++i) 
  {
    const double x = iaSolution->x_solution[i];
    if (!is_integer(x))
    {
      const double d = IAWeights::rand_excluded_middle();
      (*weights)[i] *= 1. + rand_factor * d;
    }
  }
  return true;
}

} // namespace MeshKit
