// IPData.cpp
// Interval Assignment Data for Meshkit
// Integer Program - making intervals into integer values 
//

#include "IPData.hpp"
#include <vector>
#include <math.h>
#include <cstdio>
#include <assert.h>

void IPData::initialize(const std::vector<double> &relaxed_solution)
{
  const std::size_t num_variables = relaxed_solution.size();
  // fill varIntegerBound and oldBound with zeros
  varIntegerBound.clear();
  varIntegerBound.resize(num_variables,0);
  oldBound.clear();
  oldBound.resize(num_variables,0);
  relaxedSolution = relaxed_solution; // vector copy
}

 
void IPData::constrain_integer(const int i_nonint, const int x_bound)
{
  oldBound[i_nonint] = varIntegerBound[i_nonint];;
  varIntegerBound[i_nonint] = (double) x_bound;
}

void IPData::round_solution(std::vector<double> &solution)
{
  for (int i = 0; i < solution.size(); ++i )
  {
    solution[i] = floor( solution[i] + 0.5 );
  }
}



bool IPData::is_integer(double x)
{ 
  // beware sides with more than 100 curves = 1 / 1.e-2
  return (fabs( x - floor(x+0.5) ) < 1.e-2);
}

bool IPData::solution_is_integer(const std::vector<double> &solution)
{
  bool all_int = true;
  for (int i = 0; i<solution.size(); ++i)
  {
    if (!is_integer(solution[i]))
    {
      all_int = false;
      if (1)
      {
        printf(" x %d is %f NON-INTEGER\n", i, solution[i]);
      }
      // return false; 
    }
  }
  return all_int;
}

void IPData::get_frac( double x, int &integer, double &frac) 
{
  assert(x>=0.);
  integer = floor(x); 
  frac = x - integer;
  assert( frac < 1. );
  assert( frac >= 0. );
}

double IPData::get_xl( int i ) const
{
  return floor(relaxedSolution[i]);
}

double IPData::get_xh( int i ) const
{
  return floor(relaxedSolution[i])+1;
}

double IPData::get_km( int i, double x ) const
{
  const double y =  get_xl(i) - x;
  return y>0. ? y : 0.;
}

void IPData::get_km( int i, double x, int &km_integer, double &km_frac ) const
{
  get_frac( get_km(i,x), km_integer, km_frac);
}

double IPData::get_kp( int i, double x ) const
{
  const double y = x - get_xh(i);
  return y > 0 ? y : 0.;
}

void IPData::get_kp( int i, double x, int &kp_integer, double &kp_frac ) const
{
  get_frac( get_kp(i,x), kp_integer, kp_frac);
}

