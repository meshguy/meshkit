// IARoundingFarNlp.cpp
// Interval Assignment for Meshkit
//
// Adapted from
// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_nlp.cpp 1864 2010-12-22 19:21:02Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IARoundingFarNlp.hpp"
#include "IARoundingNlp.hpp"
#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"
#include "IANlp.hpp"

#include <math.h>
#include <algorithm>

// for printf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

/* Idea: a form of IARoundingNlp with larger variable bounds, but still with a natural integer solution.
 x in [1..inf]
 xr = x optimal relaxed solution with objective function fnlp, see IANlp.xpp 
 f is piecewise linear, with corners at integer values. f slopes are unique (we hope)
 Slope definitions
 for x between xl = floor xr and xh = ceil xr, we use the difference in fnlp between xl and xh
 
 case A. xr > ceil g, g is goal I[i]
 for x above xh, 
 let h+ be fnlp ( xh+1 ) - fnlp ( xh )
 let kp be the number of intervals x is above xh 
 then slope = floor kp * h+
 for x below xl, h- = sqrt(11) / 5 h+, and slope = floor km * h-
 all this is weighted by some unique weight
 
 case B. xr in [ floor g, ceil g] 
 h+ = fnlp ( xh+1 ) - fnlp ( xh )
 h- = fnlp ( xl-1 ) - fnlp ( xl ), not used if xl == 1
 
 case C. xr < floor g
 h- = fnlp ( xl-1 ) - fnlp ( xl )
 h+ = sqrt(10)/5 h-
 If g < 2, then h- is unused, and h+ = as in case B

 // representation:
 h0 is weights 0..n-1
 h+ is weights n..2n-1
 h- is weights 2n..
 
*/


// to do
// set h
// uniquify h initially
// randomize h if solution is non-integer
// compute f
// compute f'
// compute f'' == 0

// clean up
// move the re-used part of IANlp, dealing with constraints but not obj, to some underlying class, without it being a TNLP


// the objective function we should use for deriving the weights
double IARoundingFarNlp::f_x_value( double I_i, double x_i ) const
{
  // to do revise
  return x_i > I_i ? 
  (x_i - I_i) / I_i :
  1.103402234045 * (I_i - x_i) / x_i; 
}

// destructor
IARoundingFarNlp::~IARoundingFarNlp() 
{
  data = NULL; 
  ipData = NULL;
  solution = NULL;
}

// constructor
IARoundingFarNlp::IARoundingFarNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr): 
data(data_ptr), ipData(ip_data_ptr), solution(solution_ptr), baseNlp(data_ptr, solution_ptr),
debugging(true), verbose(true) // true
{
  printf("\nIARoundingFarNlp Problem:\n");
  
  // weights for function value
  h0.resize(data->num_variables());
  hp.resize(data->num_variables());
  hm.resize(data->num_variables());
  for (int i = 0; i < data->num_variables(); ++i)
  {
    const double x = ipData->relaxedSolution[i];
    const double xl = floor(x);
    const double xh = xl + 1;
    const double g = data->I[i]; // goal
    
    const double fl = f_x_value(g, xl); 
    const double fh = f_x_value(g, xh);
    const double w = fh - fl;
    h0[i] = w;

    // use sqrt(3) to ensure hp and hm is larger than h0
    const double fp = sqrt(3) * (f_x_value(g,xh+1) - f_x_value(g,xh));
    double fm = 0.;
    
    if (xl > 1)
    {
      fm = sqrt(5) * (f_x_value(g,xl) - f_x_value(g,xl-1));
    }
    else // if (xl == 1)
    {
      // unused anyway, assign it something that won't screw up the distribution
      fm = - sqrt(5) * fabs(f_x_value(g,2) - f_x_value(g,1));
    }

    // case A
    if ( x > ceil(g) )
    {
      hp[i] = fp;
      hm[i] = -fp * sqrt(11) / 5.;      
      assert( fabs(hm[i]) < fabs(hp[i]) );
    }
    // case B
    else if ( x >= floor(g) )
    {
      hp[i] = fp;
      hm[i] = fm;
    }
    // case C
    else // x < floor(g)
    {
      hm[i] = fm;
      hp[i] = -fm * sqrt(7) / 3.;
      assert( fabs(hm[i]) > fabs(hp[i]) );
    }
        
    // invariants for all cases, to produce a super-linear preference for xl or xh
    assert( hp[i] > 0 );
    assert( hp[i] > h0[i] );
    assert( hm[i] < 0 );
    assert( hm[i] < h0[i] );

    /*
    //zzyk experiment with limited f' discontinuity
     // this makes hp < 0 and hm > 0 sometimes, so watch the other asserts
    if ( h0[i] > 0 )
    {
      hp[i] = h0[i]*2.;
      hm[i] = h0[i]*0.5;
    }
    else
    {
      hp[i] = h0[i]*0.5;
      hm[i] = h0[i]*2.;      
    }
     */
    
    if (0 && debugging)
      printf(": %f %f %f", h0[i], hp[i], hm[i]); // 0th+/i, first should be positive, second negative 
  }
  
  // uniquify and scale the weights = slopes
  const double h0_lo = 1;
  const double h0_hi = h0_lo * (10 + 2 * data->num_variables());
  // uniquify preserving h0 < hm, hp
  if (0)
  {
    // todo experiment with ranges
    // want the largest h0 to be less then the smallest hp or hm, in absolute value
    // so that a simple rounding from the relaxed solution is optimal, if feasible
    // 1.e4? 1.e3 makes the jumps in values wash out for the big chain problem, need some jumps?
    // idea, dynamically scale this based on the problem size, h0_hi = 2*num_variables;
    const double hpm_lo = h0_hi + 1.;
    const double hpm_hi = 1.e5;
    h0.uniquify(h0_lo, h0_hi);
    
    // uniquify and scale hp hm to come after h0
    IAWeights hpm();
    hpm = hp;                 // hpm = hp
    hpm.insert( hpm.end(), hm.begin(), hm.end() ); // hmp += hm
    hpm.uniquify(hpm_lo, hpm_hi); 
    hp.assign(hpm.begin(), hpm.begin()+hp.size()); // extract hp
    hm.assign(hpm.begin()+hp.size(), hpm.end());   // extract hm
  }
  // uniquify lumping h0 with hm and hp
  else {
    h0.uniquify(h0, h0_lo, h0_hi);
    
    // uniquify and scale hp hm to come after h0
    IAWeights h0pm;
    hpm = h0;                  // h0pm = h0
    h0pm.insert( h0pm.end(), hp.begin(), hp.end() ); // h0pm += hp
    h0pm.insert( h0pm.end(), hm.begin(), hm.end() ); // h0pm += hm
    h0pm.uniquify(h0_lo, h0_hi); 
    const size_t sz = data->num_variables();
    h0.assign(h0pm.begin(),      h0pm.begin()+sz);   // extract h0
    hp.assign(h0pm.begin()+sz,   h0pm.begin()+2*sz); // extract hp
    hm.assign(h0pm.begin()+2*sz, h0pm.end());        // extract hm
  
  }
  
  if (debugging)
  {
    printf("\n");
    for (int i = 0; i < data->num_variables(); ++i)
    {
      const double x = ipData->relaxedSolution[i];
      printf("x %d (relaxed %f, goal %e), h0 %10.4f hp %10.4f hm %10.4f\n", i, x, data->I[i], h0[i], hp[i], hm[i]);
    }
  }
}


// returns the size of the problem
bool IARoundingFarNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  bool base_ok = baseNlp.get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  nnz_h_lag = 0; // obj is piecewise linear. We could define h near the corners, but lets try zero first. nele_hess
  return true && base_ok;
  // need to change this if there are more variables, such as delta-minuses
}

// returns the variable bounds
bool IARoundingFarNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  const bool base_ok = baseNlp.get_bounds_info(n, x_l, x_u, m, g_l, g_u);

  // test various ranges  
  /*
  for (Index i=0; i<n; ++i) 
  {
  */
  
    /*
    // debug: x in [xl-1, xh+1]
    const double x = ipData->relaxedSolution[i];
    const double xl = floor(x);
    assert(xl >= 1.);
    if (xl > 1)
      x_l[i] = xl-1;
    else
      x_l[i] = 1;
    x_u[i] = xl + 2;
    */
    
    /*
    // debug: x in [xl,xh]
    x_l[i] = xl;
    x_u[i] = xl + 1.;
    */
    
    //if (debugging)
    //    printf("x %d (%f) in [%d,%d]\n", i, x, (int) x_l[i], (int) x_u[i]);
    
  /*
  }
  */
  
  return base_ok && true;
}

// returns the initial point for the problem
bool IARoundingFarNlp::get_starting_point(Index n, bool init_x, Number* x_init,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Minimal info is starting values for x, x_init
  // Improvement: we can provide starting values for the dual variables if we wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize x to the relaxed solution
  for (Index i=0; i<n; ++i) 
  {
    x_init[i] = ipData->relaxedSolution[i];
  }

  return true;
}

inline 
double IARoundingFarNlp::get_f_xl(int i) const
{
  return (h0[i] > 0.) ? 0. : -h0[i];
}

inline 
double IARoundingFarNlp::get_f_xh(int i) const
{
  return (h0[i] > 0.) ? h0[i] : 0.;
}


// returns the value of the objective function
bool IARoundingFarNlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == data->num_variables());

  double obj = 0.;
  for (Index i = 0; i<n; ++i)
  {
    // find which interval we're in
    const double xl = ipData->get_xl(i);
    const double xh = xl + 1;
    const double x_i = x[i];  
    
    double obj_i(-1.); // bad value
    if (x_i >xh)  // above xh
    {
      // the slope of the kth interval is k hp[i]
      // f = f(xh) + sum_{k=1:kp} k hp + kp_frac (kp+1) hp
      //   = f(xh) + kp (kp+1) / 2  hp + kp_frac (kp+1) hp
      //   = f(xh) + hp (kp+1) (kp_frac + kp / 2 )
      int kp_int;
      double kp_frac;
      ipData->get_kp(i, x_i, kp_int, kp_frac);
      obj_i = get_f_xh(i) + hp[i] * (kp_int + 1) * (kp_frac + ((double) kp_int) / 2.);
      assert(obj_i >= 0.);
    }
    else if (x_i < xl) // below xl
    {
      // the slope of the kth interval is k hm[i]
      // f = f(xl) + sum_{k=1:km} k hm + km_frac (km+1) hm
      //   = f(xl) + km (km+1) / 2  hm + km_frac (km+1) hm
      //   = f(xl) + hm (km+1) (k_frac + km / 2 )
      int km_int;
      double km_frac;
      ipData->get_km(i, x_i, km_int, km_frac);
      obj_i = get_f_xl(i) + fabs(hm[i]) * (km_int + 1) * (km_frac + ((double) km_int) / 2.); 
      assert(obj_i >= 0.);
    }
    else // middle
    { 
    // zzyk - testing truncating x within an interval around xl, continuity
    /*
    static double max_over = 0.;
    static double max_under = 0.;
    
    if (x_i < xl)
    {
      if ( xl - x_i > max_under )
      {
        max_under = xl - x_i;
        printf("max under %g\n", max_under);
      }
      x_i = xl;
    }
    else if (x_i > xh)
    {
      if ( x_i - xh > max_over )
      {
        max_over = x_i - xh;
        printf("max OVER %g\n", max_over);
      }
      x_i = xh;      
    } */
    
      obj_i = (h0[i] > 0) ? h0[i] * ( x_i - xl ) : -h0[i] * ( xh - x_i );
      assert(obj_i >= 0.);

    } 
    assert(obj_i >= 0.);

    obj += obj_i;
  }
  assert(obj >= 0.);
  obj_value = obj;
  printf(" %f", obj_value); // debug progress, is this going anywhere?
  
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IARoundingFarNlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == data->num_variables());
  for (Index i = 0; i<n; ++i)
  {
    // find which interval we're in
    if (x[i] > ipData->get_xh(i)) // above xh
    {
      // the slope of the kth interval is k hp[i]
      int kp_int;
      double kp_frac;
      ipData->get_kp(i, x[i], kp_int, kp_frac);
      grad_f[i] = (1+kp_int) * hp[i];
    }
    else if (x[i] < ipData->get_xl(i)) // below xl
    {
      // the slope of the kth interval is k hm[i]
      int km_int;
      double km_frac;
      ipData->get_km(i, x[i], km_int, km_frac);
      grad_f[i] = (1+km_int) * hm[i];
    }
    else // between xl and xh
    { 
      // slope is h0
      grad_f[i] = h0[i];
    }
  }
  return true;
}

// return the value of the constraints: g(x)
bool IARoundingFarNlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  return baseNlp.eval_g(n, x, new_x, m, g);
}

// return the structure or values of the jacobian
bool IARoundingFarNlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  return baseNlp.eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
}

//return the structure or values of the hessian
bool IARoundingFarNlp::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
	// because the constraints are linear
  // and the objective function is piecewise linear
  // the hessian for this problem is actually empty

  // to do: may need to fill with zeros so ipopt doesn't think it is globally linear
  assert(nele_hess == 0); 
  
  // This is a symmetric matrix, fill the lower left triangle only.
  if (values == NULL) {
    // return the structure. 
    ;
  }
  else {
    // return the values. 
    ;
  }

  return true;
}

void IARoundingFarNlp::finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
{
  baseNlp.finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
}
