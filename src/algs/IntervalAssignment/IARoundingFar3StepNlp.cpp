// IARoundingFar3StepNlp.cpp
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

#include "IARoundingFar3StepNlp.hpp"
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
 then slope = (1 + floor(kp)/2) * h+
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
 h0      is weights 0..n-1
 h+ = hp is weights n..2n-1
 h- = hm is weights 2n..
 
 // variable layout
 x01 is variables  n..2n-1  : x01_start
 xp is variables  2n..3n-1  : xp_start
 xm is varaibles  3n..4n-1  : xm_start
 
 // constrataint layout
 sum-equal  
 sum-even                : sum_even_start 
 x = xl + x01 + xp - xm  : x_constraint_start
 or g(j)= x - x01 - xp + xm = xl 
*/



// the objective function we should use for deriving the weights
double IARoundingFar3StepNlp::f_x_value( double I_i, double x_i ) const
{
  // to do revise
  return x_i > I_i ? 
  (x_i - I_i) / I_i :
  1.103402234045 * (I_i - x_i) / x_i; 
}

// destructor
IARoundingFar3StepNlp::~IARoundingFar3StepNlp() 
{
  data = NULL; 
  ipData = NULL;
  solution = NULL;
}

// constructor
IARoundingFar3StepNlp::IARoundingFar3StepNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr): 
data(data_ptr), ipData(ip_data_ptr), solution(solution_ptr), baseNlp(data_ptr, solution_ptr),
x01_start( data_ptr->num_variables() ), 
xp_start( 2 * data_ptr->num_variables() ), 
xm_start( 3 * data_ptr->num_variables() ),
base_n ( data_ptr->num_variables() ),
base_m ( (int) ( data_ptr->constraints.size() + data->sumEvenConstraints.size() ) ),
problem_n( (int) 4 * data_ptr->num_variables() ),
problem_m( (int) (data_ptr->constraints.size() + data_ptr->sumEvenConstraints.size() + data_ptr->num_variables())),
sum_even_start( (int) data_ptr->constraints.size()), 
x_constraint_start( (int) (data_ptr->constraints.size() + data_ptr->sumEvenConstraints.size())), 
hess_option(ZERO),
debugging(true), verbose(true) // true
{
  printf("\nIARoundingFar3StepNlp Problem:\n");
  
  
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
    IARoundingNlp::uniquify_weights(h0, h0_lo, h0_hi);
    
    // uniquify and scale hp hm to come after h0
    std::vector<double> hpm( hp );                 // hpm = hp
    hpm.insert( hpm.end(), hm.begin(), hm.end() ); // hmp += hm
    IARoundingNlp::uniquify_weights(hpm, hpm_lo, hpm_hi); 
    hp.assign(hpm.begin(), hpm.begin()+hp.size()); // extract hp
    hm.assign(hpm.begin()+hp.size(), hpm.end());   // extract hm
  }
  // uniquify lumping h0 with hm and hp
  else {
    IARoundingNlp::uniquify_weights(h0, h0_lo, h0_hi);
    
    // uniquify and scale hp hm to come after h0
    std::vector<double> h0pm( h0 );                  // h0pm = h0
    h0pm.insert( h0pm.end(), hp.begin(), hp.end() ); // h0pm += hp
    h0pm.insert( h0pm.end(), hm.begin(), hm.end() ); // h0pm += hm
    IARoundingNlp::uniquify_weights(h0pm, h0_lo, h0_hi); 
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

bool IARoundingFar3StepNlp::randomize_weights_of_non_int()
{
  
  for (int i=0; i<solution->x_solution.size(); ++i) 
  {
    const double x = solution->x_solution[i];
    if (!IPData::is_integer(x))
    {
      double d = IARoundingNlp::rand_excluded_middle();
      h0[i] *= 1. + 0.03124234 * d;
      d = IARoundingNlp::rand_excluded_middle();
      hm[i] *= 1. + 0.02894834 * d;
      d = IARoundingNlp::rand_excluded_middle();
      hp[i] *= 1. + 0.02745675 * d;
    }
  }
  return true;
}


// returns the size of the problem
bool IARoundingFar3StepNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  bool base_ok = baseNlp.get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  
  m = problem_m;
  n = problem_n;

  // nnz_jac_g == nele_jac in eval_jac_g
  // four coefficients for each g() = x - x01 - xp + xm = xl const constraint
  nnz_jac_g += 4 * base_n; 
  
  // nnz_h_lag == nele_hess in eval_h
  if (hess_option == ZERO)
    nnz_h_lag = 0.;
  else if (hess_option == ROUNDED)
    nnz_h_lag =  2 * data->num_variables(); // one for each xp, xm
  
  return true && base_ok;
  // need to change this if there are more variables, such as delta-minuses
}

// returns the variable bounds
bool IARoundingFar3StepNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  const bool base_ok = baseNlp.get_bounds_info(base_n, x_l, x_u, base_m, g_l, g_u);
  
  for (Index i=0; i<data->num_variables(); ++i)
  {
    // x01
    const int x01i = x01_start + i;
    x_l[x01i] = 0;
    x_u[x01i] = 1;
    
    // xp
    const int xpi = xp_start + i;
    x_l[xpi] = 0;
    x_u[xpi] = MESHKIT_IA_upperUnbound;
    
    // xm
    const int xmi = xm_start + i;
    x_l[xmi] = 0;
    x_u[xmi] = ipData->get_xl(i) - 1;
    assert( x_u[xmi] >= 0. );
  }

  // x = xl + x01 + xp - xm  : x_constraint_start
  // or g(j)= x - x01 - xp + xm = xl 
  for (Index j=0; j<data->num_variables(); ++j)
  {
    const int xcj = x_constraint_start + j;
    const int xl = ipData->get_xl(j);
    g_l[xcj] = xl; 
    g_u[xcj] = xl;
  }
  
  return base_ok && true;
}

// returns the initial point for the problem
bool IARoundingFar3StepNlp::get_starting_point(Index n, bool init_x, Number* x_init,
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
  for (Index i=0; i<base_n; ++i) 
  {
    const double xr = ipData->relaxedSolution[i];
    x_init[i] = xr;
    x_init[x01_start + i] = xr - floor(xr);
    x_init[xp_start + i] = 0.;
    x_init[xm_start + i] = 0.;
  }

  return true;
}

inline 
double IARoundingFar3StepNlp::get_f_xl(int i) const
{
  return (h0[i] > 0.) ? 0. : -h0[i];
}

inline 
double IARoundingFar3StepNlp::get_f_xh(int i) const
{
  return (h0[i] > 0.) ? h0[i] : 0.;
}


// returns the value of the objective function
bool IARoundingFar3StepNlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == problem_n);

  double obj = 0.;
  for (Index i = 0; i<base_n; ++i)
  {
    // x itself contributes nothing, just its components x01, xp, xm

    // xp
    {
      const double xpv = x[xp_start + i];
      int kp_int;
      double kp_frac;
      IPData::get_frac( xpv, kp_int, kp_frac );
      const double objp = hp[i] * ( kp_int * ( 1. + (kp_int + 1.) / 4. ) + kp_frac * ( 1. + kp_int / 2. ) );
      assert( xpv >= 0. );
      assert(objp >= 0.);    
      obj += objp;
    }

    // xm
    {
      const double xmv = x[xm_start + i];
      int km_int;
      double km_frac;
      IPData::get_frac( xmv, km_int, km_frac );
      const double objm = -hm[i] * ( km_int * ( 1. + (km_int + 1.) / 4. ) + km_frac * ( 1. + km_int / 2. ) );
      assert( xmv >= 0. );
      assert(objm >= 0.);    
      obj += objm;
    }

    // x01
    {
      const double x01v = x[x01_start + i];
      const double obj01 = (h0[i] > 0) ? h0[i] * x01v : h0[i] * ( x01v - 1. );
      assert(x01v >= 0.);
      assert(x01v <= 1.);
      assert(obj01 >= 0.);
      obj += obj01;
    }
  }
  assert(obj >= 0.);
  obj_value = obj;
  printf(" %f", obj_value); // debug progress, is this going down?
  
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IARoundingFar3StepNlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == problem_n);
  for (Index i = 0; i<base_n; ++i)
  {
    // x
    grad_f[i] = 0;
    
    // xp
    {
      const double xpv = x[xp_start + i];
      int kp_int;
      double kp_frac;
      IPData::get_frac( xpv, kp_int, kp_frac );    
      grad_f[xp_start+i] = hp[i] * ( 1. + kp_int / 2. );
    }
    
    // xm
    {
      const double xmv = x[xm_start + i];
      int km_int;
      double km_frac;
      IPData::get_frac( xmv, km_int, km_frac );
      grad_f[xm_start+i] = -hm[i] * ( 1. + km_int / 2. );
    }
    
    // x01
    grad_f[x01_start+i] = h0[i];
  }
  
  return true;
}

// return the value of the constraints: g(x)
bool IARoundingFar3StepNlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  const bool base_ok = baseNlp.eval_g(base_n, x, new_x, base_m, g);
  
  // x_constraints
  for (int i=0; i<base_n; ++i)
  {
    // g(k)= x - x01 - xp + xm = xl 
    const double xv = x[i];
    const double x01v = x[x01_start+i];
    const double xpv = x[xp_start+i];
    const double xmv = x[xm_start+i];
    const double g_k = xv - x01v - xpv + xmv;
    const int k = x_constraint_start + i;
    g[k] = g_k;    
  }
  
  return base_ok && true;
}

// return the structure or values of the jacobian
bool IARoundingFar3StepNlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  int base_nele_jac = baseNlp.get_neleJac();
  bool base_ok = baseNlp.eval_jac_g(base_n, x, new_x, base_m, base_nele_jac, iRow, jCol, values);
  base_nele_jac = baseNlp.get_neleJac();
  int k = base_nele_jac;

  if (values == NULL) 
  {
    // return the structure of the jacobian

    // x_constraints
    for (int i=0; i<base_n; ++i)
    {
      // g(k)= x - x01 - xp + xm = xl const
      int xi = i + x_constraint_start;
      iRow[k] = xi;
      jCol[k] = i;             // x
      ++k;
      iRow[k] = xi;
      jCol[k] = x01_start + i; // x01
      ++k;
      iRow[k] = xi;
      jCol[k] = xp_start + i;  // xp
      ++k;
      iRow[k] = xi;
      jCol[k] = xm_start + i;  // xm
      ++k;
    }
  }
  else
  {
    // return the values of the jacobian of the constraints

    // x_constraints
    for (int i=0; i<base_n; ++i)
    {
      // g(k)= x - x01 - xp + xm = xl const
      values[k++] =  1;  // x
      values[k++] = -1;  // x01
      values[k++] = -1;  // xp
      values[k++] =  1;  // xm
    }

  }
  return base_ok && true;
}

//return the structure or values of the hessian
bool IARoundingFar3StepNlp::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
	// because the constraints are linear
  // the hessian for them is actually empty

  // get_nlp_info specified the number of non-zeroes, nele_hess

  
  // This is a symmetric matrix, fill the lower left triangle only.
  if (values == NULL) {
    // return the structure. 
    if (hess_option == ZERO)
    {
      ; // empty
    }
    else if (hess_option == ROUNDED)
    {
      // Since the objective function is separable, only the diagonal is non-zero
      int k = 0;
      for ( int i = 0; i < base_n; ++i )
      {
        iRow[k] = xp_start + i;
        jCol[k] = xp_start + i;
        ++k;
      }
      for ( int i = 0; i < base_n; ++i )
      {
        iRow[k] = xm_start + i;
        jCol[k] = xm_start + i;
        ++k;
      }
    }
  }
  else {
    // return the values. 
    int k = 0;
    if (hess_option == ZERO)
    {
      ; // empty
    }
    else if (hess_option == ROUNDED)
    {
      for ( int i = 0; i < base_n; ++i )
      {
        values[k] = hp[i] / 2.; 
        // but linearly fade to zero towards xh
        const double xp = x[xp_start + i];
        if (xp < 0.5)
          values[k] *= 2. * xp;
        
        ++k;
      }
      for ( int i = 0; i < base_n; ++i )
      {
        values[k] = -hm[i] / 2.; 
        const double xm = x[xm_start + i];
        if (xm < 0.5)
          values[k] *= 2. * xm;

        ++k;
      }
      
    }
    assert( k == nele_hess );
  }

  
  return true;
}

void IARoundingFar3StepNlp::finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
{
  baseNlp.finalize_solution(status, base_n, x, z_L, z_U, base_m, g, lambda, obj_value, ip_data, ip_cq);
  // ignore the other constraints and variables, as they were just bookkeepings.
}
