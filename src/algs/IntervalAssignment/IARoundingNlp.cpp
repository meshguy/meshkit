// IARoundingNlp.cpp
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

#include "IARoundingNlp.hpp"
#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"
#include "IANlp.hpp"
#include "IASolverToolInt.hpp"

#include <math.h>
#include <limits>
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


namespace MeshKit 
{
  
bool IARoundingNlp::randomize_weights_of_non_int()
{
  IASolverToolInt sti(data, solution);
  return sti.randomize_weights_of_non_int(&weights, 0.023394588);
}


// the objective function we should use for weights
double IARoundingNlp::f_x_value( double I_i, double x_i )
{
  return x_i > I_i ? 
  (x_i - I_i) / I_i :
  1.103402234045 * (I_i - x_i) / x_i;
  // expected is 1/100 to 4 or so

  // the range of magnitudes of these weights is too high, the weights are too non-linear if we take them to a power
  // was
  // const double fh = IANlp::eval_R_i(data->I[i], xl+1);
}

// constructor
IARoundingNlp::IARoundingNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr,
  const bool set_silent): 
data(data_ptr), ipData(ip_data_ptr), solution(solution_ptr), baseNlp(data_ptr, solution_ptr, set_silent),
  weights(),
silent(set_silent), debugging(false), verbose(true) // true
{
  if (!silent)
  { 
    printf("\nIARoundingNLP Problem size:\n");
    printf("  number of variables: %lu\n", data->I.size());
    printf("  number of constraints: %lu\n\n", data->constraints.size());
  }
  
  weights.resize(data->I.size());
  if (0 && debugging)
    printf("raw linear +x weights: ");
  for (int i = 0; i < data->num_variables(); ++i)
  {
    const double x = ipData->relaxedSolution[i];
    const double xl = floor(x);
    // const double xl = floor( ipData.relaxedSolution[i] );
    assert(xl >= 1.);
    const double fl = f_x_value(data->I[i], xl); 
    const double fh = f_x_value(data->I[i], xl+1);
    const double w = fh - fl;
    weights[i] = w;
    
    if (0 && debugging)
      printf(" %f", w);
  }

  // ensure weights are unique
  weights.uniquify(1., 1.e4);
  
  if (debugging)
  {
    printf("\n");
    for (int i = 0; i < data->num_variables(); ++i)
    {
      const double x = ipData->relaxedSolution[i];
      const double xl = floor(x);
      printf("x %d (%f) in [%d,%d] weight %10.4f\n", i, x, (int) xl, (int) (xl + 1.), weights[i]);
    }
  }
}


// n = number of variables
// m = number of constraints (not counting variable bounds)


IARoundingNlp::~IARoundingNlp() {data = NULL;}

// returns the size of the problem
bool IARoundingNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  bool base_ok = baseNlp.get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  nnz_h_lag = 0; // nele_h
  return true && base_ok;
  // need to change this if there are more variables, such as delta-minuses
}

// returns the variable bounds
bool IARoundingNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  const bool base_ok = baseNlp.get_bounds_info(n, x_l, x_u, m, g_l, g_u);

  for (Index i=0; i<n; ++i) 
  {
    const double x = ipData->relaxedSolution[i];
    const double xl = floor(x);
    assert(xl >= 1.);
    x_l[i] = xl;
    x_u[i] = xl + 1.;
//    if (debugging)
//      printf("x %d (%f) in [%d,%d]\n", i, x, (int) x_l[i], (int) x_u[i]);
  }
  
  return true && base_ok; //means what?
}

// returns the initial point for the problem
bool IARoundingNlp::get_starting_point(Index n, bool init_x, Number* x_init,
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


// returns the value of the objective function
bool IARoundingNlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == data->num_variables());

  double obj = 0.;
  for (Index i = 0; i<n; ++i)
  {
    const double xl = floor( ipData->relaxedSolution[i] );
    obj += weights[i] * (x[i] - xl);
  }

  obj_value = obj;

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IARoundingNlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  //printf("E ");

  assert(n == data->num_variables());

  for (Index i = 0; i<n; ++i)
  {
    grad_f[i] = weights[i];
  }

  return true;
}

// return the value of the constraints: g(x)
bool IARoundingNlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  return baseNlp.eval_g(n, x, new_x, m, g);
}

// return the structure or values of the jacobian
bool IARoundingNlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  return baseNlp.eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
}

//return the structure or values of the hessian
bool IARoundingNlp::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
	// because the constraints are linear
  // and the objective function is linear
  // the hessian for this problem is actually empty

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

void IARoundingNlp::finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
{
  baseNlp.finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
}

} // namespace MeshKit
