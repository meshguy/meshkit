// IARoundingHeuristicMINLP.cpp
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

#include "meshkit/IANlp.hpp"
#include "meshkit/IARoundingHeuristicMINLP.hpp"
#include "meshkit/IAData.hpp"
#include "meshkit/IPData.hpp"
#include "meshkit/IASolution.hpp"

#include <math.h>

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

// constructor
IARoundingHeuristicMINLP::IARoundingHeuristicMINLP(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr): 
data(data_ptr), ip_data(ip_data_ptr), solution(solution_ptr),
debugging(true), verbose(true), baseNlp(data_ptr, solution_ptr) // true
{
  printf("\nIARoundingHeuristicMINLP Problem size:\n");
  printf("  number of variables: %lu\n", data->I.size());
  printf("  number of constraints: %lu\n\n", data->constraints.size());
}


// n = number of variables
// m = number of constraints (not counting variable bounds)


IARoundingHeuristicMINLP::~IARoundingHeuristicMINLP() {data = NULL;}

// returns the size of the problem
bool IARoundingHeuristicMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  return baseNlp.get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
}

// returns the variable bounds
bool IARoundingHeuristicMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{

  const bool base_ok = baseNlp.get_bounds_info(n, x_l, x_u, m, g_l, g_u);
  
  // add bounds from integer constraints
  if (ip_data->varIntegerBound.size())
  {
    for (Index i=0; i<n; ++i) 
    {
      const double b = ip_data->varIntegerBound[i];
      if (b != 0) // if 0 then it is unconstrained so far
      {
        // b is a lower bound if it is bigger than the goal
        if (b > data->I[i]) 
        {
          x_l[i] = b;
          x_u[i] = MESHKIT_IA_upperUnbound; 
        }
        // otherwise b is an upper bound
        else
        {
          x_l[i] = 1.0; 
          x_u[i] = b; 
        }
      }
    }
  }

  //printf("b ");
  return true && base_ok; //means what?
}

// returns the initial point for the problem
bool IARoundingHeuristicMINLP::get_starting_point(Index n, bool init_x, Number* x_init,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  const bool base_ok = baseNlp.get_starting_point(n, init_x, x_init, init_z, z_L, z_U, m, init_lambda, lambda);
  
  // initialize x to the prior solution
  if (ip_data->varIntegerBound.size())
    for (int i = 0; i<data->I.size(); ++i)
    {
      // todo: test if we need to modify this for x violating the variable bounds?
      x_init[i]= solution->x_solution[i]; 
    }
  
  return true && base_ok;
}


// returns the value of the objective function
bool IARoundingHeuristicMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  return baseNlp.eval_f(n, x, new_x, obj_value);
}

// return the gradient of the objective function grad_{x} f(x)
bool IARoundingHeuristicMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  return baseNlp.eval_grad_f(n, x, new_x, grad_f);
}

// return the value of the constraints: g(x)
bool IARoundingHeuristicMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  return baseNlp.eval_g(n, x, new_x, m, g);
}

// return the structure or values of the jacobian
bool IARoundingHeuristicMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  return baseNlp.eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
}

//return the structure or values of the hessian
bool IARoundingHeuristicMINLP::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  return baseNlp.eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nele_hess, iRow, jCol, values);
}

void IARoundingHeuristicMINLP::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
                                  const IpoptData* ip_data,
                                  IpoptCalculatedQuantities* ip_cq)
{
  // overwrites or fills in solution->x_solution
  baseNlp.finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
}
