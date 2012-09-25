// IAIntWaveNlp.cpp
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

#include "IAIntWaveNlp.hpp"
#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"
#include "IANlp.hpp"

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

using namespace Ipopt;

/* structure
 
 as base problem: objective function is f, same
 as base problem: sum-equal constraints
 
 add:
 constraints:
 cosine wave for each primal variable, forcing integrality
 cosine wave for each sum-even constraint, forcing integrality
 */

namespace MeshKit {
    
// constructor
IAIntWaveNlp::IAIntWaveNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr): 
data(data_ptr), ipData(ip_data_ptr), solution(solution_ptr), baseNlp(data_ptr, solution_ptr),
problem_n((int)data_ptr->I.size()), 
problem_m((int)(data_ptr->constraints.size() + 2*data_ptr->sumEvenConstraints.size() + data_ptr->num_variables())),
sum_even_constraint_start((int)(data_ptr->constraints.size() + data_ptr->sumEvenConstraints.size())),
x_int_constraint_start((int)(data_ptr->constraints.size() + 2*data_ptr->sumEvenConstraints.size())),
base_n((int)data_ptr->num_variables()),
base_m((int)(data_ptr->constraints.size() + data_ptr->sumEvenConstraints.size())),
PI( 2. * acos(0.0) ),
debugging(true), verbose(true) // true
{
  printf("\nIAIntWaveNLP Problem size:\n");
  printf("  number of variables: %lu\n", problem_n);
  printf("  number of constraints: %lu\n\n", problem_m);  
}

IAIntWaveNlp::~IAIntWaveNlp() {data = NULL; ipData = NULL;}

// returns the size of the problem
bool IAIntWaveNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  bool base_ok = baseNlp.get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);

  // constraints for sum-even
  int num_even_entries = 0;
  int num_even_h_entries = 0;
  for (std::vector<IAData::sumEvenConstraintRow>::const_iterator i=data->sumEvenConstraints.begin(); i != data->sumEvenConstraints.end(); ++i)
  {
  	num_even_entries += i->M.size();
  }
  nnz_jac_g += num_even_entries;

  // constraints for x-int
  nnz_jac_g += data->num_variables();
  
  // hessian elements for sum-even
  build_hessian();
  nnz_h_lag += hessian_vector.size();
  
  // hessian elements for x-int = main diagonal, already counted by objective function
  // nnz_h_lag += data->num_variables();
  
  return true && base_ok;
  // need to change this if there are more variables, such as delta-minuses
}

// returns the variable bounds
bool IAIntWaveNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  const bool base_ok = baseNlp.get_bounds_info(n, x_l, x_u, m, g_l, g_u);

  for (int i = 0; i < data->sumEvenConstraints.size(); ++i)
  {
    // cos( pi * sum_evens) == 1
    const int k = i+sum_even_constraint_start;
    g_l[k] = 1.; 
    g_u[k] = 1.; 
  }
  
  for (int i = 0; i < data->num_variables(); ++i)
  {
    // cos( 2 pi x) == 1
    const int k = i+x_int_constraint_start;
    g_l[k] = 1.;
    g_u[k] = 1.;
  }
  return true && base_ok; //means what?
}

// returns the initial point for the problem
bool IAIntWaveNlp::get_starting_point(Index n, bool init_x, Number* x_init,
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
bool IAIntWaveNlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  baseNlp.eval_f(base_n,x,new_x,obj_value);
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IAIntWaveNlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  return baseNlp.eval_grad_f(n, x, new_x, grad_f);
}

// return the value of the constraints: g(x)
bool IAIntWaveNlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  bool base_ok = baseNlp.eval_g(base_n, x, new_x, base_m, g);
  
  // cos( pi * sum_evens) == 1
  for (int i = 0; i < data->sumEvenConstraints.size(); ++i)
  {
    const int k = i+sum_even_constraint_start;
    double s = data->sumEvenConstraints[i].rhs;
    for (Index j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
    {
      s += x[ data->sumEvenConstraints[i].M[j].col ] * data->sumEvenConstraints[i].M[j].val;
    }
    const double gk = cos( PI * s );
    g[k] = gk;    
  }
  
  // cos( 2 pi x) == 1
  for (int i = 0; i < data->num_variables(); ++i)
  {
    const int k = i+x_int_constraint_start;
    const double gk = cos( 2. * PI * x[i] );
    g[k] = gk;
  }

  return true && base_ok;
}

// return the structure or values of the jacobian
bool IAIntWaveNlp::eval_jac_g(Index n, const Number* x, bool new_x,
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

    // g = cos( pi * sum_evens) == 1
    // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
    for (Index i = 0; i< data->sumEvenConstraints.size(); ++i)
    {
      for (Index j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
    		iRow[k] = i + (int) data->constraints.size();
        jCol[k] = data->sumEvenConstraints[i].M[j].col;
        ++k;
      }
    }
    // g = cos( 2 pi x) == 1
    // g' = - 2 pi sin ( 2 pi x ), for x_i
    for (int i=0; i<data->num_variables(); ++i)
    {
      iRow[k] = i;
      jCol[k] = i;
      ++k;
    }
  }
  else
  {
    // return the values of the jacobian of the constraints
    
    
    // g = cos( pi * sum_evens) == 1
    // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
    for (Index i = 0; i< data->sumEvenConstraints.size(); ++i)
    {
      double s = data->sumEvenConstraints[i].rhs;
      for (Index j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
        s += x[ data->sumEvenConstraints[i].M[j].col ] * data->sumEvenConstraints[i].M[j].val;
      }
      const double jac_gk = -PI * cos( PI * s );
      for (Index j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
        values[k++] = jac_gk;
      }
    }
    // g = cos( 2 pi x) == 1
    // g' = - 2 pi sin ( 2 pi x ), for x_i
    for (int i=0; i<data->num_variables(); ++i)
    {
      const double jac_gk = -2. * PI * sin( 2. * PI * x[i] );
      values[k++] = jac_gk;
    }
  }
  
  return true && base_ok;
}


void IAIntWaveNlp::add_hessian_entry( int i, int j, int &k )
{
  if ( j > i )
  {
    int i_temp = i;
    i = j;
    j = i_temp;
  }
  SparseMatrixEntry sme(i, j, k);
  if (hessian_map.insert( sme.key(), sme ))
  {
    hessian_vector.push_back( sme );
    ++k;
    assert( hessian_vector.size() == k );
    assert( hessian_map.size() == k );
  }
  // else it was already inserted, nothing to do
}

int IAIntWaveNlp::SparseMatrixEntry::n(0);

void IAIntWaveNlp::build_hessian()
{
  hessian_vector.clear();
  hessian_map.clear();
  SparseMatrixEntry::n = data->num_variables();
  int kk = 0;

  // objective function hessian - the main diagonal
  // important to add these first and in this order, for baseNlp to do the right thing
  for (int i = 0; i < data->num_variables(); ++i)
  {  
    add_hessian_entry( i, i, kk );
  }
  assert( kk == data->num_variables() );
  
  // sum_evens
    // g = cos( pi * sum_evens) == 1
    // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
    // g''= -pi^2 cos ( pi * sum_evens ), for each pair of variables contributing to the sum
    // assuming all the coefficients are 1    
  for (int c = 0; c< data->sumEvenConstraints.size(); ++c)
  {
    for (Index i = 0; i < data->sumEvenConstraints[i].M.size(); ++i)
    {
      int ii = data->sumEvenConstraints[c].M[i].col;
      for (Index j = 0; j <=i; ++j)
      {
        int jj = data->sumEvenConstraints[c].M[j].col;
        add_hessian_entry( ii, jj, kk );
      }
    }
  }
  
  // x integer
  // these are just the diagonals again, already added so skip
  // nele_hess = hessian_vector.size();
}

int IAIntWaveNlp::get_hessian_k( int i, int j ) const
{
  if ( i == j )
    return i;
  SparseMatrixEntry sme(i, j, -1 );
  int k = hessian_map.find( sme.key() ).first;
  return k;
}

//return the structure or values of the hessian
bool IAIntWaveNlp::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  baseNlp.eval_h(base_n, x, new_x, obj_factor, base_m, lambda, new_lambda, data->num_variables(), iRow, jCol, values);

  // hessian entry i,j is:
  // obj_factor fij + sum_k lambda[k] gkij
  // where fij =   d^2 f / d x_i d x_j
  //       gkij =  d^2 g[k] / d x_i d x_j
  // and d denotes partial derivative

  // first k entries are diagonal of objective function

  int k = data->num_variables(); // this is index where we start, for values and structure
    
  // This is a symmetric matrix, fill the lower left triangle only.
  if (values == NULL) {
    // return the structure. 
    for (int k = 0; k < hessian_vector.size(); ++k)
    {
      iRow[k] = hessian_vector[k].i;
      jCol[k] = hessian_vector[k].j;
    }
  }
  else {
    // return the values. 
    // g = cos( pi * sum_evens) == 1
    // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
    // g''= -pi^2 cos ( pi * sum_evens ), for each pair of variables contributing to the sum
    // assuming all the coefficients are 1
    for (Index i = 0; i< data->sumEvenConstraints.size(); ++i)
    {
      for (Index j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
        double contribution = 
        iRow[k] = i + sum_even_constraint_start;
        jCol[k] = data->sumEvenConstraints[i].M[j].col;
        k++;
      }
    }
    // g = cos( pi * sum_evens) == 1
    // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
    // g''= -pi^2 cos ( pi * sum_evens ), for each pair of variables contributing to the sum
    for (Index i = 0; i< data->sumEvenConstraints.size(); ++i)
    {
      double s = data->sumEvenConstraints[i].rhs;
      for (Index j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
        s += x[ data->sumEvenConstraints[i].M[j].col ] * data->sumEvenConstraints[i].M[j].val;
      }
    
    zzyk , take trig function of s
    zzyk add that entry times lambda into each of the hessian entries
      const double jac_gk = -PI * cos( PI * s );
      for (Index j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
        values[k++] = jac_gk;
      }
    }

    // g = cos( 2 pi x) == 1
    // g' = - 2 pi sin ( 2 pi x ), for x_i
    // g'' = - 4 pi^2 cos ( 2 pi x), for x_i only
    for (int i=0; i<data->num_variables(); ++i)
    {
      const double jac_gk = -2. * PI * sin( 2. * PI * x[i] );
      values[k++] = obj_factor * jac_gk;
    }
    
    // g = cos( 2 pi x) == 1
    // g' = - 2 pi sin ( 2 pi x ), for x_i
    // g'' = - 4 pi^2 cos ( 2 pi x), for x_i only
    for (int i=0; i<data->num_variables(); ++i)
    {
      iRow[k] = i + x_int_constraint_start;
      jCol[k] = i;
      k++;
    }
    

    ;
  }

  return true;
}

void IAIntWaveNlp::finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
{
  baseNlp.finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
  
  zzyk also report on how close the integer and sum-even constraints were satisfied!
}

} // namespace MeshKit
