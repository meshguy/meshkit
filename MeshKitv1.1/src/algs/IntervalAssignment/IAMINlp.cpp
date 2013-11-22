// IAMINlp.cpp
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

#include "IAMINlp.hpp"
#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"
#include "IANlp.hpp"

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

namespace MeshKit 
{
    
double p_norm = 3; // remove this, we use linear objectives here

// constructor
IAMINlp::IAMINlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr): 
data(data_ptr), ip_data(ip_data_ptr), solution(solution_ptr),
debugging(true), verbose(true) // true
{

  printf("\nIAMINlp Problem size:\n");
  printf("  number of variables: %lu\n", data->I.size());
  printf("  number of constraints: %lu\n\n", data->constraints.size());
}


// n = number of variables
// m = number of constraints (not counting variable bounds)

// apparent Meshkit style conventions:
// ClassFileName
// ClassName
// #ifndef MESHKIT_CLASSNAME_HPP
// #define MESHKIT_CLASSNAME_HPP
// namespace IAMeshKit
// public class_member_function()
// protected classMemberData

IAMINlp::~IAMINlp() {data = NULL;}

// returns the size of the problem
bool IAMINlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // printf("A ");
  // number of variables
  n = data->num_variables();

  // number of constraints
  m = (int) data->constraints.size() + (int) data->sumEvenConstraints.size();

  // number of non-zeroes in the Jacobian of the constraints
  size_t num_entries=0;
  for (std::vector<IAData::constraintRow>::const_iterator i=data->constraints.begin(); i != data->constraints.end(); ++i)
  {
  	num_entries += i->M.size();
  }
  for (std::vector<IAData::sumEvenConstraintRow>::const_iterator i=data->sumEvenConstraints.begin(); i != data->sumEvenConstraints.end(); ++i)
  {
  	num_entries += i->M.size();
  }
  nnz_jac_g = (int) num_entries;

  // number of non-zeroes in the Hessian 
  // diagonal entries for the objective function +
  // none for the constraints, since they are all linear
  nnz_h_lag = data->num_variables();

  // C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool IAMINlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
   //printf("B ");

  // The n and m we gave IPOPT in get_nlp_info are passed back to us.
  assert(n == data->num_variables());
  assert(m == (int)(data->constraints.size() + data->sumEvenConstraints.size()));

  // future interval upper and lower bounds:
  //   for midpoint subdivision, the lower bound may be 2 instead
  //   User may specify different bounds for some intervals
  //   Implement this by having another vector of lower bounds and one of upper bounds

  // relaxed problem
  if (ip_data->varIntegerBound.size() == 0)
  {
    // variables have lower bounds of 1 and no upper bounds
    for (Index i=0; i<n; ++i) 
    { 
      x_l[i] = 1.0; 
      x_u[i] = MESHKIT_IA_upperUnbound; 
    }
  }
  else
  {
    for (Index i=0; i<n; ++i) 
    {
      const double b = ip_data->varIntegerBound[i];
      if (b == 0) // unconstrained so far
      {
        x_l[i] = 1.0;
        x_u[i] = MESHKIT_IA_upperUnbound;
      }
      else
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
  
  // constraint bounds
  for (unsigned int i = 0; i<data->constraints.size(); ++i)
  {
    g_l[i] = data->constraints[i].lowerBound; 
    g_u[i] = data->constraints[i].upperBound; 
  }
  for (unsigned int i = 0; i<data->sumEvenConstraints.size(); ++i)
  {
    const int j = i + (int) data->constraints.size();
    g_l[j] = 4;
    g_u[j] = MESHKIT_IA_upperUnbound;
  }

  //printf("b ");
  return true; //means what?
}

// returns the initial point for the problem
bool IAMINlp::get_starting_point(Index n, bool init_x, Number* x_init,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  //printf("C ");

  // Minimal info is starting values for x, x_init
  // Improvement: we can provide starting values for the dual variables if we wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  assert(n==(int)data->I.size());
  assert(data->num_variables()==(int)data->I.size());
  
  // initialize x to the goals
  if (ip_data->varIntegerBound.size() == 0)
    for (unsigned int i = 0; i<data->I.size(); ++i)
    {
      x_init[i]=data->I[i];
    }
  // initialize x to the prior solution
  else
    for (unsigned int i = 0; i<data->I.size(); ++i)
    {
      // todo: test if we need to modify this for x violating the variable bounds?
      x_init[i]= solution->x_solution[i]; 
    }
  
  return true;
}


// r
inline Number IAMINlp::eval_r_i(const Number& I_i, const Number& x_i)
{
  //return x_i / I_i + I_i / x_i;
  return x_i > I_i ? 
    (x_i - I_i) / I_i :
    (I_i - x_i) / x_i;
}

inline Number IAMINlp::eval_grad_r_i(const Number& I_i, const Number& x_i)
{
  // return (1.0 / I_i) - I_i / (x_i * x_i);
  return x_i > I_i ? 
    1. / I_i :
    -I_i / (x_i * x_i);
}

inline Number IAMINlp::eval_hess_r_i(const Number& I_i, const Number& x_i)
{
  // return 2.0 * I_i / (x_i * x_i * x_i);
  return x_i > I_i ? 
    0 : 
    2 * I_i / (x_i * x_i * x_i);
}

inline Number IAMINlp::eval_R_i(const Number& I_i, const Number& x_i)
{
  if (x_i <= 0.)
    return MESHKIT_IA_upperUnbound;
  
  // old function, sum
  // return 2. + ( x_i / I_i ) * ( x_i / I_i ) + ( I_i / x_i ) * ( I_i / x_i );
  //  const register double r = (x_i - I_i) / I_i;
  //  const register double r2 = r*r;
  // return 2. + r2 + 1./r2;
  
  return pow( eval_r_i(I_i,x_i), p_norm );
}
  

inline Number IAMINlp::eval_grad_R_i(const Number& I_i, const Number& x_i)
{
  if (x_i <= 0.)
    return MESHKIT_IA_lowerUnbound;

  return p_norm * pow( eval_r_i(I_i,x_i), p_norm-1) * eval_grad_r_i(I_i, x_i);
}

inline Number IAMINlp::eval_hess_R_i(const Number& I_i, const Number& x_i)
{
  if (x_i <= 0.)
    return MESHKIT_IA_upperUnbound;
  
  // old function, sum
  // const register double I2 = I_i * I_i;
  // return 2. * ( 1. / I2 + 3. * I2 / (x_i * x_i * x_i * x_i) );

  if (p_norm == 2) 
    return (x_i > I_i) ? 2 / (I_i * I_i) : ( I_i / (x_i * x_i * x_i) ) * ( 6. * I_i / x_i - 4. ); 
  
  if (x_i > I_i ) 
    return p_norm * (p_norm - 1) * pow( (x_i - I_i) / I_i, p_norm -2) / (I_i * I_i); 
    
  // p (p-1) (I/x-1)^(p-2) I^2 / x^4 + p (I/x-1)^(p-1) 2 I / x^3  
  const register double f = (I_i - x_i) / x_i;
  const register double fp2 = pow( f, p_norm - 2);
  const register double x3 = x_i * x_i * x_i;
  return p_norm * ( (p_norm - 1) * fp2 * I_i * I_i / ( x3 * x_i ) + f * fp2 * 2. * I_i / x3 );
  
}

// s
inline Number IAMINlp::eval_s_i(const Number& I_i, const Number& x_i)
{
  return eval_r_i(x_i, I_i) * x_i;
}

inline Number IAMINlp::eval_grad_s_i(const Number& I_i, const Number& x_i)
{
  return eval_r_i(x_i, I_i) + x_i * eval_grad_r_i(I_i, x_i);
}

inline Number IAMINlp::eval_hess_s_i(const Number& I_i, const Number& x_i)
{
  return 2. * eval_grad_r_i(I_i, x_i) + x_i * eval_hess_r_i(I_i, x_i);
}

inline Number IAMINlp::eval_S_i(const Number& I_i, const Number& x_i)
{
  if (x_i <= 0.)
    return MESHKIT_IA_upperUnbound;
  const double s = eval_s_i(I_i,x_i);
  
  if (p_norm == 1)
    return s;
    
  if (p_norm == 2)
    return s*s;
  
  assert(p_norm > 2);
  return pow(s, p_norm );
}


inline Number IAMINlp::eval_grad_S_i(const Number& I_i, const Number& x_i)
{
  if (x_i <= 0.)
    return MESHKIT_IA_lowerUnbound;
  
  const double s = eval_s_i(I_i,x_i);
  const double sp = eval_grad_s_i(I_i, x_i);

  if (p_norm == 1)
    return s;
  
  assert(p_norm>1);
  if (p_norm == 2)
    return 2. * s * sp; 

  assert(p_norm > 2);
  return p_norm * pow( s, p_norm-1) * sp;
}

inline Number IAMINlp::eval_hess_S_i(const Number& I_i, const Number& x_i)
{
  if (x_i <= 0.)
    return MESHKIT_IA_upperUnbound;

  const double spp = eval_hess_s_i(I_i, x_i);
  
  if (p_norm == 1)
    return spp;
  
  const double s = eval_s_i(I_i, x_i);
  const double sp = eval_grad_s_i(I_i, x_i);

  if (p_norm==2)
    return 2. * ( sp * sp + s * spp );
  
  assert(p_norm>2);
  return p_norm * pow(s, p_norm - 2) * ( (p_norm - 1) * sp * sp + s * spp );
}


// experimental results: 
// best so far is f=R and p_norm = 3
// f=s doesn't work well for x < I, as it tends to push the x to 1 so the weight is small !

// returns the value of the objective function
bool IAMINlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  //printf("D ");

  assert(n == data->num_variables());

  // future: perhaps try max of the f, infinity norm
  
  double obj = 0.;
  for (Index i = 0; i<n; ++i)
  {
    obj += eval_R_i( data->I[i], x[i] );
  }

  obj_value = obj;

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IAMINlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  //printf("E ");

  assert(n == data->num_variables());

  for (Index i = 0; i<n; ++i)
  {
    grad_f[i] = eval_grad_R_i( data->I[i], x[i] ); 
  }

  return true;
}

// return the value of the constraints: g(x)
bool IAMINlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  //printf("F ");

  assert(n == data->num_variables());
  assert(m == (int)(data->constraints.size() + data->sumEvenConstraints.size()));
  
  
  for (Index i = 0; i<(int)data->constraints.size(); ++i)
  {
    //double g_i = constraints[i].rhs;
    double g_i = 0.;
    for (Index j = 0; j < (int)data->constraints[i].M.size(); ++j)
    {
      g_i += x[ data->constraints[i].M[j].col ] * data->constraints[i].M[j].val;
    }
    g[i] = g_i;
  }
  for (Index i = 0; i<(int)data->sumEvenConstraints.size(); ++i)
  {
    //double g_i = constraints[i].rhs;
    const int k = (int) data->constraints.size() + i;
    double g_k = 0.;
    for (Index j = 0; j < (int)data->sumEvenConstraints[i].M.size(); ++j)
    {
      g_k += x[ data->sumEvenConstraints[i].M[j].col ] * data->sumEvenConstraints[i].M[j].val;
    }
    g[k] = g_k;
  }

  return true;
}

// return the structure or values of the jacobian
bool IAMINlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  //printf("G ");

  assert(m == (int)(data->constraints.size() + data->sumEvenConstraints.size()));
  
  if (values == NULL) {
    // return the structure of the jacobian
    Index k=0;
    for (Index i = 0; i<(int)data->constraints.size(); ++i)
    {
      for (Index j = 0; j < (int)data->constraints[i].M.size(); ++j)
      {
    		iRow[k] = i;
        jCol[k] = data->constraints[i].M[j].col;
        ++k;
      }
    }
    for (Index i = 0; i< (int)data->sumEvenConstraints.size(); ++i)
    {
      for (Index j = 0; j < (int)data->sumEvenConstraints[i].M.size(); ++j)
      {
    		iRow[k] = i + (int) data->constraints.size();
        jCol[k] = data->sumEvenConstraints[i].M[j].col;
        ++k;
      }
    }

  }
  else {
    // return the values of the jacobian of the constraints
    Index k=0;
    for (Index i = 0; i < (int)data->constraints.size(); ++i)
    {
      for (Index j = 0; j < (int)data->constraints[i].M.size(); ++j)
      {
        values[k++] = data->constraints[i].M[j].val;
      }
    }
    for (Index i = 0; i < (int)data->sumEvenConstraints.size(); ++i)
    {
      for (Index j = 0; j < (int)data->sumEvenConstraints[i].M.size(); ++j)
      {
        values[k++] = data->sumEvenConstraints[i].M[j].val;
      }
    }
  }

  return true;
}

//return the structure or values of the hessian
bool IAMINlp::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  //printf("H ");

  // get_nlp_info specified the number of non-zeroes, should be
  assert(nele_hess == data->num_variables()); 
  
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

	// because the constraints are linear
    // the hessian for this constraints are actually empty
    
    // Since the objective function is separable, only the diagonal is non-zero
    for (Index idx=0; idx<data->num_variables(); ++idx)
    {
      iRow[idx] = idx;
      jCol[idx] = idx;
    }

  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

	// because the constraints are linear
    // the hessian for this problem is actually empty

    // fill the objective portion
    for (Index idx=0; idx<data->num_variables(); ++idx)
    {
      values[idx] = obj_factor * eval_hess_R_i(data->I[idx], x[idx]);
    }
  }

  return true;
}

void IAMINlp::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.
  //printf("I ");

  // copy solution into local storage x_solution
  solution->x_solution.clear(); // clear contents
  std::vector<double>(solution->x_solution).swap(solution->x_solution); // zero capacity
  solution->x_solution.reserve(n); // space for new solution
  for (Index i=0; i<n; i++) 
  {
    solution->x_solution.push_back( x[i] );  // values of new solution
  }
  assert( (int)solution->x_solution.size() == n );
  solution->obj_value = obj_value;
  assert(obj_value >= 0.);
  
  if (debugging)
  {
    printf("NLP solution:\n");       
    printf("x[%d] = %e\n", 0, x[0]);
    double rhs = 0., lhs = 0.;
    if (verbose)
    {
      printf("legend: coeff x_i (solution, goal, ratio; f(x), f'(x); F(x), F'(x) )\n");
      for (unsigned int j = 0; j < data->constraints.size(); ++j)
      {
        printf("constraint %d: ", j);
        const IAData::constraintRow & c = data->constraints[j];
        for (std::vector<IAData::sparseEntry>::const_iterator i = c.M.begin(); i < c.M.end(); ++i)
        {
          const double xv = x[i->col];
          const double gv = data->I[i->col];
          const double r = xv > gv ? (xv-gv) / gv : (gv - xv) / xv;
          printf(" %1.0f x_%d (%1.3f, %1.1f, %1.1f; %2.2g, %2.2g; %2.2g, %2.2g) ", 
                 i->val, i->col, xv, gv, 
                 r, eval_r_i(gv, xv), eval_grad_r_i(gv, xv), 
                 eval_R_i(gv,xv), eval_grad_R_i(gv,xv) );
          if (i->val > 0) 
            lhs += i->val * x[i->col];
          else
            rhs += i->val * x[i->col];
        }
        if (data->constraints.front().upperBound == data->constraints.front().lowerBound)
          printf(" = %1.1f", data->constraints.front().upperBound);  
        else
          printf(" in [%1.1f,%1.1f]", data->constraints.front().upperBound, data->constraints.front().lowerBound );
        printf(" <=> %f %f in ...\n", lhs, rhs );
      }
    
      printf("\n\nSolution of the primal variables, x\n");
      for (Index i=0; i<n; i++) {
        printf("x[%d] = %e\n", i, x[i]);
      }
      
      printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
      for (Index i=0; i<n; i++) {
        printf("z_L[%d] = %e\n", i, z_L[i]);
      }
      for (Index i=0; i<n; i++) {
        printf("z_U[%d] = %e\n", i, z_U[i]);
      }
      printf("\nFinal value of the constraints:\n");
      for (Index i=0; i<m ;i++) {
        printf("g(%d) = %e\n", i, g[i]);
      }
    }
    
    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj_value);
    
  }
}

} // namespace MeshKit
