// IABendNlp.cpp
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

#include "IABendNlp.hpp"
#include "IAData.hpp"
#include "IPData.hpp"
#include "IPBend.hpp"
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
  

// the objective function we should use for weights
double IABendNlp::f_x_value( double I_i, double x_i )
{
  return x_i > I_i ? 
  (x_i - I_i) / I_i :
  /* 1.133402234045 * */(I_i - x_i) / x_i;
  // I don't think this extra constant factor is necessary to break ties.
  // could just call baseNlp::eval_r_i

  // expected is 1/100 to 4 or so

  //??? is the following true?
  // the range of magnitudes of these weights is too high, the weights are too non-linear if we take them to a power
  // was
  // const double fh = IANlp::eval_R_i(data->I[i], xl+1);
}

// constructor
IABendNlp::IABendNlp(const IAData *data_ptr, const IPData *ip_data_ptr, const IPBendData *ip_bend_data_ptr,
                     IASolution *solution_ptr, IAWeights *weight_ptr, const bool set_silent): 
                     data(data_ptr), ipData(ip_data_ptr), ipBendData(ip_bend_data_ptr), solution(solution_ptr), 
                     baseNlp(data_ptr, solution_ptr, set_silent),
                     base_n((int)data_ptr->num_variables()),
                     base_m((int)(data_ptr->constraints.size() + data_ptr->sumEvenConstraints.size())),
                     problem_n(data_ptr->num_variables() + (int) weight_ptr->size()),
                     problem_m( (int)(data_ptr->constraints.size() + data_ptr->num_variables() + ( evenConstraintsActive ? 2 : 1) * data_ptr->sumEvenConstraints.size())),  
                     weights( weight_ptr ),
                     silent(set_silent), debugging(true), verbose(true) // true
{
  if (!silent)
  { 
    printf("\nIABendNlp Problem size:\n");
    printf("  number of variables: %lu\n", data->I.size());
    printf("  number of equality constraints: %lu\n", data->constraints.size());
    printf("  number of even constraints: %lu\n\n", data->sumEvenConstraints.size());
    // report on bends?
  }
}


// n = number of variables
// m = number of constraints (not counting variable bounds)


IABendNlp::~IABendNlp() {data = NULL;}

// returns the size of the problem
bool IABendNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  if (debugging) printf("IABendNlp::get_nlp_info\n");
  
  problem_m = (int) (data->constraints.size() + data->num_variables() + ( evenConstraintsActive ? 2 : 1) * data->sumEvenConstraints.size());
  wave_even_constraint_start = base_m + data->num_variables(); // base + deltas

  bool base_ok = baseNlp.get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  assert(n == base_n);
  assert(m == base_m);
  
  // add n for delta variables
  n += (Index) weights->size();
  
  // add m for delta-computing constraints
  m += (Index) data->num_variables();
  
  // add m for even-constraints
  if ( evenConstraintsActive )
    m += (Index) data->sumEvenConstraints.size();
  
  // Jacobian entries

  // data assumption checks
  if (debugging)
  {
    int num_deltas = 0;
    for (IPBendVec::const_iterator i=ipBendData->bendVec.begin(); 
         i != ipBendData->bendVec.end(); ++i )
    {
      num_deltas += i->num_deltas();
      printf("num_deltas: %d = %d + %d \n", i->num_deltas(), i->numDeltaPlus, i->numDeltaMinus);
    }
    assert( num_deltas == (int) weights->size() );
  }
  
  // base_program adds equality and even>4 constraints
  
  // delta constraints, x = sum deltas
  nnz_jac_g += weights->size() + data->num_variables();
  
  // wave even constraints
  // const Index base_nnz_jac_g = nnz_jac_g;
  if (evenConstraintsActive)
  {
    int num_even_entries = 0;
    for (std::vector<IAData::sumEvenConstraintRow>::const_iterator i=data->sumEvenConstraints.begin(); i != data->sumEvenConstraints.end(); ++i)
    {
      num_even_entries += i->M.size();
    }
    nnz_jac_g += num_even_entries;
  }

  // hessian for even constraints
  // objective functions for deltas are all linearized, so second derivatives are zero
  nnz_h_lag = 0; // nele_h // overwrite the non-linear objective
  build_hessian();
  nnz_h_lag = (Index) hessian_vector.size();
  
  if (debugging)
  {
    printf("IABendNlp::get_nlp_info: n %d, m %d, nnz_jac_g %d, nnz_h_lag, %d\n",
           n, m, nnz_jac_g, nnz_h_lag);
  }
  
  return true && base_ok;
}

// returns the variable bounds
bool IABendNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  if (debugging) printf("IABendNlp::get_bounds_info\n");
  
  const bool base_ok = baseNlp.get_bounds_info(base_n, x_l, x_u, base_m, g_l, g_u);

  // delta variable bounds
  for (IPBendVec::const_iterator i=ipBendData->bendVec.begin(); 
       i != ipBendData->bendVec.end(); ++i )
  {
    for (int j = 0; j < i->numDeltaPlus; ++j)
    {
      int k = i->deltaIStart + j;
      assert(k<n);
      x_l[k] = 0;
      x_u[k] = 1.;      
      //    if (debugging)
      //      printf("x %d (%f) delta+_%d in [%d,%d]\n", *i, x, j, (int) x_l[k], (int) x_u[k]);
    }
    // last delta is unlimitted
    const int last_plus_k = i->deltaIStart + i->numDeltaPlus - 1;
    assert(last_plus_k < n);
    x_u[ last_plus_k ] = MESHKIT_IA_upperUnbound;

    for (int j = 0; j < i->numDeltaMinus; ++j)
    {
      int k = i->deltaIStart + i->numDeltaPlus + j;
      assert(k<n);
      x_l[k] = 0;
      x_u[k] = 1.;      
      //    if (debugging)
      //      printf("x %d (%f) delta+_%d in [%d,%d]\n", *i, x, j, (int) x_l[k], (int) x_u[k]);
    }
    // last delta is only limited by x > 1.
    const int last_minus_k = i->deltaIStart + i->numDeltaPlus + i->numDeltaMinus - 1;
    assert(last_minus_k < n);
    x_u[ last_minus_k ] = i->xl - i->numDeltaMinus;
    
    // or MESHKIT_IA_upperUnbound;
  }
  
  // delta constraint bounds : they are equality constraints
  for (Index i = 0; i < (int) ipBendData->bendVec.size(); ++i)
  {
    const int k = i + base_m;
    assert(k < m);
    g_l[k] = 0.;
    g_u[k] = 0.;
  }
  
  // sum-even wave constraint bounds
  if (evenConstraintsActive)
  {
    for (unsigned int i = 0; i < data->sumEvenConstraints.size(); ++i)
    {
      // parabola(s) == 1 
      // equality makes ipopt think it is overdetermined, so use inequality :) 
      const int k = i + wave_even_constraint_start;
      assert(k < m);
      g_l[k] = 0.9999; // 1 - 1.e-4 means the sum will be within 1.e-2 ( because 2 = 4/2 ) of an integer. to do: tweak this tolerance?
      g_u[k] = MESHKIT_IA_upperUnbound; // 1.
      // the actual value of g will be below 1, so this contrains it to be 1 just as well
    }
    assert(wave_even_constraint_start + (int) data->sumEvenConstraints.size() == m);
  }
  else
  {
    assert((int) ipBendData->bendVec.size() + base_m == m);
  }
  
  if (debugging)
  {
    for (int i = 0; i < n; ++i)
    {
      printf("x[%d] in [%g,%g]\n", i, x_l[i], x_u[i]);
    }
    for (int j = 0; j < m; ++j)
    {
      printf("g[%d] in [%g,%g]\n", j, g_l[j], g_u[j]);
    }
  }
  return true && base_ok;
}

// returns the initial point for the problem
bool IABendNlp::get_starting_point(Index n, bool init_x, Number* x_init,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  if (debugging) printf("IABendNlp::get_starting_point\n");

  // idea: use last solution, not necessarily relaxed one ?? would that be better or worse?

  // Minimal info is starting values for x, x_init
  // Improvement: we can provide starting values for the dual variables if we wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize x to the relaxed solution
  {
    Index i;
    for (i=0; i<data->num_variables(); ++i) 
    {
      assert(i<n);
      x_init[i] = ipData->relaxedSolution[i];
    }
  
    // alternative:
    // deltas are all zero
    // that point is actually not feasible because the delta constraints are not satisfied
    for ( ; i < n; ++i )
    {
      assert(i<n);
      x_init[i] = 0.;
    }
  }
  
  // delta - add non-integer part of x to the delta
  {
    assert( (int) ipBendData->bendVec.size() == data->num_variables() );
    for (int i = 0; i < (int) ipBendData->bendVec.size(); ++i )
    {
      double d = ipData->relaxedSolution[i] - ipBendData->bendVec[i].xl;
      // assert( fabs(d) < 1.01 );
      if ( d > 0. )
      {
        for (int di = 0; di < ipBendData->bendVec[i].numDeltaPlus; ++di)
        {
          const int k = ipBendData->bendVec[i].deltaIStart + di;
          assert(k<n);
          if (( d <= 1. ) || ( di + 1 == ipBendData->bendVec[i].numDeltaPlus))
          {
            x_init[k] = d;
            break;
          }
          else
          {
            x_init[k] = 1.;
            d -= 1.;
          }
        }
      }
      else
      {
        d = fabs(d);
        for (int di = 0; di < ipBendData->bendVec[i].numDeltaMinus; ++di)
        {
          const int k = ipBendData->bendVec[i].deltaIStart + ipBendData->bendVec[i].numDeltaPlus + di;
          assert(k<n);
          if (( d <= 1. ) || ( di + 1 == ipBendData->bendVec[i].numDeltaMinus))
          {
            x_init[k] = fabs(d);
            d = 0;
            break;
          }
          else
          {
            x_init[k] = 1.;
            d -= 1.;
          }
        }
        
      }
    }
  }
  return true;
}


// returns the value of the objective function
bool IABendNlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  if (debugging) printf("IABendNlp::eval_f\n");
  assert(n == data->num_variables() + (int) weights->size());

  // only the deltas contribute
  double obj = 0.;
  Index i(0);
  Index k(data->num_variables());
  for ( ; k<n; ++i, ++k)
  {
    obj += (*weights)[i] * x[k];
  }
  obj_value = obj;

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IABendNlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  const bool debug_loc = false;
  if (debugging && debug_loc) printf("IABendNlp::eval_grad_f\n");

  assert(n == data->num_variables() + (int) weights->size());

  // only the deltas contribute
  Index i;
  for (i = 0; i<data->num_variables(); ++i)
  {
    grad_f[i] = 0;
  }

  // deltas contribute their weight
  Index k(data->num_variables());
  for (i = 0; i < (int) weights->size(); ++i, ++k)
  {
    grad_f[k] = (*weights)[i];
  }

  return true;
}

// return the value of the constraints: g(x)
bool IABendNlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  const bool debug_loc = false;
  if (debugging && debug_loc) printf("IABendNlp::eval_g\n");
  bool base_ok = baseNlp.eval_g(base_n, x, new_x, base_m, g);

  // sum x - deltas
  //       x = xl + delta_plusses - delta_minuses
  // g = (     xl + delta_plusses - delta_minuses - x  ), constrained to be zero
  for (Index i = 0; i < (int) ipBendData->bendVec.size(); ++i)
  {
    const int k = i + base_m;
    const IPBend &bend = ipBendData->bendVec[i];
    double gk = bend.xl - x[i];
    if (debugging && debug_loc)
      printf("constraint %d: sum of deltas for x_%d: xl (%d) - x (%f) ", k, i, bend.xl, x[i]);
    for (int j = 0; j < bend.numDeltaPlus; ++j )
    {
      int di = bend.deltaIStart + j;
      gk += x[di];
      if (debugging && debug_loc)
        printf(" + %f", x[di]);
    }
    for (int j = 0; j < bend.numDeltaMinus; ++j )
    {
      int di = bend.deltaIStart + bend.numDeltaPlus + j;
      gk -= x[di];
      if (debugging && debug_loc)
        printf(" - %f", x[di]);
    }
    g[k] = gk;
    if (debugging && debug_loc)
      printf(" = %f\n", gk);
  }
  
  // sum evens, parabola value
  if (evenConstraintsActive)
  {
    for (unsigned int i = 0; i < data->sumEvenConstraints.size(); ++i)
    {
      const int k = i + wave_even_constraint_start;
      double s = baseNlp.eval_even_sum(i,x);
      const double gk = eval_g_int_s(s); // e.g. (s - nearest_even_number)^2;
      g[k] = gk;    
      if (debugging)
      {
        printf("constraint %d: wave even %d, sum = %f, g = %f\n",
               k, i, s, gk );
      }
    }

  }
  
  return true && base_ok;
}

// return the structure or values of the jacobian
bool IABendNlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  const bool debug_loc = false;
  if (debugging && debug_loc) printf("IABendNlp::eval_jac_g\n");
  
  int base_nele_jac = baseNlp.get_neleJac();
  bool base_ok = baseNlp.eval_jac_g(base_n, x, new_x, base_m, base_nele_jac, iRow, jCol, values);
  base_nele_jac = baseNlp.get_neleJac();
  int kk = base_nele_jac;
  
  // sum x - deltas
  //       x = xl + delta_plusses - delta_minuses
  // g = (     xl + delta_plusses - delta_minuses - x  ), constrained to be zero
  // g' =       0 + 1 dp_indexes  - 1 dm_indexes - 1 x_i
  

  if (values == NULL) 
  {
    if (debugging && debug_loc)
    {
      printf("jacobian g structure:\n");
    }
    // deltas
    for (Index i = 0; i < (int) ipBendData->bendVec.size(); ++i)
    {
      const int k = i + base_m;
      const IPBend &bend = ipBendData->bendVec[i];
      if (debugging && debug_loc)
      {
        printf("constraint %d: jac for sum-deltas for x_%d, ", k, i);
      }
      iRow[kk] = k;
      jCol[kk] = i; // x[i]
      ++kk;
      for (int j = 0; j < bend.numDeltaPlus; ++j )
      {
        int di = bend.deltaIStart + j;
        iRow[kk] = k;
        jCol[kk] = di;
        if (debugging && debug_loc)
        {
          printf(" deltaplus %d (x[%d])", j, di );
        }
        ++kk;
      }
      for (int j = 0; j < bend.numDeltaMinus; ++j )
      {
        int di = bend.deltaIStart + bend.numDeltaPlus + j;
        iRow[kk] = k;
        jCol[kk] = di;
        if (debugging && debug_loc)
        {
          printf(" deltaminus %d (x[%d])", j, di );
        }
        ++kk;
      }
      if (debugging && debug_loc) printf("\n");
    }
    
    // sum evens
    if (evenConstraintsActive)
    {
      for (unsigned int i = 0; i< data->sumEvenConstraints.size(); ++i)
      {
        const int k = i + wave_even_constraint_start;
        if (debugging)
        {
          printf("wave even constraint %d: jac, for natural even constraint %d, ", k, i);
        }
        
        for (unsigned int j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
        {
          iRow[kk] = k;
          jCol[kk] = data->sumEvenConstraints[i].M[j].col;
          if (debugging)
          {
            printf(" kk %d (i %d, j %d)", kk, iRow[kk], jCol[kk]);
          }
          ++kk;
        }
      }
      if (debugging) printf("\n");
    }
    assert( kk == nele_jac );
    if (debugging && debug_loc)
    {
        printf("Jacobian g sparse structure:\n");
      for (int i = 0; i < nele_jac; ++i)
      {
        printf("entry %d: row %d col %d\n", i, iRow[i], jCol[i] );
              
      }
    }
  }
  else
  {
    if (debugging)
    {
      printf("jacobian g values:\n");
    }

    // deltas
    for (Index i = 0; i < (int) ipBendData->bendVec.size(); ++i)
    {
      if (debugging && debug_loc)
      {
        printf("constraint %d: jac for sum-deltas for x_%d, ", i + base_m, i);
      }

      const IPBend &bend = ipBendData->bendVec[i];
      values[kk++] = -1; // - x[i]
      if (debugging) printf("-1 (x_%d) ", i);

      for (int j = 0; j < bend.numDeltaPlus; ++j )
      {
        values[kk++] = 1; // + deltaplus_j
        if (debugging) printf("+1 (d+_%d) ", j);
      }
      for (int j = 0; j < bend.numDeltaMinus; ++j )
      {
        values[kk++] = -1; // - deltaplus_j
        if (debugging) printf("-1 (d-_%d) ", j);
      }
      if (debugging) printf("\n");
    }
    
    // sum evens
    if (evenConstraintsActive)
    {
      
      for (unsigned int i = 0; i< data->sumEvenConstraints.size(); ++i)
      {
        const int k = i + wave_even_constraint_start;
        if (debugging)
        {
          printf("wave even constraint %d: jac, for natural even constraint %d, ", k, i);
        }
        const double s = baseNlp.eval_even_sum(i, x);
        const double jac_gk = eval_jac_int_s(s); // e.g. -PI * cos( PI * s );
        for (unsigned int j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
        {
          const double coeff = data->sumEvenConstraints[i].M[j].val;
          values[kk++] = coeff * jac_gk;
          if (debugging)
          {
            printf("  %d: x_%d gradient %f * %f = %f\n", kk-1,
                   data->sumEvenConstraints[i].M[j].col, coeff, jac_gk, coeff * jac_gk);
          }
        }
        if (debugging) printf("\n");
      }
    }
    
    assert( kk == nele_jac );
    if (debugging && debug_loc)
    {
      printf("Jacobian g sparse values:\n");
      for (int i = 0; i < nele_jac; ++i)
      {
        printf("entry %d: %f\n", i, values[i]);        
      }
    }

  }
  
  return true && base_ok;
}

  
  
IABendNlp::SparseMatrixEntry::SparseMatrixEntry(const int iset, const int jset, const int kset)
{
  if ( jset > iset )
  {
    i = jset;
    j = iset;
  }
  else
  {
    i = iset;
    j = jset;
  }
  k = kset;
  assert(j <= i);
}
  
void IABendNlp::add_hessian_entry( int i, int j, int &k )
{
  SparseMatrixEntry sme(i, j, k);
  
  if (hessian_map.insert( std::make_pair(sme.key(), sme) ).second)
  {
    hessian_vector.push_back( sme );
    ++k;
    assert( (int) hessian_vector.size() == k );
    assert( (int) hessian_map.size() == k );
  }
  // else it was already inserted, nothing to do
}
  
int IABendNlp::SparseMatrixEntry::n(0);
  
void IABendNlp::build_hessian() 
{
  if (debugging) printf("IABendNlp::build_hessian\n");

  // only build once
  if (hessian_vector.size())
    return;
  
  hessian_vector.clear();
  hessian_map.clear();
  SparseMatrixEntry::n = data->num_variables() + (int) weights->size(); // variables + deltas
  int kk = 0;
  
  if (!evenConstraintsActive)
  {
    if (debugging)
      printf("even constraints not active => empty hessian\n");
    return;
  }
  
  // sum_evens
  // g = cos( pi * sum_evens) == 1
  // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
  // g''= -pi^2 cos ( pi * sum_evens ), for each pair of variables contributing to the sum
  // assuming all the coefficients are 1    
  for (unsigned int c = 0; c< data->sumEvenConstraints.size(); ++c)
  {
    for (unsigned int i = 0; i < data->sumEvenConstraints[c].M.size(); ++i)
    {
      int ii = data->sumEvenConstraints[c].M[i].col;
      for (unsigned int j = 0; j <=i; ++j)
      {
        int jj = data->sumEvenConstraints[c].M[j].col;
        if (debugging)
          printf("hessian adding %d %d %d\n", ii, jj, kk);
        add_hessian_entry( ii, jj, kk );
      }
    }
  }
  
  // nele_hess = hessian_vector.size();
  if (debugging)
  {
    printf("==========\nBuilt Hessian\n");
    print_hessian();
    printf("==========\n");
  }
}
  
int IABendNlp::get_hessian_k( int i, int j )
{
  if ( i == j )
    return i;
  SparseMatrixEntry sme(i, j, -1 );
  SparseMatrixEntry &entry = hessian_map[ sme.key() ]; 
  return entry.k;
}
  
void IABendNlp::print_hessian()
{
  printf("Packed Hessian:\n");
  for (unsigned int kk = 0; kk < hessian_vector.size(); ++kk)
  {
    SparseMatrixEntry &sme = hessian_vector[kk];
    printf(" %d: (%d, %d)\n", sme.k, sme.i, sme.j);
  }
  printf("\n");          
  
  printf("Random Access Hessian in sequence\n");
  {
    int k = 0;
    SparseMatrixMap::iterator iter;
    for (iter = hessian_map.begin(); iter != hessian_map.end(); ++iter, ++k)
    {
      const SparseMatrixEntry &sme = iter->second;
      printf(" %d: (%d, %d) k:%d key:%d\n", k, sme.i, sme.j, sme.k, sme.key() );
    }
    printf("\n");          
  }
  
  printf("Random Access Hessian in square:\n");
  for (int i = 0; i < data->num_variables(); ++i)
  {
    for (int j = 0; j < data->num_variables(); ++j)
    {
      SparseMatrixEntry sme_search(i, j, -1 );
      SparseMatrixMap::iterator iter = hessian_map.find(sme_search.key());
      if (iter == hessian_map.end())
        printf("   . ");
      else
        printf(" %3d ", iter->second.k);
    }
    printf("\n");          
  }
  printf("\n");          
  
}
  
//return the structure or values of the hessian
bool IABendNlp::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  if (debugging) printf("IABendNlp::eval_h\n");

  // fill values with zeros
  if (values)
  {
    for (unsigned int kk = 0; kk < hessian_vector.size(); ++kk)
    {
      values[kk] = 0.;
    }
    // debug, print x
  }
  
  // base hessian is empty, as all those constraints are linear and the obj. function is linearized 
  
  // hessian entry i,j is:
  // obj_factor fij + sum_k lambda[k] gkij
  // where fij =   d^2 f / d x_i d x_j
  //       gkij =  d^2 g[k] / d x_i d x_j
  // and d denotes partial derivative
  
  // This is a symmetric matrix, fill the lower left triangle only.
  if (values == NULL) {
    assert((int) hessian_vector.size() == nele_hess);
    // return the structure. 
    for (unsigned int kk = 0; kk < hessian_vector.size(); ++kk)
    {
      iRow[kk] = hessian_vector[kk].i;
      jCol[kk] = hessian_vector[kk].j;
    } 
  } // structure
  else {
    // return the values. 
    // g =    ( sum_evens - nearest_sum_even )^2 == 0
    // g' = 2 ( sum_evens - nearest_sum_even )  for each of the variables contributing to the sum
    // g''= 2                                   for each pair of variables contributing to the sum
    // assuming all the coefficients are 1
    if (evenConstraintsActive)
    {
      if (debugging)
      {
        printf("\nwave even non-zero hessian values:");
      }
      for (unsigned int i = 0; i< data->sumEvenConstraints.size(); ++i)
      {
        if (debugging)
        {
          printf("\n%d constraint: ", i);
        }
        const int k = i + wave_even_constraint_start; // index of the constraint in the problem, = lambda to use
        
        const double s = baseNlp.eval_even_sum(i, x); // sum of the variable values
        // second derivative of wave function, assuming coefficients of one
        const double wavepp = eval_hess_int_s(s); // e.g. ( -PI * PI * cos( PI * s ) ); 
        const double h_value = lambda[k] * wavepp;        
        // assign that value to all pairs, weighted by coeff_i * coeff_j
        for (unsigned int ii = 0; ii < data->sumEvenConstraints[i].M.size(); ++ii)
        {
          const int var_i_index = data->sumEvenConstraints[i].M[ii].col;
          const double coeff_i =  data->sumEvenConstraints[i].M[ii].val;
          for (unsigned int jj = 0; jj < data->sumEvenConstraints[i].M.size(); ++jj)
          {
            const int var_j_index = data->sumEvenConstraints[i].M[jj].col;
            const double coeff_j =  data->sumEvenConstraints[i].M[jj].val;
            
            const int kk = get_hessian_k(var_i_index, var_j_index);
            assert( kk >= 0 && kk < nele_hess );
            values[kk] += coeff_i * coeff_j * h_value;
            if (debugging)
            {
              printf("  lambda[%d] * coeff_%d * coeff_%d * d^2 wave / d x_%d d x_%d = %f * %f * %f * %f\n", 
                     k, var_i_index, var_j_index, var_i_index, var_j_index, 
                     lambda[k], coeff_i, coeff_j, wavepp );
            }
          }
        }
      }
      if (debugging)
      {
        printf("\n");
      }
      assert((int) hessian_vector.size() == nele_hess);
    }

  } // values
  
  return true;
}

void IABendNlp::finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
{
  if (debugging) printf("IABendNlp::finalize_solution\n");
  // by passing n and not base_n, m and not base_m,
  // we all get the deltas and sum-even values
  baseNlp.finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
}

} // namespace MeshKit
