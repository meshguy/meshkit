// IAIntWaveNlp.cpp
// Interval Assignment for Meshkit
//
#include "IAIntWaveNlp.hpp"

#include <math.h>
#include <limits>

#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"

/* structure
 
 constraints:
 cosine wave for each sum-even constraint, forcing integrality
 cosine wave for each primal variable, forcing integrality
 */

namespace MeshKit {
  
// constructor
IAIntWaveNlp::IAIntWaveNlp(const IAData *data_ptr, const IPData *ip_data_ptr, IASolution *solution_ptr, bool set_silent): 
  data(data_ptr), ipData(ip_data_ptr), solution(solution_ptr), baseNlp(data_ptr, solution_ptr),
  problem_n((int)data_ptr->I.size()), 
  problem_m((int)(data_ptr->constraints.size() + 2*data_ptr->sumEvenConstraints.size() + data_ptr->num_variables())),
  wave_even_constraint_start((int)(data_ptr->constraints.size() + data_ptr->sumEvenConstraints.size())),
  wave_int_constraint_start((int)(data_ptr->constraints.size() + 2*data_ptr->sumEvenConstraints.size())),
  base_n((int)data_ptr->num_variables()),
  base_m((int)(data_ptr->constraints.size() + data_ptr->sumEvenConstraints.size())),
  silent(set_silent), debugging(true), verbose(true) // true
{
  printf("\nIAIntWaveNLP Problem size:\n");
  printf("  number of variables: %d\n", problem_n);
  printf("  number of base (equal, even>4) constraints: %d\n", base_m);  
  printf("  number of wave-even constraints: %lu\n", data_ptr->sumEvenConstraints.size());
  printf("  number of wave-int constraints: %d\n", data_ptr->num_variables());
  printf("  total constraints: %d\n\n", problem_m);  
}


IAIntWaveNlp::~IAIntWaveNlp() {data = NULL; ipData = NULL;}

// returns the size of the problem
bool IAIntWaveNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                           Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  bool base_ok = baseNlp.get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  
  // n = number of variables = same as in base problem
  // m = number of constraints = equality constraints(side_1=side_2) + even_constraints(sum>4) set by base, plus
  //   wave even constraints
  //   wave integrality constraints
  m += (Index) data->sumEvenConstraints.size() + (Index) data->num_variables();
  assert( m == problem_m );

  if (debugging)
  {
    printf("IAIntWaveNlp::get_nlp_info\n");
    printf("base m=%d, wave m=%d\n", base_m, m);
  }
  // nnz_jac_g = number of non-zero entries in the jacobian of the constraints
  // equality constraints counted in the base program
  // constraints for sum-even
  
  // wave even constraints
  const Index base_nnz_jac_g = nnz_jac_g;
  int num_even_entries = 0;
  for (std::vector<IAData::sumEvenConstraintRow>::const_iterator i=data->sumEvenConstraints.begin(); i != data->sumEvenConstraints.end(); ++i)
  {
    num_even_entries += i->M.size();
  }
  nnz_jac_g += num_even_entries;
  
  // wave x-integer constraints
  nnz_jac_g += data->num_variables();
  
  if (debugging)
  {
    printf("nnz_jac_g = %d: base = %d, wave even = %d, wave int = %d\n", 
           nnz_jac_g, base_nnz_jac_g, num_even_entries, data->num_variables());
  }
  
  // hessian elements, second derivatives of objective and constraints
  // objectives are double counted, so we do = here rather than +=
  build_hessian();
  nnz_h_lag = (Index) hessian_vector.size();
  
  if (debugging)
  {
    printf("IAIntWaveNlp::get_nlp_info");
    printf(" n=%d, m=%d, nnz_jac_g=%d, num_even_entrites=%d, nnz_h_lag=%d\n", n, m, nnz_jac_g, num_even_entries, nnz_h_lag);
  }
  
  return true && base_ok;
  // need to change this if there are more variables, such as delta-minuses
}
  
// returns the variable bounds
bool IAIntWaveNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                   Index m, Number* g_l, Number* g_u)
{
  const bool base_ok = baseNlp.get_bounds_info(n, x_l, x_u, base_m, g_l, g_u);
  
  for (unsigned int i = 0; i < data->sumEvenConstraints.size(); ++i)
  {
    // cos( pi * sum_evens) == 1
    // equality makes ipopt think it is overdetermined, so use inequality :) 
    const int k = i + wave_even_constraint_start;
    g_l[k] = 1.; // 1.
    g_u[k] = MESHKIT_IA_upperUnbound; // 1.
  }
  
  for (int i = 0; i < data->num_variables(); ++i)
  {
    // cos( 2 pi x) == 1
    const int k = i+ wave_int_constraint_start;
    g_l[k] = 1.; // 1.
    g_u[k] = MESHKIT_IA_upperUnbound; // 1.
  }
  return true && base_ok;
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
  if (debugging)
    printf("IAIntWaveNlp::eval_g\n");
  
  bool base_ok = baseNlp.eval_g(base_n, x, new_x, base_m, g);
  
  // cos( pi * sum_evens) == 1
  for (unsigned int i = 0; i < data->sumEvenConstraints.size(); ++i)
  {
    const int k = i + wave_even_constraint_start;
    double s = baseNlp.eval_even_sum(i,x);
    const double gk = eval_g_int_s(s); // e.g. cos( PI * s );
    g[k] = gk;    
    if (debugging)
    {
      printf("IAIntWaveNlp::eval_g wave even constraint %d(%d), sum = %f, g = %f\n",
             i, k, s, gk );
    }
  }
  
  // cos( 2 pi x) == 1
  for (int i = 0; i < data->num_variables(); ++i)
  {
    const int k = i + wave_int_constraint_start;
    const double gk = eval_g_int_x(x[i]); // e.g. cos( 2. * PI * x[i] );
    g[k] = gk;
    if (debugging)
    {
      printf("IAIntWaveNlp::eval_g wave int constraint %d(%d), x_%d = %f, g = %f\n",
             i, k, i, x[i], gk );
    }
  }

  if (debugging)
    printf("IAIntWaveNlp::eval_g done.\n");

  return true && base_ok;
}

// return the structure or values of the jacobian
bool IAIntWaveNlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  if (debugging)
  {
    printf("IAIntWaveNlp::eval_jac_g");
    if (values)
      printf(" values\n");
    else 
      printf(" structure\n");
  }

  int base_nele_jac = baseNlp.get_neleJac(); // might be 0 if not computed yet, values == NULL
  bool base_ok = baseNlp.eval_jac_g(base_n, x, new_x, base_m, base_nele_jac, iRow, jCol, values);
  base_nele_jac = baseNlp.get_neleJac(); // set to true value
  int k = base_nele_jac;

  printf("base_nele_jac = %d, nele_jac = %d\n", base_nele_jac, nele_jac); // should be +1 for each equal and even constraint entry
  
  if (values == NULL) 
  {
    // return the structure of the jacobian
  
    // g = cos( pi * sum_evens) == 1
    // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
    if (debugging)
    {
      printf("base entries: ");
      for (int bi = 0; bi < k; ++bi)
      {
        printf(" %d (%d,%d)", bi, iRow[bi], jCol[bi]);        
      }
      printf("\nwave even non-zero entries: ");
    }
    for (unsigned int i = 0; i< data->sumEvenConstraints.size(); ++i)
    {
      for (unsigned int j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
    		iRow[k] = i + wave_even_constraint_start;
        jCol[k] = data->sumEvenConstraints[i].M[j].col;
        if (debugging)
        {
          printf(" %d (%d,%d)", k, iRow[k], jCol[k]);
        }
        ++k;
      }
    }

    // g = cos( 2 pi x) == 1
    // g' = - 2 pi sin ( 2 pi x ), for x_i
    if (debugging)
    {
      printf("\nwave int non-zero entries: ");
    }
    for (int i=0; i<data->num_variables(); ++i)
    {
      iRow[k] = i + wave_int_constraint_start;
      jCol[k] = i;
      if (debugging)
      {
        printf(" %d (%d,%d)", k, iRow[k], jCol[k]);
      }
      ++k;
    }
    if (debugging)
    {
      printf("\n");
      printf("k = %d, nele_jac = %d\n", k, nele_jac);
    }
    assert(k == nele_jac);
  }
  else
  {
    // return the values of the jacobian of the constraints    
    
    // g = cos( pi * sum_evens) == 1
    // g' = -pi * coeff_i * sin ( pi * sum_evens ), for each of the variables x_i contributing to the sum
    if (debugging)
    {
      printf("base values: ");
      for (int bi = 0; bi < k; ++bi)
      {
        printf(" %d (%f)", bi, values[bi]);
      }
      printf("\nwave even non-zero jacobian values: ");
    }
    for (unsigned int i = 0; i< data->sumEvenConstraints.size(); ++i)
    {
      const double s = baseNlp.eval_equal_sum(i, x);
      const double jac_gk = eval_jac_int_s(s); // e.g. -PI * cos( PI * s );
      if (debugging)
        printf("\n%d even wave: ", i);
      for (unsigned int j = 0; j < data->sumEvenConstraints[i].M.size(); ++j)
      {
        const double coeff = data->sumEvenConstraints[i].M[j].val;
        values[k++] = coeff * jac_gk;
        if (debugging)
        {
          printf("  %d: x_%d gradient %f * %f = %f\n", k-1,
                 data->sumEvenConstraints[i].M[j].col, coeff, jac_gk, coeff * jac_gk);
        }
      }
    }
    if (debugging)
      printf("\n");

    // g = cos( 2 pi x) == 1
    // g' = - 2 pi sin ( 2 pi x ), for x_i
    if (debugging)
    {
      printf("\nwave int non-zero jacobian values: ");
    }
    for (int i=0; i<data->num_variables(); ++i)
    {
      const double jac_gk = eval_jac_int_x(x[i]); // e.g. -2. * PI * sin( 2. * PI * x[i] );
      values[k++] = jac_gk;
      if (debugging)
      {
        printf("\n%d: x_%d (%f) gradient %f", k-1, i, x[i], jac_gk);
      }
    }
    if (debugging)
      printf("\n");
    assert(k == nele_jac);
  }
  
  return true && base_ok;
}

IAIntWaveNlp::SparseMatrixEntry::SparseMatrixEntry(const int iset, const int jset, const int kset)
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
  
void IAIntWaveNlp::add_hessian_entry( int i, int j, int &k )
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

int IAIntWaveNlp::SparseMatrixEntry::n(0);

void IAIntWaveNlp::build_hessian() 
{
  // only build once
  if (hessian_vector.size())
    return;
  
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
  for (unsigned int c = 0; c< data->sumEvenConstraints.size(); ++c)
  {
    for (unsigned int i = 0; i < data->sumEvenConstraints[i].M.size(); ++i)
    {
      int ii = data->sumEvenConstraints[c].M[i].col;
      for (unsigned int j = 0; j <=i; ++j)
      {
        int jj = data->sumEvenConstraints[c].M[j].col;
        add_hessian_entry( ii, jj, kk );
      }
    }
  }
  
  // x integer
  // these are just the diagonals again, already added so skip
  // nele_hess = hessian_vector.size();
  if (debugging)
  {
    printf("==========\nBuilt Hessian\n");
    print_hessian();
    printf("==========\n");
  }
}

int IAIntWaveNlp::get_hessian_k( int i, int j )
{
  if ( i == j )
    return i;
  SparseMatrixEntry sme(i, j, -1 );
  SparseMatrixEntry &entry = hessian_map[ sme.key() ]; 
  return entry.k;
}

void IAIntWaveNlp::print_hessian()
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
bool IAIntWaveNlp::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  // fill values with zeros
  if (values)
  {
    for (unsigned int kk = 0; kk < hessian_vector.size(); ++kk)
    {
      values[kk] = 0.;
    }
    // debug, print x
  }
    
  // get structure, or values, from objective function, which is just the diagonal entries
  baseNlp.eval_h(base_n, x, new_x, obj_factor, base_m, lambda, new_lambda, data->num_variables(), iRow, jCol, values);

  // hessian entry i,j is:
  // obj_factor fij + sum_k lambda[k] gkij
  // where fij =   d^2 f / d x_i d x_j
  //       gkij =  d^2 g[k] / d x_i d x_j
  // and d denotes partial derivative

  // first k entries are diagonal of objective function

  // This is a symmetric matrix, fill the lower left triangle only.
  if (values == NULL) {
    // return the structure. 
    for (unsigned int kk = 0; kk < hessian_vector.size(); ++kk)
    {
      iRow[kk] = hessian_vector[kk].i;
      jCol[kk] = hessian_vector[kk].j;
    } 
  } // structure
  else {
    // return the values. 
    // g = cos( pi * sum_evens) == 1
    // g' = -pi * sin ( pi * sum_evens ), for each of the variables contributing to the sum
    // g''= -pi^2 cos ( pi * sum_evens ), for each pair of variables contributing to the sum
    // assuming all the coefficients are 1
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
    }
    
    // g = cos( 2 pi x) == 1
    // g' = - 2 pi sin ( 2 pi x ), for x_i
    // g'' = - 4 pi^2 cos ( 2 pi x), for x_i only
    {
      if (debugging)
      {
        printf("\nwave int non-zero hessian values:\n");
      }
      for (int i=0; i < data->num_variables(); ++i)
      {
        const int k = i + wave_int_constraint_start;
        // diagonal entries, again
        
        const double hg_ii = eval_hess_int_x(x[i]); // e.g. -4. * PI * PI * cos( 2. * PI * x[i] );
        values[i] += lambda[k] * hg_ii;

        if (debugging)
        {
          printf("x_%d (%f) : lambda(%f) * d^2 wave(x_ii)/ dx_i^2 (%f)\n", i, x[i], lambda[k], hg_ii);
        }
      }
    }
    if (debugging)
      printf("\n");
  } // values

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
  
  // todo report on how close the integer and sum-even constraints were satisfied!
  // or do that in the caller
}

} // namespace MeshKit
