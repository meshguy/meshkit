// IAMilp.cpp
// Interval Assignment for Meshkit
//
#include "IAMilp.hpp"
#include "IAData.hpp"
#include "IASolution.hpp"
#include <assert.h>

#include <iostream>
#include <sstream>
#include <string>

//usr/local/include/glpk.h
#include "glpk.h" 
#include <stdio.h>
#include <math.h>

class GlpRepresentation
{
  public: 
  glp_prob *glp;
};


// to do: design description of sum-even constraints, represent those in the MILP

// constructor
IAMilp::IAMilp(const IAData *data_ptr, IASolution *solution_ptr): 
data(data_ptr), solution(solution_ptr), lp(NULL), naturalConstraintsSet(false), deltasSet(false),
weightsSet(false), maxesSet(false), solved(false),
xStart(1), deltaStart(1), constraint_tolerance(1.0e-4), integrality_tolerance(1.0e-4),
debugging(true) //zzyk true
{
  assert(data_ptr);
  assert(solution_ptr);
  printf("\nIA MILP Problem size:\n");
  printf("  number of variables: %lu\n", data->I.size());
  printf("  number of constraints: %lu\n\n", data->constraints.size());
}


IAMilp::~IAMilp() {data = NULL;}

// index of first delta variable
// organized by x_i - I_i = dp_i - dn_i
// where col(x_i) = i+1, col(dp_i) = i + I_i + 1, col(dn_i) = i + 2*I_i +1

inline 
int IAMilp::x_i(int i)
{
  return i+xStart;
}

inline 
int IAMilp::delta_plus_i(int i) 
{ 
  return i + deltaStart;
}

inline 
int IAMilp::delta_minus_i(int i) 
{ 
  return i + deltaStart + data->num_variables(); 
}

inline 
int IAMilp::delta_j(int j) 
{ 
  return j + deltaConstraintStart;
}


bool IAMilp::create_problem()
{
  lp = new GlpRepresentation();
  lp->glp = glp_create_prob();
  if (!lp || !lp->glp) 
    return false;
  if (debugging) {
    glp_set_prob_name(lp->glp, "IA MILP");
    glp_set_obj_name(lp->glp, "interval deviations");
  }
  glp_set_obj_dir(lp->glp, GLP_MIN); // minimize objective function
  return true;
}

bool IAMilp::destroy_problem()
{
  naturalConstraintsSet = false;
  deltasSet = false;
  weightsSet = false;
  maxesSet = false;
  solved = false;
  glp_delete_prob(lp->glp);
  delete lp;
  lp = NULL;
  return true;
}

int glpk_bound_type(double low, double high)
{
  if (low < high) 
  {
    if (low == MESHKIT_IA_lowerUnbound)
      if (high == MESHKIT_IA_upperUnbound )
        return GLP_FR;  // -inf, inf
      else
        return GLP_UP;  // -inf, h
      else
        if (high == MESHKIT_IA_upperUnbound )
          return GLP_LO;  // l, inf
        else
          return  GLP_DB; // l, h
  }
  return GLP_FX; // equality l = h
 }

bool IAMilp::set_natural_constraints()
{
  // don't call this twice
  if (lp==NULL)
    return false;
  if (naturalConstraintsSet)
    return true;
  
  
  // allocate columns for variables (curves)
  // and set upper and lower bounds on them
  // cols must be allocated before the rows reference their indices
  int nv = data->num_variables(); 
  xStart = glp_add_cols(lp->glp, nv); // set here, used by x_i(i)
  for (int i=0; i<data->num_variables(); ++i)
  {
    if (debugging)
    {
      std::stringstream ss;
      ss << "x" << x_i(i);
      glp_set_col_name(lp->glp, x_i(i), ss.str().c_str());
    }
    glp_set_col_bnds(lp->glp, x_i(i), GLP_LO, 1., MESHKIT_IA_upperUnbound); // 1, inf    
    // enforce integer solution
    // glp_set_col_kind sets (changes) the kind of j-th column (structural variable) as specified by the parameter kind:
    // GLP_CV = continuous variable; GLP_IV = integer variable; GLP_BV = binary variable.
    glp_set_col_kind(lp->glp, x_i(i), GLP_IV);
  }

  // allocate rows for constraints
  const int num_rows = (int) data->constraints.size();
  const int next_free_row = glp_add_rows(lp->glp, num_rows);
  assert(next_free_row==1);
  
  // upper and lower bounds on constraints
  std::vector<int> ind;
  std::vector<double> val;
  for (int j=0; j<data->constraints.size(); ++j)
  {
    if (debugging)
    {
      std::stringstream ss;
      ss << "constraint " << j+1;
      glp_set_row_name(lp->glp, j+1, ss.str().c_str());
    }
    const int glp_j = next_free_row+j;
    const double low = data->constraints[j].lowerBound;
    const double high = data->constraints[j].upperBound;
    const int bound_type = glpk_bound_type( low, high );
    glp_set_row_bnds(lp->glp, glp_j, bound_type, low, high);
    
    // coefficients
    const int nnz = (int) data->constraints[j].M.size();
    ind.reserve(nnz);
    val.reserve(nnz);
    for (int k = 0; k<nnz; ++k)
    {
      ind[k] = x_i(data->constraints[j].M[k].col);
      val[k] = data->constraints[j].M[k].val;
    }
    // -1 since glp indexes from 1...
    glp_set_mat_row(lp->glp, glp_j, nnz, ind.data()-1, val.data()-1);    
  }
    
  naturalConstraintsSet = true;
  return true;
}

void IAMilp::get_x_bounds(const int i, double &lobound, double &hibound)
{
  const double xbound = solution->x_solution[i]; 
  const double Ibound = data->I[i];
  if (xbound < Ibound)
  {
    lobound = xbound;
    hibound = Ibound;
  }
  else
  {
    lobound = Ibound;
    hibound = xbound;
  }
}

bool IAMilp::set_deltas()
{
  if (!naturalConstraintsSet)
    return false;
  if (deltasSet)
    return true;
  
  // allocate columns for deltas 
  // and set upper and lower bounds on them
  int nv = 2*data->num_variables(); 
  deltaStart = glp_add_cols(lp->glp, nv);
  for (int i=0; i<data->num_variables(); ++i)
  {
    if (debugging)
    {
      std::stringstream ssp;
      ssp << "delta_plus_x" << x_i(i);
      glp_set_col_name(lp->glp, delta_plus_i(i), ssp.str().c_str());

      std::stringstream ssm;
      ssm << "delta_minus_x" << x_i(i);
      glp_set_col_name(lp->glp, delta_minus_i(i), ssm.str().c_str());
    }
    glp_set_col_bnds(lp->glp, delta_plus_i(i), GLP_LO, 0., MESHKIT_IA_upperUnbound); // 0, inf
    glp_set_col_bnds(lp->glp, delta_minus_i(i), GLP_LO, 0., MESHKIT_IA_upperUnbound); // 0, inf
  }
  
  // allocate rows for delta-setting constraints
  const int num_rows = data->num_variables();
  deltaConstraintStart = glp_add_rows(lp->glp, num_rows);
  for (int j=0; j<num_rows; ++j)
  {
    if (debugging)
    {
      std::stringstream ss;
      // index is correct, contraint for the variable == column of x_j
      ss << "row" << delta_j(j) << "_delta_constraints_for_x" << x_i(j); 
      glp_set_row_name(lp->glp, delta_j(j), ss.str().c_str());
    }
    // old way:
    // x_i - I_i = p_i - m_i  <=> x_i - p_i + m_i = I_i
    // rhs bound
    // const double bound = data->I[j]; 

    // new way: 
    // base delta on passed in nlp x_solution *and* original problem goals
    // we trust and use the nlp solution value
    // lower <= x - delta_plus + delta_minus <= upper
    // and we have non-zero deltas only when the milp solution is outside the interval [x_nlp,I]
    
    // to do: not sure if the above is sufficient, I may need to have a small weight to push variables back towards the goals, 
    // to decide who gets to reap the improvement when it is possible to reduce several back 
    // towards the goal.
    double lobound, hibound;
    get_x_bounds(j,lobound,hibound);
    glp_set_row_bnds(lp->glp, delta_j(j), GLP_DB, lobound, hibound);
    
    // coefficients
    const int ind[3] = {x_i(j), delta_plus_i(j), delta_minus_i(j)};
    const double val[3] = {1., -1., 1.};
    // -1 since indexes from 1
    glp_set_mat_row(lp->glp, delta_j(j), 3, ind-1, val-1);
  }
  
  deltasSet = true;
  return true;
}

bool IAMilp::set_sum_even_constratins()
{
  return true;
}

// old paper started with an integer solution, but not us, so the bounds here will not be identical.
bool IAMilp::set_bounds_1() 
{
  // from old paper:
  // Curve intervals can increase by at most one
  if (debugging) printf("MILP bounds 1: [floor x, ceil x + 1]\n");
  
  for (int i=0; i<data->num_variables(); ++i)
  {
    const double x_nlp = solution->x_solution[i];  // nlp solution
    // const double Ibound = data->I[i];
    double lo = floor(x_nlp);
    if (lo<1.) lo = 1.;
    const double hi = ceil(x_nlp + 0.99);
    glp_set_col_bnds(lp->glp, x_i(i), GLP_DB, lo, hi);
  }
  
  return true;
}

bool IAMilp::set_bounds_2() 
{
  // from old paper:
  // Curve intervals can increase or decrease by one
  if (debugging) printf("MILP bounds 2: [floor x-1, ceil x+1]\n");
  
  for (int i=0; i<data->num_variables(); ++i)
  {
    const double x_nlp = solution->x_solution[i];  // nlp solution
    // const double Ibound = data->I[i];
    double lo = floor(x_nlp - 0.99);
    if (lo<1.) lo = 1.;
    const double hi = ceil(x_nlp + 0.99);
    glp_set_col_bnds(lp->glp, x_i(i), GLP_DB, lo, hi);
  }
  
  return true;
}

bool IAMilp::set_bounds_3() 
{
  // from old paper:
  // Curve intervals can double, but can’t decrease.
  if (debugging) printf("MILP bounds 3: [floor x, ceil 2x]\n");
  
  for (int i=0; i<data->num_variables(); ++i)
  {
    const double x_nlp = solution->x_solution[i];  // nlp solution
    // const double Ibound = data->I[i];
    double lo = floor(x_nlp);
    if (lo<1.) lo = 1.;
    const double hi = ceil(x_nlp * 2);
    glp_set_col_bnds(lp->glp, x_i(i), GLP_DB, lo, hi);
  }
  
  return true;
}

bool IAMilp::set_bounds_4() 
{
  // from old paper:
  // Curve intervals can double, and can decrease by at most one. Interval-sum variables are bounded between
  // the floor and twice the ceiling of the pseudo-relaxed solution.
  if (debugging) printf("MILP bounds 4: [floor x - 1, ceil 2x]\n");

  
  for (int i=0; i<data->num_variables(); ++i)
  {
    const double x_nlp = solution->x_solution[i];  // nlp solution
    // const double Ibound = data->I[i];
    double lo = floor(x_nlp - 0.99);
    if (lo<1.) lo = 1.;
    const double hi = ceil(x_nlp * 2);
    glp_set_col_bnds(lp->glp, x_i(i), GLP_DB, lo, hi);
  }
  
  return true;
}

bool IAMilp::set_bounds_A() 
{
  // new idea: bound the change based on the number of curves involved in the constraint.  
  if (debugging) printf("MILP bounds A: [floor x - c, ceil x +c], where c is max number of curves in a constraint involving x.\n");
  
  // count the number of edges each curves constraint list
  std::vector<int> num_same_side( data->num_variables(), 0 );
  std::vector<int> num_opposite_side( data->num_variables(), 0 );
  // allocate rows for constraints
  const int num_constraints = (int) data->constraints.size();
  for (int j=0; j<num_constraints; ++j)
  {
    // number of non-zeros in the constraint
    const int nnz = (int) data->constraints[j].M.size();
    // number of constraints on each side
    int num_plus(0), num_minus(0);
    for (int k = 0; k<nnz; ++k)
    {
      if (data->constraints[j].M[k].val > 0)
        ++num_plus;
      else
        ++num_minus;
    }
    for (int k = 0; k<nnz; ++k)
    {
      int i = data->constraints[j].M[k].col;
      int ss, os;
      if (data->constraints[j].M[k].val > 0)
      {
        ss = num_plus; 
        os = num_minus;
      }
      else
      {
        os = num_plus; 
        ss = num_minus;
      }
      if (num_same_side[i] < ss)
        num_same_side[i] = ss;
      if (num_opposite_side[i] < os)
        num_opposite_side[i] = os;
    }
  }  
  for (int i=0; i<data->num_variables(); ++i)
  {
    const double x_nlp = solution->x_solution[i];  // nlp solution
    double lo = floor(x_nlp - num_opposite_side[i] - num_same_side[i]); // to do: work out some reasonable bounds here
    if (lo<1.) lo = 1.;
    double hi = ceil(x_nlp + num_opposite_side[i] + num_same_side[i]);  // to do: work out some reasonable bounds here
    glp_set_col_bnds(lp->glp, x_i(i), GLP_DB, lo, hi);
  }  
  return true;
}


bool IAMilp::set_bounds_B() 
{
  // unbounded
  if (debugging) printf("MILP bounds B: [1, infinity]\n");

  for (int i=0; i<data->num_variables(); ++i)
  {
    glp_set_col_bnds(lp->glp, x_i(i), GLP_LO, 1., MESHKIT_IA_upperUnbound);
  }  
  return true;
}

bool IAMilp::weight_deltas_1()
{
  if (!deltasSet)
    return false;

  // attach weights based on passed in nlp x_solution *and* original problem goals
  weights_minus.clear();
  weights_plus.clear();
  weights_minus.resize(data->num_variables());
  weights_plus.resize(data->num_variables());
  
  for (int i=0; i<data->num_variables(); ++i)
  {
    double lobound, hibound;
    get_x_bounds(i,lobound,hibound);
    assert( lobound + 1e-2 > 1.0 );
    assert( lobound <= hibound );
    weights_minus[i] = 1.3  / lobound;
    weights_plus[i] = 1.0 / hibound; // or should this be 1.0 / lobound?
  }

  weightsSet = true;
  return true;
}

bool IAMilp::set_maxes()
{
  if (!weightsSet)
    return false;
  
  // Mp = max ( delta_plus )
  // Mm = max ( delta_minus )
  // M = max( Mp, Mm )
  
  // Mwp = max ( w_p * delta_plus )
  // Mwm = max ( w_m * delta_minus )
  // Mw = max( Mwp, Mwm )

  // allocate rows for max-delta-setting constraints
  const int num_rows = 2*data->num_variables()+2;
  mwp_j = glp_add_rows(lp->glp, num_rows);
  mwm_j = mwp_j + data->num_variables();
  mw_j = mwp_j + 2 * data->num_variables();
  
  mwp_i = glp_add_cols(lp->glp, 3);
  mwm_i = mwp_i + 1;
  mw_i = mwp_i + 2;
  
  glp_set_col_bnds( lp->glp, mwp_i, GLP_LO, 0., MESHKIT_IA_upperUnbound);
  glp_set_col_bnds( lp->glp, mwm_i, GLP_LO, 0., MESHKIT_IA_upperUnbound);
  glp_set_col_bnds( lp->glp, mw_i, GLP_LO, 0., MESHKIT_IA_upperUnbound);
  
  if (debugging)
  {
    std::stringstream ss;
    ss << "var max weighted deltas plus";
    glp_set_col_name(lp->glp, mwp_i, ss.str().c_str());
    ss.clear(); ss.str("");
    ss << "var max weighted deltas minus";
    glp_set_col_name(lp->glp, mwm_i, ss.str().c_str());
    ss.clear(); ss.str("");
    ss << "var max weighted deltas plus and minus";
    glp_set_col_name(lp->glp, mw_i, ss.str().c_str());

    for (int i = 0; i<data->num_variables(); ++i)
    {
      ss.clear(); ss.str("");
      int r = mwp_j + i;
      ss << "row" << r << " max weighted delta_plus" << x_i(i) << " constraint";
      glp_set_row_name(lp->glp, r, ss.str().c_str());
      ss.clear(); ss.str("");
      r = mwm_j + i;
      ss << "row" << r << " max weighted delta_minus" << x_i(i) << " constraint";
      glp_set_row_name(lp->glp, r, ss.str().c_str());
    }
    ss.clear(); ss.str("");
    int r = mw_j; 
    ss << "row" << r << " weighted_max_delta_ge_deltaplus";
    glp_set_row_name(lp->glp, r, ss.str().c_str());
    ss.clear(); ss.str(""); 
    ++r; 
    ss << "row" << r << " weighted_max_delta_ge_deltaminus";
    glp_set_row_name(lp->glp, r, ss.str().c_str());
  }

  for (int i = 0; i<data->num_variables(); ++i)
  {
    // mwp - dp_i >= 0
    glp_set_row_bnds(lp->glp, mwp_j + i, GLP_LO, 0., MESHKIT_IA_upperUnbound);

    // mwn - dp_i >= 0
    glp_set_row_bnds(lp->glp, mwm_j + i, GLP_LO, 0., MESHKIT_IA_upperUnbound);

    // coefficients
    const int indp[2] = {mwp_i, delta_plus_i(i)};
    const double valp[2] = {1., -weights_plus[i] };
    // -1 since indexes from 1
    glp_set_mat_row(lp->glp, mwp_j + i, 2, indp-1, valp-1);

    // coefficients
    const int indm[2] = {mwm_i, delta_minus_i(i)};
    const double valm[2] = {1., -weights_minus[i]};
    // -1 since indexes from 1
    glp_set_mat_row(lp->glp, mwm_j + i, 2, indm-1, valm-1);
  }

  // mw > mwp, mwm
  glp_set_row_bnds(lp->glp, mw_j, GLP_LO, 0., MESHKIT_IA_upperUnbound);
  const double val[2] = {1., -1.};
  const int indpp[2] = {mw_i, mwp_i};
  glp_set_mat_row(lp->glp, mw_j, 2, indpp-1, val-1);
  glp_set_row_bnds(lp->glp, mw_j+1, GLP_LO, 0., MESHKIT_IA_upperUnbound);
  const int indmm[2] = {mw_i, mwm_i};
  glp_set_mat_row(lp->glp, mw_j+1, 2, indmm-1, val-1);

  maxesSet = true;
  return true;
  
}

bool IAMilp::set_objectives_1()
{
  // to do: set objective function
  // to do: set bounds on x_i for some strategy
  
  /*
  // sum of deltas for testing
  for (int i = 0; i<data->num_variables(); ++i)
  {
    glp_set_obj_coef(lp->glp, delta_plus_i(i),  1.0);    
    glp_set_obj_coef(lp->glp, delta_minus_i(i), 1.0);    
  }
  */
  
  // obj = (num_variables + 1) * max_weighted_detlas + sum weighted_deltas
  glp_set_obj_coef(lp->glp, mw_i, 1 + data->num_variables() );
  // sum : this is the part that makes the MILP slow
  // to do: experiment with lex min max, as that might still be faster than branch and cut
  for (int i = 0; i<data->num_variables(); ++i)
  {
//    glp_set_obj_coef(lp->glp, delta_plus_i(i),  weights_plus[i]);    
//    glp_set_obj_coef(lp->glp, delta_minus_i(i), weights_minus[i]);    
  }
  
  return true;
}


bool IAMilp::glpk_solve(bool &optimal)
{
  assert(!solved); // already solved?
  optimal = false;
  
  // print problem we are solving
  if (debugging)
    glp_write_lp(lp->glp, 0, "zzykoutput"); // only works in command line, not within xcode

  
  // to do: identify independent sub-problems and solve them separately
  // to do: sum-even constraints
  // to do: figure out how to tell it about the feasible relaxed solution to the nlp, or to limit its time searching for the simplex solution
  
  // MILP branch and cut: limit of 40 variables or so is practical for time
  // time grows super-linearly
  // this is for the lp time limit I think
  //void lpx_set_real_parm(LPX *lp, int parm, double val); val is milliseconds?
  // see tm_lim instead
  lpx_set_real_parm( lp->glp, LPX_K_TMLIM, 10. );  // this doesn't seem to do anything
  // simplex for relaxed solution
//  bool relaxed_success = glp_simplex(lp->glp, NULL) == 0;
/*   {
   solved = true;
   return true;
   }
   */
  
  
  // find integer solution using branch and cut
  glp_iocp parm;
  glp_init_iocp(&parm);
  if (!debugging)
    parm.msg_lev = GLP_MSG_OFF;
  // parm.pp_tech = GLP_PP_ALL; 
//  parm.presolve = GLP_OFF; // not needed if do explicit solve above
  parm.presolve = GLP_ON; // solve relaxed solution first, remove redundant constraints, fixed variables, ...
  // parm.mip_gap = 0; // experiment with larger values to get good-enough integer solutions sooner
  parm.tm_lim = 500.; // time limit, in milliseconds. 1000 = 1 second. to do: modify this for which set of bounds
  int status = glp_intopt(lp->glp, &parm);
  if (status == GLP_ETMLIM) // timed out
  {
    // OK if it is non-optimal, but not if it is non-integer
    solved = solution_is_integer() && solution_satisfies_constraints();
    return solved;
  }
  else if (status == 0 || status == GLP_EMIPGAP)
  {
    optimal = true;
    solved = true;
    return true;
  }
  
  return false;
}


bool IAMilp::solution_satisfies_constraints()
{
  bool unsatisfied_found(false);
  for (int j = 0; j < data->constraints.size(); ++j)
  {
    double slack=0.;
    const IAData::constraintRow & c = data->constraints[j];
    for (std::vector<IAData::sparseEntry>::const_iterator i = c.M.begin(); i < c.M.end(); ++i)
    {
      const double xv = glp_mip_col_val( lp->glp, x_i(i->col));
      slack += xv * i->val;
    }
    if (data->constraints.front().upperBound == data->constraints.front().lowerBound)
    {
      if ( fabs(slack - data->constraints.front().upperBound) > constraint_tolerance)
      {
        unsatisfied_found = true;
        if (debugging) 
          print_constraint(j);
        else
          return false;
      }
    }
    else
    {
      printf(" in [%1.1f,%1.1f]", data->constraints.front().upperBound, data->constraints.front().lowerBound );
      if ( ( (data->constraints.front().upperBound - slack) < -constraint_tolerance ) ||
           ( (data->constraints.front().lowerBound - slack) > constraint_tolerance ) )
      {
        unsatisfied_found = true;
        if (debugging) 
        {
          printf("Unsatisfied constraint: ");
          print_constraint(j);
        }
        else
          return false;
      }
    }
  }
  return !unsatisfied_found;
}

void IAMilp::print_constraint(int j)
{
  //printf("constraint %d: ", j);
  double slack=0.;
  const IAData::constraintRow & c = data->constraints[j];
  for (std::vector<IAData::sparseEntry>::const_iterator i = c.M.begin(); i < c.M.end(); ++i)
  {
    // nlp solution const double xv = solution->x_solution[i->col];
    // const double xv = glp_get_col_prim( lp->glp, x_i(i->col) ); relaxed solution only
    const double xv = glp_mip_col_val( lp->glp, x_i(i->col));
    slack += xv * i->val;
    printf(" %1.0f x%d (%1.3f) ", 
           i->val, i->col, xv );
  }
  if (data->constraints.front().upperBound == data->constraints.front().lowerBound)
  {
    printf(" = %1.1f ", data->constraints.front().upperBound );  
    if ( fabs(slack - data->constraints.front().upperBound) > constraint_tolerance )
      printf(" UNSATISFIED ");
  }
  else
  {
    printf(" in [%1.1f,%1.1f]", data->constraints.front().upperBound, data->constraints.front().lowerBound );
    if ( (data->constraints.front().upperBound - slack) < -constraint_tolerance )
      printf(" UNSATISFIED BEYOND UPPERBOUND ");
    if ( (data->constraints.front().lowerBound - slack) > constraint_tolerance )
      printf(" UNSATISFIED BELOW LOWERBOUND ");
    
  }
  printf(" (%1.1f)\n", slack );
}

bool IAMilp::solution_is_integer()
{
  bool non_integer_found(false);
  for (int i=0; i<data->num_variables(); ++i)
  {
    const double x = glp_mip_col_val(lp->glp, x_i(i));
    if ( x < 0. || fabs(x - round(x)) > integrality_tolerance)
    {
      if (debugging)
      {
        printf("Noninteger variable: ");
        print_solution(i);
        non_integer_found = true;
      }
      else
        return false;
    }
  }
  return !non_integer_found;
}

void IAMilp::print_solution(int i)
{
  const double x = glp_mip_col_val(lp->glp, x_i(i));
  const double xp = glp_mip_col_val( lp->glp, delta_plus_i(i) );
  const double xm = glp_mip_col_val( lp->glp, delta_minus_i(i) );
  printf("%d: goal (%1.1f) x_nlp (%1.1f): x (%1.1f) plus (%1.1f) minus (%1.1f)\n", 
         x_i(i), data->I[i], solution->x_solution[i], x, xp, xm);
}


bool IAMilp::get_solution()
{
  // print solution
  if (debugging)
  {
    printf("MILP solution:\n");  
    // print delta values
    printf("x* and deltas\n");
    for (int i=0; i<data->num_variables(); ++i)
    {
      print_solution(i);
    }
    printf("constraints\n");
    for (int j = 0; j < data->constraints.size(); ++j)
    {
      print_constraint(j);
    }
    // todo: print objective value
  }
  
  // collect the solution into the vector
  for (int i = 0; i<data->num_variables(); ++i)
  {
    solution->x_solution[i] = glp_mip_col_val( lp->glp, x_i(i) );
  }

  // objective function values
  // z = glp_get_obj_val(lp);
  return true;
}

// return the solution
bool IAMilp::solve()
{
  bool success;
  success = create_problem();
  assert(success);

  success = set_natural_constraints();
  assert(success);
  
  success = set_sum_even_constratins();
  assert(success);
  
  success = set_deltas();
  assert(success);
  
  success = weight_deltas_1();
  assert(success);
  
  success = set_maxes();
  assert(success);

  success = set_objectives_1();
  assert(success); 

  // cycle through bounds
  // accept any integer solution, even if suboptimal
  for (int c=0; c<6; ++c)
  {
    // preference order: 2, A, 1, 3, 4, B
    switch (c)
    {
      case 0: 
        success = set_bounds_1(); 
        break;
      case 1:
        success = set_bounds_2(); 
        break;
      case 2:
        success = set_bounds_A();
        break;
      case 3: 
        success = set_bounds_3(); 
        break;
      case 4: 
        success = set_bounds_4(); 
        break;
      case 5: 
        success = set_bounds_B(); 
        break;
    }
    assert(success);
    
    // to do: cycle through solving sum of deltas, then add max of deltas or other heuristics to improve the solution
    // to do: grant longer time limit for set_bounds_B
  
    bool optimal(false);
    success = glpk_solve(optimal);
    if (success && !optimal)
    {
      if (debugging)
        printf("Stopping at sub-optimal but integer MILP solution.\n");
    }
    if (success)
      break;
  }
  
  success = get_solution();
  assert(success);
  
  destroy_problem();
  return true;
}