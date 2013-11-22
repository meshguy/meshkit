// IAMilp.hpp
// Interval Assignment for Meshkit
//
// This is the solver-translation from Meshkit to the GLPK library
// for solving the MILP (integer solution, no fractional intervals)
//

#ifndef MESHKIT_IA_IAMILP_HP
#define MESHKIT_IA_IAMILP_HP

/** Optimization problem for getting an integer solution near the relaxed NLP solution.
 * 
 * Solved using GLPK
 *
 * Problem structure:
 *  minimize F = max(f_i) + sum of (f_i)
            where f_i = a*(x-I) + b(I-x) : x is assigned intervals, I is goal intervals
 *  st. g = sum side_1 - sum side_2 = 0 : for a mapping constraint
 *          sum side_1 + sum side_2 - sum side_3 >= 0 : triangle inequalities for midpoint subdivision, etc
 *          x_nlp-1 <= x <= x_nlp + num_curves on surface, x_nlp is the nlp solution
 *
 */

#include <vector>

#incude "IAWeights.hpp"

// forward declarations
class IAData;
class IASolution;
class GlpRepresentation;

class IAMilp
{
public:
  /** default constructor */
  // on input, solution_ptr contains a non-integer solution to the interval assignment problem
  IAMilp(const IAData *data_ptr, IASolution *solution_ptr); 

  /** default destructor */
  virtual ~IAMilp();

  // Find an integer solution.
  // on completion, the integer solution is placed into the solution_ptr that was passed to the constructor
  bool solve();
  
private:  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IAMilp();
  IAMilp(const IAMilp&);
  IAMilp& operator=(const IAMilp&);
  //@}
  
  // input data
  const IAData *data;
  // solution data
  IASolution *solution;
  
  // scratch space
  // objective weights of deltas
  IAWeights weights_plus;
  IAWeights weights_minus;
  
  const bool debugging;

  // much roundoff is tolerated before treating it as integer or satisfied
  const double constraint_tolerance;
  const double integrality_tolerance;
  
  // glpk-implementation specific 
  GlpRepresentation *lp;
  bool naturalConstraintsSet;
  bool deltasSet;
  bool weightsSet;
  bool maxesSet;
  bool solved;
  
  // translate our indices to glpk indexing
  int xStart; // column index
  int x_i(int i); // glpk index for curve_i
  int deltaStart; // column index
  int delta_plus_i(int i); // glpk col index for curve_i's value above its goal
  int delta_minus_i(int i);// glpk col index for curve_i's value below its goal
  int deltaConstraintStart; // first row index of constraint that forces x_i - I_i = p_i - n_i
  int delta_j(int j); // continuation of above 
  int mwp_i;  // column index of maximum of weighted delta_plus variable
  int mwp_j;  // first row index of constraint that calculates mwp variable
  int mwm_i;  // column index of maximum of weighted deltas_minus variable
  int mwm_j;  // first row index of constraint that calculates mwm variable
  int mw_i;  // column index of maximum of weighted deltas_minus variable
  int mw_j;  // row index of constaint that calculates mw > mwp variable, mw_j+1 calcs mw > mwm
  
  // substeps in setting up and solving the solution
  bool create_problem();
  bool destroy_problem();
  bool set_natural_constraints();
  bool set_deltas();
  bool weight_deltas_1();
  bool set_maxes();
  bool set_objectives_1();
  bool set_sum_even_constratins();
  bool set_bounds_1();
  bool set_bounds_2();
  bool set_bounds_3();
  bool set_bounds_4();
  bool set_bounds_A();
  bool set_bounds_B();
  bool glpk_solve(bool &optimal);
  bool get_solution();
  bool solution_is_integer();
  bool solution_satisfies_constraints();
  void print_constraint(int j);
  void print_solution(int i);

  
  
  // utilities
  // fill with [I,x_solution] or [x_solution,I], depending on which is smaller and defines a proper interval
  void get_x_bounds(const int i, double &lobound, double &hibound);
};

#endif
