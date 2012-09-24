// IASolver.cpp
// Interval Assignment for Meshkit
//
#include "IASolverInt.hpp"
#include "IARoundingNlp.hpp"
#include "IASolverRelaxed.hpp"
// #include "IAMilp.hpp" // glpk based solution too slow

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "IpIpoptApplication.hpp"

IASolverInt::IASolverInt(const IAData * ia_data_ptr, const IASolution *relaxed_solution_ptr) 
: iaData(ia_data_ptr), debugging(true), myianlp(NULL) 
{ 
  ipData.initialize(relaxed_solution_ptr->x_solution);
}

/** default destructor */
IASolverInt::~IASolverInt() { iaData = NULL;}

bool IASolverInt::calculate_rounding_value(const int i, double &obj_increase, int &y_int )
{
  const double x = x_solution[i];
  y_int = floor(x + 0.5); // round to nearest integer
  double y = y_int;
  if ( fabs(x-y) > 1e-4) // not already integer within tolerance 
  {
    // non-integer found
    // round x away from I
    const double goal = iaData->I[i];
    // it was forced away from its goal, round in the forced direction further
    // else leave it rounded to the nearest integer
    if (fabs(x - goal) > 1.e-2) 
    {
      y_int = (x>goal) ? ceil(x) : floor(x);
      y = y_int;
    }
    assert(y >= 1.);
    
    // if it is already constrained to y_int, then consider it as already satisfying it.
    if ( ipData.varIntegerBound[i] == y )
    {
      obj_increase = 0.;
      return false;
    }

    // zzyk
    
    // todo: if x is very close to the goal, maybe we should retain freedom to round it the other direction, try both, etc.
 
    // is this the best non-integer we've found so far?
    // best means the one that changing it to integer increases the obj function the least
    // double f = myianlp->eval_R_i( goal, x ); 
    // compare to relaxed solution, not the latest pseudo-integer solution
    double f = myianlp->eval_R_i( goal, ipData.relaxedSolution[i] ); 
    double g = myianlp->eval_R_i( goal, y );
    obj_increase = g - f;           
    // assert( obj_increase > 0. ); this is only invariant if we use f = R(x), not f = R(relaxed).
    return true;
  }
  //else
  obj_increase = 0.;
  return false;
}

bool IASolverInt::find_one_non_integer(int &i_nonint, int &x_bound)
{
  i_nonint = -1; // none found
  double obj_increase_nonint( std::numeric_limits<double>::max() );
  // find the one variable whose rounding causes the least change to the objective function
  // to do: alt is to pick one at random, with some chance
  for (int i = 0; i<iaData->num_variables(); ++i)
  {
    double obj_increase;
    int y_int;
    if (calculate_rounding_value(i, obj_increase, y_int))
    {
      if (obj_increase < obj_increase_nonint)
      {
        i_nonint = i;
        x_bound = y_int;
        obj_increase_nonint = obj_increase;
      }
    }
  }
  if (i_nonint > -1)
    return true;
  return false;
}
 
void IASolverInt::constrain_integer(const int i_nonint, const int x_bound)
{
  ipData.oldBound[i_nonint] = ipData.varIntegerBound[i_nonint];;
  ipData.varIntegerBound[i_nonint] = (double) x_bound;
}

// find the one variable whose rounding causes the least change to the objective function
// unused at the moment
bool IASolverInt::constrain_one_non_integer(int &i, int &b)
{
  if (find_one_non_integer(i, b))
  {
    constrain_integer(i, b);
    if (debugging)
    {
      report_one_constraint(i, b);
    }
    return true;
  }
  return false;
}

/*
if (i>num_variables())
{
  assert(phase >= 1);
  const int j = i - num_variables();
  printf("  constraining sum[%d] from %e to %c %d\n", j, sum_even_value(j), '>', b);
}
else

 
 if (i>num_variables())
 {
 assert(phase >= 1);
 const int j = i - num_variables();
 printf("  new sum[%d] value is %e %c %d\n", j, sum_even_value(j), '>', b);
 }
 else
 {
*/

void IASolverInt::report_one_constraint(const int i, const int b)
{
  // assumes x_solution[i] is current non-integer solution
  printf("  constraining x[%d] from %e (%e) to %c %d\n", i, x_solution[i], iaData->I[i], b > iaData->I[i] ? '>' : '<', b);
}

void IASolverInt::report_one_non_integer(const int i, const int b)
{
  printf("  new x[%d] value is ", i );
  if (fabs( floor(x_solution[i]+0.5) - x_solution[i] ) < 1.0e-2)
    printf("%d", (int) floor(x_solution[i]+0.5));
  else
    printf("%e", x_solution[i]);
  printf(" (%e) %c %d\n", iaData->I[i], b > iaData->I[i] ? '>' : '<', b);
}

void IASolverInt::back_off(RoundingMap &rounding_map)
{
  // undo half of the constraints
  // todo: 
  RoundingMap::iterator undo_start = rounding_map.begin(); 
  int num_to_keep = (int) rounding_map.size() / 2; // roundoff, rounding to zero is OK
  std::advance( undo_start, num_to_keep );
  for(RoundingMap::const_iterator i = undo_start; i != rounding_map.end(); ++i)
  {
    const int x_i = i->second.first;
    if (debugging)
    {
      printf(" backing off:");
      report_one_constraint(x_i, ipData.oldBound[x_i]);
    }
    ipData.varIntegerBound[x_i] =  ipData.oldBound[x_i]; 
  }
  rounding_map.erase(undo_start, rounding_map.end()); 
}

bool IASolverInt::find_many_non_integer(RoundingMap &rounding_map)
{
  // to do: this heuristic is bad because it easily paints us into a corner 
  // to do: we should instead pick some constant fraction of the variable per *side* of a constraint, 
  // provided there is at least one other unconstrained variable on that side. And never constrain the last k-sides (k=2, or 3? or 4?) until the end, when we just constrain one at a time.
  
  // since we don't know how many are non-integer, put them all in the map, then discard the ones we won't constrain
  rounding_map.clear();
  for (int i = 0; i<iaData->num_variables(); ++i)
  {
    double obj_increase;
    int y_int;
    if (calculate_rounding_value(i, obj_increase, y_int))
    {
      // add to map
      // printf("inserting obj %e, i %d, x %e, y %d\n", obj_increase, i, x_solution[i], y_int);
      rounding_map.insert(std::make_pair(obj_increase,std::make_pair(i,y_int)));  
    }
  }
  if (rounding_map.empty())    
    return false;
  //const unsigned long num_to_round = 1;
  const unsigned long num_to_round = 1; // 1 + rounding_map.size() / 8; // todo: experiment with the 16 part
  if (debugging)
  {
    printf("Found %ld of %d non-integers, rounding %ld of them.",
           rounding_map.size(), iaData->num_variables(), num_to_round);
  }
  // remove everything after the "num_to_round" element
  RoundingMap::iterator r = rounding_map.begin();
  advance(r, num_to_round);
  rounding_map.erase(r, rounding_map.end()); 
  // todo: if erasing (or subsequent iterations) is slow due to tree rebalancing, then instead we should transfer the first num_to_round elements to a std::vector and return that instead, and discard the map

  if (debugging)
  {
    if (rounding_map.size() > 1)
      printf(" obj_increase range %e to %e\n",
             rounding_map.begin()->first, rounding_map.rbegin()->first);
    else if (rounding_map.size() == 1)
      printf( " obj_increase expected %e\n",
             rounding_map.begin()->first);
    else
      printf( " rounding nothing\n");
  }
  return true;
}

bool IASolverInt::constrain_many_non_integer(RoundingMap &rounding_map)
{
  if (find_many_non_integer(rounding_map))
  {
    for(RoundingMap::const_iterator i = rounding_map.begin(); i != rounding_map.end(); ++i)
    {
      const int x_i = i->second.first;
      const int b = i->second.second;
      constrain_integer(x_i, b); // index, new integer bound
      if (debugging)
      {
        report_one_constraint(x_i, b);
      }
    }
    // todo: heuristic to decrease the number of constrained variables if the solver returns infeasible, or it decreased these variables beyond the new upper bounds . save the old bounds to restore them.
    return true;
  }
  return false; 
}

void IASolverInt::report_many_non_integer(RoundingMap &rounding_map)
{
  for(RoundingMap::const_iterator i = rounding_map.begin(); i != rounding_map.end(); ++i)
  {
    report_one_non_integer(i->second.first, i->second.second);
  }
}

bool IASolverInt::solve_minlp()
{
  // set up and call the separate IARoundingNlp, which has a linear objective to get a natural integer solution 
  // the intuition is this will solve integrality for  most variables all at once

  if (debugging)
  {
    printf("IASolver::solve_minlp() - ");        
    printf("Attempting to find a naturally-integer solution to an NLP.\n");
  }

  
  // solver setup  
  using namespace Ipopt;
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  
  // convergence parameters
  // see $IPOPTDIR/Ipopt/src/Interfaces/IpIpoptApplication.cpp
  // our real criteria are: all integer, constraints satisfied. How to test the "all_integer" part?
  app->Options()->SetNumericValue("tol", 1e-2); //"converged" if NLP error<this, default is 1e-7. Obj are scaled to be >1, so e-2 is plenty
  app->Options()->SetNumericValue("max_cpu_time", sqrt( iaData->num_variables() ) ); // max time allowed in seconds
  app->Options()->SetIntegerValue("max_iter", 3 * iaData->num_variables() ); // max number of iterations
  // app->Options()->SetNumericValue("primal_inf_tol", 1e-2 ); 
  app->Options()->SetNumericValue("dual_inf_tol", 1e-2 ); // how close to infeasibility?
  app->Options()->SetNumericValue("constr_viol_tol", 1e-6 ); // by how much can constraints be violated?
  app->Options()->SetNumericValue("compl_inf_tol", 1e-2 ); // max norm of complementary condition

  // second criteria convergence parameters: quit if within this tol for many iterations
  app->Options()->SetIntegerValue("acceptable_iter", 4 + sqrt( iaData->num_variables() ) ); //as "tol"
  app->Options()->SetNumericValue("acceptable_tol", 1e-1 ); //as "tol"

  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetIntegerValue("print_level", 1);  // 0 to 12, most. Default is 5
  // uncomment next line to write the solution to an output file
  // app->Options()->SetStringValue("output_file", "IA.out");  
  // The following overwrites the default name (ipopt.opt) of the options file
  // app->Options()->SetStringValue("option_file_name", "IA.opt");
  
  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }
  
  // problem setup
  // a couple of different models, simplest to more complex
  IARoundingNlp *myianlp = new IARoundingNlp(iaData, &ipData, this);
  // IARoundingFarNlp *myianlp = new IARoundingFarNlp(iaData, &ipData, this);
  // IARoundingFar3StepNlp *myianlp = new IARoundingFar3StepNlp(iaData, &ipData, this); // haven't tested this. It compiles and runs but perhaps isn't correct
  // IAIntWaveNlp *myianlp = new IAIntWaveNlp(iaData, &ipData, this); // haven't tested this. It compiles and runs but perhaps isn't correct
  SmartPtr<TNLP> mynlp = myianlp; // Ipopt requires the use of smartptrs!

  bool try_again = true;
  int iter = 0;
  do {
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp); // the inherited IANlp
    
    if (status == Solve_Succeeded) {
      printf("\n\n*** The problem solved!\n");
    }
    else {
      printf("\n\n*** The problem FAILED!\n");
      return false;
    }
    
    // The problem should have been feasible, but it is possible that it had no integer solution.
    // figure out which variables are still integer
    
    // check solution for integrality
    
    // debug
    bool sum_even_sat;
    bool constraints_sat = IASolverRelaxed::constraints_satisfied( iaData, this, sum_even_sat, true );
    
    IASolution nlp_solution;
    nlp_solution.x_solution = x_solution; // vector copy
    IPData::round_solution(nlp_solution.x_solution);
    bool sum_even_sat2;
    if (IASolverRelaxed::constraints_satisfied(iaData, &nlp_solution, sum_even_sat2, debugging))
    {
      // rounding was a valid integer solution
      x_solution = nlp_solution.x_solution;
      return true;
    }

    // find out which vars were not integer, 
    // try rounding their weights and resolving
    // bool int_sat = solution_is_integer();
 //   myianlp->randomize_weights_of_non_int(); // try again? debug
    
    ++iter;
    try_again = iter < 4 + sqrt(iaData->num_variables());
    try_again = false; // debug
  } while (try_again);

  
  // todo: update partially-integer solution, perhaps using ipData - figure out how we're going to use it, first, for what structure makes sense.
  
  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.  
  return status == Solve_Succeeded;
  
}

void IASolverInt::cleanup()
{
  // delete myianlp;
  myianlp = NULL;
}


bool IASolverInt::solve_rounding_heuristic()
{
  // given non-integer feasible solution, round variables to obtain an integer one
  // if the input solution is (partially) integer, then try to respect those, but it may be overconstrained already and we may need to relax those.
  // set up and call the separate IARoundingHeuristicMINLP, which is rounding one at a time using the NLP cubic objective used for the relaxed solution

  using namespace Ipopt;
  
  if (debugging)
  {
    printf("Applying rounding heuristic to get integer solution\n");
  }

  // solver setup
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 1e-7); // 2 seems close enough, could do less, say .1
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetIntegerValue("print_level", 1);  // 0 to 12, most. Default is 5
  // uncomment next line to write the solution to an output file
  // app->Options()->SetStringValue("output_file", "IA.out");  
  // The following overwrites the default name (ipopt.opt) of the options file
  // app->Options()->SetStringValue("option_file_name", "IA.opt");
  
  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }
  
  // problem setup
  // match the relaxed problem, but with extra constraints
  IARoundingHeuristicMINLP *myianlp = new IARoundingHeuristicMINLP(iaData, &ipData, this);
  SmartPtr<TNLP> mynlp = myianlp; // Ipopt requires the use of smartptrs!

  for (int phase=0; phase<1; ++phase)
  {
    if (phase==0)
    {
      //set bounds variables that are already integer
      ; //todo
    }
    else
    {
      // relaxe those prior bounds
      // todo
      ;
    }
      
      // Ask Ipopt to solve the problem
      status = app->OptimizeTNLP(mynlp); // the inherited IANlp
    
    if (status == Solve_Succeeded) {
      printf("\n\n*** The initial mi-nlp problem solved!\n");
    }
    else {
      printf("\n\n*** The initial mi-nlp problem FAILED!\n");
      return false;
    }
    
    // make integer
    // Pick the variable whose rounding to the next further integer from its goal results in the lowest individual contribution (or increase in contribution from its current value) to the objective function.
    // Round it that way, setting a lower bound on its value (but leaving objective, etc., the same)
    // Re-solve the NLP. Some old vars might have their delta directions or values flip. That is OK.
    // Iterate.
    
    int iter = 0;
    printf("\n\n*** Solving integer solution.\n"); 
    
/* 
    // many at a time speedup, could skip and still get correct output
    if (debugging)
      printf("Many-at-a-time rounding\n");
    RoundingMap rounding_map;
    // phase = -1;
    while (constrain_many_non_integer(myianlp,ipData,rounding_map))
    { 
      
      if (debugging)
        printf("MINLP many-at-once iter %d\n", iter);
      else
        printf("%d ", iter);
      
      // solve for many new constraints, but reduce those constraints if the problem becomes infeasible
      bool try_again = false;
      do 
      {
        try_again = false;
        status = app->OptimizeTNLP(mynlp); // double-check that resolving works, w/out initialize, etc 
        if (status != Solve_Succeeded && status != Solved_To_Acceptable_Level)
        {
          back_off(ipData,rounding_map);
          try_again = rounding_map.size();
        }
      }
      while (try_again);
      if (rounding_map.size() == 0) // no progress was made
        break;
      
      // to do: what to do if status != Solve_Succeeded ?
      // assert(b >= iaData->I[i] && x_solution[i] >= b || b < iaData->I[i] && x_solution[i] <= b); // assert the new constraint was satisfied, one-at-a-time
      if (debugging)
      {
        // report_one_non_integer(i, b);
        report_many_non_integer(rounding_map);
      }
      ++iter;
      // todo: quit at an iteration limit? but a single variable could be increased several times...
    }
*/
    
    // one at a time
    if (debugging)
      printf("\nOne-at-a-time rounding.\n");
    int i; // non-integer variable to make integer
    int b; // value to make it integer
    iter = 0;
    
    while (constrain_one_non_integer(i,b))
    {
      
      if (debugging)
        printf("MINLP one-at-a-time iter %d\n", iter);
      else
        printf("%d ", iter);
      
      status = app->OptimizeTNLP(mynlp); // double-check that resolving works, w/out initialize, etc 
      // to do: what to do if status != Solve_Succeeded ?
      const double g = iaData->I[i];
      const double x = x_solution[i];
      assert(b >= g && x + 1.e-2 >= b || b < g && x <= b + 1.e-2); // assert the new constraint was satisfied, one-at-a-time
      if ((status != Solve_Succeeded && status != Solved_To_Acceptable_Level) || 
          ! (b >= g && x + 1.e-2 >= b || b < g && x <= b + 1.e-2))
      {
        // see ipopt/CoinIpopt/Ipopt/src/Interfaces/IpReturnCodes_inc.h for the enum of possible values
        // status = Infeasible_Problem_Detected; 
        if (status == Infeasible_Problem_Detected )
          printf("IPOPT Infeasible_Problem_Detected\n");
        else
          printf("IPOPT bad solution return code %d\n", status);
        break;
      }
      if (debugging)
      {
        report_one_non_integer(i, b);
      }
      ++iter;
      // todo: quit at an iteration limit? but a single variable could be increased several times...
    }
    
    if (status != Solve_Succeeded)
      if (phase == 0)
      {
        if (debugging)
          printf("phase 0 failed, trying again with relaxed integer constraints.\n");
      }
      else
        printf("***Integer Fail!***\n");
      else
      {
        if (debugging)
        {
          // print_solution(ipData);
        }
        // successful !
        break;
      }
  }
  
  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.  
  return status == Solve_Succeeded;
  
}

/*
bool IASolverInt::solve_milp()
{
  IAMilp milp(this,this);
  bool succeeded=milp.solve();
  return succeeded;
}  
*/

bool IASolverInt::solve()
{
  solve_minlp();
  bool success = solution_is_integer();
  if (success)
  {
    printf("minlp produced integer solution\n");
  }
  else
  {
    // todo: rather than applying the rounding heuristic, implement a form of IARoundingNlp with larger variable bounds, but still with a natural integer solution, by extending x to also depend on a delta_plus and delta_minus extending x beyond xl and xl+1, i.e. x = xl (const) + xh (0-1 var) + delta_plus - delta_minus. With linear objective with weight for xh as before, but weight for delta_plus to be f( xl + 2 ) - f (xl + 1), delta_minus f( xl - 2) - f(xl -1)
    return false; // debug;
    solve_rounding_heuristic();
    success = solution_is_integer();
  }
  return success;
}

