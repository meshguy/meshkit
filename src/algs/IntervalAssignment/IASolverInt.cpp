// IASolver.cpp
// Interval Assignment for Meshkit
//
#include "IASolverInt.hpp"
#include "IARoundingNlp.hpp"
#include "IASolverRelaxed.hpp"
#include "IAMINlp.hpp"
// #include "IAMilp.hpp" // glpk based solution too slow

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "IpIpoptApplication.hpp"

IASolverInt::IASolverInt(const IAData * ia_data_ptr, const IASolution *relaxed_solution_ptr) 
: iaData(ia_data_ptr), debugging(true)
{ 
  ipData.initialize(relaxed_solution_ptr->x_solution);
}

/** default destructor */
IASolverInt::~IASolverInt() { iaData = NULL;}


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
  ;
}

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
    // solve_rounding_heuristic();
    success = solution_is_integer();
  }
  return success;
}

