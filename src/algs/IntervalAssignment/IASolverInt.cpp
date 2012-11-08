// IASolverInt.cpp
// Interval Assignment for Meshkit
//
#include "IASolverInt.hpp"
#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"
#include "IARoundingNlp.hpp"
#include "IASolverRelaxed.hpp"
#include "IAMINlp.hpp"
#include "IAIntWaveNlp.hpp"

// #include "IAMilp.hpp" // glpk based solution too slow

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "IpIpoptApplication.hpp"

namespace MeshKit 
{
    
IASolverInt::IASolverInt(const IAData * ia_data_ptr, IASolution *relaxed_solution_ptr, 
  const bool set_silent) 
  : IASolverToolInt(ia_data_ptr, relaxed_solution_ptr, true), 
  silent(set_silent), debugging(true)
{ 
  ip_data(new IPData);
  // initialize copies relaxed solution, then we can overwrite relaxed_solution_pointer
  ip_data()->initialize(relaxed_solution_ptr->x_solution); 
}

/** default destructor */
IASolverInt::~IASolverInt() 
{
  delete ip_data();
}
    
bool IASolverInt::solve_intwave()
{
  if (debugging)
  {
    printf("IASolverInt::solve_intwave() - ");        
    printf("Attempting to enforce an integer and even solution to the relaxed NLP by adding sine-wave constraints.\n");
  }
  
  // solver setup  
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
/* try leaving defaults
  // convergence parameters
  // see $IPOPTDIR/Ipopt/src/Interfaces/IpIpoptApplication.cpp
  // our real criteria are: all integer, constraints satisfied. How to test the "all_integer" part?
  app->Options()->SetNumericValue("tol", 1e-6); //"converged" if NLP error<this, default is 1e-7. Obj are scaled to be >1, so e-2 is plenty // was 1e-2
  app->Options()->SetNumericValue("max_cpu_time", sqrt( iaData->num_variables() ) ); // max time allowed in seconds
  app->Options()->SetIntegerValue("max_iter", 3 * (10 + iaData->num_variables() ) ); // max number of iterations
  // app->Options()->SetNumericValue("primal_inf_tol", 1e-2 ); 
  app->Options()->SetNumericValue("dual_inf_tol", 1e-2 ); // how close to infeasibility? // was 1e-2
  app->Options()->SetNumericValue("constr_viol_tol", 1e-2 ); // by how much can constraints be violated?
  app->Options()->SetNumericValue("compl_inf_tol", 1e-6 ); // max norm of complementary condition // was 1e-2
  
  // second criteria convergence parameters: quit if within this tol for many iterations
  // was  app->Options()->SetIntegerValue("acceptable_iter", 4 + sqrt( iaData->num_variables() ) ); //as "tol"
  app->Options()->SetNumericValue("acceptable_tol", 1e-6 ); //as "tol" was 1e-1
  
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  // print level 0 to 12, Ipopt default is 5
  const int print_level = (silent) ? 0 : 1;  // simple info is 1, debug at other values
  app->Options()->SetIntegerValue("print_level", print_level);  
  // uncomment next line to write the solution to an output file
  // app->Options()->SetStringValue("output_file", "IA.out");  
  // The following overwrites the default name (ipopt.opt) of the options file
  // app->Options()->SetStringValue("option_file_name", "IA.opt");
  
  */
  
  // Intialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Ipopt::Solve_Succeeded) {
    if (!silent)
      printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }
  
  IAIntWaveNlp *myianlp = new IAIntWaveNlp(iaData, ipData, iaSolution, silent);
  Ipopt::SmartPtr<Ipopt::TNLP> mynlp = myianlp; // Ipopt requires the use of smartptrs!

  bool try_again = true;
  int iter = 0;
  
  // print();
  bool solution_ok = false;
  
  do {
    if (debugging)
    {
      print();
      printf("%d IntWave iteration\n", iter );
      // build the hessian, evaluate it and f at the current solution?
    }
      
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp); // the inherited IANlp
    
    // see /CoinIpopt/build/include/coin/IpReturnCodes_inc.h
    /*
    Solve_Succeeded=0,
    Solved_To_Acceptable_Level=1,
    Infeasible_Problem_Detected=2,
    Search_Direction_Becomes_Too_Small=3,
    Diverging_Iterates=4,
    User_Requested_Stop=5,
    Feasible_Point_Found=6,
    
    Maximum_Iterations_Exceeded=-1,
    Restoration_Failed=-2,
    Error_In_Step_Computation=-3,
    Maximum_CpuTime_Exceeded=-4,
    Not_Enough_Degrees_Of_Freedom=-10,
    Invalid_Problem_Definition=-11,
    Invalid_Option=-12,
    Invalid_Number_Detected=-13,
    
    Unrecoverable_Exception=-100,
    NonIpopt_Exception_Thrown=-101,
    Insufficient_Memory=-102,
    Internal_Error=-199
     */

    bool solved_full = false;
    bool solved_partial = false;
    bool solver_failed = false;
    bool bad_problem = false;

    switch (status) {
      case Ipopt::Solve_Succeeded:
      case Ipopt::Solved_To_Acceptable_Level:
      case Ipopt::Feasible_Point_Found:
        solved_full = true;
        break;
      case Ipopt::Maximum_Iterations_Exceeded:
      case Ipopt::User_Requested_Stop:
      case Ipopt::Maximum_CpuTime_Exceeded:
        solved_partial = true;
        break;
      case Ipopt::Infeasible_Problem_Detected:
      case Ipopt::Not_Enough_Degrees_Of_Freedom:
      case Ipopt::Invalid_Problem_Definition:
      case Ipopt::Invalid_Option:
      case Ipopt::Invalid_Number_Detected:
        bad_problem = true;
        break;
      case Ipopt::Search_Direction_Becomes_Too_Small:
      case Ipopt::Restoration_Failed:
      case Ipopt::Diverging_Iterates:
      case Ipopt::Error_In_Step_Computation:
      case Ipopt::Unrecoverable_Exception:
      case Ipopt::NonIpopt_Exception_Thrown:
      case Ipopt::Insufficient_Memory:
      case Ipopt::Internal_Error:        
        solver_failed = true;
        break;
        
      default:
        break;
    }
  
    if (!silent)
    {
      if (solved_full) {
        printf("\n\n*** IntWave solved!\n");
      }
      else {
        printf("\n\n*** IntWave FAILED!\n");
      }
    }
    
    if (debugging)
    {
      printf("\nChecking solution.\n");
      bool integer_sat = solution_is_integer(true);
      bool even_sat = even_constraints( false, true);
      bool equal_sat = equal_constraints( false, true );
      printf("IntWave solution summary, %s, equal-constraints %s, even-constraints %s.\n", 
             integer_sat ? "integer" : "NON-INTEGER",
             equal_sat ? "satisfied" : "VIOLATED", 
             even_sat ? "satisfied" : "VIOLATED" );
      if (!integer_sat)
        printf("investigate integer neighborhood\n");
      if (!even_sat)
        printf("investigate even neighborhood\n");
      if (!equal_sat)
        printf("investigate equal neighborhood\n");
    }
    
    
    IASolution nlp_solution;
    nlp_solution.x_solution = ia_solution()->x_solution; // vector copy
    IASolverToolInt sti( ia_data(), &nlp_solution );
    sti.round_solution();
    if (debugging)
      printf("Checking rounded solution, ignoring even constraints.\n");
    if (sti.equal_constraints(false, debugging))
    {
      // also even constraints
      if (debugging)
        printf("Rounding worked.\n");
      
      // rounding was a valid integer solution
      ia_solution()->x_solution.swap( nlp_solution.x_solution );
      // ia_solution()->obj_value is no longer accurate, as it was for the non-rounded solution
      return true;
    }
    
    // todo: detect and act
    // may have converged to a locally optimal, but non-feasible solution
    // if so, try a new starting point
    
    // check solution feasibility, even when not debugging
    
    if ( solved_full || solved_partial )
    {
      bool integer_sat = solution_is_integer(false);
      bool even_sat = even_constraints( false, false);
      bool equal_sat = equal_constraints( false, false );
      if ( integer_sat && even_sat && equal_sat )
        return true;
    }

    // find out which vars were not integer, 
    // try moving to a farther starting point resolving
    
 
    try_again = false; 
  } while (try_again);
  
  
  // now 
  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.  
  return solution_ok;
  
}


bool IASolverInt::solve_minlp()
{
  // set up and call the separate IARoundingNlp, which has a linear objective to get a natural integer solution 
  // the intuition is this will solve integrality for  most variables all at once

  if (debugging)
  {
    printf("IASolverInt::solve_minlp() - ");        
    printf("Attempting to find a naturally-integer solution to an NLP.\n");
  }

  
  // solver setup  
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  
  // convergence parameters
  // see $IPOPTDIR/Ipopt/src/Interfaces/IpIpoptApplication.cpp
  // our real criteria are: all integer, constraints satisfied. How to test the "all_integer" part?
  app->Options()->SetNumericValue("tol", 1e-6); //"converged" if NLP error<this, default is 1e-7. Obj are scaled to be >1, so e-2 is plenty // was 1e-2
  app->Options()->SetNumericValue("max_cpu_time", sqrt( iaData->num_variables() ) ); // max time allowed in seconds
  app->Options()->SetIntegerValue("max_iter", 3 * iaData->num_variables() ); // max number of iterations
  // app->Options()->SetNumericValue("primal_inf_tol", 1e-2 ); 
  app->Options()->SetNumericValue("dual_inf_tol", 1e-6 ); // how close to infeasibility? // was 1e-2
  app->Options()->SetNumericValue("constr_viol_tol", 1e-6 ); // by how much can constraints be violated?
  app->Options()->SetNumericValue("compl_inf_tol", 1e-6 ); // max norm of complementary condition // was 1e-2

  // second criteria convergence parameters: quit if within this tol for many iterations
// was  app->Options()->SetIntegerValue("acceptable_iter", 4 + sqrt( iaData->num_variables() ) ); //as "tol"
  app->Options()->SetNumericValue("acceptable_tol", 1e-6 ); //as "tol" was 1e-1

  app->Options()->SetStringValue("mu_strategy", "adaptive");
  // print level 0 to 12, Ipopt default is 5
  const int print_level = (silent) ? 0 : 1; 
  app->Options()->SetIntegerValue("print_level", print_level);  
  // uncomment next line to write the solution to an output file
  // app->Options()->SetStringValue("output_file", "IA.out");  
  // The following overwrites the default name (ipopt.opt) of the options file
  // app->Options()->SetStringValue("option_file_name", "IA.opt");
  
  // Intialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Ipopt::Solve_Succeeded) {
    if (!silent)
      printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }
  
  // problem setup
  // a couple of different models, simplest to more complex
  IARoundingNlp *myianlp = new IARoundingNlp(iaData, ipData, iaSolution, silent);
  // IARoundingFarNlp *myianlp = new IARoundingFarNlp(iaData, ipData, this);
  // IARoundingFar3StepNlp *myianlp = new IARoundingFar3StepNlp(iaData, ipData, this); // haven't tested this. It compiles and runs but perhaps isn't correct
  // IAIntWaveNlp *myianlp = new IAIntWaveNlp(iaData, ipData, this); // haven't tested this. It compiles and runs but perhaps isn't correct
  Ipopt::SmartPtr<Ipopt::TNLP> mynlp = myianlp; // Ipopt requires the use of smartptrs!

  bool try_again = true;
  int iter = 0;
  do {
    printf("%d rounding iteration\n", iter );
    
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp); // the inherited IANlp
    
    if (!silent)
    {
      if (status == Ipopt::Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
      }
      else {
        printf("\n\n*** The problem FAILED!\n");
      }
    }
    
    // The problem should have been feasible, but it is possible that it had no integer solution.
    // figure out which variables are still integer
    
    // check solution for integrality and constraint satified    
    if (debugging)
    {
      printf("\nChecking Natural (non-rounded) solution.\n");
      bool integer_sat = solution_is_integer(true);
      bool even_sat = even_constraints( false, true);
            bool equal_sat = equal_constraints( false, true );
      printf("Natural solution summary, %s, equal-constraints %s, even-constraints %s.\n", 
             integer_sat ? "integer" : "NON-INTEGER",
             equal_sat ? "satisfied" : "VIOLATED", 
             even_sat ? "satisfied" : "VIOLATED" );
      if (!integer_sat)
        printf("investigate this\n");
    }
    
    IASolution nlp_solution;
    nlp_solution.x_solution = ia_solution()->x_solution; // vector copy
    IASolverToolInt sti( ia_data(), &nlp_solution );
    sti.round_solution();
    if (debugging)
      printf("Checking rounded solution, ignoring even constraints.\n");
    if (sti.equal_constraints(false, debugging))
    {
      // also even constraints
      if (debugging)
        printf("Rounding worked.\n");

      // rounding was a valid integer solution
      ia_solution()->x_solution.swap( nlp_solution.x_solution );
      // ia_solution()->obj_value is no longer accurate, as it was for the non-rounded solution
      return true;
    }

    // find out which vars were not integer, 
    // try rounding their weights and resolving
    // bool int_sat = solution_is_integer();
    ++iter;
    try_again = iter < 4 + sqrt(iaData->num_variables());
    if (try_again)
    {
      if (debugging)
        printf("rounding failed, randomizing weights\n");
    
      myianlp->randomize_weights_of_non_int(); // try again? debug
    }
    else if (debugging)
      printf("giving up on rounding to non-integer solution\n");

    // try_again = false; // debug
  } while (try_again);

  
  // todo: update partially-integer solution, perhaps using ipData - figure out how we're going to use it, first, for what structure makes sense.
  
  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.  
  return status == Ipopt::Solve_Succeeded;
  
}

void IASolverInt::cleanup()
{
  ;
}

bool IASolverInt::solve()
{
  
  // debug, try solve_intwave instead 
  // longer term, use intwave as a backup when the faster and simpler milp doesn't work.
  // unfortunately, it appears to find local minima that are far from optimal, even when starting in a well
  // return solve_intwave();
  
  // todo: add a simple rounding term for the even-constraints to minlp  
  solve_minlp();
  bool success = solution_is_integer();
  if (success)
  {
    if (!silent)
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

} // namespace MeshKit
