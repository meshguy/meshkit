// IASolverRelaxed.cpp
// Interval Assignment for Meshkit
//
#include "IASolverRelaxed.hpp"
#include "IANlp.hpp"

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "IpIpoptApplication.hpp"

namespace MeshKit {
    
IASolverRelaxed::IASolverRelaxed(const IAData *ia_data, IASolution *relaxed_solution,   
                                 const bool set_silent) : 
  IASolverTool( ia_data, relaxed_solution, true ),  // debug
  p_norm(3), 
//silent(set_silent), debugging(false) {}
silent(false), debugging(true) {}

/** default destructor */
IASolverRelaxed::~IASolverRelaxed() {}

bool IASolverRelaxed::solve()
{
  
  // solve the nlp to get a non-integral solution, which we hope is close to a good integer solution    
  // adapted from HS071 ipopt example

  // p_norm set in constructor. 3 seems to work well, comes close to lex-max-min
  // smaller p has the effect of valuing the fidelity of shorter curves over longer curves more
  // larger p approaches min max
  IANlp *myianlp = new IANlp(iaData, iaSolution, silent);
  Ipopt::SmartPtr<TNLP> mynlp = myianlp; // Ipopt requires the use of smartptrs!
  
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 1e-7); // 2 seems close enough, could do less, say .1
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  // print level 0 to 12, most. Ipopt Default is 5
  int print_level = (silent) ? 0 : 1; // 1, 5
  // int print_level = 5;
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
      printf("\n\n*** Error during ipopt initialization!\n");
    return (int) status;
  }
  
  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp); // the inherited IANlp
  // todo: also check for a valid solution even if ! Solve_Succeeded, such as a sub-optimal time-out
  bool is_solved = (status == Ipopt::Solve_Succeeded);
  bool is_satisfied = is_solved && equal_constraints( false, debugging );
  // don't check even-ness, as those are like the integrality constraints and are not solved here
  
  if (!silent)
  {
    if (is_solved) {
      printf("\n\n*** The relaxed problem solved!");
      if (!is_satisfied)
        printf(" But equality-constraints were VIOLATED!");
      printf("\n");
    }
    else {
      printf("\n\n*** The relaxed problem FAILED!\n");
    }
  }
  return is_satisfied;
}  

} // namespace MeshKit
