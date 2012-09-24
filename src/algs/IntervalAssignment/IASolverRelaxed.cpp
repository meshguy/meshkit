// IASolverRelaxed.cpp
// Interval Assignment for Meshkit
//
#include "IASolverRelaxed.hpp"
#include "IANlp.hpp"

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "IpIpoptApplication.hpp"

IASolverRelaxed::IASolverRelaxed(const IAData *ia_data, IASolution *relaxed_solution) : 
iaData(ia_data), relaxedSolution(relaxed_solution), p_norm(3), debugging(true) {}

/** default destructor */
IASolverRelaxed::~IASolverRelaxed() { iaData = NULL; relaxedSolution=NULL;}


double IASolverRelaxed::check_constraint(const int i, const IAData * ia_data, const IASolution *solution, bool sum_even, bool print_me )
{
  double g_i = 0.;
  if (!sum_even)
  {
    if (solution)
    {
      for (int j = 0; j < ia_data->constraints[i].M.size(); ++j)
      {
        g_i += solution->x_solution[ ia_data->constraints[i].M[j].col ] * ia_data->constraints[i].M[j].val;
      }
    }
    if (print_me)
    {
      printf("constraint %d:", i);
      for (int j = 0; j < ia_data->constraints[i].M.size(); ++j)
      {
        printf( " %f x_%d", ia_data->constraints[i].M[j].val, ia_data->constraints[i].M[j].col );
        if (solution)
          printf(" (%f)", 
                 solution->x_solution[ ia_data->constraints[i].M[j].col ] * ia_data->constraints[i].M[j].val );
      }
      printf("\n");
    }
  }
  else
  {
    if (solution)
    {
      for (int j = 0; j < ia_data->sumEvenConstraints[i].M.size(); ++j)
      {
        g_i += solution->x_solution[ ia_data->sumEvenConstraints[i].M[j].col ] * ia_data->sumEvenConstraints[i].M[j].val;
      }
      if (print_me)
      {
        printf("sum-constraint %d:", i);
        for (int j = 0; j < ia_data->sumEvenConstraints[i].M.size(); ++j)
        {
          printf( " %f x_%d", ia_data->sumEvenConstraints[i].M[j].val, ia_data->sumEvenConstraints[i].M[j].col );
          if (solution) 
            printf(" (%f) ", 
                   solution->x_solution[ ia_data->sumEvenConstraints[i].M[j].col ] * ia_data->sumEvenConstraints[i].M[j].val );
        }
        printf("\n");
      }
    }
  }
  
  return g_i; 
}


bool IASolverRelaxed::constraints_satisfied( const IAData * ia_data, const IASolution *solution, bool &sum_even_satisfied, bool print_me )
{
  // compare to 
  // bool IANlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  // does the same thing with the data in ipopt format

  bool satisfied = true;
  for (int i = 0; i<ia_data->constraints.size(); ++i)
  {
    double g_i = check_constraint(i, ia_data, solution, false, false);
    if (g_i + 1.e-2 < ia_data->constraints[i].lowerBound || g_i - 1.e-2 > ia_data->constraints[i].upperBound)
    {
      satisfied = false;
      if (print_me)
      {
        printf( "VIOLATED ");
        check_constraint(i, ia_data, solution, false, true);
      }
    }
  }
  
  sum_even_satisfied = true;
  for (Index i = 0; i<ia_data->sumEvenConstraints.size(); ++i)
  {
    double g_k = check_constraint(i, ia_data, solution, true, false);
    int g = floor(g_k + 1.e-2);
    if (g % 2)
    {    
      sum_even_satisfied = false;
      if (print_me)
      {
        printf( "VIOLATED ");
        printf("  2k <> %d = %f = ", g, g_k);
        check_constraint(i, ia_data, solution, true, true);
      }
    }
  }
  
  return satisfied; 
}


bool IASolverRelaxed::solve()
{
  using namespace Ipopt;
  
  // solve the nlp to get a non-integral solution, which we hope is close to a good integer solution    
  // adapted from HS071 ipopt example

  // p_norm set in constructor. 3 seems to work well, comes close to lex-max-min
  // smaller p has the effect of valuing the fidelity of shorter curves over longer curves more
  // larger p approaches min max
  IANlp *myianlp = new IANlp(iaData, relaxedSolution);
  SmartPtr<TNLP> mynlp = myianlp; // Ipopt requires the use of smartptrs!
  
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
    printf("\n\n*** Error during ipopt initialization!\n");
    return (int) status;
  }
  
  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp); // the inherited IANlp
  
  if (status == Solve_Succeeded) {
    printf("\n\n*** The relaxed problem solved!\n");
  }
  else {
    printf("\n\n*** The relaxed problem FAILED!\n");
    return false;
  }
  return true;
}  

