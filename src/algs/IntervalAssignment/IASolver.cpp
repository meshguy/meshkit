// IASolver.cpp
// Interval Assignment for Meshkit
//
#include "IASolver.hpp"
#include "IPData.hpp"
// #include "IAMilp.hpp"
#include "IASolverRelaxed.hpp"
#include "IASolverInt.hpp"

#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <limits.h>

namespace MeshKit 
{
    
IASolver::IASolver() : debugging(false) {}

/** default destructor */
IASolver::~IASolver() {}

void IASolver::print_solution() // IPData &ip_data
{
  printf("\nIA solution:\n");
  for (int i=0; i<num_variables(); ++i)
  {
    printf("%d x (goal %e) ",i,I[i]);
    //printf(" relaxed %e solution ", ip_data.relaxedSolution[i]);
    int x_int = floor( x_solution[i] + 0.5);
    if (fabs( x_int - x_solution[i])<1.0e-2)
      printf("%d\n",x_int);
    else
      printf("%e NON-INTEGER\n",x_solution[i]);
  }
  for (unsigned int i=0; i< sumEvenConstraints.size(); ++i)
  {
    double d = IASolver::sum_even_value(i, this, this);
    printf("%d sum-even = %d (%e)",i, (int) floor(d+0.5), d);
    if (!IASolver::is_even(d))
      printf(" NON-EVEN\n");
    else
      printf("\n");
  }
  printf("objective function value %e\n", obj_value);
  printf("\n");
}

bool IASolver::solve_relaxed()
{
  IASolverRelaxed relaxed(this,this,!debugging); 
  bool succeeded=relaxed.solve();
  return succeeded;
}  


/* milp solution is too slow
bool IASolver::solve_int()
{
  IAMilp milp(this,this);
  bool succeeded=milp.solve();
  return succeeded;
}  
*/

bool IASolver::solve_int()
{
  // set relaxed solution to this's solution (from solve_relaxed)
  IASolverInt solver_int(this,this,!debugging);
  bool succeeded=solver_int.solve();
  // replace the solution with the int solution
  x_solution = solver_int.x_solution; // vector copy
  return succeeded;
}  

bool IASolver::solve_even()
{
  bool succeeded=true;
  return succeeded;
}  


bool IASolver::solve()
{
  
  // todo: subdivide problem into independent sub-problems for speed
  
  bool relaxed_succeeded = solve_relaxed();
  bool int_succeeded = relaxed_succeeded && solve_int();
  bool even_succeeded = int_succeeded && solve_even();
  
  if (debugging)
  {
    printf("==========IA summary:\n"); 
    if (relaxed_succeeded)
    {
      printf(" relaxed succeeded\n");
      if (int_succeeded)
      {
        printf(" integer succeeded\n");
        if (even_succeeded)
        {
          printf(" even succeeded\n");
        }
        else
        {
          printf(" even failed\n");        
        }
      }
      else
      {
        printf(" integer failed\n");        
      }
    }
    else
    {
      printf(" relaxed failed\n");
    }  
    print_solution();
  }
	  
	  
	  
  return even_succeeded;
}

double IASolver::sum_even_value(int i, const IAData *ia_data, const IASolution *current_solution)
{
  // number of non-zeros (coefficients) in the constraint
  const int nnz = (int) ia_data->sumEvenConstraints[i].M.size();
  double sum = - ia_data->sumEvenConstraints[i].rhs;
  for (int j = 0; j<nnz; ++j)
  {
    // if x is near integer, force it to be exactly integer to avoid an accumulation of roundoff that 
    // would cause us to not recognize that we already have an (even) integer solution.
    double x = current_solution->x_solution[ ia_data->sumEvenConstraints[i].M[j].col ];
    if (fabs(x - floor( x + 0.5 )) < 1.e-2)
      x = floor(x + 0.5);
    sum += x * ia_data->sumEvenConstraints[i].M[j].val;
  }
  return sum;
}

double IASolver::sum_even_value(int i)
{
  return sum_even_value(i, this, this);
}


bool IASolver::is_even(double y)
{
  int e = floor( y + 0.5 );
  if (e % 2) 
    return false;
  if ( fabs( y - e ) < 1.0e-4 )
    return true;
  return false;
}

} // namespace MeshKit
