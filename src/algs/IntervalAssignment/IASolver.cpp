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
    
IASolver::IASolver() : debugging(true) {}

/** default destructor */
IASolver::~IASolver() {}

void IASolver::set_test_problem()
{
    // test problem 0 
  // trivial 2-sided mapping problem
  printf("constructing decoupled test problem.");
  int num_pairs = 1; // each pair of sides is decoupled from all other pairs.
  // test scalability, relaxed nlp: 100,000 constraints in 1 second. milp: 40 variables in 1 second, grows exponentially!
  for (int i = 0; i<num_pairs; ++i)
  { 
    // goals x_{2i} = 2, x_{2i+1} = 2
    I.push_back( i + 1); // x_{2i}
    I.push_back( i + 3); // x_{2i+1}
    
    // constrain x0 - x1 = 0
    constraints.push_back( constraintRow() );
    constraints.back().upperBound = constraints.back().lowerBound = 0.; // equality
    constraints.back().M.push_back( sparseEntry(2*i,1.) );
    constraints.back().M.push_back( sparseEntry(2*i+1,-1.) );
  }
   
  
  /*
   // test problem 1
   printf("constructing decoupled test problem.");
   int num_sides = 160; // test scalability, 100,000 constraints in 1 second
   for (int i = 0; i<num_sides; ++i)
   { 
     // goals x_{2i} = 1, x_{2i+1} = 4
     // 1 and 4 result in a natural integer solution of 2, which MILP verifies very quickly
     I.push_back(1.1);
     I.push_back(12);
   
     // constrain x0 - x1 = 0
     constraints.push_back( constraintRow() );
     constraints.back().upperBound = constraints.back().lowerBound = 0.; // equality
     constraints.back().M.push_back( sparseEntry(2*i,1.) );
     constraints.back().M.push_back( sparseEntry(2*i+1,-1.) );
   }
    */
  
  
  /*
   // test problem 2
   printf("constructing coupled test problem - long chain");
   int num_sides = 16; // test scalability: results is 100,000 constraints in 1 second, for relaxed nlp
   for (int i = 0; i<num_sides; ++i)
   { 
     // goals x_{i} = i+1;
     I.push_back(i+1);
   
     // constrain x0 - x1 = 0
     if (i>0)
     {
       constraints.push_back( constraintRow() );
       constraints.back().upperBound = constraints.back().lowerBound = 0.; // equality
       constraints.back().M.push_back( sparseEntry(i-1,  1.) );
       constraints.back().M.push_back( sparseEntry(i  , -1.) );
     }
   } 
   */
  
    /*  
  // test problem 3, sides with more than one variable, with random goals
  printf("constructing coupled test problem - long chain\n");
  srand(10234);
  int num_sides = 16000; // test scalability: 20000 gives 20,000 constraints, 100,000 variables in 1 second relaxed solution
  int curve_id =0;
  std::vector<int> test_sum_even;
  for (int i = 0; i<num_sides; ++i)
  { 
    bool is_prior_constraint = constraints.size() > 0;
    // add a new forward constraint
    bool is_current_constraint = (i<num_sides-1);
    if (is_current_constraint)
    {
      constraints.push_back( IAData::constraintRow() );
      constraints.back().upperBound = constraints.back().lowerBound = 0.; // equality
    }
    std::vector<IAData::constraintRow>::iterator current_constraint = constraints.end();
    std::vector<IAData::constraintRow>::iterator prior_constraint = constraints.end();
    if (is_current_constraint)
      current_constraint = constraints.end() - 1;
    if (is_prior_constraint) 
      prior_constraint = current_constraint -1;
    
    int num_curves = 2 + (rand() % 9); // 
    if (i==0)
    {
      printf("sum-even side: ");
    }
    for (int j = 0; j < num_curves; j++)
    {
      // add a sum-even constraint for the first side
      if (i==0)
      {
        test_sum_even.push_back(curve_id);
        printf(" %d", curve_id);
      }
      if (is_prior_constraint)
        prior_constraint->M.push_back( IAData::sparseEntry(curve_id, -1.) );
      if (is_current_constraint)
        current_constraint->M.push_back( IAData::sparseEntry(curve_id, 1.) );
      int goal_intervals = (1 + (rand() % 3)) * (1 + (rand() % 15)); // 1 -- 32, was 4, 8
      I.push_back(goal_intervals);
      ++curve_id;
    }
    if (i==0)
    {
      constrain_sum_even(test_sum_even,1);
      test_sum_even.clear();
      printf("\n");
    }
  }
  //add also some sum-even constraints
  */  
  
  // sum-even constraints test problems
/*
  // test problem 4, a simple sum-even constraint
  int num_surfaces = 12; // 12
  int num_curves_per_surface = 4; // 4
  int num_shared_curves = 1; // 2
  
  int num_curves = 0;
  for (int i = 0; i < num_surfaces; ++i)
  {
    // gather the indices for the sum-even constraint
    int start_curve = num_curves - num_shared_curves;
    if (start_curve < 0)
      start_curve = 0;
    std::vector<int>curve_indices;
    if (debugging)
      printf("%d sum-even:",i);
    for (int j = 0; j < num_curves_per_surface; ++j)
    {
      curve_indices.push_back(start_curve+j);
      if (debugging)
        printf(" %d",start_curve+j);
    }
    num_curves = start_curve + num_curves_per_surface;
    const int rhs = 0; // test 0, -1
    constrain_sum_even(curve_indices,rhs);
    if (debugging)
      printf(" =%d\n",rhs);
  }
  // assign random goals to the curves
  for (int i = (int) I.size(); i < num_curves; ++i )
  {
    double goal = 1 + ((double) (rand() % 59)) / 10.; // 1 to 6.9
    // force an odd sum for testing purposes
    //if (i==0)
    //  goal += 1.;
    I.push_back(goal);
  }
 */
  
} 

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
  IASolverRelaxed relaxed(this,this);
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
  IASolverInt solver_int(this,this);
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
