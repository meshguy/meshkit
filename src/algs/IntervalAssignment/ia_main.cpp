// ia_main.cpp
// test IntervalMatching interface IAInterface for MeshKit

#include "TestUtil.hpp"
#include "meshkit/MKCore.hpp"

#include "meshkit/IAInterface.hpp"
#include "meshkit/IAVariable.hpp"

#include <stdio.h>
#include <iostream>

void set_decoupled_pairs(MeshKit::IAInterface &ia_interface, int num_pairs, int goal1, int goal2)
{
  // trivial 2-sided mapping problem
  // we can make multiple pairs, each pair is independent, 
  // and pair i (in 0..num_pairs-1) has sides with one curve each with goals
  // i+goal1 and i+goal2, 
  //
  // test scalability, relaxed nlp: 100,000 constraints in 1 second. milp: 40 variables in 1 second, grows exponentially!
  for (int i = 0; i<num_pairs; ++i)
  { 
    // goals x_{2i} = 2, x_{2i+1} = 2
    // x_{2i}, goal: i + goal1
    MeshKit::IAVariable *v1 = ia_interface.create_variable( NULL, MeshKit::IAVariable::SOFT, i + goal1);
    MeshKit::IAVariable *v2 = ia_interface.create_variable( NULL, MeshKit::IAVariable::SOFT, i + goal2);

    // constrain x_{2i} - x_{2i+1} = 0
    MeshKit::IAInterface::IAVariableVec side1, side2;
    side1.push_back(v1);
    side2.push_back(v2);
    ia_interface.constrain_sum_equal(side1, side2);
  }
}

/*
void IASolver::set_test_problem()
{
*/
  
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
  }
 */
  

void test_one_pair()
{
  MeshKit::MKCore *mk = new MeshKit::MKCore();
  MeshKit::IAInterface *ia_interface =
    (MeshKit::IAInterface*) mk->construct_meshop("IntervalAssignment");

  set_decoupled_pairs(*ia_interface, 1, 1, 3);
  ia_interface->execute_this(); 
  //todo: also check that the solution value for each variable is 2
  // CHECK(success); 
}

int main(int argv, char* argc[])
{
  int one_pair = RUN_TEST(test_one_pair);
//  int abrt = RUN_TEST(test_abort);
//  int expt = RUN_TEST(test_exception);
//  int succ = RUN_TEST(test_success);
  
  if (one_pair) 
    std::cerr << "IntervalAssignment failed for one_pair!" << std::endl;
  
  int success = one_pair; // + !abrt + !expt + succ;
  return success;
}
