// main file for stubbed version of interval assignment
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <cmath>

#include "IAInterface.hpp"  // from stubs
#include "IAVariable.hpp" // from stubs

MeshKit::IAInterface *new_ia_interface()
{
  return new MeshKit::IAInterface;
}

void delete_ia_interface(MeshKit::IAInterface * iainterface)
{
  delete iainterface;
}


bool check_solution_correctness( MeshKit::IAInterface *ia_interface, 
                                 std::vector< std::pair<int,int> > &correct_solution)
{
  const bool verbose_output = true;
  const bool debug = true;
  bool all_correct = true;
  MeshKit::IAInterface::VariableVec::const_iterator b = ia_interface->variables_begin();
  MeshKit::IAInterface::VariableVec::const_iterator e = ia_interface->variables_end();
  MeshKit::IAInterface::VariableVec::const_iterator i = b;
  unsigned int c = 0;
  if (debug)
    std::cout << "Checking Solution Correctness" << std::endl;
  for ( ; i != e; ++i, ++c )
  {
    const MeshKit::IAVariable *v = *i;
    assert(v);
    const int x = v->get_solution();
    assert(c < correct_solution.size() );
    const int lo = correct_solution[c].first;
    const int hi = correct_solution[c].second;
    if (debug)
      std::cout << "Checking variable " << c << " solution " << x << " in " 
              << "[" << lo << "," << hi << "]?" << std::endl;
    if (x < lo)
    {
      if (verbose_output)
        std::cout << "ERROR: Variable " << c << " solution " << x << " BELOW " 
             << "[" << lo << "," << hi << "]" << std::endl;
      all_correct = false;
    }
    if (x > hi)
    {
      if (verbose_output)
        std::cout << "ERROR: Variable " << c << " solution " << x << " ABOVE " 
             << "[" << lo << "," << hi << "]" << std::endl;
      all_correct = false;
    }
  }
  if (debug)
    std::cout << "done checking solution correctness." << std::endl;
  return all_correct;
}

void set_decoupled_pairs(MeshKit::IAInterface *ia_interface, 
				         int num_pairs, double goal1, double goal2,
				         std::vector< std::pair<int,int> > &correct_solution)
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
    const double g1 = i + goal1;
    const double g2 = i + goal2;
    MeshKit::IAVariable *v1 = ia_interface->create_variable( NULL, MeshKit::SOFT, g1);
    MeshKit::IAVariable *v2 = ia_interface->create_variable( NULL, MeshKit::SOFT, g2);
    const double compromise = sqrt( g1 * g2 );
    double lo = floor(compromise); 
    if ( ( compromise - lo ) < 0.1 )
      --lo;
    if ( lo < 1. )
      lo = 1.;
    double hi = ceil(compromise);
    if ( (hi - compromise) < 0.1 )
      ++hi;
    correct_solution.push_back( std::make_pair( lo, hi ) );
    correct_solution.push_back( std::make_pair( lo, hi ) );

    // constrain x_{2i} - x_{2i+1} = 0
    MeshKit::IAInterface::IAVariableVec side1, side2;
    side1.push_back(v1);
    side2.push_back(v2);
    ia_interface->constrain_sum_equal(side1, side2);
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
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();

 
  std::vector< std::pair<int,int> > correct_solution;
  set_decoupled_pairs(ia_interface, 1, 1, 3, correct_solution);
//  set_decoupled_pairs(ia_interface, 1, 3.2, 12.1, correct_solution);
  ia_interface->execute_this(); 
  bool solution_correct = check_solution_correctness( ia_interface, correct_solution );
  assert( solution_correct );
}

void test_many_pairs()
{
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();

  std::vector< std::pair<int,int> > correct_solution;
  set_decoupled_pairs(ia_interface, 8, 3.2, 12.1, correct_solution);
  set_decoupled_pairs(ia_interface, 1, 3.2, 12.1, correct_solution);
  set_decoupled_pairs(ia_interface, 8, 7.7, 4.2, correct_solution);
  set_decoupled_pairs(ia_interface, 40, 1.1, 5.2, correct_solution);
  set_decoupled_pairs(ia_interface, 40, 1.6, 4.5, correct_solution);
  set_decoupled_pairs(ia_interface, 4, 1.5, 1.5, correct_solution);
  set_decoupled_pairs(ia_interface, 4, 1, 1, correct_solution);
  
  ia_interface->execute_this(); 
  
  bool solution_correct = check_solution_correctness( ia_interface, correct_solution );
  assert( solution_correct );
  
  delete_ia_interface( ia_interface );
}

void mapping_test() 
{
  
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();

  MeshKit::ModelEnt* me_curves[4] = {0,0,0,0};
  // stubbed
  me_curves[0] = new MeshKit::ModelEnt();
  me_curves[1] = new MeshKit::ModelEnt();
  me_curves[2] = new MeshKit::ModelEnt();
  me_curves[3] = new MeshKit::ModelEnt();
  // convert moab entity handles to ModelEnt* and place in me_curves... ask Tim how
  
  MeshKit::MEntVector side1, side2;
  side1.push_back(me_curves[0]); side2.push_back(me_curves[2]);
  ia_interface->constrain_sum_equal(ia_interface->make_constraint_group(side1), 
                                    ia_interface->make_constraint_group(side2));
  side1.clear(); side2.clear();
  side1.push_back(me_curves[1]); side2.push_back(me_curves[3]);
  ia_interface->constrain_sum_equal(ia_interface->make_constraint_group(side1), 
                                    ia_interface->make_constraint_group(side2));

  // if there are loops, and the loops have strictly less than 4 curves, then
  // ia_interface->constrain_sum_even( ia_interface->make_constraint_group(curves in loop) );
  
  delete_ia_interface( ia_interface );
}

int main(int argv, char* argc[])
{
  // stubbed  
  test_one_pair();
  test_many_pairs();
  mapping_test();

  return 0;
}
