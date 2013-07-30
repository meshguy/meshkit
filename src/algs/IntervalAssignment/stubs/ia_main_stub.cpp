

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
  const bool debug = false;
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


void set_mapping_chain( MeshKit::IAInterface *ia_interface, const int num_sides, 
                       const bool grow_goal_by_i,
                       const int goal_m1, const int goal_m2, 
                       const int num_curve_min, const int num_curve_max )
{
  // test problem 3, sides with more than one variable, with random goals
  printf("constructing coupled test problem - mapping chain\n");
//  srand(10234); // for scaling by curves
// srand(6893498); // for scaling by faces when the other results in a bend
  MeshKit::IAInterface::IAVariableVec side1, side2;
  int num_vars = 0;
  for (int i = 0; i<num_sides; ++i)
  { 
    // move side2 to side1
    side2.swap( side1 );
    
    // create new side2
    side2.clear();
    assert( num_curve_min > 0 );
    int num_curves = num_curve_min;
    if ( num_curve_max > num_curve_min )
      num_curves += (rand() % (1 + num_curve_max - num_curve_min) ); 
    for (int j = 0; j < num_curves; j++)
    {
      int goal_intervals = (1 + (rand() % goal_m1)) * (1 + (rand() % goal_m2)); 
      if (grow_goal_by_i)
        goal_intervals += num_vars;
      MeshKit::IAVariable *v = ia_interface->create_variable( NULL, MeshKit::SOFT, goal_intervals);
      side2.push_back(v);
    }
    num_vars += num_curves;
    
    // if we have two non-trivial opposite sides, then constrain them to be equal
    if (side1.size() && side2.size())
    {
      ia_interface->constrain_sum_equal(side1, side2);
    }
    
    // add a sum-even constraint
    if (0 && i==0)
    {
      printf("sum-even side: %d", i);
      assert( side2.size() );
      ia_interface->constrain_sum_even(side2);
    }
    
    // todo: try some hard-sets and non-trivial rhs
  }
  printf("problem size: %d variables, %d equal constraints\n", num_vars, num_sides-1);

}

void set_mapping_face( MeshKit::IAInterface *ia_interface )
{
  // test problem 3, sides with more than one variable, with random goals
  printf("constructing mapping face with several curves per side\n");
  srand(10234);
  MeshKit::IAInterface::IAVariableVec side1, side2; 

  int num_vars[] = {1, 2};
  int goals[] = {10, 100, 50};

  int num_curves = num_vars[0];
  for (int j = 0; j < num_curves; j++)
  {
    int goal_intervals = goals[j];
    MeshKit::IAVariable *v = ia_interface->create_variable( NULL, MeshKit::SOFT, goal_intervals);
    side1.push_back(v);
  }

  num_curves = num_vars[1];
  for (int j = 0; j < num_curves; j++)
  {
    int goal_intervals = goals[num_vars[0] + j];
    MeshKit::IAVariable *v = ia_interface->create_variable( NULL, MeshKit::SOFT, goal_intervals);
    side2.push_back(v);
  }

  
    // if we have two non-trivial opposite sides, then constrain them to be equal
    if (side1.size() && side2.size())
    {
      ia_interface->constrain_sum_equal(side1, side2);
    }
    
}


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

void set_mapping_side( MeshKit::IAInterface *ia_interface, 
                      unsigned int num_curves, MeshKit::IAInterface::IAVariableVec &side, std::vector<double> &goals )
{
  // see set_mapping_face below for usage  
  int num_new_curves = (int) num_curves - (int) side.size(); 
  if (num_new_curves < 0)
    num_new_curves = 0;
  for (int j = 0; j < num_new_curves; j++)
  {
    double goal_intervals = 1.; // default
    if (j < (int) goals.size())
      goal_intervals = goals[j]; // use the passed in goals
    else
      goals.push_back(goal_intervals); // pass back the assigned goals
    MeshKit::IAVariable *v = ia_interface->create_variable( NULL, MeshKit::SOFT, goal_intervals);
    side.push_back(v);
  }
}


void set_mapping_face( MeshKit::IAInterface *ia_interface, 
                      unsigned int num_curves_1, MeshKit::IAInterface::IAVariableVec &side_1, std::vector<double> &goals_1,
                      unsigned int num_curves_2, MeshKit::IAInterface::IAVariableVec &side_2, std::vector<double> &goals_2) 
{
  // create a mapping face, with num_curves_? one each side.
  // two uses, independent for each side:
  // 1. Use the passed-in goals to create and fill in the ia_variables (so they can be used in the other form)
  // 2. Use the passed-in ia_variables, so to get two faces sharing a side or part of a side
  // a third use is a blend: re-use some variables and create some new ones. 
  // Goals default to 1 if unspecified or not enough are passed in.
  
//  printf("set mapping face with %u and %u opposite curves.\n", num_curves_1, num_curves_2);
//  printf("side_1: re-using %lu variables and creating %lu new ones\n", side_1.size(), num_curves_1 - side_1.size());
//  printf("side_2: re-using %lu variables and creating %lu new ones\n", side_2.size(), num_curves_2 - side_2.size());
  
  // MeshKit::IAInterface::IAVariableVec side1, side2; 
  
  set_mapping_side( ia_interface, num_curves_1, side_1, goals_1 );
  set_mapping_side( ia_interface, num_curves_2, side_2, goals_2 );

  // if we have two non-trivial opposite sides, then constrain them to be equal
  if (side_1.size() && side_2.size())
  {
    ia_interface->constrain_sum_equal(side_1, side_2);
  }
  
}


void set_half_integer(  MeshKit::IAInterface *ia_interface, int version)
{
  MeshKit::IAInterface::IAVariableVec side_1, side_2;
  std::vector<double> goals_1, goals_2;

  switch (version) {
    // loop, pointed
    case 0:
    {
      const double g1 = 107; // the one-curve side
      const double g2 =  53.5; // the multi-curve looping side
      // gives 1/2-integer solution
      
      goals_1.push_back(g1);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      set_mapping_face(ia_interface, 1, side_1, goals_1, 2, side_2, goals_2);

      MeshKit::IAInterface::IAVariableVec first_side_1 = side_1;
      MeshKit::IAVariable *common_0 = side_2[0];
      MeshKit::IAVariable *common_1 = side_2[1];
   
      side_1.clear();
      side_2.clear();
      side_1.push_back(common_0);
      goals_2.clear();
      goals_2.push_back(g2);
      set_mapping_face(ia_interface, 1, side_1, goals_1, 1, side_2, goals_2);
      
      side_1.clear();
      goals_1.clear();
      goals_1.push_back(g2);
      set_mapping_face(ia_interface, 1, side_1, goals_1, 1, side_2, goals_2);

      side_2.clear();
      side_2.push_back(common_1);
      set_mapping_face(ia_interface, 1, side_1, goals_1, 1, side_2, goals_2);
      
      // need to have many replicates of the single side, in order to have its integer weight overcome the opposite half-integer weights
      // MeshKit::IAVariable *weight_me = side_1[0]; 
      for (int r = 0; r < 7; r++)
      {
        side_2.clear();
        goals_2.clear();
        goals_2.push_back(g1);
        set_mapping_face(ia_interface, 1, first_side_1, goals_1, 1, side_2, goals_2);
      }

    }
      break;

    // loop, flat with extra curve for slack
    case 1:
    {
      const double g1 = 107; // the one-curve side
      const double g2 =  52.5; // the multi-curve looping side
      const double g3 =  2; // the multi-curve looping side

      goals_1.push_back(g1);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      goals_2.push_back(g3); // the slack curve
      set_mapping_face(ia_interface, 1, side_1, goals_1, 3, side_2, goals_2);
      
      MeshKit::IAInterface::IAVariableVec first_side_1 = side_1;
      MeshKit::IAVariable *common_0 = side_2[0];
      MeshKit::IAVariable *common_1 = side_2[1];
      
      side_1.clear();
      side_2.clear();
      side_1.push_back(common_0);
      goals_2.clear();
      goals_2.push_back(g2); // 8, for a different test problem testing flattening leading to an unbounded solution
      set_mapping_face(ia_interface, 1, side_1, goals_1, 1, side_2, goals_2);
      
      side_1.clear();
      goals_1.clear();
      goals_1.push_back(g2); // 8
      set_mapping_face(ia_interface, 1, side_1, goals_1, 1, side_2, goals_2);
      
      side_2.clear();
      side_2.push_back(common_1);
      set_mapping_face(ia_interface, 1, side_1, goals_1, 1, side_2, goals_2);
      
      for (int r = 0; r < 7; r++)
      {
        side_2.clear();
        goals_2.clear();
        goals_2.push_back(g1);
        set_mapping_face(ia_interface, 1, first_side_1, goals_1, 1, side_2, goals_2);
      }

    }
      break;
      
      // loop, pointed, but two curves per side around the loop
    case 2:
    {
      const double g1 = 107; // the one-curve side
      const double g2 =  26.75; // the multi-curve looping side 
      // above set to 26.3, 70 replicates makes some 26 and some 27
      // using 700 replicates with 26.75 gives a 1/2-integer (26.6) solution

      goals_1.push_back(g1);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      set_mapping_face(ia_interface, 1, side_1, goals_1, 4, side_2, goals_2);
      
      MeshKit::IAInterface::IAVariableVec first_side_1 = side_1;
      MeshKit::IAVariable *common_0  = side_2[0];
      MeshKit::IAVariable *common_0a = side_2[1];
      MeshKit::IAVariable *common_1  = side_2[2];
      MeshKit::IAVariable *common_1a = side_2[3];
      
      side_1.clear();
      side_1.push_back(common_0);
      side_1.push_back(common_0a);
      side_2.clear();
      goals_2.clear();
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      set_mapping_face(ia_interface, 2, side_1, goals_1, 2, side_2, goals_2);
      
      side_1.clear();
      goals_1.clear();
      goals_1.push_back(g2);
      goals_1.push_back(g2);
      set_mapping_face(ia_interface, 2, side_1, goals_1, 2, side_2, goals_2);
      
      side_2.clear();
      side_2.push_back(common_1);
      side_2.push_back(common_1a);
      set_mapping_face(ia_interface, 2, side_1, goals_1, 2, side_2, goals_2);
      
      for (int r = 0; r < 700; r++)
      {
        side_2.clear();
        goals_2.clear();
        goals_2.push_back(g1);
        set_mapping_face(ia_interface, 1, first_side_1, goals_1, 1, side_2, goals_2);
      }

    }
      break;

      // loop, flat with extra curve for slack, and two curves per side around the loop
    case 3:
    {
      const double g1 = 107; // the one-curve side
      const double g2 =  26.5; // the multi-curve looping side 
      const double g3 =  2; // the multi-curve looping side 

      goals_1.push_back(g1);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      goals_2.push_back(g3); // the slack curve
      set_mapping_face(ia_interface, 1, side_1, goals_1, 5, side_2, goals_2);
      // todo, fiddle with weights to get 1/2-integer solution
      
      MeshKit::IAInterface::IAVariableVec first_side_1 = side_1;
      MeshKit::IAVariable *common_0  = side_2[0];
      MeshKit::IAVariable *common_0a = side_2[1];
      MeshKit::IAVariable *common_1  = side_2[2];
      MeshKit::IAVariable *common_1a = side_2[3];
      
      side_1.clear();
      side_1.push_back(common_0);
      side_1.push_back(common_0a);
      side_2.clear();
      goals_2.clear();
      goals_2.push_back(g2);
      goals_2.push_back(g2);
      set_mapping_face(ia_interface, 2, side_1, goals_1, 2, side_2, goals_2);
      
      side_1.clear();
      goals_1.clear();
      goals_1.push_back(g2);
      goals_1.push_back(g2);
      set_mapping_face(ia_interface, 2, side_1, goals_1, 2, side_2, goals_2);
      
      side_2.clear();
      side_2.push_back(common_1);
      side_2.push_back(common_1a);
      set_mapping_face(ia_interface, 2, side_1, goals_1, 2, side_2, goals_2);

      side_1 = first_side_1;
      for (int r = 0; r < 0; r++) // no replicates needed to get 1/2 integer solution
      {
        side_2.clear();
        goals_2.clear();
        goals_2.push_back(g1);
        set_mapping_face(ia_interface, 1, first_side_1, goals_1, 1, side_2, goals_2);
      }

    }
      break;

    default:
      break;
  }
}

void test_half_integer()
{
  // create interface
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  
  for (int v=0; v<4; ++v) // 0 <= v < 4
  {
    ia_interface->destroy_data();
  
    // set up model
    set_half_integer(ia_interface, v);
    printf("test problem %d\n", v);
  
    // solve ia
    ia_interface->execute_this(); 
  }
  
  delete_ia_interface( ia_interface );
}

void test_scaling_by_curves()
{
 
  int num_tests = 6;
  int num_curves = 1000;
  double factor = sqrt(10); // 10./3.;

  for (int i = 0; i < num_tests; ++i)
  {
    MeshKit::IAInterface *ia_interface = new_ia_interface();
    ia_interface->destroy_data();

//    void set_mapping_chain( MeshKit::IAInterface *ia_interface, const int num_sides, 
//                           const bool grow_goal_by_i,
//                           const int goal_m1, const int goal_m2, 
//                           const int num_curve_min, const int num_curve_max )

    set_mapping_chain( ia_interface, 2, false, 3, 15, num_curves, num_curves);
    ia_interface->execute_this(); 
    
    delete_ia_interface( ia_interface );
    num_curves *= factor;
  }
  
}

void test_map_skew()
{
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();
  
  std::vector< std::pair<int,int> > correct_solution;
  set_mapping_face(ia_interface);
  //  set_decoupled_pairs(ia_interface, 1, 3.2, 12.1, correct_solution);
  ia_interface->execute_this(); 
}

void test_one_pair()
{
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();

  std::vector< std::pair<int,int> > correct_solution;
  set_decoupled_pairs(ia_interface, 1, 10., 1000., correct_solution);
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

void test_long_chain()
{
  
  // test scalability: 20000 gives 20,000 constraints, 100,000 variables in 1 second relaxed solution
//  set_mapping_chain(ia_interface, 16000, false, 3, 15, 2, 11);
  // IANLP paper: 160, 160*3.3333, 1600, 1600*3.3333, 
  int num_curves = 160;
  int num_factors = 5; // 0, 1, 2, 3, 4, 5, 6
  int srandseed[] = {10234, 6893498, 10234, 102, 102, 72346};
  //                 0       1       2      3      4      5
  for (int i = 0; i < num_factors; ++i)
  {
    MeshKit::IAInterface *ia_interface = new_ia_interface();
    ia_interface->destroy_data();

    //  srand(10234); // for scaling by curves
    // srand(6893498); // for scaling by faces when the other results in a bend
    srand( srandseed[i] );
    std::cout << " i = " << i << " seed = " << srandseed[i] << std::endl;

    set_mapping_chain(ia_interface, num_curves, false, 3, 15, 2, 11);
    // goal distribution is gaussian in [1, 32]

    ia_interface->execute_this(); 
  
    // bool solution_defined = check_solution( ia_interface );

    delete_ia_interface( ia_interface );
    
    num_curves *= sqrt(10.);
  }
}


void test_growing_chain()
{
  // test problem 2
  // printf("constructing growing chain, coupled test problem\n");
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();
  
  // goals are 1, 2, 3, 4, ... 16
  // one curve per side
  set_mapping_chain(ia_interface, 16, true, 1, 1, 1, 1);

  ia_interface->execute_this(); 
  
  // bool solution_defined = check_solution( ia_interface );
  
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
  
  ia_interface->execute_this(); 
  
  delete_ia_interface( ia_interface );
}


void paving_test() 
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
  for (int i = 0; i < 4; ++i)
  {
    me_curves[i]->meshIntervalSize = 1.; 
    me_curves[i]->entMeasure = 6. + (i==0 ? 1 : 0); //force odd
    me_curves[i]->meshIntervals = me_curves[i]->entMeasure / me_curves[i]->meshIntervalSize;
  }
  
  MeshKit::MEntVector loop1;
  loop1.push_back(me_curves[0]); 
  loop1.push_back(me_curves[1]);
  loop1.push_back(me_curves[2]); 
  loop1.push_back(me_curves[3]);
  ia_interface->constrain_sum_even(ia_interface->make_constraint_group(loop1));
  
  ia_interface->execute_this(); 
  
  delete_ia_interface( ia_interface );
}

int main(int argv, char* argc[])
{
  // stubbed  
  
  // 2013 solution for implicit non-one coefficients 
  test_half_integer();
  
  // 2013 solution for explicit sum-even's
  
  // For NLIA paper
//  test_map_skew();
//  test_scaling_by_curves();
//     test_long_chain(); // ma27 has memory issues over 50,000 constraints

//  test_one_pair();
  // test_many_pairs();
//  test_growing_chain();
  
  //paving_test();

/* 100 long chain test
  srand(9384757); // debug
  for (int i = 0; i < 100; ++i)
    test_long_chain();
*/
  
 // mapping_test();

  return 0;
}
