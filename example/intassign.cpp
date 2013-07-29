/*!
\example intassign.cpp

\section intassign_cpp_title <pretty-name-of-this-file>

\subsection intassign_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection intassign_cpp_goal Goal

\subsection intassign_cpp_cw Code Walkthrough

\subsection intassign_cpp_in Input
\image html intassign.in.jpg
There is no input.

\subsection intassign_cpp_out Output
\image html intassign.out.jpg

\subsection intassign_cpp_src Source Code
*/
/*! \file test_intassign.cpp \test
 *
 * test_intassign.cpp
 *
 */

// ia_main.cpp 
// test IntervalMatching interface IAInterface for MeshKit


#include "meshkit/MKCore.hpp"

#include "meshkit/IAInterface.hpp"
#include "meshkit/IAVariable.hpp"
#include "meshkit/TFIMapping.hpp"

#include <stdio.h>
#include <iostream>
#ifdef HAVE_ACIS
#define INTASSIGN_TEST_FILE "quadface.sat"
#elif defined(HAVE_OCC)
#define INTASSIGN_TEST_FILE "quadface.stp"
#endif

MeshKit::MKCore *mk;

MeshKit::IAInterface *new_ia_interface()
{
  return
      (MeshKit::IAInterface*) mk->construct_meshop("IntervalAssignment");
}

void delete_ia_interface(MeshKit::IAInterface *)
{
  //nada, mk takes care of it
  ;
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
  srand(10234);
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

      // if we have two non-trivial opposite sides, then constrain them to be equal
      if (side1.size() && side2.size())
        {
          ia_interface->constrain_sum_equal(side1, side2);
        }

      // add a sum-even constraint
      if (i==0)
        {
          // printf("sum-even side: %d", i);
          assert( side2.size() );
          ia_interface->constrain_sum_even(side2);
        }

      // todo: try some hard-sets and non-trivial rhs
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

void test_one_pair()
{
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();

  std::vector< std::pair<int,int> > correct_solution;
  set_decoupled_pairs(ia_interface, 1, 1, 3, correct_solution);
  //  set_decoupled_pairs(ia_interface, 1, 3.2, 12.1, correct_solution);
  ia_interface->execute_this();
  bool solution_correct = check_solution_correctness( ia_interface, correct_solution );
  CHECK( solution_correct );
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
  CHECK( solution_correct );
  
  delete_ia_interface( ia_interface );
}

void test_long_chain()
{
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();
  
  // test scalability: 20000 gives 20,000 constraints, 100,000 variables in 1 second relaxed solution
  set_mapping_chain(ia_interface, 16000, false, 3, 15, 2, 11);
  // goal distribution is gaussian in [1, 32]

  ia_interface->execute_this();
  
  // bool solution_defined = check_solution( ia_interface );

  delete_ia_interface( ia_interface );
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

  std::string file_name = TestDir + "/" + INTASSIGN_TEST_FILE;
  printf("opening %s\n", file_name.c_str());
  mk->load_geometry_mesh(file_name.c_str(), NULL);

  //check the number of geometrical edges
  MeshKit::MEntVector surfs, loops;
  mk->get_entities_by_dimension(2, surfs);
  MeshKit::ModelEnt *this_surf = (*surfs.rbegin());

  // request a specific size
  //mk->sizing_function(0.1, true);
  MeshKit::MEntVector curves;
  std::vector<int> senses, loop_sizes;
  this_surf->boundary(1, curves, &senses, &loop_sizes);
  CHECK_EQUAL(4, (int)curves.size());

  MeshKit::SizingFunction esize(mk, -1, -1);
  surfs[0]->sizing_function_index(esize.core_index());

  MeshKit::MEntVector side1, side2;
  side1.push_back(curves[0]); side2.push_back(curves[2]);
  ia_interface->constrain_sum_equal(ia_interface->make_constraint_group(side1),
                                    ia_interface->make_constraint_group(side2));
  side1.clear(); side2.clear();
  side1.push_back(curves[1]); side2.push_back(curves[3]);
  ia_interface->constrain_sum_equal(ia_interface->make_constraint_group(side1),
                                    ia_interface->make_constraint_group(side2));

  // if there are loops, and the loops have strictly less than 4 curves, then
  // ia_interface->constrain_sum_even( ia_interface->make_constraint_group(curves in loop) );
  
  //now, do the TFIMapping
  (MeshKit::TFIMapping*) mk->construct_meshop("TFIMapping", surfs);
  mk->setup_and_execute();
  mk->save_mesh("intassign.exo");
  //
  delete_ia_interface( ia_interface );
}

int main(int argv, char* argc[])
{
  // currently unable to create more than one mk called IntervalAssignment
  mk = new MeshKit::MKCore();

  test_one_pair();
  // run same test twice to check data clearing integrity
  test_one_pair();
  test_many_pairs();
  mapping_test();
  test_growing_chain();
  test_long_chain();

  test_exception();
  test_success();

  delete mk;
  int success = one_pair + one_pair2 + many_pairs + map_test + growing_chain; // + !abrt + !expt + succ;
  return success;
}
