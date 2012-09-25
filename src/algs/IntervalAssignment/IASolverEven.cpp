// IASolverEven.cpp
// Interval Assignment for Meshkit
//
#include "IASolverEven.hpp"
#include "IAMINlp.hpp"
#include "IPData.hpp"
// #include "IAMilp.hpp"

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "IpIpoptApplication.hpp"

IASolverEven::IASolverEven(const IAData *ia_data, const IASolution *int_solution) :
iaData(ia_data), debugging(true) 
{ 
  ipData.initialize(int_solution->x_solution);
}

/** default destructor */
IASolverEven::~IASolverEven() { iaData=NULL; }


double IASolverEven::sum_even_value(int i, const IAData *ia_data, const IASolution *current_solution)
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

double IASolverEven::sum_even_value(int i)
{
  return sum_even_value(i, iaData, this);
}


bool IASolverEven::is_even(double y)
{
  int e = floor( y + 0.5 );
  if (e % 2) 
    return false;
  if ( fabs( y - e ) < 1.0e-4 )
    return true;
  return false;
}

double IASolverEven::distance_to_even(int i)
{
  double sum = sum_even_value(i);
  if (is_even(sum))
    return 0.;
  return sum;
  /*
  double lower_even = 2. * floor( sum / 2.);
  if (lower_even < 4.) 
    lower_even = 4.;
  const double upper_even = 2. * ceil( sum / 2. );
  const double u = upper_even - sum;
  const double l = sum - lower_even;
  if ( u > l )
    return l;
  else
    return u;  
   */
}

bool IASolverEven::calculate_sumeven_value(const int i, double &obj_increase, int &y_int )
{
  // to do: make a version of this where x can flip which side of the goal it is on.
  // the issue is e.g. goal = 1.1, x gets rounded to 1, then x can never flip to the other side of 1.1! but if goal was 1, then it could
  const double x = x_solution[i];
  y_int = floor(x + 0.5); // round to nearest integer - this version assumes we are already integer within roundoff.
  // round x away from I
  const double goal = iaData->I[i]; // shorthand reminder
  if (x>goal)
    ++y_int;
  else 
    --y_int;
  if (y_int >= 1)
  {
    // compare to relaxed solution, not the latest pseudo-integer solution
    double f = myianlp->eval_R_i( goal, ipData.relaxedSolution[i] ); 
    double g = myianlp->eval_R_i( goal, y_int );
    obj_increase = g - f;           
    // assert( obj_increase > 0. ); this is only invariant if we use f = R(x), not f = R(relaxed).
    return true;
  }
  //else
  obj_increase = 0.;
  return false;

}

bool IASolverEven::calculate_rounding_value(const int i, double &obj_increase, int &y_int )
{
  const double x = x_solution[i];
  y_int = floor(x + 0.5); // round to nearest integer
  double y = y_int;
  if ( fabs(x-y) > 1e-4) // not already integer within tolerance 
  {
    // non-integer found
    // round x away from I
    const double goal = iaData->I[i];
    // it was forced away from its goal, round in the forced direction further
    // else leave it rounded to the nearest integer
    if (fabs(x - goal) > 1.e-2) 
    {
      y_int = (x>goal) ? ceil(x) : floor(x);
      y = y_int;
    }
    assert(y >= 1.);
    
    // if it is already constrained to y_int, then consider it as already satisfying it.
    if ( ipData.varIntegerBound[i] == y )
    {
      obj_increase = 0.;
      return false;
    }

    // todo: if x is very close to the goal, maybe we should retain freedom to round it the other direction, try both, etc.
 
    // is this the best non-integer we've found so far?
    // best means the one that changing it to integer increases the obj function the least
    // double f = myianlp->eval_R_i( goal, x ); 
    // compare to relaxed solution, not the latest pseudo-integer solution
    double f = myianlp->eval_R_i( goal, ipData.relaxedSolution[i] ); 
    double g = myianlp->eval_R_i( goal, y );
    obj_increase = g - f;           
    // assert( obj_increase > 0. ); this is only invariant if we use f = R(x), not f = R(relaxed).
    return true;
  }
  //else
  obj_increase = 0.;
  return false;
}

bool IASolverEven::find_one_non_integer(int &i_nonint, int &x_bound)
{
  i_nonint = -1; // none found
  double obj_increase_nonint( std::numeric_limits<double>::max() );
  // find the one variable whose rounding causes the least change to the objective function
  // to do: alt is to pick one at random, with some chance
  for (int i = 0; i<iaData->num_variables(); ++i)
  {
    double obj_increase;
    int y_int;
    if (calculate_rounding_value(i, obj_increase, y_int))
    {
      if (obj_increase < obj_increase_nonint)
      {
        i_nonint = i;
        x_bound = y_int;
        obj_increase_nonint = obj_increase;
      }
    }
  }
  if (i_nonint > -1)
    return true;
  // everybody is already integer, enforce some sum-even constraint
  // pick an un-satisfied sum-even constraint
  int num_noneven_found = 0;
  // algorithm parameter
  // sum_strategy = 0 pick biggest, 1: uniformly at random 
  // to do: test these further. Are there bad cases where randomness helps? Some random sequences give poor-quality solutions, rounding a 1.5 up to 3 or 4.
  const int sum_strategy = 1; // rand() % 2; // 0
  printf("sum_strategy=%d ", sum_strategy);
  double biggest_non_even_value = 0.;
  int biggest_non_even_i = -1;
  for (unsigned int i = 0; i < iaData->sumEvenConstraints.size(); ++i)
  {
    const double d = distance_to_even(i);
    if (d > 0.)
    {
      ++num_noneven_found;
      if (sum_strategy == 0 && d > biggest_non_even_value || 
          sum_strategy == 1 && (rand() % num_noneven_found == 0))
      {
        biggest_non_even_value = d;
        biggest_non_even_i = i;
      }
    }
  }
  // now pick an curve of the sum-even-constraint to round
  if (biggest_non_even_i > -1)
  {
    const int c = biggest_non_even_i;
    // algorithm parameter
    // 0: pick biggest non-even , 1: uniformly at random 
    // to do: test these further
    const int curve_strategy = 1; // rand() % 2; // 0
    printf("curve_strategy=%d\n", curve_strategy);
    int num_found = 0;
    double smallest_increase = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < iaData->sumEvenConstraints[c].M.size(); ++i)
    {
      int j = iaData->sumEvenConstraints[c].M[i].col;
      assert( iaData->sumEvenConstraints[c].M[i].val == 1. );
      const double x = x_solution[j];
      if (!is_even(x))
      {
        ++num_found;
        int y_int;
        double obj_increase;
        if (calculate_sumeven_value(j, obj_increase, y_int))
        {
          if (curve_strategy == 0 && obj_increase < smallest_increase || 
              curve_strategy == 1 && (rand() % num_found == 0))
          {
            i_nonint = j;
            x_bound = y_int;
            smallest_increase = obj_increase;
          }
        }
      }
    }
    // if no non-even curves were found, pick an odd one
    if (i_nonint == -1)
    {
      for (unsigned int i = 0; i < iaData->sumEvenConstraints[c].M.size(); ++i)
      {
        int j = iaData->sumEvenConstraints[c].M[i].col;
        assert( iaData->sumEvenConstraints[c].M[i].val == 1. );
        ++num_found;
        int y_int;
        double obj_increase;
        if (calculate_sumeven_value(j, obj_increase, y_int))
        {
          if (curve_strategy == 0 && obj_increase < smallest_increase || 
              curve_strategy == 1 && (rand() % num_found == 0))
          {
            i_nonint = j;
            x_bound = y_int;
            smallest_increase = obj_increase;
          }
        }
      }
    }
    if (i_nonint == -1)
      // failed to achieve sum-even solution
      return false;
    return true;
  }
  return false;
}
 
void IASolverEven::constrain_integer(const int i_nonint, const int x_bound)
{
  ipData.oldBound[i_nonint] = ipData.varIntegerBound[i_nonint];;
  ipData.varIntegerBound[i_nonint] = (double) x_bound;
}

// find the one variable whose rounding causes the least change to the objective function
// unused at the moment
bool IASolverEven::constrain_one_non_integer( int &i, int &b)
{
  if (find_one_non_integer(i, b))
  {
    constrain_integer(i, b);
    if (debugging)
    {
      report_one_constraint(i, b);
    }
    return true;
  }
  return false;
}

/*
if (i>num_variables())
{
  assert(phase >= 1);
  const int j = i - num_variables();
  printf("  constraining sum[%d] from %e to %c %d\n", j, sum_even_value(j), '>', b);
}
else

 
 if (i>num_variables())
 {
 assert(phase >= 1);
 const int j = i - num_variables();
 printf("  new sum[%d] value is %e %c %d\n", j, sum_even_value(j), '>', b);
 }
 else
 {
*/

void IASolverEven::report_one_constraint(const int i, const int b)
{
  // assumes x_solution[i] is current non-integer solution
  printf("  constraining x[%d] from %e (%e) to %c %d\n", i, x_solution[i], iaData->I[i], b > iaData->I[i] ? '>' : '<', b);
}

void IASolverEven::report_one_non_integer(const int i, const int b)
{
  printf("  new x[%d] value is ", i );
  if (fabs( floor(x_solution[i]+0.5) - x_solution[i] ) < 1.0e-2)
    printf("%d", (int) floor(x_solution[i]+0.5));
  else
    printf("%e", x_solution[i]);
  printf(" (%e) %c %d\n", iaData->I[i], b > iaData->I[i] ? '>' : '<', b);
}

void IASolverEven::back_off(RoundingMap &rounding_map)
{
  // undo half of the constraints
  // todo: 
  RoundingMap::iterator undo_start = rounding_map.begin(); 
  int num_to_keep = (int) rounding_map.size() / 2; // roundoff, rounding to zero is OK
  std::advance( undo_start, num_to_keep );
  for(RoundingMap::const_iterator i = undo_start; i != rounding_map.end(); ++i)
  {
    const int x_i = i->second.first;
    if (debugging)
    {
      printf(" backing off:");
      report_one_constraint(x_i, ipData.oldBound[x_i]);
    }
    ipData.varIntegerBound[x_i] =  ipData.oldBound[x_i]; 
  }
  rounding_map.erase(undo_start, rounding_map.end()); 
}

bool IASolverEven::find_many_non_integer(RoundingMap &rounding_map)
{
  // to do: this heuristic is bad because it easily paints us into a corner 
  // to do: we should instead pick some constant fraction of the variable per *side* of a constraint, 
  // provided there is at least one other unconstrained variable on that side. And never constrain the last k-sides (k=2, or 3? or 4?) until the end, when we just constrain one at a time.
  
  // since we don't know how many are non-integer, put them all in the map, then discard the ones we won't constrain
  rounding_map.clear();
  for (int i = 0; i<iaData->num_variables(); ++i)
  {
    double obj_increase;
    int y_int;
    if (calculate_rounding_value(i, obj_increase, y_int))
    {
      // add to map
      // printf("inserting obj %e, i %d, x %e, y %d\n", obj_increase, i, x_solution[i], y_int);
      rounding_map.insert(std::make_pair(obj_increase,std::make_pair(i,y_int)));  
    }
  }
  if (rounding_map.empty())    
    return false;
  //const unsigned long num_to_round = 1;
  const unsigned long num_to_round = 1; // 1 + rounding_map.size() / 8; // todo: experiment with the 16 part
  if (debugging)
  {
    printf("Found %ld of %d non-integers, rounding %ld of them.",
           rounding_map.size(), iaData->num_variables(), num_to_round);
  }
  // remove everything after the "num_to_round" element
  RoundingMap::iterator r = rounding_map.begin();
  advance(r, num_to_round);
  rounding_map.erase(r, rounding_map.end()); 
  // todo: if erasing (or subsequent iterations) is slow due to tree rebalancing, then instead we should transfer the first num_to_round elements to a std::vector and return that instead, and discard the map

  if (debugging)
  {
    if (rounding_map.size() > 1)
      printf(" obj_increase range %e to %e\n",
             rounding_map.begin()->first, rounding_map.rbegin()->first);
    else if (rounding_map.size() == 1)
      printf( " obj_increase expected %e\n",
             rounding_map.begin()->first);
    else
      printf( " rounding nothing\n");
  }
  return true;
}

bool IASolverEven::constrain_many_non_integer(RoundingMap &rounding_map)
{
  if (find_many_non_integer(rounding_map))
  {
    for(RoundingMap::const_iterator i = rounding_map.begin(); i != rounding_map.end(); ++i)
    {
      const int x_i = i->second.first;
      const int b = i->second.second;
      constrain_integer(x_i, b); // index, new integer bound
      if (debugging)
      {
        report_one_constraint(x_i, b);
      }
    }
    // todo: heuristic to decrease the number of constrained variables if the solver returns infeasible, or it decreased these variables beyond the new upper bounds . save the old bounds to restore them.
    return true;
  }
  return false; 
}

void IASolverEven::report_many_non_integer(RoundingMap &rounding_map)
{
  for(RoundingMap::const_iterator i = rounding_map.begin(); i != rounding_map.end(); ++i)
  {
    report_one_non_integer(i->second.first, i->second.second);
  }
}

void IASolverEven::cleanup()
{
  myianlp = NULL; // smart pointers should take care of deleting this.
}

bool IASolverEven::solve()
{
  using namespace Ipopt;
  
  // solve the nlp to get a non-integral solution, which we hope is close to a good integer solution    
  // =======================================
  // from HS071 ipopt example

  // smaller p has the effect of valuing the fidelity of shorter curves over longer curves more
  // larger p approaches min max
  myianlp = new IAMINlp(iaData, &ipData, this);
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
    printf("\n\n*** Error during initialization!\n");
    cleanup();
    return (int) status;
  }
  
  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp); // the inherited IANlp
  
  if (status == Solve_Succeeded) {
    printf("\n\n*** The relaxed problem solved!\n");
  }
  else {
    printf("\n\n*** The relaxed problem FAILED!\n");
    cleanup();
    return false;
  }
  
  // make integer
  // Pick the variable whose rounding to the next further integer from its goal results in the lowest individual contribution (or increase in contribution from its current value) to the objective function.
  // Round it that way, setting a lower bound on its value (but leaving objective, etc., the same)
  // Re-solve the NLP. Some old vars might have their delta directions or values flip. That is OK.
  // Iterate.

  ipData.initialize(x_solution);
  int iter = 0;
  printf("\n\n*** Solving integer solution.\n"); 

    
  /* */
  // many at a time speedup, could skip and still get correct output
  if (debugging)
    printf("Many-at-a-time rounding\n");
  RoundingMap rounding_map;
  // phase = -1;
  while (constrain_many_non_integer(rounding_map))
  /* */
  {
    
    if (debugging)
      printf("MINLP many-at-once iter %d\n", iter);
    else
      printf("%d ", iter);
    
    // solve for many new constraints, but reduce those constraints if the problem becomes infeasible
    bool try_again = false;
    do 
    {
      try_again = false;
      status = app->OptimizeTNLP(mynlp); // double-check that resolving works, w/out initialize, etc 
      if (status != Solve_Succeeded && status != Solved_To_Acceptable_Level)
      {
        back_off(rounding_map);
        try_again = rounding_map.size();
      }
    }
    while (try_again);
    if (rounding_map.size() == 0) // no progress was made
      break;
      
    // to do: what to do if status != Solve_Succeeded ?
    // assert(b >= I[i] && x_solution[i] >= b || b < I[i] && x_solution[i] <= b); // assert the new constraint was satisfied, one-at-a-time
    if (debugging)
    {
      // report_one_non_integer(i, b);
      report_many_non_integer(rounding_map);
    }
    ++iter;
    // todo: quit at an iteration limit? but a single variable could be increased several times...
  }

  
  // one at a time
  if (debugging)
    printf("\nOne-at-a-time rounding, including sum-even constraints.\n");
  int i; // non-integer variable to make integer
  int b; // value to make it integer
  iter = 0;

  while (constrain_one_non_integer(i,b))
  {
      
    if (debugging)
      printf("MINLP one-at-a-time iter %d\n", iter);
    else
      printf("%d ", iter);
      
    status = app->OptimizeTNLP(mynlp); // double-check that resolving works, w/out initialize, etc 
    // to do: what to do if status != Solve_Succeeded ?
    const double g = iaData->I[i];
    const double x = x_solution[i];
    assert(b >= g && x + 1.e-2 >= b || b < g && x <= b + 1.e-2); // assert the new constraint was satisfied, one-at-a-time
    if ((status != Solve_Succeeded && status != Solved_To_Acceptable_Level) || 
        ! (b >= g && x + 1.e-2 >= b || b < g && x <= b + 1.e-2))
    {
      // see ipopt/CoinIpopt/Ipopt/src/Interfaces/IpReturnCodes_inc.h for the enum of possible values
      // status = Infeasible_Problem_Detected; 
      if (status == Infeasible_Problem_Detected )
        printf("IPOPT Infeasible_Problem_Detected\n");
      else
        printf("IPOPT bad solution return code %d\n", status);
      break;
    }
    if (debugging)
    {
      report_one_non_integer(i, b);
    }
    ++iter;
    // todo: quit at an iteration limit? but a single variable could be increased several times...
  }
  
  if (status != Solve_Succeeded)
    printf("***Integer Fail!***\n");
  else
  {
    if (debugging)
    {
      // print_solution();
    }
  }
  
  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted. 
  cleanup();
  return status == Solve_Succeeded;
  
}

