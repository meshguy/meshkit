// IASolverTool.cpp
// Interval Assignment for Meshkit
//

#include "IASolverTool.hpp"
#include "IAData.hpp"
#include "IASolution.hpp"
#include "IPData.hpp" 

#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <limits.h>

namespace MeshKit 
{

IASolverTool::~IASolverTool()
{
  iaData = NULL; 
  iaSolution=NULL;
}
  
bool IASolverTool::is_integer(const double x, double &x_int_double) const
{
  x_int_double = ( floor( x + 0.5 ) );
  // beware cases where this epsilon is too large, e.g. sides with more than 100 curves = 1 / 1.e-2
  return fabs( x - x_int_double ) < 1.e-2; 
}

bool IASolverTool::is_integer(const double x, int &x_int) const
{
  double x_int_double;
  bool is_int = is_integer( x, x_int_double );
  x_int = (int) ( x_int_double );
  return is_int; 
}

  
bool IASolverTool::is_integer(const double x) const
{
  double x_int_double;
  return is_integer( x, x_int_double );
}

  
bool IASolverTool::valid_solution() const
{
  return (iaSolution && iaData && (iaSolution->x_solution.size() >= iaData->I.size())); // could be bigger if deltas are retained
}
  
    
void IASolverTool::print(const bool do_print_solution, const bool do_print_nonint,
                         const bool do_print_equal_constraints, const bool do_print_nonequal,
                         const bool do_print_even_constraints, const bool do_print_noneven ) const
{
  
  int info_case;
  /*  
   four cases
   0 nothing
   1 data_and_solution
   2 data_only
   3 solution_only 
   */
  if (iaData)
  {
    if (valid_solution())
      info_case = 1; //data_and_solution = true;
    else
      info_case = 2; // data_only = true;    
  }
  else
  {
    if (iaSolution && iaSolution->x_solution.size())
      info_case = 3; // solution_only = true;
    else
      info_case = 0; // nothing = true;
  }
  
  // header
  switch (info_case) {
    case 0:
      printf("no data, no solution\n");
      return;
      break;

    case 1:
      printf("\nIA data and solution:\n");
      break;
      
    case 2:
      printf("\nIA data:\n");
      break;

    case 3:
      printf("\nIA solution:\n");
      break;

    default:
      break;
  }

  // variables
  if ( info_case == 1 || info_case == 2 ) // data exists
  {
    printf("%d vars\n", iaData->num_variables());
    for (int i=0; i<iaData->num_variables(); ++i)
    {
      printf("%d x (goal %e) ",i, iaData->I[i]);
      if (do_print_solution && info_case == 1)
      {
        //printf(" relaxed %e solution ", ip_data.relaxedSolution[i]);
        const double x = iaSolution->x_solution[i];
        int x_int;
        if (is_integer(x,x_int))
          printf("%d\n",x_int);
        else
          printf("%e  NON-INTEGER\n",x);
      }
      else
      {
        printf("\n");
      }
    }
  }
  else if ( info_case == 3 && do_print_solution)
  {
    printf("%d x solution values\n", iaData->num_variables());
    for (unsigned int i=0; i<iaSolution->x_solution.size(); ++i)
    {
      printf("x_%d ",i);
      //printf(" relaxed %e solution ", ip_data.relaxedSolution[i]);
      const double x = iaSolution->x_solution[i];
      int x_int;
      if (is_integer(x,x_int))
        printf("%d\n",x_int);
      else
        printf("%e  NON-INTEGER\n",x);
    }
  }
  
  if ( info_case == 1 || info_case == 2 ) // data exists
  {
    if (do_print_equal_constraints || do_print_nonequal)
    {
      printf("%lu equality constraints:\n", iaData->constraints.size());
      equal_constraints(do_print_equal_constraints, do_print_nonequal );
    }
    if (do_print_even_constraints || do_print_noneven)
    {
      printf("%lu even constraints:\n", iaData->sumEvenConstraints.size());
      even_constraints(do_print_even_constraints, do_print_noneven);
    }
  }
  
  if (do_print_solution && valid_solution())
    printf("objective function value %e\n", iaSolution->obj_value);
  printf("\n");
}
  
void IASolverTool::print_solution() const
{
  print( true, false, false, false, false, false );
}
  
void IASolverTool::print_problem() const
{
  print( false, true, false, false, false, false );
}
  
bool IASolverTool::is_even(double y, int &y_even) const
{
  // nearest even value, including 0 and 2
  y_even = (int) floor( 2. * floor( y / 2. + 0.5 ) + 1.e-2 );
  if ( fabs( y - y_even ) < 1.0e-2 )
    return true;
  return false;
}
  
bool IASolverTool::is_even(double y) const
{
  int y_even;
  return is_even(y, y_even);
}
  
double IASolverTool::equal_value(const int i) const
{
  double g_i = 0.; // no rhs, it is placed into the bounds instead
  if (valid_solution())
  {
    for (unsigned int j = 0; j < iaData->constraints[i].M.size(); ++j)
    {
      g_i += iaSolution->x_solution[ iaData->constraints[i].M[j].col ] * iaData->constraints[i].M[j].val;
    }
  }
  return g_i;
}
  
  
double IASolverTool::even_value(const int i) const
{
  double g_i = - iaData->sumEvenConstraints[i].rhs;
  if (valid_solution())
  {
    for (unsigned int j = 0; j < iaData->sumEvenConstraints[i].M.size(); ++j)
    {
      double x = iaSolution->x_solution[ iaData->sumEvenConstraints[i].M[j].col ];
      double x_int_double;
      if (is_integer(x, x_int_double))
        x = x_int_double;          
      g_i += x * iaData->sumEvenConstraints[i].M[j].val;
    }
  }       
  return g_i;
}

void IASolverTool::even_floor_ceil(double s, double &s_floor, double &s_ceil) const
{
  s_floor = 2. * floor(s / 2);
  s_ceil = s_floor + 2.;
}
  
void IASolverTool::int_floor_ceil(double x, double &x_floor, double &x_ceil) const
{
  x_floor = floor(x);
  x_ceil = x_floor + 1.;
}
  
  // compare to 
  // bool IANlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  // does the same thing with the data in ipopt format

std::pair< bool, double> IASolverTool::equal_constraint(const int i, const bool print_me, const bool print_unsatisfied ) const
{
  const double g_i = equal_value(i);
  // epsilons values
  const double epsilon_lower = 1.e-2;
  const double epsilon_upper = 1.e-2;
  bool satisfied = 
    (g_i + epsilon_lower > iaData->constraints[i].lowerBound) || 
    (g_i < iaData->constraints[i].upperBound + epsilon_upper);
  if (print_me || (print_unsatisfied && !satisfied))
  {
    printf("equal constraint %d:", i);
    for (unsigned int j = 0; j < iaData->constraints[i].M.size(); ++j)
    {
      printf( " %f x_%d", iaData->constraints[i].M[j].val, iaData->constraints[i].M[j].col );
      if (valid_solution())
        printf(" (%f)", 
               iaSolution->x_solution[ iaData->constraints[i].M[j].col ] * iaData->constraints[i].M[j].val );
    }
    if (valid_solution())
      printf(" = %f", g_i);
    printf(" in [%f, %f]", iaData->constraints[i].lowerBound, iaData->constraints[i].upperBound);
    if (!satisfied)
      printf("  VIOLATED");
    printf("\n");
  }
  return std::make_pair(satisfied, g_i);
}
  
std::pair< bool, double> IASolverTool::even_constraint(const int i, const bool print_me, const bool print_unsatisfied ) const
{
  const double g_i = even_value(i);
  // epsilon values
  int g_even;
  bool satisfied = is_even(g_i, g_even);
  if (print_me || (print_unsatisfied && !satisfied))
  {
    printf("even constraint %d:", i);
    for (unsigned int j = 0; j < iaData->sumEvenConstraints[i].M.size(); ++j)
    {
      printf( " %f x_%d", iaData->sumEvenConstraints[i].M[j].val, iaData->sumEvenConstraints[i].M[j].col );
      if (valid_solution())
        printf(" (%f) ", 
               iaSolution->x_solution[ iaData->sumEvenConstraints[i].M[j].col ] 
               * iaData->sumEvenConstraints[i].M[j].val );
    }
    if (valid_solution())
    {
      if (satisfied)
        printf(" = %f = %d", g_i, g_even);
      else
        printf(" = %f != %d NOT-EVEN", g_i, g_even);
    }
    printf("\n");
  }
  return std::make_pair(satisfied, g_i);
}
  
  
bool IASolverTool::equal_constraints( bool print_me, bool print_unsatisfied ) const
{
  bool equal_satisfied = true;
  for (unsigned int i = 0; i<iaData->constraints.size(); ++i)
  {
    if (! equal_constraint(i, print_me, print_unsatisfied).first )
      equal_satisfied = false;
  }
  return equal_satisfied;
}
  
bool IASolverTool::even_constraints( bool print_me, bool print_unsatisfied ) const
{
  bool even_satisfied = true;
  for (unsigned int i = 0; i<iaData->sumEvenConstraints.size(); ++i)
  {
    if (! even_constraint(i, print_me, print_unsatisfied).first )
      even_satisfied = false;
  }
  return even_satisfied;
}
  
bool IASolverTool::all_constraints( bool print_me, bool print_unsatisfied ) const
{
  bool equal_satisfied = equal_constraints( print_me, print_unsatisfied );
  bool even_satisfied = even_constraints( print_me, print_unsatisfied );
  return equal_satisfied && even_satisfied;
}


} // namespace MeshKit
