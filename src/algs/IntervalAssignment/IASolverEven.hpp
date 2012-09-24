// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVEREVEN_HP
#define MESHKIT_IA_IASOLVEREVEN_HP

#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"

#include <map> //todo move this data to an auxillary class

class IANlp;
class IAMINlp;
class IPData;

class IASolverEven: public IASolution 
{
public:
  /** default constructor */
  IASolverEven(const IAData *ia_data, const IASolution *int_solution);

  /** default destructor */
  virtual ~IASolverEven();
  
  bool solve();
  // return true if solved; false if not solved (e.g. infeasible)

  // some utility functions
  static bool is_even(double y);
  static double sum_even_value(int i, const IAData *ia_data_ptr, const IASolution *solution);

private:  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IASolverEven(const IASolverEven&);
  IASolverEven& operator=(const IASolverEven&);
  //@}
  
 
  // print x that is currently non-integer
  void report_one_constraint(const int i, const int b); 
  // print x that was non-integer, but was constrained to be integer and resolved
  void report_one_non_integer(const int i, const int b); 
  
  bool calculate_rounding_value(const int i, double &obj_increase, int &y_int );
  // y_int = what x[i] should be rounded to
  // obj_increase is the increase in the objective function value if x[i] changes to y_int
  
  bool calculate_sumeven_value(const int i, double &obj_increase, int &y_int );
  // assumes x[i] is already an integer
  // y_int = what x[i] should be rounded to, in order to help satisfy a sum-even constraint
  // obj_increase is the increase in the objective function value if x[i] changes to y_int


  void constrain_integer(const int i_nonint, const int x_bound);
  bool find_one_non_integer(int &i_nonint, int &x_bound);
  bool constrain_one_non_integer(int &i, int &b);
  
  // value of objective function increase if x is rounded to the next farther integer, and index of x
  typedef std::pair<int,int> RoundingX;  // x index, x bound
  typedef std::multimap<double, RoundingX> RoundingMap; // obj value increase, roundingX
  bool find_many_non_integer(RoundingMap &rounding_map);
  bool constrain_many_non_integer(RoundingMap &rounding_map);
  void back_off(RoundingMap &rounding_map);
  void report_many_non_integer(RoundingMap &rounding_map);

  double sum_even_value(int i);
  double distance_to_even(int i);
  
  // data
  const IAData *iaData;
  IPData ipData; // holds integer solution
  IAMINlp *myianlp;

  void cleanup(); // clean up memory management of the above pointers.
  
  // debug
  const bool debugging;
  void print_solution(IPData &ip_data);
  
};

#endif
