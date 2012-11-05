// IASolverTool.hpp
// Interval Assignment for Meshkit
// common solver utility functions for checking and computing data

#ifndef MESHKIT_IA_IASOLVERTOOL_HP
#define MESHKIT_IA_IASOLVERTOOL_HP

#include <cstddef>
#include <utility> // for pair

namespace MeshKit {

  class IAData;
  class IASolution;

class IASolverTool
{
  // no data or inheritance, except debugging flags
public:
  /** default constructor */
  IASolverTool() 
  : iaData(NULL), iaSolution(NULL), debuggingTool(false) 
  {}
  IASolverTool( const IAData *ia_data, IASolution *ia_solution, bool debugging_set = false) 
  : iaData(ia_data), iaSolution(ia_solution), debuggingTool( debugging_set ) 
  {}
  IASolverTool(const bool debugging_set) 
    : iaData(NULL), iaSolution(NULL), debuggingTool(debugging_set) 
    {}
  
  /** default destructor */
  virtual ~IASolverTool();

  // member get/set
  const IAData *ia_data() const {return iaData;}
  IASolution *ia_solution() const {return iaSolution;}
  void ia_data(const IAData *set_data) {iaData = set_data;}
  void ia_solution(IASolution *set_solution) {iaSolution = set_solution;}
  virtual void set_debug(const bool debugging_set) {debuggingTool = debugging_set;}

  
  // data checking
  // data interogation

  // true if a solution has been defined, and it is the right size to match the data
  bool valid_solution() const;

  // true if x is within standardized epsilon of an integer
  bool is_integer(const double x) const; 
  // also rounds x to nearest integer (up or down)
  bool is_integer(const double x, int &x_int) const; 
  bool is_integer(const double x, double &x_int_double) const;

  // data computation
  bool is_even(double y) const;
  bool is_even(double y, int &y_even) const;
  // see also even_constraint

  // compute the value of the ith constraint
  // the even one rounds nearly-integer x to integer x
  double equal_value(int i) const;
  double even_value(int i) const;
  
  // first is true if constraint is satisfied
  // second contains the value of the constraint
  // ith constraint
  std::pair< bool, double> equal_constraint(const int i, 
                                            const bool print_me=false, 
                                            const bool print_unsatisfied=false ) const;
  std::pair< bool, double> even_constraint(const int i, 
                                           const bool print_me=false, 
                                           const bool print_unsatisfied=false ) const;
  // all constraints satisfied?
  bool equal_constraints( const bool print_me=false, 
                          const bool print_unsatisfied=false ) const;
  bool even_constraints( const bool print_me=false, 
                        const bool print_unsatisfied=false ) const;
  // all types
  bool all_constraints( const bool print_me=false, 
                        const bool print_unsatisfied=false ) const;
  
  // data printing for debugging
  void print_solution() const;
  void print_problem() const;
  void print( const bool do_print_solution, const bool do_print_nonint,
             const bool do_print_equal_constraints, const bool do_print_nonequal,
             const bool do_print_even_constraints, const bool do_print_noneven
             ) const;
  void print() const
    { print (true, true, true, true, true, true); }
  void print_violations() const
    { print(false, true, false, true, false, true); }

protected:  

  // problem definition and solution
  const IAData *iaData;
  IASolution *iaSolution;
  
  // debugging
  bool debuggingTool; 

};

} // namespace MeshKit 
#endif
