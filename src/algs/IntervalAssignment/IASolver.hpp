// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVER_HP
#define MESHKIT_IA_IASOLVER_HP

#include "IAData.hpp"
#include "IASolution.hpp"

namespace MeshKit {

class IANlp;
class IPData;

class IASolver: public IAData, public IASolution 
{
public:
  /** default constructor */
  IASolver();

  /** default destructor */
  virtual ~IASolver();
  
  bool solve();
  // return true if solved; false if not solved (e.g. infeasible)
  
  // utility functions, todo move
  static bool is_even(double y);
  static double sum_even_value(int i, const IAData *ia_data_ptr, const IASolution *solution);

private:  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IASolver(const IASolver&);
  IASolver& operator=(const IASolver&);
  //@}
  
  // workhorse
  bool solve_relaxed();
  bool solve_int();
  bool solve_even();
  double sum_even_value(int i);
  
  // debugging
  const bool debugging;
  void print_solution() const;
  void print_problem() const;
  void print(const bool do_print_solution, const bool do_print_constraints) const;

  
};

} // namespace MeshKit 
#endif
