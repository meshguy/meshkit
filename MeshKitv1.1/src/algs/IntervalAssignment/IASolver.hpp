// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVER_HP
#define MESHKIT_IA_IASOLVER_HP

#include "IASolverToolInt.hpp"

namespace MeshKit {

class IAData;
class IASolution;

class IASolver: public IASolverToolInt
{
public:
  /** default constructor */
  IASolver(IAData *ia_data, IASolution *ia_solution);

  /** default destructor */
  virtual ~IASolver();
  
  bool solve();
  // return true if solved; false if not solved (e.g. infeasible)
  
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
 
};

} // namespace MeshKit 
#endif
