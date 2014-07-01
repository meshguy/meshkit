// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVER_RELAXED_HP
#define MESHKIT_IA_IASOLVER_RELAXED_HP

#include "meshkit/IASolverTool.hpp"

namespace MeshKit {

class IASolverRelaxed : public IASolverTool
{
public:
  /** default constructor */
  // ia_data contains the problem specification
  // relaxed_solution contains the relaxed solution after solve()
  IASolverRelaxed(const IAData *ia_data, IASolution *relaxed_solution, const bool set_silent = true);  

  /** default destructor */
  virtual ~IASolverRelaxed();
  
  bool solve();
  // return true if solved; false if not solved (e.g. infeasible)

private:  
 
  // data
  int p_norm;
  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IASolverRelaxed(const IASolverRelaxed&);
  IASolverRelaxed& operator=(const IASolverRelaxed&);
  //@}

  bool constraints_satisfied();

  // debug
  const bool silent;
  const bool debugging;
  
};

} // namespace MeshKit
#endif
