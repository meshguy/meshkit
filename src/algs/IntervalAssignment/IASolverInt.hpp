// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVERINT_HP
#define MESHKIT_IA_IASOLVERINT_HP

#include <map>

#include "IASolverToolInt.hpp"
#include "IPData.hpp"

namespace MeshKit {
  
class IAData;
class IASolution;
class IAIntWaveNlp;

class IASolverInt : public IASolverToolInt
{
public:
  /** default constructor */
  IASolverInt(const IAData * ia_data_ptr, IASolution *relaxed_solution_ptr, const bool set_silent = true);

  /** default destructor */
  virtual ~IASolverInt();
  
  bool solve();
  // return true if solved; false if not solved (e.g. infeasible)
    
private:  
  // hide untrusted default methods
  //@{
  IASolverInt(const IASolverInt&);
  IASolverInt& operator=(const IASolverInt&);
  //@}
  
  // debug
  const bool silent;
  const bool debugging;  


  void cleanup();
  
  // top level
  // several different problem formulations and their solutions are possible
  enum SolverType {ROUNDING, BEND, COS, PARABOLA};
  bool solve_round();
  bool solve_wave(const SolverType solver_type);
  // workhorse, lower level routines
  bool solve_wave_workhorse(IAIntWaveNlp *mynlp);

};

} // namespace MeshKit 

#endif
