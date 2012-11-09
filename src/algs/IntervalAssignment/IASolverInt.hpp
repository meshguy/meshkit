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
  bool solve_minlp();
  bool solve_cos();
  bool solve_parabola();

  // workhorse
  bool solve_intwave(IAIntWaveNlp *mynlp);

  // mixed-integer solution stuff. 
  // MI NLP - solve a mixed integer nlp that we hope has an integer solution near the relaxed solution.
  // gets most of the variables integer
  // rounding heuristic
  // find a varaible that isn't integer yet, and whether we should increase or decrease its value
  // takes care of the remaining problem

  // value of objective function increase if x is rounded to the next farther integer, and index of x
  typedef std::pair<int,int> RoundingX;  // x index, x bound
  typedef std::multimap<double, RoundingX> RoundingMap; // obj value increase, roundingX
  bool find_many_non_integer(RoundingMap &rounding_map);
  bool constrain_many_non_integer(RoundingMap &rounding_map);
  void back_off(RoundingMap &rounding_map);
  void report_many_non_integer(RoundingMap &rounding_map);


};

} // namespace MeshKit 

#endif
