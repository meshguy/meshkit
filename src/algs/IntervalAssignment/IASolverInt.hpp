// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVERINT_HP
#define MESHKIT_IA_IASOLVERINT_HP

#include "IAData.hpp"
#include "IPData.hpp"
#include "IASolution.hpp"

#include <map> 

namespace MeshKit {

class IASolverInt: public IASolution
{
public:
  /** default constructor */
  IASolverInt(const IAData * ia_data_ptr, const IASolution *relaxed_solution_ptr);

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
  
  IPData ipData; // stores relaxed solution
  const IAData *iaData;
    // debug
  const bool debugging;  


  void cleanup();
  
  // top level
  bool solve_minlp();
  bool solution_is_integer() {return IPData::solution_is_integer(this->x_solution);}
  
  // mixed-integer solution stuff. 
  // MI NLP - solve a mixed integer nlp that we hope has an integer solution near the relaxed solution.
  // gets most of the variables integer
  // rounding heuristic
  // find a varaible that isn't integer yet, and whether we should increase or decrease its value
  // takes care of the remaining problem
  
  
  // to do, remove some of this stuff
  
  bool calculate_sumeven_value(const int i, double &obj_increase, int &y_int );
  // assumes x[i] is already an integer
  // y_int = what x[i] should be rounded to, in order to help satisfy a sum-even constraint
  // obj_increase is the increase in the objective function value if x[i] changes to y_int


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
