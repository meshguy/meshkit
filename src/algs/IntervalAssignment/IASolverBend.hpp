// IASolverBend.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVERBEND_HP
#define MESHKIT_IA_IASOLVERBEND_HP

#include <map>

#include "IASolverToolInt.hpp"
#include "IPBend.hpp"
#include "IAWeights.hpp"
// #include "IPData.hpp"
// #include "IABendNlp.hpp"

namespace MeshKit {
  
class IAData;
class IASolution;
class IABendNlp;
class IPData;

class IASolverBend : public IASolverToolInt
{
public:
  /** default constructor */
  IASolverBend(const IAData * ia_data_ptr, IASolution *relaxed_solution_ptr, const bool set_silent = true);

  /** default destructor */
  virtual ~IASolverBend();
  
  bool solve();
  // return true if solved; false if not solved (e.g. infeasible)

private:  
  // hide untrusted default methods
  //@{
  IASolverBend(const IASolverBend&);
  IASolverBend& operator=(const IASolverBend&);
  //@}
  
  // if true, then sum-even constraints are enforced, otherwise ignored.
  // must remain the same during a call to solve() = ipopt as it changes the 
  // problem size and structure.
  bool evenConstraintsActive;
  
  // debug
  const bool silent;
  const bool debugging;    
  
  // top level

  // set initial ip bends from relaxed solution
  void initialize_ip_bends();

  // update them based on fractional solutions, large deltas
  bool update_ip_bends();
  
  // call the nlp solver
  bool solve_nlp();

  void cleanup();

  // round solution to nearest integer
  // if its feasible, then replace current solution and return true
  bool round_solution();
  
  IABendNlp *myianlp;

  IPBendData bendData;

  // weight index for deltas start at delta_i_start - num_variables   
  IAWeights weights;

  // convert IPBends to IAWeights
  void add_bend_weights(unsigned int i);
  // void add_bend_sum_weights(unsigned int i, const double factor);

  // utility
  double f_x_value( double I_i, double x_i ) const;
  double fpow(double f) const;
  double get_f_xl(int i) const; // obj function value at xl
  double get_f_xh(int i) const; // obj function value at xh

  
};

} // namespace MeshKit 

#endif
