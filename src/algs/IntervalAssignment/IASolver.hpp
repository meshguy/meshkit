// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVER_HP
#define MESHKIT_IA_IASOLVER_HP

#include "IAData.hpp"
#include "IASolution.hpp"

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
  
  // debugging
  const bool debugging;
  void set_test_problem();
  void print_solution();
  
};

#endif
