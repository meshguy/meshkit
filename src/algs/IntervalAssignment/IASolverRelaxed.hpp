// IASolver.hpp
// Interval Assignment for Meshkit
// Meshkit calls this class function to find the intervals
// this class calls underlying solvers
//
#ifndef MESHKIT_IA_IASOLVER_RELAXED_HP
#define MESHKIT_IA_IASOLVER_RELAXED_HP

#include "IAData.hpp"
#include "IASolution.hpp"

#include <map> 

class IANlp;
class IPData;

class IASolverRelaxed
{
public:
  /** default constructor */
  // ia_data contains the problem specification
  // relaxed_solution contains the relaxed solution after solve()
  IASolverRelaxed(const IAData *ia_data, IASolution *relaxed_solution);  

  /** default destructor */
  virtual ~IASolverRelaxed();
  
  bool solve();
  // return true if solved; false if not solved (e.g. infeasible)
  
  static double check_constraint(const int i, const IAData * ia_data, const IASolution *solution, bool sum_even, bool print_me );
  static bool constraints_satisfied( const IAData * ia_data, const IASolution *solution, bool &sum_even_satisfied, bool print_me );
  
private:  
  const IAData *iaData;
  IASolution *relaxedSolution;
  
  // data
  int p_norm;
  
  // hide untrusted default methods
  //@{
  //  IA_NLP();
  IASolverRelaxed(const IASolverRelaxed&);
  IASolverRelaxed& operator=(const IASolverRelaxed&);
  //@}

  bool constraints_satisfied() { bool junk; return constraints_satisfied(iaData, relaxedSolution, junk, debugging); }

  // debug
  const bool debugging;
  
};

#endif
