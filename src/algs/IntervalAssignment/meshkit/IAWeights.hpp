// IAWeights.hpp
// Interval Assignment for Meshkit
//
// Weights or coefficients of 
// 
// Tools for making them unique.

#ifndef MESHKIT_IA_IAWEIGHTS_HP
#define MESHKIT_IA_IAWEIGHTS_HP

#include <vector>

namespace MeshKit {

  class IAWeights : public std::vector<double>
{

public:
  /** default constructor */
  IAWeights();

  /** default destructor */
  virtual ~IAWeights();
 
  // data
  
   // debug
  const bool debugging;
  
  // algorithms
  
  // return a number between -1 and 1, but not close to zero
  // generate double in [-1,-0.5] U [.5,1]
  static double rand_excluded_middle();
  
  // rescale and randomize so fabs of weights are in [lo, hi], 
  // and different from each other
  void uniquify(const double lo, const double hi);

  // debug
  void print() const;

};

} // namespace MeshKit 

#endif
