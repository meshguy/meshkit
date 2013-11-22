// IPData.cpp
// Interval Assignment Data for Meshkit
// Integer Program - making intervals into integer values 
//

#include "IPData.hpp"

namespace MeshKit {

void IPData::initialize(const std::vector<double> &relaxed_solution)
{
  const std::size_t num_variables = relaxed_solution.size();
  // fill varIntegerBound and oldBound with zeros
  varIntegerBound.clear();
  varIntegerBound.resize(num_variables,0);
  oldBound.clear();
  oldBound.resize(num_variables,0);
  relaxedSolution = relaxed_solution; // vector copy
}

 
void IPData::constrain_integer(const int i_nonint, const int x_bound)
{
  oldBound[i_nonint] = varIntegerBound[i_nonint];;
  varIntegerBound[i_nonint] = (double) x_bound;
}

} // namespace MeshKit
