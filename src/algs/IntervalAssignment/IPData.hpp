// IPData.hpp
// Interval Assignment Data for Meshkit
// Integer Program - making intervals into integer values 
//
// Only the IA family of classes should be using this.
//
// Ipopt appears to use Number for double.
//

#ifndef MESHKIT_IA_IPDATA_HP
#define MESHKIT_IA_IPDATA_HP

#include <vector>

namespace MeshKit 
{

class IPData
{
public:
  /** default constructor */
  IPData() {}
  
  /** default destructor */
  virtual ~IPData() {}
  
  // lower or upper integer bounds that try to force an integer solution 
  // compare to the goal to tell if it is an upper or lower bound
  std::vector<int> varIntegerBound; 
  std::vector<int> oldBound; 
  std::vector<double> relaxedSolution;

  void initialize(const std::vector<double> &relaxed_solution);
  
  void constrain_integer(const int i_nonint, const int x_bound);
  
};

// default constructors for object and its members OK
// inline IPData::IPData() {} 

// default destructor OK

} // namespace MeshKit 

#endif
