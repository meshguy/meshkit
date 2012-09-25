// IASolution.hpp
// Interval Assignment Solution for Meshkit
//

#ifndef MESHKIT_IA_IASOLUTION_HP
#define MESHKIT_IA_IASOLUTION_HP
 
#include <vector>

namespace MeshKit {
    
class IASolution
{
public:
  /** default constructor */
  IASolution() {};

  /** default destructor */
  virtual ~IASolution() {};

  std::vector<double> x_solution; // intervalSolution, the thing we're solving for.
  double obj_value; // relative measure of the goodness of the solution. Smaller is better.
};

// default constructors and destructors are fine

} // namespace MeshKit
#endif
