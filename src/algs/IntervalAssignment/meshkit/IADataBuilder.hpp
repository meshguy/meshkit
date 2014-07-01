// IADataBuilder.hpp
// Interval Assignment Data for Meshkit
//
// Tool for building IAData, the definition of an interval assignment problem.
// Only the IA family of classes should be using this.
//

#ifndef MESHKIT_IA_IADATABUILDER_HP
#define MESHKIT_IA_IADATABUILDER_HP

class IAData;
#include "meshkit/IAData.hpp"

namespace MeshKit 
{
    
class IADataBuilder
{
public:
  /** default constructor */
  IADataBuilder(IAData* ia_data_ptr) : ia_data( ia_data_ptr ) {}

  /** default destructor */
  virtual ~IADataBuilder() {}

  // add a variable 
  // return its id = index in I
  int add_variable( double goal );
                   
  // index of variables that must sum to even
  void constrain_sum_even(const std::vector<int> &curve_indices, const int rhs=0); 
  // index of variables that must add to the same values
  // here side_1 - side_2 = rhs, so fixed-size curves of side_1 subtract from rhs. rhs may be negative.
  void constrain_opposite_side_equal(const std::vector<int> &side_1, const std::vector<int> &side_2, 
                                     const int rhs=0 ); 
  
private:
  IAData *ia_data;
  
};

} // namespace MeshKit 

#endif
