// IADataBulder.cpp
// Interval Assignment for Meshkit

#include "IAData.hpp"
#include "IADataBuilder.hpp"

namespace MeshKit 
{

int IADataBuilder::add_variable( const double goal )
{
  assert(ia_data);
  int variable_index = (int) ia_data->I.size();
  ia_data->I.push_back(goal);
  return variable_index;
}
  

void IADataBuilder::constrain_sum_even(const std::vector<int> &curve_indices, const int rhs)
{
  assert(ia_data);
  if (!curve_indices.size())
    return;
  ia_data->sumEvenConstraints.push_back( IAData::sumEvenConstraintRow() );
  ia_data->sumEvenConstraints.back().M.reserve(curve_indices.size());
  ia_data->sumEvenConstraints.back().rhs = rhs;
  for (unsigned int j = 0; j < curve_indices.size(); ++j)
  {
    ia_data->sumEvenConstraints.back().M.push_back( IAData::sparseEntry( curve_indices[j], 1. ) );
  }
  assert( ia_data->sumEvenConstraints.back().M.size() == curve_indices.size() );
}

void IADataBuilder::constrain_opposite_side_equal(const std::vector<int> &side_1, const std::vector<int> &side_2, 
                                           const int rhs )
{
  assert(ia_data);
  const size_t num_curves = side_1.size() + side_2.size();
  if (!num_curves)
    return;
  ia_data->constraints.push_back( IAData::constraintRow() );
  ia_data->constraints.back().M.reserve( num_curves );
  ia_data->constraints.back().upperBound = rhs;
  ia_data->constraints.back().lowerBound = rhs;
  for (unsigned int j = 0; j < side_1.size(); ++j)
  {
    ia_data->constraints.back().M.push_back( IAData::sparseEntry( side_1[j], 1. ) );
  }
  for (unsigned int j = 0; j < side_2.size(); ++j)
  {
    ia_data->constraints.back().M.push_back( IAData::sparseEntry( side_2[j], -1. ) );
  }
  assert( ia_data->constraints.back().M.size() == num_curves );
}

} // namespace MeshKit 

