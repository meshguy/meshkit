// IAData.hpp
// Interval Assignment Data for Meshkit
//
// This is the underlying data representation for the constraints needed by Meshkit 
// and the solver interface classes. 
// Only the IA family of classes should be using this.
//
// Ipopt appears to use Number for double.
//

#ifndef MESHKIT_IA_IADATA_HP
#define MESHKIT_IA_IADATA_HP

#include <vector>
#include <assert.h>

const double MESHKIT_IA_upperUnbound = 2e19; // above 1e19, ipopt specific
const double MESHKIT_IA_lowerUnbound = -2e19; // below -1e19, ipopt specific 


class IAData
{
public:
  /** default constructor */
  IAData() {};

  /** default destructor */
  virtual ~IAData() {};

  // ipopt solver forces the use of smartpointers which causes uncontrollable memory management and 
  // hard to track memory bugs, so we have IANlp point to our data instead of inheriting this as a base class

  // number of goals
  // number of constraints
  // goals
  // variables
  // constraints 
  //   vector of constraints
  //   each constraint are the coefficients of the linear program, in sparse format

  // choice of numeric type and index type matters. 
  // Ipopt uses typedef Number double and typdef Index int, so 
  // lets use double and int for now
  
  // size_t numVariables, numConstraints; implicit, use x.size()

  std::vector<double> I; // intervalGoals;
  int num_variables() const { return (int) I.size(); }
    
  struct sparseEntry
  {
    // int row; // defined implicitly by index of contraint in constraints vector
    int col;
    double val; 
    sparseEntry(const int c, const double v) : col(c), val(v) {}
    sparseEntry() : col(-1), val(0.) {} // bad values if unspecified
  };

  struct constraintRow
  {
    std::vector<sparseEntry> M; // non-zeros in the constraint row
    // double rhs; // righthand side constant incorporated into upper and lower bound
    double upperBound, lowerBound;  // upper and lower bounds. equality for equal
    constraintRow() : upperBound(0.), lowerBound(0.) {};
    // or upperBound( MESHKIT_IA_upperUnbound), lowerBound(MESHKIT_IA_lowerUnbound) {}
  };
    
  struct sumEvenConstraintRow
  {
    std::vector<sparseEntry> M;
    int rhs; // right-hand-side constant
    sumEvenConstraintRow(): rhs(0) {}
  };
  
  std::vector<constraintRow> constraints;
  
  std::vector<sumEvenConstraintRow> sumEvenConstraints;
  
  // index of variables that must sum to even
  // to do: associate curves with some unique tag, and add an interface function that takes those tags.
  // then at solve-time associate an index in sequence with that tag.
  void constrain_sum_even(const std::vector<int> &curve_indices, const int rhs=0); 
  // here side_1 - side_2 = rhs, so fixed-size curves of side_1 subtract from rhs. rhs may be negative.
  void constrain_opposite_side_equal(const std::vector<int> &side_1, const std::vector<int> &side_2, 
                                     const int rhs=0 ); 
  
};

// default constructors for object and its members OK
// inline IAData::IAData() {} 

// default destructor OK
// to do, make a .cc file and add this there
inline
void IAData::constrain_sum_even(const std::vector<int> &curve_indices, const int rhs)
{
  if (!curve_indices.size())
    return;
  sumEvenConstraints.push_back( sumEvenConstraintRow() );
  sumEvenConstraints.back().M.reserve(curve_indices.size());
  sumEvenConstraints.back().rhs = rhs;
  for (int j = 0; j < curve_indices.size(); ++j)
  {
    sumEvenConstraints.back().M.push_back( sparseEntry( curve_indices[j], 1. ) );
  }
  assert( sumEvenConstraints.back().M.size() == curve_indices.size() );
}

inline
void IAData::constrain_opposite_side_equal(const std::vector<int> &side_1, const std::vector<int> &side_2, 
                                           const int rhs )
{
  const size_t num_curves = side_1.size() + side_2.size();
  if (!num_curves)
    return;
  constraints.push_back( constraintRow() );
  constraints.back().M.reserve( num_curves );
  constraints.back().upperBound = rhs;
  constraints.back().lowerBound = rhs;
  for (int j = 0; j < side_1.size(); ++j)
  {
    constraints.back().M.push_back( sparseEntry( side_1[j], 1. ) );
  }
  for (int j = 0; j < side_2.size(); ++j)
  {
    constraints.back().M.push_back( sparseEntry( side_2[j], -1. ) );
  }
  assert( constraints.back().M.size() == num_curves );
}


#endif
