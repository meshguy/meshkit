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
#include <cstddef>
#include <assert.h>

namespace MeshKit 
{
    
const double MESHKIT_IA_upperUnbound = 2e19; // above 1e19, ipopt specific
const double MESHKIT_IA_lowerUnbound = -2e19; // below -1e19, ipopt specific 


class IAData
{
public:
  /** default constructor */
  IAData() {}

  /** default destructor */
  virtual ~IAData() {}

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
    constraintRow() : upperBound(0.), lowerBound(0.) {}
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
  
};

// default constructors for object and its members OK
// inline IAData::IAData() {} 

} // namespace MeshKit 

#endif
