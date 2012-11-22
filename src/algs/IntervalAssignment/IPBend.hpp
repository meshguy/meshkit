// IPBend.hpp
// Interval Assignment Data for Meshkit
//
// Data for defining an objective function that is piecewise linear and bends at integer points
//

#ifndef MESHKIT_IA_IPBEND_HP
#define MESHKIT_IA_IPBEND_HP

#include <vector>

namespace MeshKit 
{

class IPBend
{
public:
  /** constructor */
  IPBend() :
    deltaIStart(0), numDeltaPlus(1), numDeltaMinus(1)
  {}
  
  /** default destructor */
  virtual ~IPBend() {}
  
  // index i of x[i] or sum-even constraint i is implicitly defined by position in IPBendVec

  int xl; // x value if all deltas are zero
  
  int deltaIStart; // index of the first delta variables in the Nlp

  // delta0 is always defined, one of them
  int numDeltaPlus; // number of plus deltas defined
  int numDeltaMinus; // number of minus deltas defined
  int num_deltas() const 
    {return numDeltaPlus + numDeltaMinus; } 
  
  // x = xl + sum ( deltas_plus ) - sum (deltas_minus );
  // delta_plus[k] for k < numDeltaPlus + 1 is in [0,1], last is [0,infinity]
  
  
};

typedef std::vector< IPBend > IPBendVec;
  
class IPBendData
{
public:
  IPBendVec bendVec;
  int numSumVars;
  int sumVarStart;
  double maxActiveVarWeight;
};
  
// default constructors for object and its members OK
// inline IPData::IPData() {} 

// default destructor OK

} // namespace MeshKit 

#endif
