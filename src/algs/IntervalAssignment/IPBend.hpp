// IPBend.hpp
// Interval Assignment Data for Meshkit
//
// Data for defining an objective function that is piecewise linear and bends at integer points
//

#ifndef MESHKIT_IA_IPBEND_HP
#define MESHKIT_IA_IPBEND_HP

#include <vector>
#include <limits>

namespace MeshKit 
{

class IPBend
{
public:
  // index of delta, magnitude of tilt = weight multiplier
  // note tilt is applied to that delta and all subsequent ones
  typedef std::pair<unsigned int, double> IPTilt;
  
  /** constructor */
  IPBend() :
    deltaIStart(0), numDeltaPlus(1), numDeltaMinus(1), 
    plusTilts(), minusTilts()  // usually these vecs are empty
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
  
  static int num_deltas_max() 
  { return std::numeric_limits<int>::max()/2 - 2; }
  
  // sorted tilts, for x > g and x < g
  std::vector< IPTilt > plusTilts, minusTilts; 
  
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
