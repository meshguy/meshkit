// IAWeights.cpp
// Interval Assignment for Meshkit
//
// Weights

#include "IAWeights.hpp"

#include <math.h>
#include <limits>
#include <algorithm>
#include <assert.h>
#include <cstdio>

namespace MeshKit 
{


// generate integers from 0..n-1
struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return current++;}
} UniqueNumber;


class FabsWeightComparer
{
public: 
  std::vector<double> *w;
  FabsWeightComparer(std::vector<double> *weightvec) {w = weightvec;}
  bool operator() (const int lhs, const int rhs)
  {
    return fabs((*w)[lhs]) < fabs((*w)[rhs]);
  }
};

IAWeights::IAWeights() : std::vector<double>(), 
debugging(true)  
{}

IAWeights::~IAWeights() 
{}
  
  // static
double IAWeights::rand_excluded_middle()
{
  // generate double in [-1,-0.5] U [.5,1]
  double d = ((double) rand() / RAND_MAX) - 0.5;
  if (d<0.)
    d -= 0.5;
  else
    d += 0.5;
  assert( d >= -1. );
  assert( d <= 1. );
  assert( d <= -0.5 || d >= 0.5 );
  return d;
}

void IAWeights::uniquify(const double lo, const double hi)
{
  assert( hi >= lo );
  assert( lo >= 0. );

  // find min an max of input
  double fabs_min_weight = std::numeric_limits<double>::max();
  double fabs_max_weight = 0.;
  for (unsigned int i = 0; i < size(); ++i)
  {
    const double w = (*this)[i];
    const double fabsw = fabs(w);
    if (fabsw < fabs_min_weight)
      fabs_min_weight = fabsw;
    if (fabsw > fabs_max_weight)
      fabs_max_weight = fabsw;
  }
  
  // relative range of input and output 
  const double input_range = fabs_max_weight - fabs_min_weight; 
  assert( input_range >= 0.);
  const double output_range = hi - lo;
  assert( output_range >= 0.);
  
  // scale the weightvec so | max | is 1.e4, and min is 1
  // the range should be well below the ipopt solver tolerance, which is 1.0e-7
  // "typically" the raw weights are between 1.e-4 and 10
  if (fabs_max_weight < 1.)
    fabs_max_weight = 1.;
  //was const double s = 1.e4 / fabs_max_weight;
  double s = output_range / input_range; // could be nan, so limit in next line
  if ( s > 1.e8 )
    s = 1.e8;
  for (unsigned int i = 0; i < size(); ++i)
  {
    const double fabsw = lo + ( (fabs((*this)[i]) - fabs_min_weight) * s );
    (*this)[i] = (*this)[i] > 0. ? fabsw : -fabsw;

    /* was
    weightvec[i] *= s;
    if (fabs(weightvec[i]) < 1.)
      if (weightvec[i] < 0.)
        weightvec[i] = -1.;
      else
        weightvec[i] = 1.;
     */
  }
  
  // uniquify the weights. 
  // ensure a random minimum ratio between consecutive weights
  // we'd really like no weight to be the sum of other weights, but randomization should catch most of the cases.
  // We randomize rather than make a deterministict fraction.
  
  // get the indices of the sorted order of the weights
  std::vector<int> sorted_fabs_weights(size());
  std::generate(sorted_fabs_weights.begin(), sorted_fabs_weights.end(), UniqueNumber );
  std::sort( sorted_fabs_weights.begin(), sorted_fabs_weights.end(), FabsWeightComparer(this) );
  
  
  srand(9384757); 
  double prior_fw = 0.;
  for (unsigned int i = 0; i < size(); ++i)
  {
    // detect consecutive identical weights and modify the later one
    const int j = sorted_fabs_weights[i]; // index of weight in weights
    double w = (*this)[j];
    double fw = fabs(w);
    if (fw - prior_fw < lo * 0.01) // relative tolerance
    {
      const double eps = lo * (0.01 + 0.02 * ((double) rand() / RAND_MAX));
      fw = prior_fw + eps; // use prior_fw rather than fw to ensure a min gap, uniform in [0.01, 0.03]
      (*this)[j] = w = (w<0.) ? -fw : fw;
    }
    prior_fw = fw;
    //printf("%d: w_%d %10.4f\n", i, j, weightvec[j]); 
  }
  
  // scale again if max was exceeded
  if (prior_fw > hi)
  {
    // assume minimum is lo, ignoring the random bit we added
    const double ss = output_range / ( prior_fw - lo );
    for (unsigned int i = 0; i < size(); ++i)
    {
      double w = (fabs((*this)[i]) - lo) * ss + lo;
      assert( w >= 0. );
      // with roundoff, could be slightly above hi, force it
      if ( w > hi )
        w = hi;
      if ( w < lo )
        w = lo;
      (*this)[i] = ((*this)[i] < 0.) ? -w : w; 
      assert( w <= hi );
      assert( w >= lo );
    }
  }
  
  if (0) // debugging
  {
    printf("unique weights with fabs in [%e, %e]\n", lo, hi);
    for (unsigned int i = 0; i < size(); ++i)
    {
      const int j = sorted_fabs_weights[i]; // index of weight in weightvec
      const double w = (*this)[j];
      printf("%d: w_%d %10.4f\n", i, j, w); 
      assert( fabs(w) <= hi );
      assert( fabs(w) >= lo );
    }
  }
  // exit(1); //zzyk
}

// debug
void IAWeights::print() const
{
  printf("weights:\n");
  for (unsigned int i = 0; i < size(); ++i)
  {
    const double w = (*this)[i];
    printf("w[%u] = %f\n", i, w);
  }
  printf("\n");

}

} // namespace MeshKit
