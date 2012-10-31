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
  
  static void round_solution(std::vector<double> &solution);

  // standardize epsilon for determining if x is an integer solution
  static bool solution_is_integer(const std::vector<double> &solution);
  static bool is_integer(double x); 

  // measures how far above or below floor( relaxed solution ) x is, in whole integers and fractional
  static void get_frac( double x, int &integer, double &frac);
  double get_xl( int i ) const; // floor of relaxed solution
  double get_xh( int i ) const; // xl + 1

  double get_km( int i, double x ) const;  // how much bigger x is than xh, or zero if less
  double get_kp( int i, double x ) const;  // xl - x, or zero if negative
  void get_km( int i, double x, int &km_integer, double &km_frac ) const; // split into integer and fractional part
  void get_kp( int i, double x, int &kp_integer, double &kp_frac ) const;


};

// default constructors for object and its members OK
// inline IPData::IPData() {} 

// default destructor OK

} // namespace MeshKit 

#endif
