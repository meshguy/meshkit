// IASolverToolInt.hpp
// Interval Assignment for Meshkit
// common integer problem tools

#ifndef MESHKIT_IA_IASOLVERTOOLINT_HP
#define MESHKIT_IA_IASOLVERTOOLINT_HP

#include "meshkit/IASolverTool.hpp"

namespace MeshKit {
  
class IPData;
class IAWeights;

class IASolverToolInt : public IASolverTool 
{
  // no data or inheritance, except debugging flags
public:
  /** default constructor */
  IASolverToolInt();
  IASolverToolInt(const bool debugging_set);
  IASolverToolInt(const IAData *ia_data, IASolution *ia_solution, bool debugging_set = false);

  /** default destructor */
  virtual ~IASolverToolInt();
  
  // member get/set
  IPData *ip_data() {return ipData;}
  void ip_data(IPData *set_data) {ipData = set_data;}
  virtual void set_debug(const bool debugging_set) {debuggingToolInt = debugging_set; IASolverTool::set_debug(debugging_set);}

  // tools
  // round solution in place
  void round_solution();
  
  // standardize epsilon for determining if x is an integer solution
  bool solution_is_integer(const bool print_non_integer = false);
  
  // measures how far above or below floor( relaxed solution ) x is, in whole integers and fractional
  void get_frac( double x, int &integer, double &frac) const;
  double get_xl( int i ) const; // floor of relaxed solution
  double get_xh( int i ) const; // xl + 1
  
  double get_km( int i, double x ) const;  // how much bigger x is than xh, or zero if less
  double get_kp( int i, double x ) const;  // xl - x, or zero if negative
  void get_km( int i, double x, int &km_integer, double &km_frac ) const; // split into integer and fractional part
  void get_kp( int i, double x, int &kp_integer, double &kp_frac ) const;
  
  // randomize the value of weights for the non-integer values
  // return true if any weights were changed
  bool randomize_weights_of_non_int(IAWeights* weights, const double rand_factor);
  
protected:  
  
  IPData *ipData;

  // debugging
  bool debuggingToolInt; 
  
};

} // namespace MeshKit 
#endif
