#ifndef SIZINGFUNCTION
#define SIZINGFUNCTION

#include <cfloat>

#include "meshkit/Types.h"

namespace MeshKit {

class MKCore;
    

/** \class SizingFunction SizingFunction.hpp "meshkit/SizingFunction.hpp"
 * \brief A sizing function used with meshing algorithms
 */
class SizingFunction
{
public:

    //! Copy constructor
  SizingFunction(const SizingFunction &sf) throw(Error);
  
    //! Bare constructor
  SizingFunction(MKCore *mkcore) throw(Error);

    //! Constructor
  SizingFunction(MKCore *mkcore, int num_int) throw(Error);

    //! Destructor
  virtual ~SizingFunction() throw(Error);

    //! Get core instance
  MKCore *mk_core() const;
  
    //! Get size with optional location
  virtual double size(double *xyz = NULL) const;
  
    //! Get the number of intervals (assuming they were set somehow somewhere else)
  virtual int intervals() const;
  
    //! Set the number of intervals
  virtual void intervals(int num_int);  
  
protected:
    //! MKCore associated with this SizingFunction
  MKCore *mkCore;

    //! Interval setting
  int thisInterval;

    //! Size setting
  double thisSize;
  
private:

};

  //! Copy constructor
SizingFunction::SizingFunction(const SizingFunction &sf) throw(Error) 
        : mkCore(sf.mk_core()), thisInterval(sf.intervals()), thisSize(sf.size())
{}

  //! Bare constructor
SizingFunction::SizingFunction(MKCore *mkcore) throw(Error) 
        : mkCore(mkcore), thisInterval(-1), thisSize(-DBL_MAX) 
{}

  //! Constructor
SizingFunction::SizingFunction(MKCore *mkcore, int num_int) throw(Error) 
        : mkCore(mkcore), thisInterval(num_int), thisSize(-DBL_MAX)
{}
    
  //! Destructor
SizingFunction::~SizingFunction() throw(Error) 
{}
  
  //! Get core instance
MKCore *SizingFunction::mk_core() const 
{
  return mkCore;
}
  
  //! Get size with optional location
double SizingFunction::size(double *) const
{
  return thisSize;
}

  //! Get the number of intervals (assuming they were set somehow somewhere else)
int SizingFunction::intervals() const
{
  return thisInterval;
}

  //! Set the number of intervals
void SizingFunction::intervals(int num_int) 
{
  thisInterval = num_int;
}

#endif

  
