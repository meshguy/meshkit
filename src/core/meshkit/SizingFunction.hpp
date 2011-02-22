#ifndef MESHKIT_SIZINGFUNCTION_HPP
#define MESHKIT_SIZINGFUNCTION_HPP

#include "meshkit/Types.hpp"
#include "meshkit/MKCore.hpp"

namespace MeshKit {

class MKCore;
    
/** \class SizingFunction SizingFunction.hpp "meshkit/SizingFunction.hpp"
 * \brief A sizing function used with meshing algorithms
 *
 * SizingFunction stores a \em requested size or # intervals; the actual # intervals, computed or
 * just requested, is stored on the ModelEnt.  The Firmness is also stored on the ModelEnt.
 * This class manages its own adding to MBCore, so no need for caller to do that.
 */
class SizingFunction
{
public:

    //! Copy constructor
  SizingFunction(const SizingFunction &sf);
  
    //! Constructor
  SizingFunction(MKCore *mkcore, int num_int = -1, double int_size = -1.0);

    //! Destructor
  virtual ~SizingFunction();

    //! Get core instance
  MKCore *mk_core() const;
  
    /** \brief Get size with optional location
     * \param xyz Location where size is requested
     * \return Size at requested location
     */
  virtual double size(double *xyz = NULL) const;
  
    //! Get the number of intervals (assuming they were set somehow somewhere else)
  virtual int intervals() const;
  
    /** \brief Set the number of intervals
     * \param num_int Intervals to be set
     */
  virtual void intervals(int num_int);  

    //! Return index of this sf in MKCore
  virtual unsigned int core_index() const;
  
protected:
    //! MKCore associated with this SizingFunction
  MKCore *mkCore;

    //! Interval setting
  int thisInterval;

    //! Size setting
  double thisSize;

    //! This SizingFunction's index in MKCore
  unsigned int coreIndex;
  
private:

};

  //! Copy constructor
inline SizingFunction::SizingFunction(const SizingFunction &sf) 
        : mkCore(sf.mk_core()), thisInterval(sf.intervals()), thisSize(sf.size())
{
  coreIndex = mkCore->add_sizing_function(this);
}

  //! Constructor
inline SizingFunction::SizingFunction(MKCore *mkcore, int num_int, double int_size) 
        : mkCore(mkcore), thisInterval(num_int), thisSize(int_size)
{
  coreIndex = mkCore->add_sizing_function(this);
}
    
  //! Destructor
inline SizingFunction::~SizingFunction() 
{
  mkCore->remove_sizing_function(coreIndex, false);
}
  
  //! Get core instance
inline MKCore *SizingFunction::mk_core() const 
{
  return mkCore;
}
  
  //! Get size with optional location
inline double SizingFunction::size(double *) const
{
  return thisSize;
}

  //! Get the number of intervals (assuming they were set somehow somewhere else)
inline int SizingFunction::intervals() const
{
  return thisInterval;
}

  //! Set the number of intervals
inline void SizingFunction::intervals(int num_int) 
{
  thisInterval = num_int;
}

    //! Return index of this sf in MKCore
inline unsigned int SizingFunction::core_index() const 
{
  return coreIndex;
}

} // namespace MeshKit

#endif

  
