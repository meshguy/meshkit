#ifndef ERROR
#define ERROR

/** \file Error.hpp
 */
#include <string>
#include <typeinfo>

namespace MeshKit {

#define MKERRCHK(err, descr) \
    {if (MK_SUCCESS != err.error_code()) throw err;}
    
#define MBERRCHK(err, descr) \
    {if (moab::MB_SUCCESS != err) throw Error(err, descr);}
    
#define IBERRCHK(err, descr) \
    {if (iBase_SUCCESS != err) throw Error(err, descr);}
    
    enum ErrorCode {
        MK_SUCCESS = 0, 
        MK_FAILURE,
        MK_NOT_FOUND,
        MK_MULTIPLE_FOUND,
        MK_MESHOP_NOT_FOUND,
        MK_NOT_IMPLEMENTED,
        MK_WRONG_DIMENSION,
        MK_ALREADY_DEFINED,
        MK_BAD_INPUT,
        MK_BAD_GEOMETRIC_EVALUATION,
        MK_INCOMPLETE_MESH_SPECIFICATION
    };
    
/** \class Error Error.hpp "meshkit/Error.hpp"
 * \brief The Error object returned from or thrown by many MeshKit functions    
 *
 * This class is derived from std::exception.
 */
class Error : public std::exception
{
public:

    /** \brief Constructor
     * \param err %Error code to which this Error should be initialized
     * \param descr String describing the type of error represented by this object
     */
  Error(int err);
  
  Error(int err, const char* format, ...)
#ifdef __GNUC__
    __attribute__((format(printf,3,4)))
#endif
    ;

  Error() {};

    //! Destructor
  virtual ~Error() throw();

    //! Return the error code
  virtual ErrorCode error_code() const;
  
    //! Standard function for retrieving a description of this Error
  virtual const char *what() const throw();

private:
    //! Locally-stored error code, from ErrorCode enumeration
  ErrorCode errorCode;

    //! String-based description of error
  std::string errDescription;
};

inline Error::Error(int err) : errorCode((ErrorCode)err) {}

inline Error::~Error() throw () {}

inline ErrorCode Error::error_code() const
{
  return errorCode;
}

inline const char *Error::what() const throw ()
{
  return errDescription.c_str();
}
      
}

#endif
