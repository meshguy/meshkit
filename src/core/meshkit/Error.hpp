#ifndef ERROR
#define ERROR

/** \file Error.hpp
 */
#include <string>
#include <typeinfo>

typedef int ErrorCode;

#define MKERRCHK(err, descr) \
    {if (MK_SUCCESS != err) throw Error(err, descr);}
    
namespace MeshKit {

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
  Error(int err, const char *descr = NULL);

    //! Copy constructor
  Error(const Error &err);

    //! Operator=
  Error &operator=(const Error &err);

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

inline Error::Error(int err, const char *descr) : errorCode(err), errDescription(descr) {}

inline Error::Error(const Error &err) : errorCode(err.error_code()), errDescription(err.what()) {}

inline Error::~Error() throw () {}
    
inline Error &Error::operator=(const Error &err) 
{
  errorCode = err.error_code();
  errDescription.assign(err.what());
  return *this;
}

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
