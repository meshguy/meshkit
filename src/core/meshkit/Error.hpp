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

    /** \brief Destructor
     *
     * Must have throw() to match std::exception.
     */
  virtual ~Error() throw();

    //! Return the error code
  virtual ErrorCode error_code() const;
  
    /** \brief Standard function for retrieving a description of this Error
     *
     * Must have throw() to match std::exception.
     * \return Text description of error.
     */
  virtual const char *what() const throw();

    /** \brief Return an error string for the specified code
     * \param err Error code whose string is being queried
     * \return String corresponding to the error code
     */
  static const char* error_str(ErrorCode err);

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
      
inline const char* Error::error_str(ErrorCode err)
{
  switch (err) {
    case MK_SUCCESS: return "Success";
    case MK_FAILURE: return "Failure";
    case MK_NOT_FOUND: return "Not found";
    case MK_MULTIPLE_FOUND: return "Multiple entities found";
    case MK_MESHOP_NOT_FOUND: return "MeshOp not found";
    case MK_NOT_IMPLEMENTED: return "Not implemented";
    case MK_WRONG_DIMENSION: return "Wrong dimension";
    case MK_ALREADY_DEFINED: return "Already defined";
    case MK_BAD_INPUT: return "Bad input";
    case MK_BAD_GEOMETRIC_EVALUATION: return "Bad geometric evaluation";
    case MK_INCOMPLETE_MESH_SPECIFICATION: return "Incomplete mesh specification";
  };
}

} // namespace MeshKit

#endif
