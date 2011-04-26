#ifndef MESHKIT_MKEXCEPTION_HPP
#define MESHKIT_MKEXCEPTION_HPP

#include <exception>
#include <string>
#include <iMesh.h>
#include <iGeom.h>

class MKException : public std::exception
{
public:
  MKException(int code, const std::string &what) : code_(code), what_(what) {}
  virtual ~MKException() throw() {}

  int code() const { return code_; }
  virtual const char * what() const throw() { return what_.c_str(); }
private:
  int code_;
  std::string what_;
};

inline void check_error(iMesh_Instance instance, int err)
{
  if (err) {
    char descr[120];
    iMesh_getDescription(instance, descr, sizeof(descr)-1);
    throw MKException(err, descr);
  }
}

inline void check_error(iGeom_Instance instance, int err)
{
  if (err) {
    char descr[120];
    iGeom_getDescription(instance, descr, sizeof(descr)-1);
    throw MKException(err, descr);
  }
}

#endif
