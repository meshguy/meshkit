#ifndef LOCALTAG_HPP
#define LOCALTAG_HPP

#include <iMesh.h>
#include <string>
#include <sstream>

#include "MKException.hpp"

class LocalTag
{
public:
  explicit LocalTag(iMesh_Instance impl, int size = 1,
                    int type = iBase_ENTITY_HANDLE) : impl_(impl)
  {
#ifdef MOAB
    init("", size, type);
#else
    init("local_tag", size, type);
#endif
  }

  LocalTag(iMesh_Instance impl, const std::string &name, int size = 1,
           int type = iBase_ENTITY_HANDLE) : impl_(impl)
  {
    init(name, size, type);
  }

  ~LocalTag()
  {
    int err;
    iMesh_destroyTag(impl_, tag_handle_, true, &err);
  }

  operator iBase_TagHandle()
  {
    return tag_handle_;
  }

private:
  void init(const std::string &name, int size, int type)
  {
    int err;
    iMesh_createTag(impl_, name.c_str(), 1, iBase_ENTITY_HANDLE,
                    &tag_handle_, &err, name.size());

	for (int i=0; err == iBase_TAG_ALREADY_EXISTS; i++) {
      std::ostringstream os;
      os << name << i;
      std::string tag_name = os.str();
      iMesh_createTag(impl_, tag_name.c_str(), 1, iBase_ENTITY_HANDLE,
                      &tag_handle_, &err, tag_name.size());
    }

    check_error(impl_, err);
  }

  iMesh_Instance impl_;
  iBase_TagHandle tag_handle_;
};

#endif
