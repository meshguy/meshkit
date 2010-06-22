#ifndef LOCALTAG_HPP
#define LOCALTAG_HPP

#include <iMesh.h>
#include <string>
#include <sstream>
#include <iostream>

// Stuff like this belongs somewhere common, but it's not a big deal right now.
#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}

class LocalTag
{
public:
  explicit LocalTag(iMesh_Instance impl,const std::string &name = "local_tag")
      : impl_(impl)
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

    ERROR("Failed to create local tag");
    // TODO: throw an exception here
  }

  ~LocalTag()
  {
    int err;
    iMesh_destroyTag(impl_, tag_handle_, true, &err);
    ERROR("Failed to force-destroy local tag.");
  }

  operator iBase_TagHandle()
  {
    return tag_handle_;
  }

  iMesh_Instance impl_;
  iBase_TagHandle tag_handle_;
};

#endif
