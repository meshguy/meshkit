#ifndef LOCALTAG_HPP
#define LOCALTAG_HPP

#include "iMesh.hh"
#include "meshkit/MKCore.hpp"
#include "meshkit/Error.hpp"
#include <string>
#include <sstream>

#include "MKException.hpp"

namespace MeshKit
{
  class LocalTag
  {
  public:
    explicit LocalTag(MKCore *mkCore, int size = 1,
                       iMesh::TagValueType type = iBase_ENTITY_HANDLE)
            : imesh_(mkCore->imesh_instance()->instance())
    {
#ifdef MOAB
      init("", size, type);
#else
      init("local_tag", size, type);
#endif
    }

    LocalTag(MKCore *mkCore, const std::string &name, int size = 1,
             iMesh::TagValueType type = iBase_ENTITY_HANDLE)
      : imesh_(mkCore->imesh_instance()->instance())
    {
      init(name, size, type);
    }

    ~LocalTag()
    {
      imesh_.destroyTag(tag_handle_, true);
    }

    operator iMesh::TagHandle()
    {
      return tag_handle_;
    }

  private:
    void init(const std::string &name, int size, iMesh::TagValueType type)
    {
      iMesh::Error err;
      err = imesh_.createTag(name.c_str(), size, type, tag_handle_);

      for (int i=0; err == iBase_TAG_ALREADY_EXISTS; i++) {
        std::ostringstream os;
        os << name << i;
        std::string tag_name = os.str();
        err = imesh_.createTag(tag_name.c_str(), size, type, tag_handle_);
      }

      IBERRCHK(err, "FIXME");
    }

    iMesh imesh_;
    iMesh::TagHandle tag_handle_;
  };
}

#endif
