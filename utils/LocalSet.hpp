#ifndef LOCALSET_HPP
#define LOCALSET_HPP

#include <iMesh.h>

#include "MKException.hpp"

class LocalSet
{
public:
  explicit LocalSet(iMesh_Instance impl, bool isList = false)
      : impl_(impl)
  {
    int err;
    iMesh_createEntSet(impl_, isList, &set_, &err);
    check_error(impl_, err);
  }

  ~LocalSet()
  {
    int err;
    iMesh_destroyEntSet(impl_, set_, &err);
  }

  operator iBase_EntitySetHandle()
  {
    return set_;
  }

  iMesh_Instance impl_;
  iBase_EntitySetHandle set_;
};

#endif
