#ifndef MESHKIT_LOCAL_SET_HPP
#define MESHKIT_LOCAL_SET_HPP

#include "meshkit/iMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/Error.hpp"

namespace MeshKit
{
  class LocalSet
  {
  public:
    explicit LocalSet(MKCore *mkCore, bool isList = false)
      : imesh_(mkCore->imesh_instance()->instance())
    {
      IBERRCHK(imesh_.createEntSet(isList, set_), "");
    }

    ~LocalSet()
    {
      imesh_.destroyEntSet(set_);
    }

    operator iBase_EntitySetHandle()
    {
      return set_;
    }

    iMesh imesh_;
    iMesh::EntitySetHandle set_;
  };
}

#endif
