#ifndef MOVEMESH_HPP
#define MOVEMESH_HPP

#include "iMesh_extensions.h"
#include "Transform.hpp"

class MoveMesh
{
public:
  MoveMesh(iMesh_Instance impl);
  virtual ~MoveMesh();

  iMesh_Instance impl() const { return impl_; }

  void move(iBase_EntityHandle *ent_handles, int num_ents,
            const copy::Transform &trans);
  void move(iBase_EntitySetHandle set_handle, const copy::Transform &trans);
private:
  iMesh_Instance impl_;
};

#endif
