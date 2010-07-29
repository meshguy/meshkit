#include "MoveMesh.hpp"

#include "CopyUtils.hpp"
#include "LocalSet.hpp"
#include "SimpleArray.hpp"

MoveMesh::MoveMesh(iMesh_Instance impl) : impl_(impl)
{}

MoveMesh::~MoveMesh()
{}

void MoveMesh::move(iBase_EntityHandle *ent_handles, int num_ents,
                    const copy::Transform &trans)
{
  int err;

  LocalSet set(impl_);
  
  iMesh_addEntArrToSet(impl_, ent_handles, num_ents, set, &err);
  check_error(impl_, err);

  move(set, trans);
}

void MoveMesh::move(iBase_EntitySetHandle set_handle,
                    const copy::Transform &trans)
{
  int err;

  SimpleArray<iBase_EntityHandle> ents;
  SimpleArray<iBase_EntityHandle> verts;
  SimpleArray<int> indices;
  SimpleArray<int> offsets;
  
  iMesh_getStructure(impl_, set_handle, ARRAY_INOUT(ents),
                     ARRAY_INOUT(verts), ARRAY_INOUT(indices),
                     ARRAY_INOUT(offsets), &err);
  check_error(impl_, err);

  trans(impl_, ARRAY_IN(verts), ARRAY_INOUT(verts));
}
