#include "MeshOpSet.hpp"
#include "meshkit/MeshOpProxy.hpp"
#include "meshkit/Error.hpp"
#include <string>

namespace MeshKit {

MeshOpSet::MeshOpSet() {}

MeshOpSet& MeshOpSet::instance()
{
  static MeshOpSet singleton;
  return singleton;
}

void MeshOpSet::register_mesh_op( MeshOpProxy* proxy )
{
    // Check if a MeshOp with the same name is already registered
  iterator other = mesh_op_no_throw( proxy->name() );
  if (other != allMeshOps.end()) {
    if (*other == proxy)
      return;
    else
      throw Error(MK_MULTIPLE_FOUND,"Conflicting MeshOp name: \"%s\"", proxy->name());
  }
  
    // Try to detect unterminated output entity lists
  const moab::EntityType* types = proxy->output_types();
  for (int i = 0; types[i] != moab::MBMAXTYPE; ++i) 
    if (i >= moab::MBMAXTYPE)
      throw Error(MK_BAD_INPUT,"Unterminated output type list for MeshOp: \"%s\"", proxy->name());
  
    // Add to list of all MeshOps
  allMeshOps.push_back(proxy);
  
    // Add to list for each dimension that it can mesh
  for (int i = 0; i < 4; ++i)
    if (proxy->can_mesh(static_cast<iBase_EntityType>(i)))
      dimMeshOps[i].push_back(proxy);
}

MeshOpSet::iterator MeshOpSet::mesh_op_no_throw( const char* op_name ) const
{
  std::string search(op_name);
  for (iterator i = allMeshOps.begin(); i != allMeshOps.end(); ++i) 
    if (search == (*i)->name())
      return i;
  return allMeshOps.end();
}

MeshOpProxy* MeshOpSet::mesh_op( const char* op_name ) const
{
  iterator i = mesh_op_no_throw( op_name );
  if (i == allMeshOps.end())
    throw Error(MK_NOT_FOUND,"Invalid MeshOp name: \"%s\"", op_name);
  return *i;
}

unsigned MeshOpSet::index( const char* op_name ) const
{
  iterator i = mesh_op_no_throw( op_name );
  if (i == allMeshOps.end())
    throw Error(MK_NOT_FOUND,"Invalid MeshOp name: \"%s\"", op_name);
  return i - allMeshOps.begin();
}

MeshOpProxy* MeshOpSet::mesh_op( unsigned index ) const
{
  if (index >= allMeshOps.size())
    new Error(MK_NOT_FOUND,"Invalid MeshOp index: %u", index);
    
  return allMeshOps[index];
}

} // namespace MeshKit
