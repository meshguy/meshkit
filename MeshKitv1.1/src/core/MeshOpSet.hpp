#ifndef MESHKIT_MESH_OP_SET_HPP
#define MESHKIT_MESH_OP_SET_HPP

#include <vector>

namespace MeshKit
{

class MeshOpProxy;

/**\brief Maintain global, singletin list of registered MeshOps
 *
 * This class implements the list of registered MeshOps.  It
 * uses the singleton pattern to provide a single global list 
 * while avoiding issues with order of initialization of static
 * objects.
 *
 * This class is intended only for internal use in \c MKCore.  Access
 * to this data maintained by this class should be done through
 * static methods in the \c MKCore class.
 */
class MeshOpSet
{
public:
  /**\brief Get singleton instance */
  static MeshOpSet& instance();
  
  /**\brief Register a mesh op. 
   *
   * Register a new MeshOp.  Will fail upon duplicate
   * names unless proxy pointer is also the same (duplicate
   * registration).
   *
   *\param proxy Proxy for MeshOp sub-class.
   */
  void register_mesh_op( MeshOpProxy* proxy );
  
  /**\brief Type of list returned by reference from member methods */
  typedef std::vector<MeshOpProxy*> OpList;
  
  /**\brief Type of iterator to use with /c OpList */
  typedef OpList::const_iterator iterator;
  
  /**\brief Get list of all MeshOps that can be used to
   *        generate mesh entities of the specified dimension.
   *\param dimension Dimension of entity to be meshed
   *\return List of MeshOps that can mesh entities of the specified dimesion.
   */
  const OpList& mesh_ops( unsigned dimension ) const
    { return dimMeshOps[dimension]; }
  
  /**\brief Get list of all mesh ops
   *
   *\return List of all registered MeshOps
   */
  const OpList& mesh_ops() const
    { return allMeshOps; }
  
  /**\brief Get MeshOpProxy by name
   * 
   * Get MeshOpProxy by name.  Throws exception if not found.
   *
   *\param name MeshOp class name.
   *\return MeshOpProxy for MeshOp with the passed name
   */
  MeshOpProxy* mesh_op( const char* name ) const;
  
  /**\brief Get index of MeshOpProxy
   * 
   * Get index of MeshOpProxy.  Throws exception if not found.
   *
   *\param name MeshOp class name.
   *\return index of MeshOpProxy for MeshOp with the passed name
   */
  unsigned index( const char* name ) const;
  
  /**\brief Get MeshOpProxy by index
   * 
   * Get MeshOpProxy by index.  Throws exception if not found.
   *
   *\param index MeshOp index.
   *\return the MeshOpProxy
   */
  MeshOpProxy* mesh_op( unsigned index ) const;

private:

  /**\brief Private constructor for singleton pattern */
  MeshOpSet();

  /**\brief Get MeshOpProxy by name
   *
   * Returns \c allMeshOps.end() if not found.
   *
   *\param name MeshOp class name.
   *\return iterator into allMeshOps
   */
  iterator mesh_op_no_throw( const char* name ) const;
   

  /**\brief List of all registered MeshOps */
  OpList allMeshOps;
  /**\brief Lists of all registered indexed by dimension of generated entities */
  OpList dimMeshOps[4];
};

} // namespace MeshKit

#endif // #ifdef MESHKIT_MESH_OP_SET_HPP
