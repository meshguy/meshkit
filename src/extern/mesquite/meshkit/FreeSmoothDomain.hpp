#ifndef MESHKIT_FREE_SMOOTH_DOMAIN_HPP
#define MESHKIT_FREE_SMOOTH_DOMAIN_HPP

/**\file FreeSmoothDomain.hpp
 *\brief Implement MeshKit::FreeSmoothDomain
 */

#include "MsqIGeom.hpp"
#include "meshkit/Types.hpp"
#include "moab/Types.hpp"

namespace MeshKit {

class MKCore;

/**\brief Implement Mesquite::MeshDomain subclass for free smooth of
 *        mesh on geometry.
 */
class FreeSmoothDomain : public Mesquite::MsqCommonIGeom
{
  public:
    
    /**\param entities Top-level ModelEnts to operate on */
    FreeSmoothDomain( MKCore* core, const MEntVector& entities );
    
    virtual ~FreeSmoothDomain();
    
    void snap_to( Mesquite::Mesh::VertexHandle entity_handle,
                  Mesquite::Vector3D& coordinat ) const;

    void vertex_normal_at( Mesquite::Mesh::VertexHandle entity_handle,
                           Mesquite::Vector3D& coordinate ) const;

    void element_normal_at( Mesquite::Mesh::ElementHandle entity_handle,
                            Mesquite::Vector3D& coordinate ) const;

    void vertex_normal_at( const Mesquite::Mesh::VertexHandle* handles,
                           Mesquite::Vector3D coordinates[],
                           unsigned count,
                           Mesquite::MsqError& err ) const;

    void closest_point( Mesquite::Mesh::VertexHandle handle,
                        const Mesquite::Vector3D& position,
                        Mesquite::Vector3D& closest,
                        Mesquite::Vector3D& normal,
                        Mesquite::MsqError& err ) const;

    void domain_DoF( const Mesquite::Mesh::VertexHandle* handle_array,
                     unsigned short* dof_array,
                     size_t num_vertices,
                     Mesquite::MsqError& err ) const;
    

    /**\brief Given a mesh entity handle, get cooresponding geometry entity */
    iBase_EntityHandle get_geometry( Mesquite::Mesh::EntityHandle mesh_ent ) const;

    void get_geometry( const Mesquite::Mesh::EntityHandle* handles_array,
                       size_t num_handles,
                       iBase_EntityHandle* geom_array,
                       Mesquite::MsqError& err ) const;

  private:
    bool haveEntGeomRelTag;
    moab::Tag entGeomRel;
    moab::Interface* moabIface;
    mutable std::vector<iBase_EntityHandle> tmpHandles;

}; // class FreeSmoothDomain

} // namespace MeshKit

#endif
