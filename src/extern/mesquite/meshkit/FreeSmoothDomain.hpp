#ifndef MESHKIT_FREE_SMOOTH_DOMAIN_HPP
#define MESHKIT_FREE_SMOOTH_DOMAIN_HPP

/**\file FreeSmoothDomain.hpp
 *\brief Implement MeshKit::FreeSmoothDomain
 */

#include "MKVersion.h"
#include "MsqIGeom.hpp"
#include "meshkit/Types.hpp"
#include "moab/Types.hpp"

namespace MeshKit {

class MKCore;

/**\brief Implement MESQUITE_NS::MeshDomain subclass for free smooth of
 *        mesh on geometry.
 */
class FreeSmoothDomain : public MESQUITE_NS::MsqCommonIGeom
{
  public:
    
    /**\param entities Top-level ModelEnts to operate on */
    FreeSmoothDomain( MKCore* core, const MEntVector& entities );
    
    virtual ~FreeSmoothDomain();
    
    void snap_to( MESQUITE_NS::Mesh::VertexHandle entity_handle,
                  MESQUITE_NS::Vector3D& coordinat ) const;

    void vertex_normal_at( MESQUITE_NS::Mesh::VertexHandle entity_handle,
                           MESQUITE_NS::Vector3D& coordinate ) const;

    void element_normal_at( MESQUITE_NS::Mesh::ElementHandle entity_handle,
                            MESQUITE_NS::Vector3D& coordinate ) const;

    void vertex_normal_at( const MESQUITE_NS::Mesh::VertexHandle* handles,
                           MESQUITE_NS::Vector3D coordinates[],
                           unsigned count,
                           MESQUITE_NS::MsqError& err ) const;

    void closest_point( MESQUITE_NS::Mesh::VertexHandle handle,
                        const MESQUITE_NS::Vector3D& position,
                        MESQUITE_NS::Vector3D& closest,
                        MESQUITE_NS::Vector3D& normal,
                        MESQUITE_NS::MsqError& err ) const;

    void domain_DoF( const MESQUITE_NS::Mesh::VertexHandle* handle_array,
                     unsigned short* dof_array,
                     size_t num_vertices,
                     MESQUITE_NS::MsqError& err ) const;
    

    /**\brief Given a mesh entity handle, get cooresponding geometry entity */
    iBase_EntityHandle get_geometry( MESQUITE_NS::Mesh::EntityHandle mesh_ent ) const;

    void get_geometry( const MESQUITE_NS::Mesh::EntityHandle* handles_array,
                       size_t num_handles,
                       iBase_EntityHandle* geom_array,
                       MESQUITE_NS::MsqError& err ) const;

  private:
    bool haveEntGeomRelTag;
    moab::Tag entGeomRel;
    moab::Interface* moabIface;
    mutable std::vector<iBase_EntityHandle> tmpHandles;

}; // class FreeSmoothDomain

} // namespace MeshKit

#endif
