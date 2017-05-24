
/*!
  \file   MsqFBiGeom.hpp
  \brief  Mesquite::MeshDomain implemented on ITAPS like FBiGeom API
  \author Iulian
  \date   2012-02-21
*/

#ifndef MSQ_FBIGEOM_HPP
#define MSQ_FBIGEOM_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "meshkit/FBiGeom.hpp"
#include "Vector3D.hpp"

using namespace MESQUITE_NS;

namespace MeshKit
{

/**\brief A Mesquite::MeshDomain implemented on top of the ITAPS FBiGeom API.
 *
 * Simple MeshDomain class implementation that queries a single iGeom
 * entity for all geometric queries.  Suitable for use when the entire
 * mesh to be smoothed lies on a single FB geometric surface.
 */
class MsqFBiGeom : public MeshDomain
{
public:

  MsqFBiGeom( FBiGeom * ifbigeom,
            iBase_EntityHandle geom_ent_handle );

  virtual ~MsqFBiGeom();

  void snap_to( Mesh::VertexHandle entity_handle,
      Vector3D& coordinat ) const;

  void vertex_normal_at( Mesh::VertexHandle entity_handle,
      Vector3D& coordinate ) const;

  void element_normal_at( MBMesquite::Mesh::ElementHandle entity_handle,
      Vector3D& coordinate ) const;
  
  void vertex_normal_at( const Mesh::VertexHandle* handles,
      Vector3D coordinates[],
                         unsigned count,
                         MsqError& err ) const;

  void closest_point( Mesh::VertexHandle handle,
                      const Vector3D& position,
                      Vector3D& closest,
                      Vector3D& normal,
                      MsqError& err ) const;

  void domain_DoF( const Mesh::VertexHandle* handle_array,
                   unsigned short* dof_array,
                   size_t num_vertices,
                   MsqError& err ) const;
private:
  
    /** A handle for the geometry entity to evaluate */
  iBase_EntityHandle geomEntHandle;
  /** FBiGeom instance used for MeshKit::iGeom-like calls */
  FBiGeom * fbigeom;
};

} // end namespace

#endif
