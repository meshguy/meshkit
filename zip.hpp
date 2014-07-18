#ifndef ZIP_HPP
#define ZIP_HPP

#include "MBCore.hpp"
#include "gen.hpp"
#include "arc.hpp"

MBInterface *MBI();
namespace zip {
  MBErrorCode t_joint( MBTag normal_tag, 
                       const MBEntityHandle vert0,              
                       const MBEntityHandle vert1,                         
                       const MBEntityHandle vert2,
                       bool debug );
/// removes the entitiy handle tri from the loaded mesh                
  MBErrorCode delete_degenerate_tris( MBEntityHandle tri );
/// checks that no triangles in the MBRange tris are degenterate. If 
/// degenerate triangles are found, they are deleted from the mesh. 
  MBErrorCode delete_degenerate_tris( MBRange tris );

  MBErrorCode delete_adj_degenerate_tris( const MBEntityHandle adj_vert );


/// merges two vertices by updating the entity handle of the deleted vert to the vert to 
/// keep in the correct arc. Uses MOAB function merge_entities to merge the vertices in 
/// the database. Also deletes the triangles adjacent to the merged vertices if one
/// becomes degenerate.
  MBErrorCode merge_verts( const MBEntityHandle keep_vert, 
                           const MBEntityHandle delete_vert,
                           std::vector<MBEntityHandle> &arc0,
                           std::vector<MBEntityHandle> &arc1 );

/// test two normal vectors to see if they point in the same direction
  MBErrorCode test_normals( const std::vector<MBCartVect> norms0, 
                            const std::vector<MBCartVect> norms1,
                            std::vector<int> &inverted_tri_indices );
  MBErrorCode test_normals( const             MBCartVect  norms0, 
                            const             MBCartVect  norms1 );

  MBErrorCode remove_inverted_tris(MBTag normal_tag, MBRange tris, const bool debug );

/// tests the watertightness of all arcs in the vector-array of moab entity handles arcs
  MBErrorCode test_zipping( const double FACET_TOL,
                            const std::vector< std::vector<MBEntityHandle> > arcs );

}

#endif
