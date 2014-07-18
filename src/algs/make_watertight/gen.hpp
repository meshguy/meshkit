#ifndef GEN_HPP
#define GEN_HPP

#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp" // for merging verts
#include "MBCartVect.hpp"

// SENSE CONVENTIONS
#define SENSE_FORWARD 1
#define SENSE_REVERSE -1
#define SENSE_UNKNOWN 0


MBInterface *MBI(); 
namespace gen {
  bool error( const bool error_has_occured, const std::string message="" );
  
  /// prints a string of the error code returned from MOAB to standard output 
  void moab_printer(MBErrorCode error_code);

  void print_vertex_cubit( const MBEntityHandle vertex );

  void print_vertex_coords( const MBEntityHandle vertex );

  void print_triangles( const MBRange tris );

  void print_triangle( const MBEntityHandle triangle, bool print_edges );

  void print_edge( const MBEntityHandle edge );

  void print_vertex_count( const MBEntityHandle input_meshset);

  void print_range( const MBRange range );

  void print_range_of_edges( const MBRange range );

  void print_arc_of_edges( const std::vector<MBEntityHandle> arc_of_edges );

  void print_arcs( const std::vector < std::vector<MBEntityHandle> > arcs );

  void print_loop( const std::vector<MBEntityHandle> loop_of_verts );

MBErrorCode find_closest_vert( const MBEntityHandle reference_vert,
                               const std::vector<MBEntityHandle> arc_of_verts,
                               unsigned &position,
                               const double dist_limit );

  MBErrorCode find_closest_vert( const double tol,
                                 const MBEntityHandle reference_vert,
                                 const std::vector<MBEntityHandle> loop_of_verts,
                                 std::vector<unsigned> &positions, 
                                 std::vector<double>   &dists);
  /*  MBErrorCode find_closest_vert( const MBEntityHandle reference_vert,
                                  const std::vector<std::vector<MBEntityHandle> > loops_of_verts,
                                  unsigned int &loop, unsigned int &position, 
                                  double &min_dist);
  */
  // Merge the range of vertices. We do not merge by edges (more
  // stringent) because we do not want to miss corner vertices.

/// finds any vertices within the MBRange vertices that are with in tol of 
/// each other and merges them
  MBErrorCode merge_vertices( MBRange vertices /* in */, 
			      const  double tol       /* in */);
			      //bool &merge_vertices_again /* out */);

/// returns the square of the distance between v0 and v1
  MBErrorCode squared_dist_between_verts( const MBEntityHandle v0, 
                                          const MBEntityHandle v1, 
                                          double &d);
/// returns the distance between v0 and v1

  double dist_between_verts( const MBCartVect v0, const MBCartVect v1 );
  MBErrorCode dist_between_verts( const MBEntityHandle v0, const MBEntityHandle v1,
                                  double &d );
  double dist_between_verts( double coords0[], double coords1[] );
  double dist_between_verts( MBEntityHandle vert0, MBEntityHandle vert1 );                             

  // Return the length of the curve defined by MBEDGEs or ordered MBVERTEXs.
  double length( std::vector<MBEntityHandle> curve );

  // Given a vertex and vector of edges, return the number of edges adjacent to the vertex.
  unsigned int n_adj_edges( MBEntityHandle vert, MBRange edges );

// Return true if the edges share a vertex. Does not check for coincident edges.
  bool edges_adjacent( MBEntityHandle edge0, MBEntityHandle edge1 );

// get the direction unit vector from one vertex to another vertex
  //MBErrorCode get_direction( MBEntityHandle from_vert, MBEntityHandle to_vert, double dir[] );
  MBErrorCode get_direction( const MBEntityHandle from_vert, const MBEntityHandle to_vert,
                           MBCartVect &dir ); 

// from http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry1
  double edge_point_dist( const MBCartVect a, const MBCartVect b, const MBCartVect c );
  double edge_point_dist( const MBEntityHandle endpt0, const MBEntityHandle endpt1, 
                          const MBEntityHandle pt );
  double edge_point_dist( const MBEntityHandle edge, const MBEntityHandle pt );

  MBErrorCode point_curve_min_dist( const std::vector<MBEntityHandle> curve, 
                                    const MBEntityHandle pt, double &min_dist,
                                    const double max_dist_along_curve );
  MBErrorCode point_curve_min_dist( const std::vector<MBEntityHandle> curve, 
                                    const MBEntityHandle pt, double &min_dist );

  double triangle_area( const MBCartVect a, const MBCartVect b, const MBCartVect c);
  MBErrorCode triangle_area( const MBEntityHandle conn[], double &area );
  MBErrorCode triangle_area( const MBEntityHandle triangle, double &area );
  double triangle_area( MBRange triangles );
  
  bool triangle_degenerate( const MBEntityHandle triangle );
  bool triangle_degenerate( const MBEntityHandle v0, const MBEntityHandle v1, const MBEntityHandle v2);

/// gets the normal vectors of all triangles in MBRange triangles and returns them as MBCartVect's in normals
  MBErrorCode triangle_normals( const MBRange triangles, std::vector<MBCartVect> &normals );
  MBErrorCode triangle_normal( const MBEntityHandle triangle, MBCartVect &normal );
  MBErrorCode triangle_normal( const MBEntityHandle v0, const MBEntityHandle v1,
                               const MBEntityHandle v2, MBCartVect &normal );
  MBErrorCode triangle_normal( const MBCartVect v0, const MBCartVect v1, 
                               const MBCartVect v2, MBCartVect &normal ); 

  // Distance between a point and line. The line is defined by two verts.
  // We are using a line and not a line segment!
  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  MBErrorCode line_point_dist( const MBEntityHandle line_pt1, const MBEntityHandle line_pt2,
                               const MBEntityHandle pt0, double &dist );

  // Project the point onto the line. Not the line segment! 
  // Change the coordinates of the pt0 to the projection.
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,     
                                     const MBEntityHandle line_pt2,             
                                     const MBEntityHandle pt0 );

  // Do not change the coords of pt0. Instead return the projected coords.              
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,     
                                     const MBEntityHandle line_pt2,             
                                     const MBEntityHandle pt0,              
                                     MBCartVect &projected_coords,
                                     double &parameter );  
  // Get the distance of pt0 from line_pt1 if pt0 is projected. No coords are
  // changed in MOAB.
  MBErrorCode point_line_projection( const MBEntityHandle line_pt1,
				     const MBEntityHandle line_pt2,
				     const MBEntityHandle pt0,
				     double &dist_along_edge  );
  
  MBErrorCode ear_clip_polygon( std::vector<MBEntityHandle> polygon_of_verts,
                                    const MBCartVect plane_normal_vector, MBRange &new_tris );

  int geom_id_by_handle( const MBEntityHandle set);

/// gets the normal vector of each triangle in tris and tags each triangle with its normal vector
  MBErrorCode save_normals( MBRange tris, MBTag normal_tag );

  MBErrorCode flip(const MBEntityHandle tri, const MBEntityHandle vert0, 
                   const MBEntityHandle vert2, const MBEntityHandle surf_set);


/// creates a set of ordered verts from the a set of ordered edges. uses commone vertex between edges to check continuity. 
  MBErrorCode ordered_verts_from_ordered_edges( const std::vector<MBEntityHandle> ordered_edges,
                                                std::vector<MBEntityHandle> &ordered_verts );
/// returns the average distance between the vertices in arc0 and arc1. assumes the vertices
/// are ordered appropriately
  MBErrorCode dist_between_arcs( bool debug,
                         const std::vector<MBEntityHandle> arc0,
                         const std::vector<MBEntityHandle> arc1,
                         double &dist );

  // skin edges are a vector of two vertex handles
    // Hold edges in an array of handles.
  struct edge {
    MBEntityHandle edge, v0, v1;
  };
  int compare_edge(const void *a, const void *b);
  MBErrorCode find_skin( MBRange tris, const int dim,                     
			 // std::vector<std::vector<MBEntityHandle> > &skin_edges,    
			 MBRange &skin_edges,                         
                         const bool );
  //MBErrorCode find_skin( const MBRange tris, const int dim, MBRange &skin_edges, const bool );
  MBErrorCode measure( const MBEntityHandle set, const MBTag geom_tag, double &size, bool debug, bool verbose );

  // Given a curve and surface set, get the relative sense.
  // From CGMA/builds/dbg/include/CubitDefines, CUBIT_UNKNOWN=-1, CUBIT_FORWARD=0, CUBIT_REVERSED=1
  MBErrorCode get_curve_surf_sense( const MBEntityHandle surf_set, const MBEntityHandle curve_set,
                                    int &sense, bool debug = false );

  MBErrorCode surface_sense( MBEntityHandle volume, int num_surfaces,const MBEntityHandle* surfaces,int* senses_out );
  MBErrorCode surface_sense( MBEntityHandle volume, MBEntityHandle surface, int& sense_out );

  MBTag get_tag( const char* name, int size, MBTagType store,MBDataType type, const void* def_value, bool create_if_missing);

  MBErrorCode delete_surface( MBEntityHandle surf, MBTag geom_tag, MBRange tris, int id, bool debug, bool verbose);

  MBErrorCode remove_surf_sense_data(MBEntityHandle del_surf, bool debug);

  MBErrorCode combine_merged_curve_senses( std::vector<MBEntityHandle> &curves, MBTag merge_tag, bool debug = false) ;

 /// used to get all mesh tags necessary for sealing a mesh
  MBErrorCode get_sealing_mesh_tags( double &facet_tol,
                             double &sme_resabs_tol,
                             MBTag &geom_tag, 
                             MBTag &id_tag, 
                             MBTag &normal_tag, 
                             MBTag &merge_tag, 
                             MBTag &faceting_tol_tag, 
                             MBTag &geometry_resabs_tag, 
                             MBTag &size_tag, 
                             MBTag &orig_curve_tag);

 /// sets the tracking and ordering options of meshsets retrieved from the mesh
  MBErrorCode get_geometry_meshsets( MBRange geometry_sets[], MBTag geom_tag, bool verbose = false);

 /// returns MB_SUCCESS if there exists geometry sets of each dimension in the model
 /// returns MB_FAILURE if there are no geometry sets of any dimension in the model
  MBErrorCode check_for_geometry_sets(MBTag geom_tag, bool verbose);

}

#endif
