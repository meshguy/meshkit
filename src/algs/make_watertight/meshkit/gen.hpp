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
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp" // for merging verts
#include "MBCartVect.hpp"

// SENSE CONVENTIONS
#define SENSE_FORWARD 1
#define SENSE_REVERSE -1
#define SENSE_UNKNOWN 0


// returns a moab interface
MBInterface *MBI(); 

namespace gen {
  bool error( const bool error_has_occured, const std::string message="" );
  
  /// prints a string of the error code returned from MOAB to standard output 
  void moab_printer(moab::ErrorCode error_code);

  void print_vertex_cubit( const moab::EntityHandle vertex );

  void print_vertex_coords( const moab::EntityHandle vertex );

  void print_triangles( const moab::Range tris );

  void print_triangle( const moab::EntityHandle triangle, bool print_edges );

  void print_edge( const moab::EntityHandle edge );

  void print_vertex_count( const moab::EntityHandle input_meshset);

  void print_range( const moab::Range range );

  void print_range_of_edges( const moab::Range range );

  void print_arc_of_edges( const std::vector<moab::EntityHandle> arc_of_edges );

  void print_arcs( const std::vector < std::vector<moab::EntityHandle> > arcs );

  void print_loop( const std::vector<moab::EntityHandle> loop_of_verts );

moab::ErrorCode find_closest_vert( const moab::EntityHandle reference_vert,
                               const std::vector<moab::EntityHandle> arc_of_verts,
                               unsigned &position,
                               const double dist_limit );

  moab::ErrorCode find_closest_vert( const double tol,
                                 const moab::EntityHandle reference_vert,
                                 const std::vector<moab::EntityHandle> loop_of_verts,
                                 std::vector<unsigned> &positions, 
                                 std::vector<double>   &dists);
  /*  moab::ErrorCode find_closest_vert( const moab::EntityHandle reference_vert,
                                  const std::vector<std::vector<moab::EntityHandle> > loops_of_verts,
                                  unsigned int &loop, unsigned int &position, 
                                  double &min_dist);
  */
  // Merge the range of vertices. We do not merge by edges (more
  // stringent) because we do not want to miss corner vertices.

/// finds any vertices within the moab::Range vertices that are with in tol of 
/// each other and merges them
  moab::ErrorCode merge_vertices( moab::Range vertices /* in */, 
			      const  double tol       /* in */);
			      //bool &merge_vertices_again /* out */);

/// returns the square of the distance between v0 and v1
  moab::ErrorCode squared_dist_between_verts( const moab::EntityHandle v0, 
                                          const moab::EntityHandle v1, 
                                          double &d);
/// returns the distance between v0 and v1

  double dist_between_verts( const MBCartVect v0, const MBCartVect v1 );
  moab::ErrorCode dist_between_verts( const moab::EntityHandle v0, const moab::EntityHandle v1,
                                  double &d );
  double dist_between_verts( double coords0[], double coords1[] );
  double dist_between_verts( moab::EntityHandle vert0, moab::EntityHandle vert1 );                             

  // Return the length of the curve defined by MBEDGEs or ordered MBVERTEXs.
  double length( std::vector<moab::EntityHandle> curve );

  // Given a vertex and vector of edges, return the number of edges adjacent to the vertex.
  unsigned int n_adj_edges( moab::EntityHandle vert, moab::Range edges );

// Return true if the edges share a vertex. Does not check for coincident edges.
  bool edges_adjacent( moab::EntityHandle edge0, moab::EntityHandle edge1 );

// get the direction unit vector from one vertex to another vertex
  //moab::ErrorCode get_direction( moab::EntityHandle from_vert, moab::EntityHandle to_vert, double dir[] );
  moab::ErrorCode get_direction( const moab::EntityHandle from_vert, const moab::EntityHandle to_vert,
                           MBCartVect &dir ); 

// from http://www.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry1
  double edge_point_dist( const MBCartVect a, const MBCartVect b, const MBCartVect c );
  double edge_point_dist( const moab::EntityHandle endpt0, const moab::EntityHandle endpt1, 
                          const moab::EntityHandle pt );
  double edge_point_dist( const moab::EntityHandle edge, const moab::EntityHandle pt );

  moab::ErrorCode point_curve_min_dist( const std::vector<moab::EntityHandle> curve, 
                                    const moab::EntityHandle pt, double &min_dist,
                                    const double max_dist_along_curve );
  moab::ErrorCode point_curve_min_dist( const std::vector<moab::EntityHandle> curve, 
                                    const moab::EntityHandle pt, double &min_dist );

  double triangle_area( const MBCartVect a, const MBCartVect b, const MBCartVect c);
  moab::ErrorCode triangle_area( const moab::EntityHandle conn[], double &area );
  moab::ErrorCode triangle_area( const moab::EntityHandle triangle, double &area );
  double triangle_area( moab::Range triangles );
  
  bool triangle_degenerate( const moab::EntityHandle triangle );
  bool triangle_degenerate( const moab::EntityHandle v0, const moab::EntityHandle v1, const moab::EntityHandle v2);

/// gets the normal vectors of all triangles in moab::Range triangles and returns them as MBCartVect's in normals
  moab::ErrorCode triangle_normals( const moab::Range triangles, std::vector<MBCartVect> &normals );
  moab::ErrorCode triangle_normal( const moab::EntityHandle triangle, MBCartVect &normal );
  moab::ErrorCode triangle_normal( const moab::EntityHandle v0, const moab::EntityHandle v1,
                               const moab::EntityHandle v2, MBCartVect &normal );
  moab::ErrorCode triangle_normal( const MBCartVect v0, const MBCartVect v1, 
                               const MBCartVect v2, MBCartVect &normal ); 

  // Distance between a point and line. The line is defined by two verts.
  // We are using a line and not a line segment!
  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  moab::ErrorCode line_point_dist( const moab::EntityHandle line_pt1, const moab::EntityHandle line_pt2,
                               const moab::EntityHandle pt0, double &dist );

  // Project the point onto the line. Not the line segment! 
  // Change the coordinates of the pt0 to the projection.
  moab::ErrorCode point_line_projection( const moab::EntityHandle line_pt1,     
                                     const moab::EntityHandle line_pt2,             
                                     const moab::EntityHandle pt0 );

  // Do not change the coords of pt0. Instead return the projected coords.              
  moab::ErrorCode point_line_projection( const moab::EntityHandle line_pt1,     
                                     const moab::EntityHandle line_pt2,             
                                     const moab::EntityHandle pt0,              
                                     MBCartVect &projected_coords,
                                     double &parameter );  
  // Get the distance of pt0 from line_pt1 if pt0 is projected. No coords are
  // changed in MOAB.
  moab::ErrorCode point_line_projection( const moab::EntityHandle line_pt1,
				     const moab::EntityHandle line_pt2,
				     const moab::EntityHandle pt0,
				     double &dist_along_edge  );
  
  moab::ErrorCode ear_clip_polygon( std::vector<moab::EntityHandle> polygon_of_verts,
                                    const MBCartVect plane_normal_vector, moab::Range &new_tris );

  int geom_id_by_handle( const moab::EntityHandle set);

/// gets the normal vector of each triangle in tris and tags each triangle with its normal vector
  moab::ErrorCode save_normals( moab::Range tris, moab::Tag normal_tag );

  moab::ErrorCode flip(const moab::EntityHandle tri, const moab::EntityHandle vert0, 
                   const moab::EntityHandle vert2, const moab::EntityHandle surf_set);


/// creates a set of ordered verts from the a set of ordered edges. uses commone vertex between edges to check continuity. 
  moab::ErrorCode ordered_verts_from_ordered_edges( const std::vector<moab::EntityHandle> ordered_edges,
                                                std::vector<moab::EntityHandle> &ordered_verts );
/// returns the average distance between the vertices in arc0 and arc1. assumes the vertices
/// are ordered appropriately
  moab::ErrorCode dist_between_arcs( bool debug,
                         const std::vector<moab::EntityHandle> arc0,
                         const std::vector<moab::EntityHandle> arc1,
                         double &dist );

  // skin edges are a vector of two vertex handles
    // Hold edges in an array of handles.
  struct edge {
    moab::EntityHandle edge, v0, v1;
  };
  int compare_edge(const void *a, const void *b);
  moab::ErrorCode find_skin( moab::Range tris, const int dim,                     
			 // std::vector<std::vector<moab::EntityHandle> > &skin_edges,    
			 moab::Range &skin_edges,                         
                         const bool );
  //moab::ErrorCode find_skin( const moab::Range tris, const int dim, moab::Range &skin_edges, const bool );
  moab::ErrorCode measure( const moab::EntityHandle set, const moab::Tag geom_tag, double &size, bool debug, bool verbose );

  // Given a curve and surface set, get the relative sense.
  // From CGMA/builds/dbg/include/CubitDefines, CUBIT_UNKNOWN=-1, CUBIT_FORWARD=0, CUBIT_REVERSED=1
  moab::ErrorCode get_curve_surf_sense( const moab::EntityHandle surf_set, const moab::EntityHandle curve_set,
                                    int &sense, bool debug = false );

  moab::ErrorCode surface_sense( moab::EntityHandle volume, int num_surfaces,const moab::EntityHandle* surfaces,int* senses_out );
  moab::ErrorCode surface_sense( moab::EntityHandle volume, moab::EntityHandle surface, int& sense_out );

  moab::Tag get_tag( const char* name, int size, moab::TagType store,MBDataType type, const void* def_value, bool create_if_missing);

  moab::ErrorCode delete_surface( moab::EntityHandle surf, moab::Tag geom_tag, moab::Range tris, int id, bool debug, bool verbose);

  moab::ErrorCode remove_surf_sense_data(moab::EntityHandle del_surf, bool debug);

  moab::ErrorCode combine_merged_curve_senses( std::vector<moab::EntityHandle> &curves, moab::Tag merge_tag, bool debug = false) ;

 /// used to get all mesh tags necessary for sealing a mesh
  moab::ErrorCode get_sealing_mesh_tags( double &facet_tol,
                             double &sme_resabs_tol,
                             moab::Tag &geom_tag, 
                             moab::Tag &id_tag, 
                             moab::Tag &normal_tag, 
                             moab::Tag &merge_tag, 
                             moab::Tag &faceting_tol_tag, 
                             moab::Tag &geometry_resabs_tag, 
                             moab::Tag &size_tag, 
                             moab::Tag &orig_curve_tag);

 /// sets the tracking and ordering options of meshsets retrieved from the mesh
  moab::ErrorCode get_geometry_meshsets( moab::Range geometry_sets[], moab::Tag geom_tag, bool verbose = false);

 /// returns moab::MB_SUCCESS if there exists geometry sets of each dimension in the model
 /// returns MB_FAILURE if there are no geometry sets of any dimension in the model
  moab::ErrorCode check_for_geometry_sets(moab::Tag geom_tag, bool verbose);

}

#endif
