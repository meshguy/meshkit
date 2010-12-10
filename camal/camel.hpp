/**
 * Copyright 2006 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */
#ifndef CAMEL_H
#define CAMEL_H

#include <vector>

#include "iGeom.h"
#include "iMesh.h"
#include "iRel.h"
#include <math.h>


#define DOT(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define LENGTH_SQ(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define DIST2(a, b) sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) \
      + (b[2]-a[2])*(b[2]-a[2]))
#define NORMALIZE(a) {double lsq = LENGTH_SQ(a); \
                      if (lsq == 0.0) return 0; \
                      lsq = 1.0 / sqrt(lsq); \
                      for (int i = 0; i < 3; i++) a[i] *= lsq;}
#define CROSS(a, b, c) {c[0] = a[1]*b[2] - a[2]*b[1]; \
                        c[1] = a[2]*b[0] - a[0]*b[2]; \
                        c[2] = a[0]*b[1] - a[1]*b[0];}
const double PI = acos(-1.0);
const double TWO_PI = 2.0 * PI;

class CMEL 
{
public:
  CMEL( iGeom_Instance geom_ptr, 
        iMesh_Instance mesh_ptr = 0, 
        iRel_Instance relate_ptr = 0,
        iRel_PairHandle relation = 0 );
  
  bool mesh_geometry(double mesh_size, int mesh_intervals,
                     const bool force_intervals,
                     const bool quadMesh = false);
  
  bool mesh_entities(iBase_EntityHandle *gentities, const int num_geom_entities, 
                     double mesh_size, int mesh_intervals, const bool force_intervals,
                     const bool quadMesh = false);
  
  bool mesh_entity(iBase_EntityHandle gentity, double mesh_size,
                   int mesh_intervals, const bool force_intervals,
                   std::vector<iBase_EntityHandle> &new_entities,
                   const bool quadMesh = false);
  
  bool mesh_boundary(iBase_EntityHandle gentity, 
                     double mesh_size, int mesh_intervals, const bool force_intervals,
                     iBase_EntityHandle **bounding_ents = NULL, int *bounding_ent_size = NULL,
                     const bool quadMesh = false);
  
  bool mesh_vertex(iBase_EntityHandle gentity, 
                   std::vector<iBase_EntityHandle> &new_entities);
  
  bool mesh_curve(iBase_EntityHandle gentity, 
                  double mesh_size,
                  int mesh_intervals, 
                  std::vector<iBase_EntityHandle> &new_entities);
  
  bool is_meshed(iBase_EntityHandle gentity);

  bool assign_mesh(iBase_EntityHandle gentity, std::vector<iBase_EntityHandle> &mesh);

  bool create_vertices_elements(iBase_EntityHandle gentity, 
                                std::vector<iBase_EntityHandle> &bdy_verts, 
                                std::vector<double> &coords, std::vector<int> &connect, 
                                int ent_type,
                                std::vector<iBase_EntityHandle> &new_entities);
  
  bool create_vertices_elements(iBase_EntityHandle gentity, 
                                std::vector<iBase_EntityHandle> &bdy_verts, 
                                double *coords, const int num_points,
                                std::vector<int> &connect, 
                                int ent_type,
                                std::vector<iBase_EntityHandle> &new_entities);
  
  bool get_attrib_intervals(iBase_EntityHandle gentity,
                            double &mesh_size,
                            int &mesh_interval);
  
  bool bdy_coords_connect(std::vector<iBase_EntityHandle> &bdy_mesh, 
                          std::vector<int> &bdy_orientations, 
                          std::vector<iBase_EntityHandle> &bdy_verts,
                          std::vector<double> &coords, 
                          std::vector<int> &connect);

  bool bdy_elements_senses_grouped(iBase_EntityHandle gentity,
                                   std::vector<iBase_EntityHandle> &elements,
                                   std::vector<int> &element_senses,
                                   std::vector<int> &group_sizes);
  
  bool bdy_elements_senses(iBase_EntityHandle gentity,
                                std::vector<iBase_EntityHandle> & elements,
                                std::vector<int> & element_senses);

  bool bdy_nodes(iBase_EntityHandle gentity,
                                 std::vector<iBase_EntityHandle> &end_nodes);
  
    // assign a temporary indexing to the vertices in bdy_verts; pass back the tag handle
    // for that indexing
  bool assign_tmp_indices(std::vector<iBase_EntityHandle> &bdy_verts,
                          iBase_TagHandle &index_tag,
                          int offset = 0);

    // return id of gentity, useful for debugging
  int get_gentity_id(iBase_EntityHandle gentity);
  
    // given an entity, return the bounding groups (loops or shells)
  bool bdy_geom_grouped(iBase_EntityHandle gentity, 
                        std::vector<iBase_EntityHandle> &group_ents,
                        std::vector<int> &group_sizes);

  enum BooleanType {INTERSECT, UNION};
  
  bool get_adjs_bool(std::vector<iBase_EntityHandle> &from_ents,
                     int to_type,
                     std::vector<iBase_EntityHandle> &to_ents,
                     CMEL::BooleanType op_type);
  
  iBase_EntityHandle shared_entity(iBase_EntityHandle ent1,
                             iBase_EntityHandle ent2,
                             int dimension);
  
  iBase_EntityHandle next_winding(iBase_EntityHandle this_edge, 
                            iBase_EntityHandle gface, 
                            int this_sense, 
                            std::vector<iBase_EntityHandle> &tmp_adjs);
  
  bool trimSurface(const char * polygon_filename, int len);

  bool grounding_line(const char * polygon_filename, int len, double widthLeft,
        double widthRight);

  iGeom_Instance geomIface;
  iMesh_Instance meshIface;
  iRel_Instance relateIface;
  iRel_PairHandle pairHandle;

private:

  void print_meshed_entity( iBase_EntityHandle geom_entity,
                            int elem_type,
                            int num_points,
                            int num_conn );
                            
  
  bool get_mesh_set(iBase_EntityHandle gentity, iBase_EntitySetHandle &mesh_set,
                    bool create_if_missing);
  
  bool bdy_connect_coords(iBase_EntityHandle gentity,
                          std::vector<int> &connect_grps,
                          std::vector<double> &coords) {return false;}
                          
  bool get_both_sense(iBase_EntityHandle this_gent, 
                      iBase_EntityHandle other_gent,
                      bool prec_ent, int &this_sense);
  
  std::vector<double> newBoundary;

  // this would correspond to the grounding line; it will be meshed separately, and it will
  // divide the surface in 2 faces. it will have to be broken in an even number of edges, conforming
  // the internal boundary will start and end at the closest points to the new boundary. Or maybe the
  // new boundary will conform to the grounding line?
  // yes, maybe start from the grounding line.

  std::vector<double> internalBoundary;

  double _widthLeft, _widthRight; // width of the grounding line, left and right of grounding line

};

#endif

//#ifdef _cplusplus
//}
//#endif
