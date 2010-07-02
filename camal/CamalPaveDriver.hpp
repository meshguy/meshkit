
#ifndef CamalPaveDriver_H
#define CamalPaveDriver_H

#include <vector>

#include "iMesh.h"


class CamalPaveDriver
{
public:
  CamalPaveDriver( iMesh_Instance input_mesh, iBase_EntitySetHandle set, 
        iMesh_Instance output_mesh);
  
  bool remesh(double mesh_size, int mesh_intervals,
                     const bool force_intervals);
  
  iMesh_Instance _meshIface;
  iMesh_Instance _meshOutface;
  iBase_EntitySetHandle _set;
#if 0
  bool mesh_entities(iBase_EntityHandle *gentities, const int num_geom_entities, 
                     double mesh_size, int mesh_intervals, const bool force_intervals);
  
  bool mesh_entity(iBase_EntityHandle gentity, double mesh_size,
                   int mesh_intervals, const bool force_intervals,
                   std::vector<iBase_EntityHandle> &new_entities);
  
  bool mesh_boundary(iBase_EntityHandle gentity, 
                     double mesh_size, int mesh_intervals, const bool force_intervals,
                     iBase_EntityHandle **bounding_ents = NULL, int *bounding_ent_size = NULL);
  
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
                           std::vector<iBase_EntityHandle> &elements,
                           std::vector<int> &element_senses);
  
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
  
  iGeom_Instance geomIface;
  iMesh_Instance meshIface;
  iRel_Instance relateIface;
  iRel_RelationHandle relationHandle;

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
#endif
  
};

#endif  // CamalPaveDriver_H

