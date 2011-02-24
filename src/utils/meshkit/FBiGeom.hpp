#ifndef FBIGEOM_HH
#define FBIGEOM_HH

/** \file FBiGeom.hpp;
 *
 */

#include "iGeom.hpp"
#include "moab/FBEngine.hpp"


/** \class FBiGeom
 * \brief C++ interface for ITAPS iGeom interface
 *
 * This class is a simple wrapper for the ITAPS iGeom interface,
 * derived from our iGeom wrapper class
 * The primary benefit to using this class is that it will use mesh - based geometry defined in
 * FBEngine class.  This file includes both declaration and definition of all iGeom class functions, i.e.
 * all functions are inlined; but they are also virtual, so not much benefit there.
 * The class can be constructed and destructed in the standard C++ way.
 *
 * For complete documentation of these functions, see the iGeom header in the CGM source
 * (http://trac.mcs.anl.gov/projects/ITAPS/browser/cgm/trunk/itaps/iGeom.h for now).
 */

namespace MeshKit {

class MKCore;

class FBiGeom : public iGeom {
  public:

    //inline FBiGeom( iGeom_Instance instance );
    
    FBiGeom( MKCore * mk, bool smooth = false);

    ~FBiGeom();
    
    Error Init();

    Error load( const char* file_name,
                       const char* options = 0 );

    int meshkit_index() { return _index;}

    // methods from iBaseVirtual that are overloaded

    inline iGeomBase::EntitySetHandle getRootSet();

    // methods from iGeom
    inline Error getEntities( EntitySetHandle set,
                                  EntityType type,
                                  std::vector<EntityHandle>& entities_out );


    // methods from iBase
    iGeom::Error getTagHandle( const char* name, TagHandle& handle_out );

    iGeom::Error createTag( const char* tag_name,
                      int tag_num_type_values,
                      TagValueType tag_type,
                      TagHandle& tag_handle_out );

    iGeom::Error getData( EntityHandle entity_handle,
                    TagHandle tag_handle,
                    void* tag_value_out );

    iGeom::Error getArrData( const EntityHandle* entity_handles,
                       int entity_handles_size,
                       TagHandle tag_handle,
                       void* tag_values_out );

    iGeom::Error setData( EntityHandle entity_handle,
                    TagHandle tag_handle,
                    const void* tag_value );

    iGeom::Error setArrData( const EntityHandle* entity_handles,
                       int entity_handles_size,
                       TagHandle tag_handle,
                       const void* tag_values );

    iGeom::Error getEntType( EntityHandle handle, EntityType& type_out );

    iGeom::Error getEntAdj( EntityHandle handle,
                                EntityType type_requested,
                                std::vector<EntityHandle>& adj_entities_out );

    iGeom::Error getEgFcSense( EntityHandle edge, EntityHandle face, int& sense );

    iGeom::Error measure( const EntityHandle* entities,
                              int entities_size,
                              double* measures );
    iGeom::Error getEntNrmlSense( EntityHandle face, EntityHandle region,
        int& sense );

    iGeom::Error getEgVtxSense( EntityHandle edge, EntityHandle vtx1, EntityHandle vtx2, int& sense );

    iGeom::Error getEntNrmlXYZ( EntityHandle entity,
                                    double x, double y, double z,
                                    double& i, double& j, double& k );

    iGeom::Error getVtxCoord( EntityHandle vertex,
                                  double& x, double& y, double& z );

    iGeom::Error getEntClosestPt( EntityHandle entity,
                                      double near_x, double near_y, double near_z,
                                      double& on_x, double& on_y, double& on_z );
    iGeom::Error getEgEvalXYZ( EntityHandle edge,
                                   double x, double y, double z,
                                   double& on_x, double& on_y, double& on_z,
                                   double& tngt_i, double& tngt_j, double& tngt_k,
                                   double& cvtr_i, double& cvtr_j, double& cvtr_k );
    iGeom::Error getFcEvalXYZ( EntityHandle face,
                                   double x, double y, double z,
                                   double& on_x, double& on_y, double& on_z,
                                   double& nrml_i, double& nrml_j, double& nrml_k,
                                   double& cvtr1_i, double& cvtr1_j, double& cvtr1_k,
                                   double& cvtr2_i, double& cvtr2_j, double& cvtr2_k );

    iGeom::Error getEntURange( EntityHandle edge,
                                   double& u_min, double& u_max );
    iGeom::Error getEntUtoXYZ( EntityHandle edge, double u,
                                   double& x, double& y, double& z );
    iGeom::Error isEntAdj( EntityHandle entity1, EntityHandle entity2,
        bool& adjacent_out );

    iGeom::Error getFaceType( EntityHandle face, std::string& type );

    iGeom::Error isEntParametric( EntityHandle entity, bool& parametric );

    iGeom::Error getEntBoundBox( EntityHandle entity,
                           double& min_x, double& min_y, double& min_z,
                           double& max_x, double& max_y, double& max_z );

#if 0
    inline Error save( const char* file_name,
                       const char* options = 0 );

    inline Error getBoundBox( double& min_x, double& min_y, double& min_z,
                              double& max_x, double& max_y, double& max_z );
    
    inline int getParametric();
    
    inline Error getNumOfType( EntitySetHandle set, EntityType type, int& count_out  );
    

     
    inline Error getEntType( EntityHandle handle, EntityType& type_out );
    inline Error getArrType( const EntityHandle* entity_handles,
                             int entity_handles_Size,
                             EntityType* types_out );
    
    inline Error getEntAdj( EntityHandle handle, 
                            EntityType type_requested,
                            std::vector<EntityHandle>& adj_entities_out );                   
    inline Error getArrAdj( const EntityHandle* entity_handles,
                            int entity_handles_size,
                            EntityType type_requested,
                            std::vector<EntityHandle>& adjacent_entity_handles_out,
                            int* offsets_out );
    inline Error getEnt2ndAdj( EntityHandle handle, 
                               EntityType bridge_dimension,
                               EntityType type_requested,
                               std::vector<EntityHandle>& adj_entities_out );                   
    inline Error getArr2ndAdj( const EntityHandle* entity_handles,
                               int entity_handles_size,
                               EntityType order_adjacent_key,
                               EntityType type_requested,
                               std::vector<EntityHandle>& adjacent_entity_handles_out,
                               int* offsets_out );
    inline Error isEntAdj( EntityHandle entity1, EntityHandle entity2, bool& adjacent_out );
    
    inline Error isArrAdj( const EntityHandle* entities1, 
                           const EntityHandle* entities2,
                           int num_entity_pairs,
                           int* is_adj_out );

    inline Error getEntClosestPt( EntityHandle entity,
                                  double near_x, double near_y, double near_z,
                                  double& on_x, double& on_y, double& on_z );
    inline Error getArrClosestPt( const EntityHandle* handles, 
                                  int handles_size,
                                  StorageOrder order,
                                  const double* near_coordinates,
                                  int near_coordinates_size,
                                  double* on_coordinates );
    
    inline Error getEntNrmlXYZ( EntityHandle entity,
                                double x, double y, double z,
                                double& i, double& j, double& k );
    inline Error getArrNrmlXYZ( const EntityHandle* entities,
                                int entities_size,
                                StorageOrder order,
                                const double* xyz,
                                int xyz_size,
                                double* ijk );
    
    inline Error getEntNrmlPlXYZ( EntityHandle entity,
                                  double x, double y, double z,
                                  double& on_x, double& on_y, double& on_z,
                                  double& i, double& j, double& k );
    inline Error getArrNrmlPlXYZ( const EntityHandle* entities,
                                  int entities_size,
                                  StorageOrder order,
                                  const double* near_xyz,
                                  int near_xyz_size,
                                  double* on_xyz,
                                  double* nrml_ijk );
    
    inline Error getEntTgntXYZ( EntityHandle entity,
                                double x, double y, double z,
                                double& i, double& j, double& k );
    inline Error getArrTgntXYZ( const EntityHandle* entities,
                                int entities_size,
                                StorageOrder order,
                                const double* xyz,
                                int xyz_size,
                                double* ijk );
    
    inline Error getFcCvtrXYZ( EntityHandle face,
                               double x, double y, double z,
                               double& i1, double& j1, double& k1,
                               double& i2, double& j2, double& k2 );
    inline Error getEgCvtrXYZ( EntityHandle edge,
                               double x, double y, double z,
                               double& i, double& j, double& k );
    inline Error getEntArrCvtrXYZ( const EntityHandle* entities,
                                   int entities_size,
                                   StorageOrder order,
                                   const double* xyz,
                                   int xyz_size,
                                   double* cvtr_1,
                                   double* cvtr_2 );
    
    inline Error getEgEvalXYZ( EntityHandle edge,
                               double x, double y, double z,
                               double& on_x, double& on_y, double& on_z,
                               double& tngt_i, double& tngt_j, double& tngt_k,
                               double& cvtr_i, double& cvtr_j, double& cvtr_k );
    inline Error getFcEvalXYZ( EntityHandle face,
                               double x, double y, double z,
                               double& on_x, double& on_y, double& on_z,
                               double& nrml_i, double& nrml_j, double& nrml_k,
                               double& cvtr1_i, double& cvtr1_j, double& cvtr1_k,
                               double& cvtr2_i, double& cvtr2_j, double& cvtr2_k );
    inline Error getArrEgEvalXYZ( const EntityHandle* edges,
                                  int edges_size,
                                  StorageOrder order,
                                  const double* near_coords,
                                  int near_coords_size,
                                  double* on_coords,
                                  double* tangent,
                                  double* curvature );
    inline Error getArrFcEvalXYZ( const EntityHandle* faces,
                                  int faces_size,
                                  StorageOrder order,
                                  const double* near_coords,
                                  int near_coords_size,
                                  double* on_coords,
                                  double* normal,
                                  double* curvature1,
                                  double* curvature2 );

    inline Error getEntBoundBox( EntityHandle entity,
                                 double& min_x, double& min_y, double& min_z,
                                 double& max_x, double& max_y, double& max_z );
    inline Error getArrBoundBox( const EntityHandle* entities,
                                 int entities_size,
                                 StorageOrder order,
                                 double* min_corners,
                                 double* max_corners );
    
    inline Error getVtxCoord( EntityHandle vertex,
                              double& x, double& y, double& z );
    inline Error getVtxArrCoords( const EntityHandle* vertices,
                                  int vertices_size,
                                  StorageOrder order,
                                  double* coords );
    
    inline Error getPntRayIntsct( double x, double y, double z,
                                  double i, double j, double k,
                                  StorageOrder order,
                                  std::vector<EntityHandle>& entities_out,
                                  std::vector<double>& points_out,
                                  std::vector<double>& params_out );
    
    inline Error getPntClsf( double x, double y, double z, 
                             EntityHandle& handle_out ) ;
    inline Error getPntArrClsf( StorageOrder order,
                                const double* coords,
                                int coords_size,
                                EntityHandle* entities_out );
    
    inline Error getEntNrmlSense( EntityHandle face, EntityHandle region, int& sense );
    inline Error getEgFcSense( EntityHandle edge, EntityHandle face, int& sense );
    inline Error getEgVtxSense( EntityHandle edge, EntityHandle vtx1, EntityHandle vtx2, int& sense );
    
    inline Error getArrNrmlSense( const EntityHandle* faces, int faces_size,
                                  const EntityHandle* vols,  int vols_size,
                                  int* senses_out );
    inline Error getEgFcArrSense( const EntityHandle* edges, int edges_size,
                                  const EntityHandle* faces, int faces_size,
                                  int* senses_out );
    inline Error getEgVtxArrSense( const EntityHandle* edges, int edges_size,
                                   const EntityHandle* vertices1, int vertices1_size,
                                   const EntityHandle* vertices2, int vertices2_size,
                                   int* senses_out );
  inline Error getSense(EntityHandle ent, EntityHandle wrt_ent, int &sense);
  inline Error getArrSense(const EntityHandle *ent, int num_ents, EntityHandle wrt_ent, int *sense);
  
    inline Error measure( const EntityHandle* entities,
                          int entities_size,
                          double* measures );
    
    inline Error getFaceType( EntityHandle face, std::string& type );
    inline Error isEntParametric( EntityHandle entity, bool& parametric );
    inline Error isArrParametric( const EntityHandle* entities,
                                  int entities_size,
                                  int* is_parametric );
    
    inline Error getEntUVtoXYZ( EntityHandle face,
                                double u, double v,
                                double& x, double& y, double& z );
    inline Error getEntUtoXYZ( EntityHandle edge, double u,
                               double& x, double& y, double& z );
    inline Error getArrUVtoXYZ( const EntityHandle* faces,
                                int faces_size,
                                StorageOrder order,
                                const double* uv,
                                int uv_size,
                                double* xyz );
    inline Error getArrUtoXYZ( const EntityHandle* edges,
                               int edges_size,
                               const double* u,
                               int u_size,
                               StorageOrder order,
                               double* xyz );
    
    inline Error getEntXYZtoUV( EntityHandle face,
                                double x, double y, double z,
                                double& u, double& v );
    inline Error getEntXYZtoU( EntityHandle edge,
                               double x, double y, double z,
                               double& u );
    inline Error getArrXYZtoUV( const EntityHandle* faces,
                                int faces_size,
                                StorageOrder order,
                                const double* coords,
                                int coords_size,
                                double* uv );
    inline Error getArrXYZtoU( const EntityHandle* edges,
                               int edges_size,
                               StorageOrder order,
                               const double* coords,
                               int coords_size,
                               double* u );
    inline Error getEntXYZtoUVHint( EntityHandle face,
                                    double x, double y, double z,
                                    double& u, double& v );
    inline Error getArrXYZtoUVHint( const EntityHandle* faces,
                                    int faces_size,
                                    StorageOrder order,
                                    const double* coords,
                                    int coords_size,
                                    double* uv );
    
    inline Error getEntNrmlUV( EntityHandle face,
                               double u, double v,
                               double& i, double& j, double& k );
    inline Error getArrNrmlUV( const EntityHandle* faces,
                               int faces_size,
                               StorageOrder order,
                               const double* uv,
                               int uv_size,
                               double* normals );
    inline Error getEntTgntU( EntityHandle edge,
                              double u, 
                              double& i, double& j, double& k );
    inline Error getArrTgntU( const EntityHandle* edges,
                              int edges_size,
                              StorageOrder order,
                              const double* u,
                              int u_size,
                              double* normals );

    inline Error getEnt1stDrvt( EntityHandle handle,
                                double u, double v,
                                double& du_i, double& du_j, double& du_k,
                                double& dv_i, double& dv_j, double& dv_k );
    inline Error getEnt2ndDrvt( EntityHandle handle,
                                double u, double v,
                                double& duu_i, double& duu_j, double& duu_k,
                                double& dvv_i, double& dvv_j, double& dvv_k,
                                double& duv_i, double& duv_j, double& duv_k );
    inline Error getArr1stDrvt( const EntityHandle* entities,
                                int entities_size,
                                StorageOrder order,
                                const double* uv,
                                int uv_size,
                                double* dvtr_u,
                                double* dvtr_v );
    inline Error getArr2ndDrvt( const EntityHandle* entities,
                                int entities_size,
                                StorageOrder order,
                                const double* uv,
                                int uv_size,
                                double* dvtr_uu,
                                double* dvtr_vv,
                                double* dvtr_uv );
    

    inline Error getFcCvtrUV( EntityHandle face, double u, double v,
                              double& i1, double& j1, double& k1,
                              double& i2, double& j2, double& k2 );
    inline Error getFcArrCvtrUV( const EntityHandle* faces, int faces_size,
                                 StorageOrder order,
                                 const double* uv, int uv_size,
                                 double* cvtr1, double* cvtr2 );

    inline Error isEntPeriodic( EntityHandle entity,
                                bool& in_u, bool& in_v );
    inline Error isArrPeriodic( const EntityHandle* entities,
                                int entities_size,
                                int* in_uv );

    inline Error isFcDegenerate( EntityHandle face, bool& is_degenerate );
    inline Error isFcArrDegenerate( const EntityHandle* faces,
                                    int faces_size,
                                    int* degenerate );
    
    inline Error getTolerance( int& type_out, double& tolerance_out );
    inline Error getEntTolerance( EntityHandle entity, double& tolerance );
    inline Error getArrTolerance( const EntityHandle* entities,
                                  int entities_size,
                                  double* tolerances );

    inline Error getEntUVRange( EntityHandle face,
                                double& u_min, double& v_min,
                                double& u_max, double& v_max );
    inline Error getEntURange( EntityHandle edge,
                               double& u_min, double& u_max );
    inline Error getArrUVRange( const EntityHandle* faces,
                                int faces_size,
                                StorageOrder order,
                                double* uv_min,
                                double* uv_max );
    inline Error getArrURange( const EntityHandle* edges,
                               int edges_size,
                               double* u_min,
                               double* u_max );

    inline Error getEntUtoUV( EntityHandle edge, 
                              EntityHandle face,
                              double edge_u,
                              double& face_u, double& face_v );
    inline Error getVtxToUV( EntityHandle vertex,
                             EntityHandle face,
                             double& u, double& v );
    inline Error getVtxToU( EntityHandle vertex,
                            EntityHandle edge,
                            double& u );
    inline Error getArrUtoUV( const EntityHandle* edges, int edges_size,
                              const EntityHandle* faces, int faces_size,
                              const double* edge_u, int edge_u_size,
                              StorageOrder order,
                              double* face_uv );
    inline Error getVtxArrToUV( const EntityHandle* vertices, int vertices_size,
                                const EntityHandle* faces, int faces_size,
                                StorageOrder order,
                                double* face_uv );
    inline Error getVtxArrToU( const EntityHandle* vertices, int vertices_size,
                               const EntityHandle* edges, int edges_size,
                               double* edge_u );
    
    
    inline Error deleteAll();
    
    inline Error deleteEnt( EntityHandle entity );

    inline Error copyEnt( EntityHandle source, EntityHandle& copy );
    
    inline Error createSphere( double radius, EntityHandle& sphere );
    inline Error createPrism( double height, int num_sides,
                              double maj_radius, double min_radius,
                              EntityHandle& prism );
    inline Error createBrick( double x, double y, double z, EntityHandle& brick );
    inline Error createCylinder( double height, double maj_rad, double min_rad,
                                 EntityHandle& cylinder );
    inline Error createTorus( double maj_rad, double min_rad, EntityHandle& torus );
    
    inline Error moveEnt( EntityHandle entity, double x, double y, double z );
    inline Error rotateEnt( EntityHandle entity, double angle,
                            double axis_x, double axis_y, double axis_z );
    inline Error reflectEnt( EntityHandle entity, 
                             double norm_x, double norm_y, double norm_z );
    inline Error scaleEnt( EntityHandle entity, 
                           double x_factor, double y_factor, double z_factor );
  
    inline Error uniteEnts( const EntityHandle* entities,
                            int entities_size,
                            EntityHandle& result_entity );
    inline Error subtractEnts( EntityHandle blank, EntityHandle tool, 
                               EntityHandle& result );
    inline Error intersectEnts( EntityHandle entity1, EntityHandle entity2,
                                EntityHandle& result );
    
    inline Error sectionEnt( EntityHandle entity,
                             double plane_x, double plane_y, double plane_z,
                             double offset,
                             bool reverse,
                             EntityHandle& result );
    
    inline Error sweepEntAboutAxis( EntityHandle entity,
                                    double angle,
                                    double axis_x, 
                                    double axis_y,
                                    double axis_z,
                                    EntityHandle& swept_entity );
    
    inline Error imprintEnts( const EntityHandle* entities, int entities_size );
    inline Error mergeEnts( const EntityHandle* entities, int entities_size, double tolerance );

#if 0
    // iterators are already defined in iGeom.hpp
/** \class EntArrIter iGeom.hpp "iGeom.hpp"
 * \brief Class for iterating over %iGeom entity arrays.
 */
    class EntArrIter {
      private:
        friend class iGeom;
        iBase_EntityArrIterator mHandle;
        iGeom_Instance mInstance;
        int mSize;
      public:
        EntArrIter() : mHandle(0), mInstance(0), mSize(0) {}
        inline ~EntArrIter();
        inline Error getNext( EntityHandle* entity_handles_out,
                              int& size_out,
                              bool& has_more_data_out );
        inline Error reset();
    };
    
/** \class EntIter iGeom.hpp "iGeom.hpp"
 * \brief Class for iterating over %iGeom entities.
 */
    class EntIter {
      private:
        friend class iGeom;
        iBase_EntityIterator mHandle;
        iGeom_Instance mInstance;
      public:
        EntIter() : mHandle(0), mInstance(0) {}
        inline ~EntIter();
        inline Error getNext( EntityHandle& entity_handle_out,
                              bool& has_more_data_out );
        inline Error reset();
    };
#endif
    inline Error initEntIter( EntitySetHandle set,
                              EntityType requested_type,
                              EntIter& iter );
    inline Error initEntArrIter( EntitySetHandle set,
                                 EntityType requested_type,
                                 int requested_array_size,
                                 EntArrIter& iter );
#endif // comment out the bulk
  private:
    
    moab::FBEngine * _fbEngine;
    bool _smooth;
    MKCore * _mk;
    unsigned int _index;
      // prohibit copying
    FBiGeom( const FBiGeom& ) {}
    void operator=(const FBiGeom&) {}
};

} // namespace MeshKit
#endif
