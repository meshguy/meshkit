#ifndef MESHKIT_ITAPS_GEOM_HPP
#define MESHKIT_ITAPS_GEOM_HPP

/** \file iGeom.hpp
 */

#include "iGeom.h"
#include "iRel.hpp"

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <string>

/** \class iGeom
 * \brief C++ interface to ITAPS iGeom interface
 *
 * This class is a simple wrapper for the ITAPS iGeom interface.  The primary benefit to using this class
 * instead of iGeom directly is that lists of handles are passed as std::vectors instead of pointers to
 * handle arrays.  This file includes both declaration and definition of all iGeom class functions, i.e.
 * all functions are inlined.  The class can be constructed and destructed in the standard C++ way; the
 * implementation of those functions call into the standard iGeom C functions newGeom and dtor.
 *
 * For complete documentation of these functions, see the iGeom header in the CGM source
 * (http://trac.mcs.anl.gov/projects/ITAPS/browser/cgm/trunk/itaps/iGeom.h for now).
 */

class iGeom {
	protected:
		iGeom_Instance mInstance;
	public:
	virtual inline iRel::IfaceType iface_type() const {
	    return iRel::IGEOM_IFACE;
	}

		typedef iBase_EntitySetHandle EntitySetHandle;
		typedef iBase_EntityHandle EntityHandle;
		typedef iBase_TagHandle TagHandle;
		typedef iBase_ErrorType Error;
		typedef iBase_EntityType EntityType;
		typedef iBase_StorageOrder StorageOrder;
		typedef iBase_TagValueType TagValueType;

		inline iGeom(const char* options = 0);
		inline iGeom(iGeom_Instance instance);

		inline iGeom(bool meshBased);

		virtual inline ~iGeom();

		virtual inline Error load(const char* file_name, const char* options = 0);

		virtual inline Error save(const char* file_name, const char* options = 0);

		virtual inline Error getBoundBox(double& min_x, double& min_y, double& min_z,
				double& max_x, double& max_y, double& max_z) const;

		virtual inline int getParametric();

		virtual inline Error getNumOfType(EntitySetHandle set, EntityType type,
				int& count_out) const;

		virtual inline Error getEntities(EntitySetHandle set, EntityType type,
				std::vector<EntityHandle>& entities_out) const;

		virtual inline Error getEntType(EntityHandle handle, EntityType& type_out) const;

		virtual inline Error getArrType(const EntityHandle* entity_handles,
				int entity_handles_Size, EntityType* types_out) const;

		virtual inline Error getEntAdj(EntityHandle handle,
				EntityType type_requested, std::vector<EntityHandle>& adj_entities_out) const;

		virtual inline Error getArrAdj(const EntityHandle* entity_handles,
				int entity_handles_size, EntityType type_requested, std::vector<
				EntityHandle>& adjacent_entity_handles_out, int* offsets_out) const;

		virtual inline Error getEnt2ndAdj(EntityHandle handle,
				EntityType bridge_dimension, EntityType type_requested, std::vector<
				EntityHandle>& adj_entities_out) const;

		virtual inline Error getArr2ndAdj(const EntityHandle* entity_handles,
				int entity_handles_size, EntityType order_adjacent_key,
				EntityType type_requested,
				std::vector<EntityHandle>& adjacent_entity_handles_out, int* offsets_out) const;

		virtual inline Error isEntAdj(EntityHandle entity1, EntityHandle entity2,
				bool& adjacent_out) const;

		virtual inline Error isArrAdj(const EntityHandle* entities1,
				const EntityHandle* entities2, int num_entity_pairs, int* is_adj_out) const;

		virtual inline Error getEntClosestPt(EntityHandle entity, double near_x,
				double near_y, double near_z, double& on_x, double& on_y, double& on_z) const;

		virtual inline Error getEntClosestPtTrimmed(EntityHandle entity, double near_x,
				double near_y, double near_z, double& on_x, double& on_y, double& on_z) const;

		virtual inline Error getArrClosestPt(const EntityHandle* handles,
				int handles_size, StorageOrder order, const double* near_coordinates,
				int near_coordinates_size, double* on_coordinates) const;

		virtual inline Error getEntNrmlXYZ(EntityHandle entity, double x, double y,
				double z, double& i, double& j, double& k) const;

		virtual inline Error getArrNrmlXYZ(const EntityHandle* entities,
				int entities_size, StorageOrder order, const double* xyz, int xyz_size,
				double* ijk) const;

		virtual inline Error getEntNrmlPlXYZ(EntityHandle entity, double x, double y,
				double z, double& on_x, double& on_y, double& on_z, double& i, double& j,
				double& k) const;

		virtual inline Error getArrNrmlPlXYZ(const EntityHandle* entities,
				int entities_size, StorageOrder order, const double* near_xyz,
				int near_xyz_size, double* on_xyz, double* nrml_ijk) const;

		virtual inline Error getEntTgntXYZ(EntityHandle entity, double x, double y,
				double z, double& i, double& j, double& k) const;

		virtual inline Error getArrTgntXYZ(const EntityHandle* entities,
				int entities_size, StorageOrder order, const double* xyz, int xyz_size,
				double* ijk) const;

		virtual inline Error getFcCvtrXYZ(EntityHandle face, double x, double y,
				double z, double& i1, double& j1, double& k1, double& i2, double& j2,
				double& k2) const;

		virtual inline Error getEgCvtrXYZ(EntityHandle edge, double x, double y,
				double z, double& i, double& j, double& k) const;

		virtual inline Error getEntArrCvtrXYZ(const EntityHandle* entities,
				int entities_size, StorageOrder order, const double* xyz, int xyz_size,
				double* cvtr_1, double* cvtr_2) const;

		virtual inline Error getEgEvalXYZ(EntityHandle edge, double x, double y,
				double z, double& on_x, double& on_y, double& on_z, double& tngt_i,
				double& tngt_j, double& tngt_k, double& cvtr_i, double& cvtr_j,
				double& cvtr_k) const;

		virtual inline Error getFcEvalXYZ(EntityHandle face, double x, double y,
				double z, double& on_x, double& on_y, double& on_z, double& nrml_i,
				double& nrml_j, double& nrml_k, double& cvtr1_i, double& cvtr1_j,
				double& cvtr1_k, double& cvtr2_i, double& cvtr2_j, double& cvtr2_k) const;

		virtual inline Error getArrEgEvalXYZ(const EntityHandle* edges,
				int edges_size, StorageOrder order, const double* near_coords,
				int near_coords_size, double* on_coords, double* tangent,
				double* curvature) const;

		virtual inline Error getArrFcEvalXYZ(const EntityHandle* faces,
				int faces_size, StorageOrder order, const double* near_coords,
				int near_coords_size, double* on_coords, double* normal,
				double* curvature1, double* curvature2) const;

		virtual inline Error
			getEntBoundBox(EntityHandle entity, double& min_x, double& min_y,
					double& min_z, double& max_x, double& max_y, double& max_z) const;

		virtual inline Error getArrBoundBox(const EntityHandle* entities,
				int entities_size, StorageOrder order, double* min_corners,
				double* max_corners) const;

		virtual inline Error getVtxCoord(EntityHandle vertex, double& x, double& y,
				double& z) const;

		virtual inline Error getVtxArrCoords(const EntityHandle* vertices,
				int vertices_size, StorageOrder order, double* coords) const;

		virtual inline Error getPntRayIntsct(double x, double y, double z, double i,
				double j, double k, StorageOrder order,
				std::vector<EntityHandle>& entities_out, std::vector<double>& points_out,
				std::vector<double>& params_out) const;

		virtual inline Error getPntClsf(double x, double y, double z,
				EntityHandle& handle_out)const;

		virtual inline Error getPntArrClsf(StorageOrder order, const double* coords,
				int coords_size, EntityHandle* entities_out) const;

		virtual inline Error getEntNrmlSense(EntityHandle face, EntityHandle region,
				int& sense) const;

		virtual inline Error getEgFcSense(EntityHandle edge, EntityHandle face,
				int& sense) const;

		virtual inline Error getEgVtxSense(EntityHandle edge, EntityHandle vtx1,
				EntityHandle vtx2, int& sense) const;

		virtual inline Error getArrNrmlSense(const EntityHandle* faces,
				int faces_size, const EntityHandle* vols, int vols_size, int* senses_out)const;

		virtual inline Error getEgFcArrSense(const EntityHandle* edges,
				int edges_size, const EntityHandle* faces, int faces_size,
				int* senses_out) const;

		virtual inline Error getEgVtxArrSense(const EntityHandle* edges,
				int edges_size, const EntityHandle* vertices1, int vertices1_size,
				const EntityHandle* vertices2, int vertices2_size, int* senses_out) const;

		virtual inline Error getSense(EntityHandle ent, EntityHandle wrt_ent,
				int &sense) const;

		virtual inline Error getArrSense(const EntityHandle *ent, int num_ents,
				EntityHandle wrt_ent, int *sense) const;

		virtual inline Error measure(const EntityHandle* entities, int entities_size,
				double* measures) const;

		virtual inline Error getFaceType(EntityHandle face, std::string& type) const;

		virtual inline Error isEntParametric(EntityHandle entity, bool& parametric) const;

		virtual inline Error isArrParametric(const EntityHandle* entities,
				int entities_size, int* is_parametric) const;

		virtual inline Error getEntUVtoXYZ(EntityHandle face, double u, double v,
				double& x, double& y, double& z) const;

		virtual inline Error getEntUtoXYZ(EntityHandle edge, double u, double& x,
				double& y, double& z) const;

		virtual inline Error getArrUVtoXYZ(const EntityHandle* faces, int faces_size,
				StorageOrder order, const double* uv, int uv_size, double* xyz) const;

		virtual inline Error getArrUtoXYZ(const EntityHandle* edges, int edges_size,
				const double* u, int u_size, StorageOrder order, double* xyz) const;

		virtual inline Error getEntXYZtoUV(EntityHandle face, double x, double y,
				double z, double& u, double& v) const;

		virtual inline Error getEntXYZtoU(EntityHandle edge, double x, double y,
				double z, double& u) const;

		virtual inline Error getArrXYZtoUV(const EntityHandle* faces, int faces_size,
				StorageOrder order, const double* coords, int coords_size, double* uv) const;

		virtual inline Error getArrXYZtoU(const EntityHandle* edges, int edges_size,
				StorageOrder order, const double* coords, int coords_size, double* u) const;

		virtual inline Error getEntXYZtoUVHint(EntityHandle face, double x, double y,
				double z, double& u, double& v) const;

		virtual inline Error getArrXYZtoUVHint(const EntityHandle* faces,
				int faces_size, StorageOrder order, const double* coords,
				int coords_size, double* uv) const;

		virtual inline Error getEntNrmlUV(EntityHandle face, double u, double v,
				double& i, double& j, double& k) const;

		virtual inline Error getArrNrmlUV(const EntityHandle* faces, int faces_size,
				StorageOrder order, const double* uv, int uv_size, double* normals) const;

		virtual inline Error getEntTgntU(EntityHandle edge, double u, double& i,
				double& j, double& k) const;

		virtual inline Error getArrTgntU(const EntityHandle* edges, int edges_size,
				StorageOrder order, const double* u, int u_size, double* normals) const;

		virtual inline Error getEnt1stDrvt(EntityHandle handle, double u, double v,
				double& du_i, double& du_j, double& du_k, double& dv_i, double& dv_j,
				double& dv_k) const;
		virtual inline Error
			getEnt2ndDrvt(EntityHandle handle, double u, double v, double& duu_i,
					double& duu_j, double& duu_k, double& dvv_i, double& dvv_j,
					double& dvv_k, double& duv_i, double& duv_j, double& duv_k) const;

		virtual inline Error getArr1stDrvt(const EntityHandle* entities,
				int entities_size, StorageOrder order, const double* uv, int uv_size,
				double* dvtr_u, double* dvtr_v) const;

		virtual inline Error getArr2ndDrvt(const EntityHandle* entities,
				int entities_size, StorageOrder order, const double* uv, int uv_size,
				double* dvtr_uu, double* dvtr_vv, double* dvtr_uv) const;

		virtual inline Error getFcCvtrUV(EntityHandle face, double u, double v,
				double& i1, double& j1, double& k1, double& i2, double& j2, double& k2) const;

		virtual inline Error getFcArrCvtrUV(const EntityHandle* faces,
				int faces_size, StorageOrder order, const double* uv, int uv_size,
				double* cvtr1, double* cvtr2) const;

		virtual inline Error isEntPeriodic(EntityHandle entity, bool& in_u,
				bool& in_v) const;

		virtual inline Error isArrPeriodic(const EntityHandle* entities,
				int entities_size, int* in_uv) const;

		virtual inline Error isFcDegenerate(EntityHandle face, bool& is_degenerate) const;

		virtual inline Error isFcArrDegenerate(const EntityHandle* faces,
				int faces_size, int* degenerate) const;

		virtual inline Error getTolerance(int& type_out, double& tolerance_out) const;
		virtual inline Error getEntTolerance(EntityHandle entity, double& tolerance) const;
		virtual inline Error getArrTolerance(const EntityHandle* entities,
				int entities_size, double* tolerances) const;

		virtual inline Error getEntUVRange(EntityHandle face, double& u_min,
				double& v_min, double& u_max, double& v_max) const;
		virtual inline Error getEntURange(EntityHandle edge, double& u_min,
				double& u_max) const;
		virtual inline Error getArrUVRange(const EntityHandle* faces, int faces_size,
				StorageOrder order, double* uv_min, double* uv_max) const;
		virtual inline Error getArrURange(const EntityHandle* edges, int edges_size,
				double* u_min, double* u_max) const;

		virtual inline Error getEntUtoUV(EntityHandle edge, EntityHandle face,
				double edge_u, double& face_u, double& face_v) const;
		virtual inline Error getVtxToUV(EntityHandle vertex, EntityHandle face,
				double& u, double& v) const;
		virtual inline Error getVtxToU(EntityHandle vertex, EntityHandle edge,
				double& u) const;
		virtual inline Error getArrUtoUV(const EntityHandle* edges, int edges_size,
				const EntityHandle* faces, int faces_size, const double* edge_u,
				int edge_u_size, StorageOrder order, double* face_uv) const;
		virtual inline Error getVtxArrToUV(const EntityHandle* vertices,
				int vertices_size, const EntityHandle* faces, int faces_size,
				StorageOrder order, double* face_uv) const;
		virtual inline Error getVtxArrToU(const EntityHandle* vertices,
				int vertices_size, const EntityHandle* edges, int edges_size,
				double* edge_u) const;

		virtual inline Error deleteAll();

		virtual inline Error deleteEnt(EntityHandle entity);

		virtual inline Error copyEnt(EntityHandle source, EntityHandle& copy);

		virtual inline Error createSphere(double radius, EntityHandle& sphere);
		virtual inline Error createPrism(double height, int num_sides,
				double maj_radius, double min_radius, EntityHandle& prism);
		virtual inline Error createBrick(double x, double y, double z,
				EntityHandle& brick);
		virtual inline Error createCylinder(double height, double maj_rad,
				double min_rad, EntityHandle& cylinder);
		virtual inline Error createCone(double height, double maj_rad,
				double min_rad, double top_rad, EntityHandle& cone);
		virtual inline Error createTorus(double maj_rad, double min_rad,
				EntityHandle& torus);

		virtual inline Error moveEnt(EntityHandle entity, double x, double y,
				double z);
		virtual inline Error rotateEnt(EntityHandle entity, double angle,
				double axis_x, double axis_y, double axis_z);
		virtual inline Error reflectEnt(EntityHandle entity, double x,
				double y, double z, double norm_x,
				double norm_y, double norm_z);
		virtual inline Error scaleEnt(EntityHandle entity, double x,
				double y, double z, double x_factor,
				double y_factor, double z_factor);

		virtual inline Error uniteEnts(const EntityHandle* entities,
				int entities_size, EntityHandle& result_entity);
		virtual inline Error subtractEnts(EntityHandle blank, EntityHandle tool,
				EntityHandle& result);
		virtual inline Error intersectEnts(EntityHandle entity1,
				EntityHandle entity2, EntityHandle& result);

		virtual inline Error sectionEnt(EntityHandle entity, double plane_x,
				double plane_y, double plane_z, double offset, bool reverse,
				EntityHandle& result);

		virtual inline Error sweepEntAboutAxis(EntityHandle entity, double angle,
				double axis_x, double axis_y, double axis_z, EntityHandle& swept_entity);

		virtual inline Error imprintEnts(const EntityHandle* entities,
				int entities_size);
		virtual inline Error mergeEnts(const EntityHandle* entities,
				int entities_size, double tolerance);

		// copied from iBaseVirtual.hpp
		iGeom_Instance instance() {
			return mInstance;
		}

		virtual inline std::string getDescription() const;

		virtual inline Error getErrorType() const;

		virtual inline EntitySetHandle getRootSet() const;

		virtual inline Error createEntSet(bool is_list, EntitySetHandle& handle_out);
		virtual inline Error destroyEntSet(EntitySetHandle handle);
		virtual inline Error isList(EntitySetHandle handle, bool& is_list);

		virtual inline Error getNumEntSets(EntitySetHandle set, int num_hops,
				int& num_sets_out) const;
		virtual inline Error getEntSets(EntitySetHandle set, int num_hops, std::vector<
				EntitySetHandle>& contained_sets_out) const;

		virtual inline Error addEntToSet(EntityHandle entity, EntitySetHandle set);
		virtual inline Error rmvEntFromSet(EntityHandle entity, EntitySetHandle set);

		virtual inline Error addEntArrToSet(const EntityHandle* entity_handles,
				int entity_handles_size, EntitySetHandle entity_set);
		virtual inline Error rmvEntArrFromSet(const EntityHandle* entity_handles,
				int entity_handles_size, EntitySetHandle entity_set);

		virtual inline Error addEntSet(EntitySetHandle to_add, EntitySetHandle add_to);
		virtual inline Error rmvEntSet(EntitySetHandle to_rmv, EntitySetHandle rmv_from);

		virtual inline Error isEntContained(EntitySetHandle set, EntityHandle ent,
				bool& contained_out) const;
		virtual inline Error isEntArrContained(EntitySetHandle containing_set,
				const EntityHandle* entity_handles, int num_entity_handles,
				bool* is_contained_out) const;
		virtual inline Error isEntSetContained(EntitySetHandle containing_set,
				EntitySetHandle contained_set, bool& contained_out) const;

		virtual inline Error addPrntChld(EntitySetHandle parent, EntitySetHandle child);
		virtual inline Error rmvPrntChld(EntitySetHandle parent, EntitySetHandle child);
		virtual inline Error isChildOf(EntitySetHandle parent, EntitySetHandle child,
				bool& is_child_out) const;
		virtual inline Error getNumChld(EntitySetHandle parent, int num_hops,
				int& num_child_out) const;
		virtual inline Error getNumPrnt(EntitySetHandle child, int num_hops,
				int& num_parent_out) const;
		virtual inline Error getChldn(EntitySetHandle parent, int num_hops, std::vector<
				EntitySetHandle>& children_out) const;
		virtual inline Error getPrnts(EntitySetHandle child, int num_hops, std::vector<
				EntitySetHandle>& parents_out) const;

		virtual inline Error subtract(EntitySetHandle set1, EntitySetHandle set2,
				EntitySetHandle& result_set_out);
		virtual inline Error intersect(EntitySetHandle set1, EntitySetHandle set2,
				EntitySetHandle& result_set_out);
		virtual inline Error unite(EntitySetHandle set1, EntitySetHandle set2,
				EntitySetHandle& result_set_out);

		virtual inline Error createTag(const char* tag_name, int tag_num_type_values,
				TagValueType tag_type, TagHandle& tag_handle_out);

		virtual inline Error destroyTag(TagHandle tag_handle, bool forced);
		virtual inline Error getTagName(TagHandle tag_handle, std::string& name_out) const;
		virtual inline Error getTagSizeValues(TagHandle tag_handle, int& size_out) const;
		virtual inline Error getTagSizeBytes(TagHandle tag_handle, int& size_out) const;
		virtual inline Error getTagHandle(const char* name, TagHandle& handle_out) const;
		virtual inline Error getTagType(TagHandle tag_handle, TagValueType& type_out) const;

		virtual inline Error setEntSetData(EntitySetHandle set_handle, TagHandle tag_handle,
				const void* tag_value);
		virtual inline Error setEntSetIntData(EntitySetHandle set_handle,
				TagHandle tag_handle, int value);
		virtual inline Error setEntSetDblData(EntitySetHandle set_handle,
				TagHandle tag_handle, double value);
		virtual inline Error setEntSetEHData(EntitySetHandle set_handle,
				TagHandle tag_handle, EntityHandle value);
		virtual inline Error setEntSetESHData(EntitySetHandle set_handle,
				TagHandle tag_handle, EntitySetHandle value);

		virtual inline Error getEntSetData(EntitySetHandle set_handle, TagHandle tag_handle,
				void* tag_value_out) const;
		virtual inline Error getEntSetIntData(EntitySetHandle set_handle,
				TagHandle tag_handle, int& value_out) const;
		virtual inline Error getEntSetDblData(EntitySetHandle set_handle,
				TagHandle tag_handle, double& value_out) const;
		virtual inline Error getEntSetEHData(EntitySetHandle set_handle,
				TagHandle tag_handle, EntityHandle& value_out) const;
		virtual inline Error getEntSetESHData(EntitySetHandle set_handle,
				TagHandle tag_handle, EntitySetHandle& value_out) const;

		virtual inline Error getAllEntSetTags(EntitySetHandle set,
				std::vector<TagHandle>& tags_out) const;
		virtual inline Error getAllTags(EntityHandle entity,
				std::vector<TagHandle>& tags_out) const;

		virtual inline Error rmvEntSetTag(EntitySetHandle set, TagHandle tag);
		virtual inline Error rmvTag(EntityHandle entity, TagHandle tag);
		virtual inline Error rmvArrTag(const EntityHandle* handles, int size, TagHandle tag);

		virtual inline Error getArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle, void* tag_values_out) const;
		virtual inline Error getIntArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle, int* tag_values_out) const;
		virtual inline Error getDblArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle, double* tag_values_out) const;
		virtual inline Error getEHArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle,
				EntityHandle* tag_values_out) const;
		virtual inline Error getESHArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle,
				EntitySetHandle* tag_values_out) const;

		virtual inline Error setArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle, const void* tag_values);
		virtual inline Error setIntArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle, const int* tag_values);
		virtual inline Error setDblArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle, const double* tag_values);
		virtual inline Error setEHArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle,
				const EntityHandle* tag_values);
		virtual inline Error setESHArrData(const EntityHandle* entity_handles,
				int entity_handles_size, TagHandle tag_handle,
				const EntitySetHandle* tag_values);

		virtual inline Error setData(EntityHandle entity_handle, TagHandle tag_handle,
				const void* tag_value);
		virtual inline Error setIntData(EntityHandle entity_handle, TagHandle tag_handle,
				int value);
		virtual inline Error setDblData(EntityHandle entity_handle, TagHandle tag_handle,
				double value);
		virtual inline Error setEHData(EntityHandle entity_handle, TagHandle tag_handle,
				EntityHandle value);
		virtual inline Error setESHData(EntityHandle entity_handle, TagHandle tag_handle,
				EntitySetHandle value);

		virtual inline Error getData(EntityHandle entity_handle, TagHandle tag_handle,
				void* tag_value_out) const;
		virtual inline Error getFacets(EntityHandle entity_handle, double dist_tolerance,
						   std::vector<double> &point, std::vector<int> &facets) const;
		virtual inline Error getIntData(EntityHandle entity_handle, TagHandle tag_handle,
				int& value_out) const;
		virtual inline Error getDblData(EntityHandle entity_handle, TagHandle tag_handle,
				double& value_out) const;
		virtual inline Error getEHData(EntityHandle entity_handle, TagHandle tag_handle,
				EntityHandle& value_out) const;
		virtual inline Error getESHData(EntityHandle entity_handle, TagHandle tag_handle,
				EntitySetHandle& value_out) const;

		virtual bool isFBiGeom() {return false;}
		// end copy from iBaseVirtual.hpp
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
				EntArrIter() :
					mHandle(0), mInstance(0), mSize(0) {
					}
				inline ~EntArrIter();
				inline Error getNext(EntityHandle* entity_handles_out, int& size_out,
						bool& has_more_data_out);
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
				EntIter() :
					mHandle(0), mInstance(0) {
					}
				inline ~EntIter();
				inline Error getNext(EntityHandle& entity_handle_out,
						bool& has_more_data_out);
				inline Error reset();
		};

		virtual inline Error initEntIter(EntitySetHandle set,
				EntityType requested_type, EntIter& iter);
		virtual inline Error initEntArrIter(EntitySetHandle set,
				EntityType requested_type, int requested_array_size, EntArrIter& iter);
	private:
		bool iGeomInstanceOwner;

		// prohibit copying
		iGeom(const iGeom&) {
		}
		void operator=(const iGeom&) {
		}
};

inline iGeom::iGeom(const char* options) :
	iGeomInstanceOwner(true)
{
	int err, len = options ? strlen(options) : 0;
	iGeom_newGeom(options, &mInstance, &err, len);
	if (iBase_SUCCESS != err) {
		mInstance = 0;
		iGeomInstanceOwner = false;
	}
}

inline iGeom::iGeom(iGeom_Instance instance) :
	iGeomInstanceOwner(false)
{
	mInstance = instance;
}
inline iGeom::iGeom(bool meshBased)
{
	mInstance = 0;
	iGeomInstanceOwner = false;
}
inline iGeom::~iGeom()
{
	if (iGeomInstanceOwner) {
		int err;
		iGeom_dtor(mInstance, &err);
	}
}

inline iGeom::Error iGeom::load(const char* file_name, const char* options)
{
	int err, len = options ? strlen(options) : 0;
	iGeom_load(mInstance, file_name, options, &err, strlen(file_name), len);
	return (Error) err;
}

inline iGeom::Error iGeom::save(const char* file_name, const char* options)
{
	int err, len = options ? strlen(options) : 0;
	iGeom_save(mInstance, file_name, options, &err, strlen(file_name), len);
	return (Error) err;
}

//inline iGeom::StorageOrder
//iGeom::getDfltStorage()
//{
//  int err, order;
//  iGeom_getDfltStorage( mInstance, &order, &err );
//  return (iBase_SUCCESS == err) ? (StorageOrder)order : iBase_UNDETERMINED;
//}

inline iGeom::Error iGeom::getNumOfType(EntitySetHandle set, EntityType type,
		int& count_out) const
{
	int err;
	iGeom_getNumOfType(mInstance, set, type, &count_out, &err);
	return (Error) err;
}

inline iGeom::Error iGeom::getEntities(EntitySetHandle set, EntityType type,
		std::vector<EntityHandle>& entities_out) const
{
	// if input vect has no allocated space, allocate some so
	// we don't accidentally ask the impl to allocate an array
	if (entities_out.capacity() == 0) {
		int count;
		Error err2 = getNumOfType(set, iBase_ALL_TYPES, count);
		if (err2 != iBase_SUCCESS)
			return err2;
		entities_out.resize(count);
	}

	// try getting results using whatever space input vector has allocated
	int err, size = 0, alloc = entities_out.capacity();
	entities_out.resize(entities_out.capacity());
	EntityHandle* ptr = &entities_out[0];
	iGeom_getEntities(mInstance, set, type, &ptr, &alloc, &size, &err);
	entities_out.resize(size);

	// if input vector was too small, try again with increased size
	if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
		alloc = entities_out.size();
		ptr = &entities_out[0];
		iGeom_getEntities(mInstance, set, type, &ptr, &alloc, &size, &err);
	}

	return (Error) err;
}

inline iGeom::Error iGeom::deleteEnt(EntityHandle handle)
{
	int err;
	iGeom_deleteEnt(mInstance, handle, &err);
	return (Error) err;
}

//inline iGeom::Error
//iGeom::deleteEntArr( const EntityHandle* entity_handles, int num_handles )
//{
//  int err;
//  iGeom_deleteEntArr( mInstance, entity_handles, num_handles, &err );
//  return (Error)err;
//}

//inline iGeom::Error
//iGeom::getAdjEntities( EntitySetHandle set,
//                       EntityType type_requestor,
//                       EntityType type_requested,
//                       std::vector<EntityHandle>& adj_entity_handles,
//                       std::vector<int>& offset )
//{
//  std::vector<EntityHandle> entities;
//  Error err = getEntities( set, type_requestor, entities );
//  if (iBase_SUCCESS != err)
//    return err;
//
//  offset.resize( entities.size() + 1 );
//  return getArrAdj( &entities[0], entities.size(), type_requested,
//                    adj_entity_handles, &offset[0] );
//}

inline iGeom::Error iGeom::initEntIter(EntitySetHandle set,
		EntityType requested_type, iGeom::EntIter& iter)
{
	int err;
	iter.mInstance = mInstance;
	iGeom_initEntIter(mInstance, set, requested_type, &iter.mHandle, &err);
	return (Error) err;
}

inline iGeom::Error iGeom::initEntArrIter(EntitySetHandle set,
		EntityType requested_type, int requested_array_size,
		iGeom::EntArrIter& iter)
{
	int err;
	iter.mInstance = mInstance;
	iter.mSize = requested_array_size;
	iGeom_initEntArrIter(mInstance, set, requested_type, requested_array_size,
			&iter.mHandle, &err);
	return (Error) err;
}

inline iGeom::EntArrIter::~EntArrIter()
{
	int err;
	if (mHandle != 0) {
		iGeom_endEntArrIter(mInstance, mHandle, &err);
		mHandle = 0;
	}
}

inline iGeom::EntIter::~EntIter()
{
	int err;
	if (mHandle != 0) {
		iGeom_endEntIter(mInstance, mHandle, &err);
		mHandle = 0;
	}
}

inline iGeom::Error iGeom::EntArrIter::getNext(EntityHandle* entity_handles,
		int& size_out, bool& has_more_data_out)
{
	int err, alloc = mSize, has_data;
	iGeom_getNextEntArrIter(mInstance, mHandle, &entity_handles, &alloc,
			&size_out, &has_data, &err);
	has_more_data_out = (has_data != 0);
	return (Error) err;
}

inline iGeom::Error iGeom::EntIter::getNext(EntityHandle& handle_out,
		bool& has_more_data_out)
{
	int err, has_data;
	iGeom_getNextEntIter(mInstance, mHandle, &handle_out, &has_data, &err);
	has_more_data_out = (has_data != 0);
	return (Error) err;
}

inline iGeom::Error iGeom::EntArrIter::reset()
{
	int err;
	iGeom_resetEntArrIter(mInstance, mHandle, &err);
	return (Error) err;
}

inline iGeom::Error iGeom::EntIter::reset()
{
	int err;
	iGeom_resetEntIter(mInstance, mHandle, &err);
	return (Error) err;
}

inline iGeom::Error iGeom::getEntType(EntityHandle handle, EntityType& type_out) const
{
	int err, result;
	iGeom_getEntType(mInstance, handle, &result, &err);
	type_out = (EntityType) result;
	return (Error) err;
}

inline iGeom::Error iGeom::getArrType(const EntityHandle* entity_handles,
		int entity_handles_size, EntityType* types_out) const
{
	int err, alloc = entity_handles_size, junk, *ptr;
	std::vector<int> storage;
	if (sizeof(EntityType) == sizeof(int))
		ptr = reinterpret_cast<int*> (types_out);
	else {
		storage.resize(entity_handles_size);
		ptr = &storage[0];
	}

	iGeom_getArrType(mInstance, entity_handles, entity_handles_size, &ptr,
			&alloc, &junk, &err);

	if (sizeof(EntityType) != sizeof(int))
		for (int i = 0; i < entity_handles_size; ++i)
			types_out[i] = (EntityType) storage[i];

	return (Error) err;
}

inline iGeom::Error iGeom::getEntAdj(EntityHandle handle,
		EntityType type_requested, std::vector<EntityHandle>& adj_entities_out) const
{
	if (adj_entities_out.capacity() == 0)
		adj_entities_out.resize(12);
	else
		adj_entities_out.resize(adj_entities_out.capacity());

	int err, alloc = adj_entities_out.size(), size = 0;
	EntityHandle* ptr = &adj_entities_out[0];
	iGeom_getEntAdj(mInstance, handle, type_requested, &ptr, &alloc, &size, &err);
	adj_entities_out.resize(size);

	if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
		alloc = adj_entities_out.size();
		ptr = &adj_entities_out[0];
		iGeom_getEntAdj(mInstance, handle, type_requested, &ptr, &alloc, &size,
				&err);
	}

	return (Error) err;
}

inline iGeom::Error iGeom::getArrAdj(const EntityHandle* entity_handles,
		int entity_handles_size, EntityType type_requested, std::vector<
		EntityHandle>& adj_entities_out, int* offsets_out) const
{
	if (adj_entities_out.capacity() == 0)
		adj_entities_out.resize(12 * entity_handles_size);
	else
		adj_entities_out.resize(adj_entities_out.capacity());

	int err, alloc = adj_entities_out.size(), size = 0;
	int off_alloc = entity_handles_size + 1, junk;
	EntityHandle* ptr = &adj_entities_out[0];
	iGeom_getArrAdj(mInstance, entity_handles, entity_handles_size,
			type_requested, &ptr, &alloc, &size, &offsets_out, &off_alloc, &junk,
			&err);
	adj_entities_out.resize(size);

	if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
		alloc = adj_entities_out.size();
		ptr = &adj_entities_out[0];
		iGeom_getArrAdj(mInstance, entity_handles, entity_handles_size,
				type_requested, &ptr, &alloc, &size, &offsets_out, &off_alloc, &junk,
				&err);
	}

	return (Error) err;
}

inline iGeom::Error iGeom::getEnt2ndAdj(EntityHandle handle,
		EntityType bridge_entity_type, EntityType type_requested, std::vector<
		EntityHandle>& adj_entities_out) const
{
	if (adj_entities_out.capacity() == 0)
		adj_entities_out.resize(12);
	else
		adj_entities_out.resize(adj_entities_out.capacity());

	int err, alloc = adj_entities_out.size(), size = 0;
	EntityHandle* ptr = &adj_entities_out[0];
	iGeom_getEnt2ndAdj(mInstance, handle, bridge_entity_type, type_requested,
			&ptr, &alloc, &size, &err);
	adj_entities_out.resize(size);

	if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
		alloc = adj_entities_out.size();
		ptr = &adj_entities_out[0];
		iGeom_getEnt2ndAdj(mInstance, handle, bridge_entity_type, type_requested,
				&ptr, &alloc, &size, &err);
	}

	return (Error) err;
}

inline iGeom::Error iGeom::getArr2ndAdj(const EntityHandle* entity_handles,
		int entity_handles_size, EntityType order_key, EntityType type_requested,
		std::vector<EntityHandle>& adj_entities_out, int* offsets_out) const
{
	if (adj_entities_out.capacity() == 0)
		adj_entities_out.resize(12 * entity_handles_size);
	else
		adj_entities_out.resize(adj_entities_out.capacity());

	int err, alloc = adj_entities_out.size(), size = 0;
	int off_alloc = entity_handles_size + 1, junk;
	EntityHandle* ptr = &adj_entities_out[0];
	iGeom_getArr2ndAdj(mInstance, entity_handles, entity_handles_size, order_key,
			type_requested, &ptr, &alloc, &size, &offsets_out, &off_alloc, &junk,
			&err);
	adj_entities_out.resize(size);

	if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
		alloc = adj_entities_out.size();
		ptr = &adj_entities_out[0];
		iGeom_getArr2ndAdj(mInstance, entity_handles, entity_handles_size,
				order_key, type_requested, &ptr, &alloc, &size, &offsets_out,
				&off_alloc, &junk, &err);
	}

	return (Error) err;
}

inline iGeom::Error iGeom::getBoundBox(double& min_x, double& min_y,
		double& min_z, double& max_x, double& max_y, double& max_z) const
{
	int err;
	iGeom_getBoundBox(mInstance, &min_x, &min_y, &min_z, &max_x, &max_y, &max_z,
			&err);
	return (Error) err;
}

inline int iGeom::getParametric()
{
	int err, result;
	iGeom_getParametric(mInstance, &result, &err);
	return result;
}

inline iGeom::Error iGeom::isEntAdj(EntityHandle entity1, EntityHandle entity2,
		bool& adjacent_out) const
{
	int err, result;
	iGeom_isEntAdj(mInstance, entity1, entity2, &result, &err);
	adjacent_out = (result != 0);
	return (Error) err;
}

inline iGeom::Error iGeom::isArrAdj(const EntityHandle* entities1,
		const EntityHandle* entities2, int num_entity_pairs, int* is_adj_out) const
{
	int err, alloc = num_entity_pairs, size = 0;
	iGeom_isArrAdj(mInstance, entities1, num_entity_pairs, entities2,
			num_entity_pairs, &is_adj_out, &alloc, &size, &err);
	return (Error) err;
}

inline iGeom::Error iGeom::getEntClosestPt(EntityHandle entity, double near_x,
		double near_y, double near_z, double& on_x, double& on_y, double& on_z) const
{
	int err;
	iGeom_getEntClosestPt(mInstance, entity, near_x, near_y, near_z, &on_x,
			&on_y, &on_z, &err);
	return (Error) err;
}

inline iGeom::Error iGeom::getEntClosestPtTrimmed(EntityHandle entity, double near_x,
                                          double near_y, double near_z, double& on_x, double& on_y, double& on_z) const
{
    int err;
    iGeom_getEntClosestPtTrimmed(mInstance, entity, near_x, near_y, near_z, &on_x, &on_y, &on_z, &err);
    return (Error) err;
}

inline iGeom::Error iGeom::getArrClosestPt(const EntityHandle* handles,
          int handles_size, StorageOrder order, const double* near_coordinates,
          int near_coordinates_size, double* on_coordinates) const
{
     int err, alloc = std::max(near_coordinates_size, 3 * handles_size), size = 0;
     iGeom_getArrClosestPt(mInstance, handles, handles_size, order,
                           near_coordinates, near_coordinates_size, &on_coordinates, &alloc, &size,
                           &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntNrmlXYZ(EntityHandle entity, double x,
          double y, double z, double& i, double& j, double& k) const
{
     int err;
     iGeom_getEntNrmlXYZ(mInstance, entity, x, y, z, &i, &j, &k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrNrmlXYZ(const EntityHandle* entities,
          int entities_size, StorageOrder order, const double* xyz, int xyz_size,
          double* ijk) const
{
     int err, alloc = std::max(xyz_size, 3 * entities_size), size = 0;
     iGeom_getArrNrmlXYZ(mInstance, entities, entities_size, order, xyz, xyz_size,
                         &ijk, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntNrmlPlXYZ(EntityHandle entity, double x,
          double y, double z, double& on_x, double& on_y, double& on_z, double& i,
          double& j, double& k) const
{
     int err;
     iGeom_getEntNrmlPlXYZ(mInstance, entity, x, y, z, &on_x, &on_y, &on_z, &i,
                           &j, &k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrNrmlPlXYZ(const EntityHandle* entities,
          int entities_size, StorageOrder order, const double* near_xyz,
          int near_xyz_size, double* on_xyz, double* nrml_ijk) const
{
     int err, alloc = std::max(near_xyz_size, 3 * entities_size), size = 0;
     iGeom_getArrNrmlPlXYZ(mInstance, entities, entities_size, order, near_xyz,
                           near_xyz_size, &on_xyz, &alloc, &size, &nrml_ijk, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntTgntXYZ(EntityHandle entity, double x,
          double y, double z, double& i, double& j, double& k) const
{
     int err;
     iGeom_getEntTgntXYZ(mInstance, entity, x, y, z, &i, &j, &k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrTgntXYZ(const EntityHandle* entities,
          int entities_size, StorageOrder order, const double* xyz, int xyz_size,
          double* ijk) const
{
     int err, alloc = std::max(xyz_size, 3 * entities_size), size = 0;
     iGeom_getArrTgntXYZ(mInstance, entities, entities_size, order, xyz, xyz_size,
                         &ijk, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getFcCvtrXYZ(EntityHandle face, double x, double y,
                                        double z, double& i1, double& j1, double& k1, double& i2, double& j2,
                                        double& k2) const
{
     int err;
     iGeom_getFcCvtrXYZ(mInstance, face, x, y, z, &i1, &j1, &k1, &i2, &j2, &k2,
                        &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEgCvtrXYZ(EntityHandle edge, double x, double y,
                                        double z, double& i, double& j, double& k) const
{
     int err;
     iGeom_getEgCvtrXYZ(mInstance, edge, x, y, z, &i, &j, &k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntArrCvtrXYZ(const EntityHandle* entities,
          int entities_size, StorageOrder order, const double* xyz, int xyz_size,
          double* cvtr_1, double* cvtr_2) const
{
     int err, alloc = std::max(xyz_size, 3 * entities_size), size = 0;
     iGeom_getEntArrCvtrXYZ(mInstance, entities, entities_size, order, xyz,
                            xyz_size, &cvtr_1, &alloc, &size, &cvtr_2, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEgEvalXYZ(EntityHandle edge, double x, double y,
                                        double z, double& on_x, double& on_y, double& on_z, double& tngt_i,
                                        double& tngt_j, double& tngt_k, double& cvtr_i, double& cvtr_j,
                                        double& cvtr_k) const
{
     int err;
     iGeom_getEgEvalXYZ(mInstance, edge, x, y, z, &on_x, &on_y, &on_z, &tngt_i,
                        &tngt_j, &tngt_k, &cvtr_i, &cvtr_j, &cvtr_k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getFcEvalXYZ(EntityHandle face, double x, double y,
                                        double z, double& on_x, double& on_y, double& on_z, double& tngt_i,
                                        double& tngt_j, double& tngt_k, double& cvtr1_i, double& cvtr1_j,
                                        double& cvtr1_k, double& cvtr2_i, double& cvtr2_j, double& cvtr2_k) const
{
     int err;
     iGeom_getFcEvalXYZ(mInstance, face, x, y, z, &on_x, &on_y, &on_z, &tngt_i,
                        &tngt_j, &tngt_k, &cvtr1_i, &cvtr1_j, &cvtr1_k, &cvtr2_i, &cvtr2_j,
                        &cvtr2_k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrEgEvalXYZ(const EntityHandle* edges,
          int edges_size, StorageOrder order, const double* near_coords,
          int near_coords_size, double* on_coords, double* tangent, double* curvature) const
{
     int err, alloc = std::max(near_coords_size, 3 * edges_size), size = 0;
     iGeom_getArrEgEvalXYZ(mInstance, edges, edges_size, order, near_coords,
                           near_coords_size, &on_coords, &alloc, &size, &tangent, &alloc, &size,
                           &curvature, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrFcEvalXYZ(const EntityHandle* faces,
          int faces_size, StorageOrder order, const double* near_coords,
          int near_coords_size, double* on_coords, double* tangent,
          double* curvature1, double* curvature2) const
{
     int err, alloc = std::max(near_coords_size, 3 * faces_size), size = 0;
     iGeom_getArrFcEvalXYZ(mInstance, faces, faces_size, order, near_coords,
                           near_coords_size, &on_coords, &alloc, &size, &tangent, &alloc, &size,
                           &curvature1, &alloc, &size, &curvature2, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntBoundBox(EntityHandle entity, double& min_x,
          double& min_y, double& min_z, double& max_x, double& max_y, double& max_z) const
{
     int err;
     iGeom_getEntBoundBox(mInstance, entity, &min_x, &min_y, &min_z, &max_x,
                          &max_y, &max_z, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrBoundBox(const EntityHandle* entities,
          int entities_size, StorageOrder order, double* min_corners,
          double* max_corners) const
{
     int err, alloc = 3 * entities_size, size = 0, order_int = order;
     iGeom_getArrBoundBox(mInstance, entities, entities_size, order_int,
                          &min_corners, &alloc, &size, &max_corners, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getVtxCoord(EntityHandle vertex, double& x,
                                       double& y, double& z) const
{
     int err;
     iGeom_getVtxCoord(mInstance, vertex, &x, &y, &z, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getVtxArrCoords(const EntityHandle* vertices,
          int vertices_size, StorageOrder order, double* coords) const
{
     int err, alloc = vertices_size, size = 0;
     iGeom_getVtxArrCoords(mInstance, vertices, vertices_size, order, &coords,
                           &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getPntRayIntsct(double x, double y, double z,
          double i, double j, double k, StorageOrder order,
          std::vector<EntityHandle>& entities_out, std::vector<double>& points_out,
          std::vector<double>& params_out) const
{
     int err, count;
     Error err2 = getNumOfType(getRootSet(), iBase_ALL_TYPES, count);
     if (err2 != iBase_SUCCESS)
          return err2;

     entities_out.resize(count);
     points_out.resize(3 * count);
     params_out.resize(2 * count);
     int entities_alloc = entities_out.size(), entities_size = 0;
     int points_alloc = points_out.size(), points_size = 0;
     int params_alloc = params_out.size(), params_size = 0;
     EntityHandle* entities_ptr = &entities_out[0];
     double * points_ptr = &points_out[0];
     double * params_ptr = &params_out[0];

     iGeom_getPntRayIntsct(mInstance, x, y, z, i, j, k, &entities_ptr,
                           &entities_alloc, &entities_size, order, &points_ptr, &points_alloc,
                           &points_size, &params_ptr, &params_alloc, &params_size, &err);
     entities_out.resize(entities_size);
     points_out.resize(points_size);
     params_out.resize(params_size);
     if (err == iBase_BAD_ARRAY_SIZE || err == iBase_BAD_ARRAY_DIMENSION) {
          entities_alloc = entities_out.size();
          points_alloc = points_out.size();
          params_alloc = params_out.size();
          entities_ptr = &entities_out[0];
          points_ptr = &points_out[0];
          params_ptr = &params_out[0];
          iGeom_getPntRayIntsct(mInstance, x, y, z, i, j, k, &entities_ptr,
                                &entities_alloc, &entities_size, order, &points_ptr, &points_alloc,
                                &points_size, &params_ptr, &params_alloc, &params_size, &err);
     }

     return (Error) err;
}

inline iGeom::Error iGeom::getPntClsf(double x, double y, double z,
                                      EntityHandle& handle_out) const
{
     int err;
     iGeom_getPntClsf(mInstance, x, y, z, &handle_out, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getPntArrClsf(StorageOrder order,
          const double* coords, int coords_size, EntityHandle* entities_out) const
{
     int err, alloc = coords_size / 3, size = 0;
     iGeom_getPntArrClsf(mInstance, order, coords, coords_size, &entities_out,
                         &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getSense(EntityHandle ent, EntityHandle wrt_ent,
                                    int &sense) const
{
     EntityType tp_wrt, tp_ent;
     Error err = getEntType(wrt_ent, tp_wrt);
     if (iBase_SUCCESS != err)
          return err;
     err = getEntType(ent, tp_ent);
     if (iBase_SUCCESS != err)
          return err;
     if (tp_wrt - tp_ent != 1)
          return iBase_FAILURE;
     switch (tp_wrt) {
     case iBase_REGION:
          return getEntNrmlSense(ent, wrt_ent, sense);
          break;
     case iBase_FACE:
          return getEgFcSense(ent, wrt_ent, sense);
     case iBase_EDGE:
          return getEgVtxSense(wrt_ent, ent, ent, sense);
     case iBase_VERTEX:
     case iBase_ALL_TYPES:
          return iBase_FAILURE;
     }
     return iBase_FAILURE;
}

inline iGeom::Error iGeom::getArrSense(const EntityHandle *ent, int num_ents,
                                       EntityHandle wrt_ent, int *sense) const
{
     EntityType tp_wrt, tp_ent;
     Error err = getEntType(wrt_ent, tp_wrt);
     if (iBase_SUCCESS != err)
          return err;
     err = getEntType(ent[0], tp_ent);
     if (iBase_SUCCESS != err)
          return err;
     if (tp_wrt - tp_ent != 1)
          return iBase_FAILURE;
     std::vector<EntityHandle> dum_wrts(num_ents, wrt_ent);
     switch (tp_wrt) {
     case iBase_REGION:
          return getArrNrmlSense(ent, num_ents, &dum_wrts[0], num_ents, sense);
          break;
     case iBase_FACE:
          return getEgFcArrSense(ent, num_ents, &dum_wrts[0], num_ents, sense);
          break;
     case iBase_EDGE:
          return getEgVtxArrSense(&dum_wrts[0], num_ents, ent, num_ents, ent,
                                  num_ents, sense);
          break;
     case iBase_VERTEX:
     case iBase_ALL_TYPES:
          return iBase_FAILURE;
          break;
     }
     return iBase_FAILURE;
}

inline iGeom::Error iGeom::getEntNrmlSense(EntityHandle face,
          EntityHandle region, int& sense) const
{
     int err;
     iGeom_getEntNrmlSense(mInstance, face, region, &sense, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEgFcSense(EntityHandle edge, EntityHandle face,
                                        int& sense) const
{
     int err;
     iGeom_getEgFcSense(mInstance, edge, face, &sense, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEgVtxSense(EntityHandle edge, EntityHandle vtx1,
          EntityHandle vtx2, int& sense) const
{
     int err;
     iGeom_getEgVtxSense(mInstance, edge, vtx1, vtx2, &sense, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrNrmlSense(const EntityHandle* faces,
          int faces_size, const EntityHandle* vols, int vols_size, int* senses_out) const
{
     int err, alloc = std::max(vols_size, faces_size), size = 0;
     iGeom_getArrNrmlSense(mInstance, faces, faces_size, vols, vols_size,
                           &senses_out, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEgFcArrSense(const EntityHandle* edges,
          int edges_size, const EntityHandle* faces, int faces_size, int* senses_out) const
{
     int err, alloc = std::max(edges_size, faces_size), size = 0;
     iGeom_getEgFcArrSense(mInstance, edges, edges_size, faces, faces_size,
                           &senses_out, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEgVtxArrSense(const EntityHandle* edges,
          int edges_size, const EntityHandle* vertices1, int vertices1_size,
          const EntityHandle* vertices2, int vertices2_size, int* senses_out) const
{
     int err, alloc = std::max(vertices1_size,
                               std::max(vertices2_size, edges_size)), size = 0;
     iGeom_getEgVtxArrSense(mInstance, edges, edges_size, vertices1,
                            vertices1_size, vertices2, vertices2_size, &senses_out, &alloc, &size,
                            &err);
     return (Error) err;
}

inline iGeom::Error iGeom::measure(const EntityHandle* entities,
                                   int entities_size, double* measures) const
{
     int err, alloc = entities_size, size = 0;
     iGeom_measure(mInstance, entities, entities_size, &measures, &alloc, &size,
                   &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getFaceType(EntityHandle face, std::string& type) const
{
     char buffer[1024];
     int err, len = sizeof(buffer);
     iGeom_getFaceType(mInstance, face, buffer, &err, &len);
     type = std::string(buffer, len);
     return (Error) err;
}

inline iGeom::Error iGeom::isEntParametric(EntityHandle entity,
          bool& parametric) const
{
     int err, result;
     iGeom_isEntParametric(mInstance, entity, &result, &err);
     parametric = (result != 0);
     return (Error) err;
}

inline iGeom::Error iGeom::isArrParametric(const EntityHandle* entities,
          int entities_size, int* is_parametric) const
{
     int err, alloc = entities_size, size = 1;
     iGeom_isArrParametric(mInstance, entities, entities_size, &is_parametric,
                           &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntUVtoXYZ(EntityHandle face, double u, double v,
          double& x, double& y, double& z) const
{
     int err;
     iGeom_getEntUVtoXYZ(mInstance, face, u, v, &x, &y, &z, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntUtoXYZ(EntityHandle edge, double u, double& x,
                                        double& y, double& z) const
{
     int err;
     iGeom_getEntUtoXYZ(mInstance, edge, u, &x, &y, &z, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrUVtoXYZ(const EntityHandle* faces,
          int faces_size, StorageOrder order, const double* uv, int uv_size,
          double* xyz) const
{
     int err, alloc = std::max(3 * uv_size / 2, 3 * faces_size), size = 0;
     iGeom_getArrUVtoXYZ(mInstance, faces, faces_size, order, uv, uv_size, &xyz,
                         &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrUtoXYZ(const EntityHandle* edges,
                                        int edges_size, const double* u, int u_size, StorageOrder order,
                                        double* xyz) const
{
     int err, alloc = std::max(3 * u_size, 3 * edges_size), size = 0;
     iGeom_getArrUtoXYZ(mInstance, edges, edges_size, u, u_size, order, &xyz,
                        &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntXYZtoUV(EntityHandle face, double x, double y,
          double z, double& u, double& v) const
{
     int err;
     iGeom_getEntXYZtoUV(mInstance, face, x, y, z, &u, &v, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntXYZtoU(EntityHandle edge, double x, double y,
                                        double z, double& u) const
{
     int err;
     iGeom_getEntXYZtoU(mInstance, edge, x, y, z, &u, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrXYZtoUV(const EntityHandle* faces,
          int faces_size, StorageOrder order, const double* coords, int coords_size,
          double* uv) const
{
     int err, alloc = std::max(2 * coords_size / 3, 2 * faces_size), size = 0;
     iGeom_getArrXYZtoUV(mInstance, faces, faces_size, order, coords, coords_size,
                         &uv, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrXYZtoU(const EntityHandle* edges,
                                        int edges_size, StorageOrder order, const double* coords, int coords_size,
                                        double* u) const
{
     int err, alloc = std::max(coords_size / 3, edges_size), size = 0;
     iGeom_getArrXYZtoU(mInstance, edges, edges_size, order, coords, coords_size,
                        &u, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntXYZtoUVHint(EntityHandle face, double x,
          double y, double z, double& u, double& v) const
{
     int err;
     iGeom_getEntXYZtoUVHint(mInstance, face, x, y, z, &u, &v, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrXYZtoUVHint(const EntityHandle* faces,
          int faces_size, StorageOrder order, const double* coords, int coords_size,
          double* uv) const
{
     int err, alloc = std::max(2 * coords_size / 3, 2 * faces_size), size = 0;
     iGeom_getArrXYZtoUVHint(mInstance, faces, faces_size, order, coords,
                             coords_size, &uv, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntNrmlUV(EntityHandle face, double u, double v,
                                        double& i, double& j, double& k) const
{
     int err;
     iGeom_getEntNrmlUV(mInstance, face, u, v, &i, &j, &k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrNrmlUV(const EntityHandle* faces,
                                        int faces_size, StorageOrder order, const double* uv, int uv_size,
                                        double* normals) const
{
     int err, alloc = std::max(3 * uv_size / 2, 3 * faces_size), size = 0;
     iGeom_getArrNrmlUV(mInstance, faces, faces_size, order, uv, uv_size,
                        &normals, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntTgntU(EntityHandle edge, double u, double& i,
                                       double& j, double& k) const
{
     int err;
     iGeom_getEntTgntU(mInstance, edge, u, &i, &j, &k, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrTgntU(const EntityHandle* edges,
                                       int edges_size, StorageOrder order, const double* u, int u_size,
                                       double* normals) const
{
     int err, alloc = std::max(3 * u_size, 3 * edges_size), size = 0;
     iGeom_getArrTgntU(mInstance, edges, edges_size, order, u, u_size, &normals,
                       &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEnt1stDrvt(EntityHandle handle, double u,
          double v, double& du_i, double& du_j, double& du_k, double& dv_i,
          double& dv_j, double& dv_k) const
{
     int err, du_alloc = 3, dv_alloc = 3, du_size = 0, dv_size = 0;
     double du[3], dv[3], *du_ptr = du, *dv_ptr = dv;
     iGeom_getEnt1stDrvt(mInstance, handle, u, v, &du_ptr, &du_alloc, &du_size,
                         &dv_ptr, &dv_alloc, &dv_size, &err);
     du_i = du[0];
     du_j = du[1];
     du_k = du[2];
     dv_i = dv[0];
     dv_j = dv[1];
     dv_k = dv[2];
     return (Error) err;
}

inline iGeom::Error iGeom::getEnt2ndDrvt(EntityHandle handle, double u,
          double v, double& duu_i, double& duu_j, double& duu_k, double& dvv_i,
          double& dvv_j, double& dvv_k, double& duv_i, double& duv_j, double& duv_k) const
{
     int err, uu_alloc = 3, uv_alloc = 3, vv_alloc = 3, uu_size = 0, uv_size = 0,
                                       vv_size = 0;
     double uu[3], uv[3], vv[3], *uu_ptr = uu, *vv_ptr = vv, *uv_ptr = uv;
     iGeom_getEnt2ndDrvt(mInstance, handle, u, v, &uu_ptr, &uu_alloc, &uu_size,
                         &vv_ptr, &vv_alloc, &vv_size, &uv_ptr, &uv_alloc, &uv_size, &err);
     duu_i = uu[0];
     duu_j = uu[1];
     duu_k = uu[2];
     dvv_i = vv[0];
     dvv_j = vv[1];
     dvv_k = vv[2];
     duv_i = uv[0];
     duv_j = uv[1];
     duv_k = uv[2];
     return (Error) err;
}

inline iGeom::Error iGeom::getArr1stDrvt(const EntityHandle* entities,
          int entities_size, StorageOrder order, const double* uv, int uv_size,
          double* dvtr_u, double* dvtr_v) const
{
     int err, allocu = std::max(3 * uv_size / 2, 3 * entities_size), sizeu = 0,
                       allocv = allocu, sizev = 0;
     std::vector<int> offset1(std::max(uv_size / 2, entities_size) + 1);
     std::vector<int> offset2(std::max(uv_size / 2, entities_size) + 1);
     int alloc1 = offset1.size(), size1 = 0, *ptr1 = &offset1[0];
     int alloc2 = offset2.size(), size2 = 0, *ptr2 = &offset2[0];
     iGeom_getArr1stDrvt(mInstance, entities, entities_size, order, uv, uv_size,
                         &dvtr_u, &allocu, &sizeu, &ptr1, &alloc1, &size1, &dvtr_v, &allocv,
                         &sizev, &ptr2, &alloc2, &size2, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArr2ndDrvt(const EntityHandle* entities,
          int entities_size, StorageOrder order, const double* uv, int uv_size,
          double* dvtr_uu, double* dvtr_vv, double* dvtr_uv) const
{
     int err, allocuu = std::max(3 * uv_size / 2, 3 * entities_size), sizeuu = 0,
                        allocvv = allocuu, sizevv = 0, allocuv = allocuu, sizeuv = 0;
     std::vector<int> offset1(std::max(uv_size / 2, entities_size) + 1);
     std::vector<int> offset2(std::max(uv_size / 2, entities_size) + 1);
     std::vector<int> offset3(std::max(uv_size / 2, entities_size) + 1);
     int alloc1 = offset1.size(), size1 = 0, *ptr1 = &offset1[0];
     int alloc2 = offset2.size(), size2 = 0, *ptr2 = &offset2[0];
     int alloc3 = offset3.size(), size3 = 0, *ptr3 = &offset3[0];
     iGeom_getArr2ndDrvt(mInstance, entities, entities_size, order, uv, uv_size,
                         &dvtr_uu, &allocuu, &sizeuu, &ptr1, &alloc1, &size1, &dvtr_vv, &allocvv,
                         &sizevv, &ptr2, &alloc2, &size2, &dvtr_uv, &allocuv, &sizeuv, &ptr3,
                         &alloc3, &size3, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getFcCvtrUV(EntityHandle face, double u, double v,
                                       double& i1, double& j1, double& k1, double& i2, double& j2, double& k2) const
{
     int err;
     iGeom_getFcCvtrUV(mInstance, face, u, v, &i1, &j1, &k1, &i2, &j2, &k2, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getFcArrCvtrUV(const EntityHandle* faces,
          int faces_size, StorageOrder order, const double* uv, int uv_size,
          double* cvtr1, double* cvtr2) const
{
     int err, alloc = std::max(3 * uv_size / 2, 3 * faces_size), size = 0;
     iGeom_getFcArrCvtrUV(mInstance, faces, faces_size, order, uv, uv_size,
                          &cvtr1, &alloc, &size, &cvtr2, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::isEntPeriodic(EntityHandle entity, bool& in_u,
          bool& in_v) const
{
     int err, u, v;
     iGeom_isEntPeriodic(mInstance, entity, &u, &v, &err);
     in_u = (u != 0);
     in_v = (v != 0);
     return (Error) err;
}

inline iGeom::Error iGeom::isArrPeriodic(const EntityHandle* entities,
          int entities_size, int* in_uv) const
{
     int err, alloc = 2 * entities_size, size = 0;
     iGeom_isArrPeriodic(mInstance, entities, entities_size, &in_uv, &alloc,
                         &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::isFcDegenerate(EntityHandle face,
          bool& is_degenerate) const
{
     int err, result;
     iGeom_isFcDegenerate(mInstance, face, &result, &err);
     is_degenerate = (result != 0);
     return (Error) err;
}

inline iGeom::Error iGeom::isFcArrDegenerate(const EntityHandle* faces,
          int faces_size, int* degenerate) const
{
     int err, alloc = faces_size, size = 0;
     iGeom_isFcArrDegenerate(mInstance, faces, faces_size, &degenerate, &alloc,
                             &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getTolerance(int& type_out, double& tolerance_out) const
{
     int err;
     iGeom_getTolerance(mInstance, &type_out, &tolerance_out, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntTolerance(EntityHandle entity,
          double& tolerance) const
{
     int err;
     iGeom_getEntTolerance(mInstance, entity, &tolerance, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrTolerance(const EntityHandle* entities,
          int entities_size, double* tolerances) const
{
     int err, alloc = entities_size, size = 0;
     iGeom_getArrTolerance(mInstance, entities, entities_size, &tolerances,
                           &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntUVRange(EntityHandle face, double& u_min,
          double& v_min, double& u_max, double& v_max) const
{
     int err;
     iGeom_getEntUVRange(mInstance, face, &u_min, &v_min, &u_max, &v_max, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntURange(EntityHandle edge, double& u_min,
                                        double& u_max) const
{
     int err;
     iGeom_getEntURange(mInstance, edge, &u_min, &u_max, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrUVRange(const EntityHandle* faces,
          int faces_size, StorageOrder order, double* uv_min, double* uv_max) const
{
     int err, alloc = faces_size, size = 0;
     iGeom_getArrUVRange(mInstance, faces, faces_size, order, &uv_min, &alloc,
                         &size, &uv_max, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrURange(const EntityHandle* edges,
                                        int edges_size, double* u_min, double* u_max) const
{
     int err, alloc = edges_size, size = 0;
     iGeom_getArrURange(mInstance, edges, edges_size, &u_min, &alloc, &size,
                        &u_max, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getEntUtoUV(EntityHandle edge, EntityHandle face,
                                       double edge_u, double& face_u, double& face_v) const
{
     int err;
     iGeom_getEntUtoUV(mInstance, edge, face, edge_u, &face_u, &face_v, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getVtxToUV(EntityHandle vertex, EntityHandle face,
                                      double& u, double& v) const
{
     int err;
     iGeom_getVtxToUV(mInstance, vertex, face, &u, &v, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getVtxToU(EntityHandle vertex, EntityHandle edge,
                                     double& u) const
{
     int err;
     iGeom_getVtxToU(mInstance, vertex, edge, &u, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getArrUtoUV(const EntityHandle* edges,
                                       int edges_size, const EntityHandle* faces, int faces_size,
                                       const double* edge_u, int edge_u_size, StorageOrder order, double* face_uv) const
{
     int err, alloc = std::max(edge_u_size, std::max(edges_size, faces_size));
     int size = 0;
     iGeom_getArrUtoUV(mInstance, edges, edges_size, faces, faces_size, edge_u,
                       edge_u_size, order, &face_uv, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getVtxArrToUV(const EntityHandle* vertices,
          int vertices_size, const EntityHandle* faces, int faces_size,
          StorageOrder order, double* face_uv) const
{
     int err, alloc = std::max(vertices_size, faces_size), size = 0;
     iGeom_getVtxArrToUV(mInstance, vertices, vertices_size, faces, faces_size,
                         order, &face_uv, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::getVtxArrToU(const EntityHandle* vertices,
                                        int vertices_size, const EntityHandle* edges, int edges_size,
                                        double* edge_u) const
{
     int err, alloc = std::max(vertices_size, edges_size), size = 0;
     iGeom_getVtxArrToU(mInstance, vertices, vertices_size, edges, edges_size,
                        &edge_u, &alloc, &size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::deleteAll()
{
     int err;
     iGeom_deleteAll(mInstance, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::copyEnt(EntityHandle source, EntityHandle& copy)
{
     int err;
     iGeom_copyEnt(mInstance, source, &copy, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::createSphere(double radius, EntityHandle& sphere)
{
     int err;
     iGeom_createSphere(mInstance, radius, &sphere, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::createPrism(double height, int num_sides,
                                       double maj_radius, double min_radius, EntityHandle& prism)
{
     int err;
     iGeom_createPrism(mInstance, height, num_sides, maj_radius, min_radius,
                       &prism, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::createBrick(double x, double y, double z,
                                       EntityHandle& brick)
{
     int err;
     iGeom_createBrick(mInstance, x, y, z, &brick, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::createCylinder(double height, double maj_rad,
          double min_rad, EntityHandle& cylinder)
{
     int err;
     iGeom_createCylinder(mInstance, height, maj_rad, min_rad, &cylinder, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::createCone(double height, double maj_rad,
          double min_rad, double rad_top, EntityHandle& cone)
{
     int err;
     iGeom_createCone(mInstance, height, maj_rad, min_rad, rad_top, &cone, &err);
     return (Error) err;
}


inline iGeom::Error iGeom::createTorus(double maj_rad, double min_rad,
                                       EntityHandle& torus)
{
     int err;
     iGeom_createTorus(mInstance, maj_rad, min_rad, &torus, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::moveEnt(EntityHandle entity, double x, double y,
                                   double z)
{
     int err;
     iGeom_moveEnt(mInstance, entity, x, y, z, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::rotateEnt(EntityHandle entity, double angle,
                                     double axis_x, double axis_y, double axis_z)
{
     int err;
     iGeom_rotateEnt(mInstance, entity, angle, axis_x, axis_y, axis_z, &err);
     return (Error) err;
}
inline iGeom::Error iGeom::reflectEnt(EntityHandle entity, double x,
                                      double y, double z, double norm_x,
                                      double norm_y, double norm_z)
{
     int err;
     iGeom_reflectEnt(mInstance, entity, x, y, z, norm_x, norm_y, norm_z, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::scaleEnt(EntityHandle entity, double x,
                                      double y, double z, double x_factor,
                                      double y_factor, double z_factor)
{
     int err;
     iGeom_scaleEnt(mInstance, entity, x, y, z, x_factor, y_factor, z_factor, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::uniteEnts(const EntityHandle* entities,
                                     int entities_size, EntityHandle& result_entity)
{
     int err;
     iGeom_uniteEnts(mInstance, entities, entities_size, &result_entity, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::subtractEnts(EntityHandle blank, EntityHandle tool,
                                        EntityHandle& result)
{
     int err;
     iGeom_subtractEnts(mInstance, blank, tool, &result, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::intersectEnts(EntityHandle entity1,
          EntityHandle entity2, EntityHandle& result)
{
     int err;
     iGeom_intersectEnts(mInstance, entity1, entity2, &result, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::sectionEnt(EntityHandle entity, double plane_x,
                                      double plane_y, double plane_z, double offset, bool reverse,
                                      EntityHandle& result)
{
     int err;
     iGeom_sectionEnt(mInstance, entity, plane_x, plane_y, plane_z, offset,
                      reverse, &result, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::sweepEntAboutAxis(EntityHandle entity, double angle,
          double axis_x, double axis_y, double axis_z, EntityHandle& swept_entity)
{
     int err;
     iGeom_sweepEntAboutAxis(mInstance, entity, angle, axis_x, axis_y, axis_z,
                             &swept_entity, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::imprintEnts(const EntityHandle* entities,
                                       int entities_size)
{
     int err;
     iGeom_imprintEnts(mInstance, entities, entities_size, &err);
     return (Error) err;
}

inline iGeom::Error iGeom::mergeEnts(const EntityHandle* entities,
                                     int entities_size, double tolerance)
{
     int err;
     iGeom_mergeEnts(mInstance, entities, entities_size, tolerance, &err);
     return (Error) err;
}


iGeom::Error
iGeom::getErrorType() const
{
     int err;
     iGeom_getErrorType( mInstance, &err );
     return (Error)err;
}

std::string
iGeom::getDescription() const
{
     std::vector<char> buffer(1024);
     iGeom_getDescription( mInstance, &buffer[0], buffer.size() );
     return std::string(&buffer[0]);
}




iGeom::EntitySetHandle
iGeom::getRootSet() const
{
     int err;
     EntitySetHandle result;
     iGeom_getRootSet( mInstance, &result, &err );
     return iBase_SUCCESS == err ? result : 0;
}



iGeom::Error
iGeom::createEntSet( bool is_list, EntitySetHandle& handle_out )
{
     int err;
     iGeom_createEntSet( mInstance, is_list, &handle_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::destroyEntSet( EntitySetHandle handle )
{
     int err;
     iGeom_destroyEntSet( mInstance, handle, &err );
     return (Error)err;
}

iGeom::Error
iGeom::isList( EntitySetHandle handle, bool& is_list )
{
     int err, result;
     iGeom_isList( mInstance, handle, &result, &err );
     is_list = (result != 0);
     return (Error)err;
}

iGeom::Error
iGeom::getNumEntSets( EntitySetHandle set, int num_hops, int& num_sets_out ) const
{
     int err;
     iGeom_getNumEntSets( mInstance, set, num_hops, &num_sets_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getEntSets( EntitySetHandle set, int num_hops,
                   std::vector<EntitySetHandle>& contained_sets_out ) const
{
     int err, count;
     iGeom_getNumEntSets( mInstance, set, num_hops, &count, &err );
     if (iBase_SUCCESS != err)
          return (Error)err;
     contained_sets_out.resize(count);
     int alloc = contained_sets_out.size(), size;
     EntitySetHandle* ptr = &contained_sets_out[0];
     iGeom_getEntSets( mInstance, set, num_hops, &ptr, &alloc, &size, &err );
     return (Error)err;
}

iGeom::Error
iGeom::addEntToSet( EntityHandle entity, EntitySetHandle set )
{
     int err;
     iGeom_addEntToSet( mInstance, entity,set, &err );
     return (Error)err;
}

iGeom::Error
iGeom::rmvEntFromSet( EntityHandle entity, EntitySetHandle set )
{
     int err;
     iGeom_rmvEntFromSet( mInstance, entity,set, &err );
     return (Error)err;
}

iGeom::Error
iGeom::addEntArrToSet( const EntityHandle* entity_handles,
                       int entity_handles_size,
                       EntitySetHandle entity_set )
{
     int err;
     iGeom_addEntArrToSet( mInstance, entity_handles, entity_handles_size, entity_set, &err );
     return (Error)err;
}

iGeom::Error
iGeom::rmvEntArrFromSet( const EntityHandle* entity_handles,
                         int entity_handles_size,
                         EntitySetHandle entity_set )
{
     int err;
     iGeom_rmvEntArrFromSet( mInstance, entity_handles, entity_handles_size, entity_set, &err );
     return (Error)err;
}

iGeom::Error
iGeom::addEntSet( EntitySetHandle to_add, EntitySetHandle add_to )
{
     int err;
     iGeom_addEntSet( mInstance, to_add, add_to, &err );
     return (Error)err;
}

iGeom::Error
iGeom::rmvEntSet( EntitySetHandle to_rmv, EntitySetHandle rmv_from )
{
     int err;
     iGeom_rmvEntSet( mInstance, to_rmv, rmv_from, &err );
     return (Error)err;
}

iGeom::Error
iGeom::isEntContained( EntitySetHandle set, EntityHandle ent, bool& contained_out ) const
{
     int err, result;
     iGeom_isEntContained( mInstance, set, ent, &result, &err );
     contained_out = (result != 0);
     return (Error)err;
}

iGeom::Error
iGeom::isEntArrContained( EntitySetHandle containing_set,
                          const EntityHandle* entity_handles,
                          int num_entity_handles,
                          bool* is_contained_out ) const
{
     int err, *ptr = 0, alloc = 0, size = 0;
     iGeom_isEntArrContained( mInstance, containing_set,
                              entity_handles, num_entity_handles,
                              &ptr, &alloc, &size, &err );
     if (iBase_SUCCESS != err)
          return (Error)err;
     for (int i = 0; i < num_entity_handles; ++i)
          is_contained_out[i] = (ptr[i] != 0);
     free(ptr);
     return iBase_SUCCESS;
}

iGeom::Error
iGeom::isEntSetContained( EntitySetHandle containing_set,
                          EntitySetHandle contained_set,
                          bool& contained_out ) const
{
     int err, result;
     iGeom_isEntSetContained( mInstance, containing_set, contained_set, &result, &err );
     contained_out = (result != 0);
     return (Error)err;
}

iGeom::Error
iGeom::addPrntChld( EntitySetHandle parent, EntitySetHandle child )
{
     int err;
     iGeom_addPrntChld( mInstance, parent, child, &err );
     return (Error)err;
}

iGeom::Error
iGeom::rmvPrntChld( EntitySetHandle parent, EntitySetHandle child )
{
     int err;
     iGeom_rmvPrntChld( mInstance, parent, child, &err );
     return (Error)err;
}

iGeom::Error
iGeom::isChildOf( EntitySetHandle parent, EntitySetHandle child, bool& is_child_out ) const
{
     int err, result;
     iGeom_isChildOf( mInstance, parent, child, &result, &err );
     is_child_out = (result != 0);
     return (Error)err;
}

iGeom::Error
iGeom::getNumChld( EntitySetHandle parent, int num_hops, int& num_child_out ) const
{
     int err;
     iGeom_getNumChld( mInstance, parent, num_hops, &num_child_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getNumPrnt( EntitySetHandle child, int num_hops, int& num_parent_out ) const
{
     int err;
     iGeom_getNumPrnt( mInstance, child, num_hops, &num_parent_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getChldn( EntitySetHandle parent, int num_hops,
                 std::vector<EntitySetHandle>& children_out ) const
{
     int err, count;
     iGeom_getNumChld( mInstance, parent, num_hops, &count, &err );
     if (iBase_SUCCESS != err)
          return (Error)err;
     children_out.resize(count);
     int alloc = children_out.size(), size;
     EntitySetHandle* ptr = &children_out[0];
     iGeom_getEntSets( mInstance, parent, num_hops, &ptr, &alloc, &size, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getPrnts( EntitySetHandle child, int num_hops,
                 std::vector<EntitySetHandle>& parents_out ) const
{
     int err, count;
     iGeom_getNumPrnt( mInstance, child, num_hops, &count, &err );
     if (iBase_SUCCESS != err)
          return (Error)err;
     parents_out.resize(count);
     int alloc = parents_out.size(), size;
     EntitySetHandle* ptr = &parents_out[0];
     iGeom_getEntSets( mInstance, child, num_hops, &ptr, &alloc, &size, &err );
     return (Error)err;
}


iGeom::Error
iGeom::subtract( EntitySetHandle set1, EntitySetHandle set2,
                 EntitySetHandle& result_set_out )
{
     int err;
     iGeom_subtract( mInstance, set1, set1, &result_set_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::intersect( EntitySetHandle set1, EntitySetHandle set2,
                  EntitySetHandle& result_set_out )
{
     int err;
     iGeom_intersect( mInstance, set1, set1, &result_set_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::unite( EntitySetHandle set1, EntitySetHandle set2,
              EntitySetHandle& result_set_out )
{
     int err;
     iGeom_unite( mInstance, set1, set1, &result_set_out, &err );
     return (Error)err;
}


iGeom::Error
iGeom::createTag( const char* tag_name,
                  int tag_num_type_values,
                  TagValueType tag_type,
                  TagHandle& tag_handle_out )
{
     int err;
     iGeom_createTag( mInstance, tag_name, tag_num_type_values, tag_type,
                      &tag_handle_out, &err, strlen(tag_name) );
     return (Error)err;
}


iGeom::Error
iGeom::destroyTag( TagHandle tag_handle, bool forced )
{
     int err;
     iGeom_destroyTag( mInstance, tag_handle, forced, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getTagName( TagHandle tag_handle, std::string& name_out ) const
{
     int err;
     char buffer[1024];
     memset( buffer, 0, sizeof(buffer) );
     iGeom_getTagName( mInstance, tag_handle, buffer, &err, sizeof(buffer) );
     name_out = buffer;
     return (Error)err;
}

iGeom::Error
iGeom::getTagSizeValues( TagHandle tag_handle, int& size_out ) const
{
     int err;
     iGeom_getTagSizeValues( mInstance, tag_handle, &size_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getTagSizeBytes( TagHandle tag_handle, int& size_out ) const
{
     int err;
     iGeom_getTagSizeBytes( mInstance, tag_handle, &size_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getTagHandle( const char* name, TagHandle& handle_out ) const
{
     int err;
     iGeom_getTagHandle( mInstance, name, &handle_out, &err, strlen(name) );
     return (Error)err;
}

iGeom::Error
iGeom::getTagType( TagHandle tag_handle, TagValueType& type_out ) const
{
     int err, result;
     iGeom_getTagType( mInstance, tag_handle, &result, &err );
     type_out = (TagValueType)result;
     return (Error)err;
}


iGeom::Error
iGeom::setEntSetData( EntitySetHandle set_handle,
                      TagHandle tag_handle,
                      const void* tag_value )
{
     int err, size = 1;
     iGeom_getTagSizeBytes( mInstance, tag_handle, &size, &err );
     iGeom_setEntSetData( mInstance, set_handle, tag_handle,
                          tag_value, size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::setEntSetIntData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         int value )
{
     int err;
     iGeom_setEntSetIntData( mInstance, set_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setEntSetDblData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         double value )
{
     int err;
     iGeom_setEntSetDblData( mInstance, set_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setEntSetEHData( EntitySetHandle set_handle,
                        TagHandle tag_handle,
                        EntityHandle value )

{
     int err;
     iGeom_setEntSetEHData( mInstance, set_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setEntSetESHData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         EntitySetHandle value )

{
     int err;
     iGeom_setEntSetESHData( mInstance, set_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getEntSetData( EntitySetHandle set_handle,
                      TagHandle tag_handle,
                      void* tag_value_out ) const
{
     int err, alloc = std::numeric_limits<int>::max(), size;
     iGeom_getEntSetData( mInstance, set_handle, tag_handle,
                          &tag_value_out, &alloc, &size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::getEntSetIntData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         int& value_out ) const
{
     int err;
     iGeom_getEntSetIntData( mInstance, set_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getEntSetDblData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         double& value_out ) const
{
     int err;
     iGeom_getEntSetDblData( mInstance, set_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getEntSetEHData( EntitySetHandle set_handle,
                        TagHandle tag_handle,
                        EntityHandle& value_out ) const

{
     int err;
     iGeom_getEntSetEHData( mInstance, set_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getEntSetESHData( EntitySetHandle set_handle,
                         TagHandle tag_handle,
                         EntitySetHandle& value_out ) const

{
     int err;
     iGeom_getEntSetESHData( mInstance, set_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getAllEntSetTags( EntitySetHandle set,
                         std::vector<TagHandle>& tags_out ) const
{
     if (tags_out.capacity() == 0)
          tags_out.resize( 32 );
     else
          tags_out.resize( tags_out.capacity() );

     int err, alloc = tags_out.size(), size = 0;
     TagHandle* ptr = &tags_out[0];
     iGeom_getAllEntSetTags( mInstance, set, &ptr, &alloc, &size, &err );
     tags_out.resize(size);

     if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
          alloc = tags_out.size();
          ptr = &tags_out[0];
          iGeom_getAllEntSetTags( mInstance, set, &ptr, &alloc, &size, &err );
     }

     return (Error)err;
}

iGeom::Error
iGeom::getAllTags( EntityHandle entity,
                   std::vector<TagHandle>& tags_out ) const

{
     if (tags_out.capacity() == 0)
          tags_out.resize( 32 );
     else
          tags_out.resize( tags_out.capacity() );

     int err, alloc = tags_out.size(), size = 0;
     TagHandle* ptr = &tags_out[0];
     iGeom_getAllTags( mInstance, entity, &ptr, &alloc, &size, &err );
     tags_out.resize(size);

     if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
          alloc = tags_out.size();
          ptr = &tags_out[0];
          iGeom_getAllTags( mInstance, entity, &ptr, &alloc, &size, &err );
     }

     return (Error)err;
}

iGeom::Error
iGeom::rmvEntSetTag( EntitySetHandle set, TagHandle tag )
{
     int err;
     iGeom_rmvEntSetTag( mInstance, set, tag, &err );
     return (Error)err;
}

iGeom::Error
iGeom::rmvTag( EntityHandle entity, TagHandle tag )
{
     int err;
     iGeom_rmvTag( mInstance, entity, tag, &err );
     return (Error)err;
}

iGeom::Error
iGeom::rmvArrTag( const EntityHandle* handles, int size, TagHandle tag )
{
     int err;
     iGeom_rmvArrTag( mInstance, handles, size, tag, &err );
     return (Error)err;
}


iGeom::Error
iGeom::getArrData( const EntityHandle* entity_handles,
                   int entity_handles_size,
                   TagHandle tag_handle,
                   void* tag_values_out ) const
{
     int err, alloc = std::numeric_limits<int>::max(), size;
     iGeom_getArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                       &tag_values_out, &alloc, &size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::getIntArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      int* tag_values_out ) const
{
     int err, alloc = std::numeric_limits<int>::max(), size;
     iGeom_getIntArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                          &tag_values_out, &alloc, &size, &err);
     return (Error)err;
}

iGeom::Error iGeom::getFacets(EntityHandle entity_handle, double dist_tolerance,
                              std::vector<double> &points, std::vector<int> &facets) const
{
  // try getting results using whatever space input vector has allocated
  bool already_allocd = false;
  int err, points_size = 0, points_alloc = points.capacity(),
    facets_size = 0, facets_alloc = facets.capacity();
  if (points_alloc && !facets_alloc) points_alloc = 0;
  else if (!points_alloc && facets_alloc) facets_alloc = 0;
  if (points_alloc && facets_alloc) already_allocd = true;
  double *points_ptr = (points_alloc ? &points[0] : NULL);
  int *facets_ptr = (facets_alloc ? &facets[0] : NULL);
  iGeom_getFacets(mInstance, entity_handle, dist_tolerance,
                  &points_ptr, &points_alloc, &points_size,
                  &facets_ptr, &facets_alloc, &facets_size, &err);

  if (iBase_SUCCESS == err && already_allocd) {
    points.resize(points_size);
    facets.resize(facets_size);
    return (Error)err;
  }
  else if (iBase_BAD_ARRAY_DIMENSION == err || iBase_BAD_ARRAY_SIZE == err) {
  // if input vector was too small, try again requesting allocation
    points_alloc = 0; points_ptr = NULL;
    facets_alloc = 0; facets_ptr = NULL;
    iGeom_getFacets(mInstance, entity_handle, dist_tolerance,
                    &points_ptr, &points_alloc, &points_size,
                    &facets_ptr, &facets_alloc, &facets_size, &err);
  }

  if (iBase_SUCCESS == err) {
    points.resize(points_size);
    std::copy(points_ptr, points_ptr+points_size, points.begin());
    facets.resize(facets_size);
    std::copy(facets_ptr, facets_ptr+facets_size, facets.begin());
    free(points_ptr); free(facets_ptr);
  }

  return (Error)err;
}

iGeom::Error
iGeom::getDblArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      double* tag_values_out ) const
{
     int err, alloc = std::numeric_limits<int>::max(), size;
     iGeom_getDblArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                          &tag_values_out, &alloc, &size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::getEHArrData( const EntityHandle* entity_handles,
                     int entity_handles_size,
                     TagHandle tag_handle,
                     EntityHandle* tag_values_out ) const

{
     int err, alloc = std::numeric_limits<int>::max(), size;
     iGeom_getEHArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                         &tag_values_out, &alloc, &size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::getESHArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      EntitySetHandle* tag_values_out ) const

{
     int err, alloc = std::numeric_limits<int>::max(), size;
     iGeom_getESHArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                          &tag_values_out, &alloc, &size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::setArrData( const EntityHandle* entity_handles,
                   int entity_handles_size,
                   TagHandle tag_handle,
                   const void* tag_values )
{
     int err, size = 1;
     iGeom_getTagSizeBytes( mInstance, tag_handle, &size, &err );
     iGeom_setArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                       tag_values, size*entity_handles_size,
                       &err );
     return (Error)err;
}

iGeom::Error
iGeom::setIntArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      const int* tag_values )
{
     int err, size = 1;
     iGeom_getTagSizeValues( mInstance, tag_handle, &size, &err );
     iGeom_setIntArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                          tag_values, size*entity_handles_size, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setDblArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      const double* tag_values )
{
     int err, size = 1;
     iGeom_getTagSizeValues( mInstance, tag_handle, &size, &err );
     iGeom_setDblArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                          tag_values, size*entity_handles_size, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setEHArrData( const EntityHandle* entity_handles,
                     int entity_handles_size,
                     TagHandle tag_handle,
                     const EntityHandle* tag_values )
{
     int err, size = 1;
     iGeom_getTagSizeValues( mInstance, tag_handle, &size, &err );
     iGeom_setEHArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                         tag_values, size*entity_handles_size, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setESHArrData( const EntityHandle* entity_handles,
                      int entity_handles_size,
                      TagHandle tag_handle,
                      const EntitySetHandle* tag_values )
{
     int err, size = 1;
     iGeom_getTagSizeValues( mInstance, tag_handle, &size, &err );
     iGeom_setESHArrData( mInstance, entity_handles, entity_handles_size, tag_handle,
                          tag_values, size*entity_handles_size, &err );
     return (Error)err;
}



iGeom::Error
iGeom::setData( EntityHandle entity_handle,
                TagHandle tag_handle,
                const void* tag_value )
{
     int err, size = 1;
     iGeom_getTagSizeBytes( mInstance, tag_handle, &size, &err );
     iGeom_setData( mInstance, entity_handle, tag_handle,
                    tag_value, size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::setIntData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   int value )
{
     int err;
     iGeom_setIntData( mInstance, entity_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setDblData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   double value )
{
     int err;
     iGeom_setDblData( mInstance, entity_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setEHData( EntityHandle entity_handle,
                  TagHandle tag_handle,
                  EntityHandle value )

{
     int err;
     iGeom_setEHData( mInstance, entity_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::setESHData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   EntitySetHandle value )

{
     int err;
     iGeom_setESHData( mInstance, entity_handle, tag_handle, value, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getData( EntityHandle entity_handle,
                TagHandle tag_handle,
                void* tag_value_out ) const
{
     int err, alloc = std::numeric_limits<int>::max(), size;
     iGeom_getData( mInstance, entity_handle, tag_handle,
                    &tag_value_out, &alloc, &size, &err);
     return (Error)err;
}

iGeom::Error
iGeom::getIntData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   int& value_out ) const
{
     int err;
     iGeom_getIntData( mInstance, entity_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getDblData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   double& value_out ) const
{
     int err;
     iGeom_getDblData( mInstance, entity_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getEHData( EntityHandle entity_handle,
                  TagHandle tag_handle,
                  EntityHandle& value_out ) const
{
     int err;
     iGeom_getEHData( mInstance, entity_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

iGeom::Error
iGeom::getESHData( EntityHandle entity_handle,
                   TagHandle tag_handle,
                   EntitySetHandle& value_out ) const
{
     int err;
     iGeom_getESHData( mInstance, entity_handle, tag_handle, &value_out, &err );
     return (Error)err;
}

#endif
