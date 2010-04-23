/*
 * QslimDecimation.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: iulian
 */

#include "QslimDecimation.h"
#include "MBRange.hpp"

// for proximity searches
#include "MBAdaptiveKDTree.hpp"

#include "Mat4.h"
#include "defs.h"
#include "quadrics.h"
#include <time.h>

// those are used in model navigation/simplification
#include "primitives.h"

// this is the global thing, that everybody will use
MBInterface * mb;
MBTag uniqIDtag; // this will be used to mark vertices MBEntityHandles
MBTag validTag;

MBRange verts; // original list of vertices, that are part of the original triangles
MBRange triangles;
MBRange edgs;

QslimOptions opts; // external

int uniqID(MBEntityHandle v) {
	int val;
	MBErrorCode rval = mb->tag_get_data(uniqIDtag, &v, 1, &val);
	assert(rval==MB_SUCCESS);
	return val;
}
// the vertices are not deleted anymore, just invalidated
// the edges are deleted, though, and triangles
int ehIsValid(MBEntityHandle v) {
	unsigned char val;
	MBErrorCode rval = mb->tag_get_data(validTag, &v, 1, &val);
	assert(rval==MB_SUCCESS);
	return (int) val;
}

// include here the main classes used for decimation

#include "Heap.h"
// prox grid is used for proximity grid only
//#include "ProxGrid.h"

class pair_info: public Heapable {
public:
	MBEntityHandle v0, v1; // Vertex *v0, *v1;

	Vec3 candidate;
	double cost;

	//pair_info ( Vertex *a,Vertex *b ) { v0=a; v1=b; cost=HUGE; }
	pair_info(MBEntityHandle a, MBEntityHandle b) {
		v0 = a;
		v1 = b;
		cost = HUGE;
	}

	bool isValid() {
		return ehIsValid(v0) && ehIsValid(v1);
	}//  :v0->isValid() && v1->isValid(); }
};

typedef buffer<pair_info *> pair_buffer;

class vert_info {
public:

	pair_buffer pairs;

	Mat4 Q;
	double norm;

	vert_info() :
		Q(Mat4::zero) {
		pairs.init(2);
		norm = 0.0;
	}
};

// these are
static Heap *heap;
static array<vert_info> vinfo;
static double proximity_limit; // distance threshold squared

int validFaceCount;
int validVertCount;

////////////////////////////////////////////////////////////////////////
//
// Low-level routines for manipulating pairs
//

static inline vert_info& vertex_info(MBEntityHandle v)//Vertex *v )
{
	//  MBEntityHandle should return an integer tag with an
	//  index in the big array of vert_info
	//   something like: return tag
	//   for the time being, we can return the simple id...
	return vinfo(uniqID(v));
}

static
bool check_for_pair(MBEntityHandle v0, MBEntityHandle v1)//Vertex *v0, Vertex *v1 )
{
	const pair_buffer& pairs = vertex_info(v0).pairs;

	for (int i = 0; i < pairs.length(); i++) {
		if (pairs(i)->v0 == v1 || pairs(i)->v1 == v1)
			return true;
	}

	return false;
}

static pair_info *new_pair(MBEntityHandle v0, MBEntityHandle v1)//  Vertex *v0, Vertex *v1 )
{
	vert_info& v0_info = vertex_info(v0);
	vert_info& v1_info = vertex_info(v1);

	pair_info *pair = new pair_info(v0, v1);
	v0_info.pairs.add(pair);
	v1_info.pairs.add(pair);

	return pair;
}

static
void delete_pair(pair_info *pair) {
	vert_info& v0_info = vertex_info(pair->v0);
	vert_info& v1_info = vertex_info(pair->v1);

	v0_info.pairs.remove(v0_info.pairs.find(pair));
	v1_info.pairs.remove(v1_info.pairs.find(pair));

	if (pair->isInHeap())
		heap->kill(pair->getHeapPos());

	delete pair;
}

////////////////////////////////////////////////////////////////////////
//
// The actual guts of the algorithm:
//
//     - pair_is_valid
//     - compute_pair_info
//     - do_contract
//

static
bool pair_is_valid(MBEntityHandle u, MBEntityHandle v)// Vertex *u, Vertex *v )
{
	//
	Vec3 vu = getVec3FromMBVertex(mb, u);
	Vec3 vv = getVec3FromMBVertex(mb, v);
	return norm2(vu - vv) < proximity_limit;
	//return  norm2 ( *u - *v ) < proximity_limit;
}

static
int predict_face(MBEntityHandle tria, MBEntityHandle v1, MBEntityHandle v2, /*Face& F, Vertex *v1, Vertex *v2,*/
Vec3& vnew, Vec3& f1, Vec3& f2, Vec3& f3) {
	int nmapped = 0;
	const MBEntityHandle * conn;
	int num_nodes;
	MBErrorCode rval = mb->get_connectivity(tria, conn, num_nodes);
	assert(3==num_nodes && rval == MB_SUCCESS);
	if (conn[0] == v1 || conn[0] == v2) {
		f1 = vnew;
		nmapped++;
	} else
		f1 = getVec3FromMBVertex(mb, conn[0]);

	if (conn[1] == v1 || conn[1] == v2) {
		f2 = vnew;
		nmapped++;
	} else
		f2 = getVec3FromMBVertex(mb, conn[1]);

	if (conn[2] == v1 || conn[2] == v2) {
		f3 = vnew;
		nmapped++;
	} else
		f3 = getVec3FromMBVertex(mb, conn[2]);

	// find vertices in face tria
	/*
	 if ( F.vertex ( 0 ) == v1 || F.vertex ( 0 ) == v2 )
	 { f1 = vnew;  nmapped++; }
	 else f1 = *F.vertex ( 0 );

	 if ( F.vertex ( 1 ) == v1 || F.vertex ( 1 ) == v2 )
	 { f2 = vnew;  nmapped++; }
	 else f2 = *F.vertex ( 1 );

	 if ( F.vertex ( 2 ) == v1 || F.vertex ( 2 ) == v2 )
	 { f3 = vnew;  nmapped++; }
	 else f3 = *F.vertex ( 2 );
	 */
	return nmapped;
}

#define MESH_INVERSION_PENALTY 1e9

static
double pair_mesh_positivity(/* Model& M,*/MBEntityHandle v1, MBEntityHandle v2, /*Vertex *v1, Vertex *v2,*/
		Vec3& vnew) {
	std::vector<MBEntityHandle> changed;

	// :  construct the list of faces influenced by the
	//   moving of vertices v1 and v2 into vnew
	//M.contractionRegion ( v1, v2, changed );
	MBErrorCode rval = contractionRegion(mb, v1, v2, changed);
	if (rval != MB_SUCCESS) {
		std::cout << "error in getting adjacency information vs: "
				<< mb->id_from_handle(v1) << " " << mb->id_from_handle(v2)
				<< "\n";
	}

	// double Nsum = 0;
	if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
		*opts.logfile << " positivity for v1, v2: " << mb->id_from_handle(v1)
				<< " " << mb->id_from_handle(v2) << std::endl;

	for (int i = 0; i < changed.size(); i++) {
		//Face& F = *changed ( i );
		MBEntityHandle F = changed[i];
		Vec3 f1, f2, f3;

		int nmapped = predict_face(F, v1, v2, vnew, f1, f2, f3);

		//
		// Only consider non-degenerate faces
		if (nmapped < 2) {
			//Plane Pnew ( f1, f2, f3 );
#if 0
			Vec3 normalNew = Pnew.normal();
			if ( normalNew[Z] < positivityMin )
			positivityMin=normalNew[Z]; // Z direction!!!
			if (opts.logfile && opts.selected_output&OUTPUT_CONTRACTIONS )
			*opts.logfile << "Triangle " << mb->id_from_handle(F)
			<< " nmapped " << nmapped << std::endl;
			if (opts.logfile && positivityMin<=0 && opts.selected_output&OUTPUT_CONTRACTIONS )
			*opts.logfile << "Triangle " << mb->id_from_handle(F)
			<< " normal Z:" << normalNew[Z] << std::endl;
#endif
			double positiv = (f2[0] - f1[0]) * (f3[1] - f1[1])
					- (f2[1] - f1[1]) * (f3[0] - f1[0]);
			if (positiv <= 0) {
				if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
					*opts.logfile << "Triangle " << mb->id_from_handle(F)
							<< " nmapped " << nmapped << " orient: " << positiv
							<< std::endl;
				return MESH_INVERSION_PENALTY * 10;
			}
		}
	}

	//return (-Nmin) * MESH_INVERSION_PENALTY;

	return 0.0;
}

static
double pair_mesh_penalty( /*Model& M, Vertex *v1, Vertex *v2,*/
MBEntityHandle v1, MBEntityHandle v2, Vec3& vnew) {
	std::vector<MBEntityHandle> changed;

	//   construct the list of faces influenced by the
	//   moving of vertices v1 and v2 into vnew
	//M.contractionRegion ( v1, v2, changed );
	MBErrorCode rval = contractionRegion(mb, v1, v2, changed);
	if (rval != MB_SUCCESS) {
		std::cout << "error in getting adjacency information vs: "
				<< mb->id_from_handle(v1) << " " << mb->id_from_handle(v2)
				<< "\n";
	}

	// double Nsum = 0;
	double Nmin = 0;

	for (int i = 0; i < changed.size(); i++) {
		//Face& F = *changed ( i );
		MBEntityHandle F = changed[i];
		Vec3 f1, f2, f3;

		int nmapped = predict_face(F, v1, v2, vnew, f1, f2, f3);
		//
		// Only consider non-degenerate faces
		if (nmapped < 2) {
			Plane Pnew(f1, f2, f3);
			Plane p = trianglePlane(mb, F);
			double delta = Pnew.normal() * p.normal(); //  Pnew.normal() * F.plane().normal();

			if (Nmin > delta)
				Nmin = delta;
		}
	}

	//return (-Nmin) * MESH_INVERSION_PENALTY;
	if (Nmin < 0.0)
		return MESH_INVERSION_PENALTY;
	else
		return 0.0;
}

static
void compute_pair_info(pair_info *pair /* Model * pM0,*/) {
	MBEntityHandle v0 = pair->v0;
	MBEntityHandle v1 = pair->v1;

	// Vertex *v0 = pair->v0;
	// Vertex *v1 = pair->v1;

	vert_info& v0_info = vertex_info(v0);
	vert_info& v1_info = vertex_info(v1);

	Mat4 Q = v0_info.Q + v1_info.Q;
	double norm = v0_info.norm + v1_info.norm;

	pair->cost = quadrix_pair_target(Q, v0, v1, pair->candidate);

	if (opts.will_weight_by_area)
		pair->cost /= norm;

	if (opts.will_preserve_mesh_quality)
		pair->cost += pair_mesh_penalty(/* *pM0,*/v0, v1, pair->candidate);

	if (opts.height_fields)
		pair->cost += pair_mesh_positivity(/* *pM0, */v0, v1, pair->candidate);

	if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)

		*opts.logfile << "pair ( v" << mb->id_from_handle(v0) << " v"
				<< mb->id_from_handle(v1) << " ) cost: " << -pair->cost
				<< std::endl;
	//
	// NOTICE:  In the heap we use the negative cost.  That's because
	//          the heap is implemented as a MAX heap.
	//
	if (pair->isInHeap()) {
		heap->update(pair, -pair->cost);
	} else {
		heap->insert(pair, -pair->cost);
	}
	if (opts.logfile && opts.selected_output & OUTPUT_COST) {
		heap_node *top = heap->top();
		pair_info *pairTop = (pair_info *) top->obj;
		*opts.logfile << " i/u pair (" << uniqID(pair->v0) << "," << uniqID(
				pair->v1) << ")=" << pair->cost << "  min : (" << uniqID(
				pairTop->v0) << "," << uniqID(pairTop->v1) << ") "
				<< pairTop->cost << std::endl;
	}
}
void recomputeChangedPairsCost(std::vector<MBEntityHandle> & changed,/* Model *pM0, Vertex *v0,*/
MBEntityHandle v0) {
	//
	for (int i = 0; i < changed.size(); i++) {

		MBEntityHandle F = changed[i];
		const MBEntityHandle * conn;
		int num_nodes;
		MBErrorCode rval = mb->get_connectivity(F, conn, num_nodes);
		//if (!F->isValid())
		//   continue;
		// recompute the pair that is not connected to vertex v0
		// loop over all the vertices of F that are not v0, and recompute the costs
		// of all the pairs associated  that do not contain v0
		// we do not have to recreate or delete any pair, we just recompute what we have
		// some will be recomputed 2 times, but it is OK
		for (int k = 0; k < 3; k++) {
			MBEntityHandle v = conn[k];
			if (v == v0)
				continue;
			vert_info & v_info = vertex_info(v);
			for (int j = 0; j < v_info.pairs.length(); j++) {
				pair_info *p = v_info.pairs(j);
				if (p->v0 == v0 || p->v1 == v0)
					continue; // do not recompute cost of pairs already computed
				if (opts.logfile && (opts.selected_output & OUTPUT_COST))
					*opts.logfile << "recompute cost of pair (v" << uniqID(
							p->v0) + 1 << " v" << uniqID(p->v1) + 1 << ")"
							<< std::endl;
				compute_pair_info(p);
			}

		}

	}
}

static
void do_contract(pair_info *pair) {

	MBEntityHandle v0 = pair->v0;
	MBEntityHandle v1 = pair->v1;
	vert_info& v0_info = vertex_info(v0);
	vert_info& v1_info = vertex_info(v1);
	Vec3 vnew = pair->candidate;
	if (opts.logfile && (opts.selected_output & OUTPUT_COST)) {
		*opts.logfile << "---- do contract v0:" << uniqID(v0) + 1 << " v1:"
				<< uniqID(v1) + 1 << std::endl;
	}
	int i;

	//
	// Make v0 be the new vertex
	v0_info.Q += v1_info.Q;
	v0_info.norm += v1_info.norm;

	//
	// Perform the actual contraction
	std::vector<MBEntityHandle> changed;
	MBErrorCode rval1 = contract(mb, v0, v1, vnew, changed);
	assert (MB_SUCCESS == rval1);

#ifdef SUPPORT_VCOLOR
	//
	// If the vertices are colored, color the new vertex
	// using the average of the old colors.
	v0->props->color += v1->props->color;
	v0->props->color /= 2;
#endif

	//
	// Remove the pair that we just contracted
	delete_pair(pair);
	//
	// Recalculate pairs associated with v0
	for (i = 0; i < v0_info.pairs.length(); i++) {
		pair_info *p = v0_info.pairs(i);
		compute_pair_info(p);
	}

	//
	// Process pairs associated with now dead vertex

	static pair_buffer condemned(6); // collect condemned pairs for execution
	condemned.reset();

	for (i = 0; i < v1_info.pairs.length(); i++) {
		pair_info *p = v1_info.pairs(i);

		MBEntityHandle u;
		if (p->v0 == v1)
			u = p->v1;
		else if (p->v1 == v1)
			u = p->v0;
		else
			std::cerr << "YOW!  This is a bogus pair." << std::endl;

		if (!check_for_pair(u, v0)) {
			p->v0 = v0;
			p->v1 = u;
			v0_info.pairs.add(p);
			compute_pair_info(p);
		} else
			condemned.add(p);
	}

	for (i = 0; i < condemned.length(); i++)
		// Do you have any last requests?
		delete_pair(condemned(i));
	// only now you can delete the vertex v1 from database
	// MBErrorCode rval = mb->delete_entities(&v1, 1);
	// no, it is better to invalidate the vertex, do not delete it
	// maybe we will delete at the end all that are invalid ??
	int invalid = 0;
	MBErrorCode rval = mb->tag_set_data(validTag, &v1, 1, &invalid);

	v1_info.pairs.reset(); // safety precaution
	recomputeChangedPairsCost(changed, v0);

}

////////////////////////////////////////////////////////////////////////
//
// External interface: setup and single step iteration
//

bool decimate_quadric(MBEntityHandle v, Mat4& Q) {
	if (vinfo.length() > 0) {
		vert_info & vinf = vertex_info(v);
		Q = vinf.Q;
		return true;
	} else
		return false;
}

// it is assumed it is mb, moab interface
void decimate_contract() {
	heap_node *top;
	pair_info *pair;

	for (;;) {
		top = heap->extract();
		if (!top)
			return;
		pair = (pair_info *) top->obj;

		//
		// This may or may not be necessary.  I'm just not quite
		// willing to assume that all the junk has been removed from the
		// heap.
		if (pair->isValid())
			break;

		delete_pair(pair);
	}

	if (opts.logfile && (opts.selected_output & OUTPUT_COST))
		*opts.logfile << "#$cost " << validFaceCount << " before contract: "
				<< pair->cost << std::endl;

	do_contract(pair);

	if (opts.logfile && (opts.selected_output & OUTPUT_COST))
		*opts.logfile << "#$cost " << validFaceCount << std::endl;

	validVertCount--; // Attempt to maintain valid vertex information
}

double decimate_error(MBEntityHandle v) {
	vert_info& info = vertex_info(v);

	double err = quadrix_evaluate_vertex(v, info.Q);

	if (opts.will_weight_by_area)
		err /= info.norm;

	return err;
}

double decimate_min_error() {
	heap_node *top;
	pair_info *pair;

	for (;;) {
		top = heap->top();
		if (!top)
			return -1.0;
		pair = (pair_info *) top->obj;

		if (pair->isValid())
			break;

		top = heap->extract();
		delete_pair(pair);
	}

	return pair->cost;
}
#if 0
double decimate_max_error ( )
{
	real max_err = 0;

	for ( int i=0; i<m.vertCount(); i++ )
	if ( m.vertex ( i )->isValid() )
	{
		max_err = MAX ( max_err, decimate_error ( m.vertex ( i ) ) );
	}

	return max_err;
}
#endif

QslimDecimation::QslimDecimation(iMesh_Instance mesh,
		iBase_EntitySetHandle root_set) {
	m_mesh = mesh;
	m_InitialSet = root_set;// it is not necessarily the root set
}

QslimDecimation::~QslimDecimation() {
}
int QslimDecimation::decimate(QslimOptions & iOpts) {
	// opts is external
	opts = iOpts;
	mb = reinterpret_cast<MBInterface *> (m_mesh);
	// need to get all the triangles from the set
	// also all the edges, and all vertices
	//
	if (NULL == mb)
		return 1;// error
	//MBEntityHandle mbSet = reinterpret_cast<MBEntityHandle>(m_InitialSet);

	clock_t start_time = clock();
	int err = Init();
	if (err) {
		std::cerr << "Error in initialization of decimation. Abort\n";
		return 1;
	}
	clock_t init_time = clock();
	std::cout << "   Initialization  " << (double) (init_time - start_time)
			/ CLOCKS_PER_SEC << " s.\n";

	int faces_reduction = validFaceCount - opts.face_target;
	int counter = 0, interval = 0;
	clock_t currentTime = init_time;
	while (validFaceCount > opts.face_target && decimate_min_error()
			< opts.error_tolerance) {
		int initf = validFaceCount;
		// big routine
		decimate_contract();
		counter += (initf - validFaceCount);
		if (counter > faces_reduction / opts.timingIntervals) {
			// print some time stats
			clock_t p10_time = clock();
			std::cout << "     " << ++interval << "/" << opts.timingIntervals
					<< " reduce to " << validFaceCount << " faces in "
					<< (double) (p10_time - currentTime) / CLOCKS_PER_SEC
					<< " s, total:" << (double) (p10_time - init_time)
					/ CLOCKS_PER_SEC << " s.\n";
			counter = 0;
			currentTime = p10_time;
		}
	}

	clock_t finish_time = clock();
	std::cout << "   Decimation: " << (double) (finish_time - init_time)
			/ CLOCKS_PER_SEC << " s.\n";

	MBRange::const_reverse_iterator rit;
	if (opts.useDelayedDeletion) {

		// delete triangles and edges that are invalid
		for (rit = triangles.rbegin(); rit != triangles.rend(); rit++) {
			MBEntityHandle v = *rit;
			// check the validity
			if (ehIsValid(v))
				continue;
			MBErrorCode rval = mb->delete_entities(&v, 1);
		}
		// maybe we should delete all edges, but for now, we will keep them
		for (rit = edgs.rbegin(); rit != edgs.rend(); rit++) {
			MBEntityHandle v = *rit;
			// check the validity
			if (ehIsValid(v))
				continue;
			MBErrorCode rval = mb->delete_entities(&v, 1);
		}

	}
	// delete them one by one
	for (rit = verts.rbegin(); rit != verts.rend(); rit++) {
		MBEntityHandle v = *rit;
		// check the validity
		if (ehIsValid(v))
			continue;
		MBErrorCode rval = mb->delete_entities(&v, 1);
	}
	clock_t delete_vTime = clock();
	std::cout << "   Delete Vertices: "
			<< (double) (delete_vTime - finish_time) / CLOCKS_PER_SEC
			<< " s.\n";
	// we need to delete the tags we created; they are artificial
	//

	return 0;
}

int QslimDecimation::Init() {
	int i, j;

	MBEntityHandle * set = reinterpret_cast<MBEntityHandle *> (&m_InitialSet);
	MBErrorCode rval = mb->get_entities_by_type(*set, MBTRI, triangles);
	validFaceCount = triangles.size();// this gets reduced every time we simplify the model

	// create all the edges if not existing
	mb->get_adjacencies(triangles, 1, true, edgs, MBInterface::UNION);

	// MBRange verts;// the vertices are always there, they do not need to be created
	mb->get_adjacencies(triangles, 0, true, verts, MBInterface::UNION);
	int numNodes = verts.size();
	validVertCount = numNodes; // this will be kept
	vinfo.init(numNodes);
	// set a unique integer tag with the position in vinfo array
	//  this will be used instead of v->uniqID in the vinfo array
	int def_data = -1;

	rval = mb->tag_create("uniqID", sizeof(int), MB_TAG_DENSE, uniqIDtag,
			&def_data);
	if (MB_SUCCESS != rval)
		return 1;

	unsigned char def_data_bit = 1;// valid by default
	rval = mb->tag_create("valid", 1, MB_TAG_BIT, validTag, &def_data_bit);
	if (MB_SUCCESS != rval)
		return 1;

	// set tag for each vertex; this will not be changed during simplification
	i = 0; // for index
	for (MBRange::iterator it = verts.begin(); it != verts.end(); it++, i++) {
		// this index i will be the index in the vinfo array
		MBEntityHandle v = *it;
		rval = mb->tag_set_data(uniqIDtag, &v, 1, &i);
		// the default valid data is 1; set it to 0 only to "mark" the vertex invalid for future deletion
		// is it really necessary if we put the default value as 1 anyway

		rval = mb->tag_set_data(validTag, &v, 1, &def_data_bit);
	}
	MBRange::iterator it;
	if (opts.logfile && opts.selected_output & OUTPUT_MODEL_DEFN) {
		for (it = verts.begin(); it != verts.end(); it++) {
			MBEntityHandle v = *it;
			double coords[3];
			rval = mb->get_coords(&v, 1, coords);
			*opts.logfile << "v: " << uniqID(v) << " " << mb->id_from_handle(v)
					<< " " << coords[0] << " " << coords[1] << " " << coords[2]
					<< std::endl;
		}
	}
	std::cout << "  Decimate:  Distributing shape constraints." << std::endl;

	if (opts.will_use_vertex_constraint)
		for (it = verts.begin(); it != verts.end(); it++) {
			MBEntityHandle v = *it;
			vertex_info(v).Q = quadrix_vertex_constraint(v);
		}

	if (opts.will_use_plane_constraint) {
		for (it = triangles.begin(); it != triangles.end(); it++) {

			MBEntityHandle tr = *it;
			Mat4 Q = quadrix_plane_constraint(tr);
			double norm = 0.0;

			if (opts.will_weight_by_area) {
				norm = 1; // triangle area : m_model->face ( i )->area();
				Q *= norm;
			}
			const MBEntityHandle * conn;
			int num_nodes;
			rval = mb->get_connectivity(tr, conn, num_nodes);
			for (j = 0; j < 3; j++) {
				vert_info& vj_info = vertex_info(conn[j]);
				vj_info.Q += Q;
				vertex_info(conn[j]).norm += norm;

			}
		}
	}

	if (opts.will_constrain_boundaries) {
		std::cout << "  Decimate:  Accumulating discontinuity constraints."
				<< std::endl;
		for (it = edgs.begin(); it != edgs.end(); it++) {
			MBEntityHandle edg = *it;
			if (is_border(edg)) {
				const MBEntityHandle * conn;
				int num_nodes;
				rval = mb->get_connectivity(edg, conn, num_nodes);
				if (MB_SUCCESS != rval)
					return 1;// fail
				Mat4 B = quadrix_discontinuity_constraint(edg);
				double norm = 0.0;

				if (opts.will_weight_by_area) {
					Vec3 ve1 = getVec3FromMBVertex(mb, conn[0]);
					Vec3 ve2 = getVec3FromMBVertex(mb, conn[1]);
					norm = norm2(ve1 - ve2);
					B *= norm;
				}

				B *= opts.boundary_constraint_weight;
				vert_info& v0_info = vertex_info(conn[0]);
				vert_info& v1_info = vertex_info(conn[1]);
				v0_info.Q += B;
				v0_info.norm += norm;
				v1_info.Q += B;
				v1_info.norm += norm;
			}
		}
	}

	std::cout << "  Decimate:  Allocating heap." << std::endl;
	heap = new Heap(edgs.size());

	int pair_count = 0;

	std::cout << "  Decimate:  Collecting pairs [edges]." << std::endl;
	for (it = edgs.begin(); it != edgs.end(); it++) {
		MBEntityHandle edg = *it;
		const MBEntityHandle * conn;
		int num_nodes;
		rval = mb->get_connectivity(edg, conn, num_nodes);
		if (MB_SUCCESS != rval)
			return 1;// fail
		pair_info *pair = new_pair(conn[0], conn[1]);
		compute_pair_info(pair);
		pair_count++;
	}

	if (opts.pair_selection_tolerance < 0) {
		opts.pair_selection_tolerance = 1;//m_model->bounds.radius * 0.05;
		std::cout << "  Decimate:  Auto-limiting at 5% of model radius."
				<< std::endl;
	}
	proximity_limit = opts.pair_selection_tolerance
			* opts.pair_selection_tolerance;
	if (proximity_limit > 0) {
		std::cout << "  Decimate:  Collecting pairs [limit="
				<< opts.pair_selection_tolerance << "]." << std::endl;
		// use adaptive kd tree to find proximity vertices
		MBEntityHandle tree_root = 0;
		MBAdaptiveKDTree kd(mb, true);
		rval = kd.build_tree(verts, tree_root);
		if (rval != MB_SUCCESS) {
			std::cout << "Can't build tree for vertices" << std::endl;
			return 1;
		}

		for (it = verts.begin(); it != verts.end(); it++) {
			MBRange closeVertices;
			closeVertices.clear();
			MBEntityHandle v = *it;
			double coords[3];
			mb->get_coords(&v, 1, coords);
			//MBCartVect v1(coords);
			std::vector<MBEntityHandle> leaves; // sets of vertices close by
			kd.leaves_within_distance(tree_root, coords,
					opts.pair_selection_tolerance, leaves);
			// add to the list of close vertices
			for (j = 0; j < leaves.size(); j++) {
				rval = mb->get_entities_by_type(leaves[j], MBVERTEX,
						closeVertices);// add to the list
			}

			for (MBRange::iterator it2 = closeVertices.begin(); it2
					!= closeVertices.end(); it2++) {
				MBEntityHandle vclose = *it2;
				if (v == vclose)
					continue;
				double coords2[3];
				mb->get_coords(&vclose, 1, coords2);

				//MBCartVect v2(coords2);
				double dd = (coords[0] - coords2[0]) * (coords[0] - coords2[0])
						+ (coords[1] - coords2[1]) * (coords[1] - coords2[1])
						+ (coords[2] - coords2[2]) * (coords[2] - coords2[2]);

				if (dd > proximity_limit)
					continue;

#ifdef SAFETY
				assert ( pair_is_valid ( v1,v2 ) );
#endif
				if (!check_for_pair(v, vclose)) {
					pair_info *pair = new_pair(v, vclose);
					compute_pair_info(pair);
					pair_count++;
				}
			}

		}
	} else
		std::cout << "  Decimate:  Ignoring non-edge pairs [limit=0]."
				<< std::endl;

	std::cout << "  Decimate:  Designated " << pair_count << " pairs."
			<< std::endl;

	return 0;// no error
}
