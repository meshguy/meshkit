/**
 this file will prepare the paver call; basically, create the
 smooth surface, and pass it to the Pave routine call
 SmoothFaceEval will be derived from CAMAL eval, and will call the smooth faceter that I will
 extract from CGM. The arithmetic is explained in Face-Based Surfaces
   for 3D Mesh generation
 // input
 */
//
// mesh something with paver from CAMAL
//
#include "CamalPaveDriver.hpp"

#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"

#ifdef USE_CGM
// copy from facets.cpp
#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "GeometryQueryTool.hpp"
#include "FacetModifyEngine.hpp"
#include "CubitObserver.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Body.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "Loop.hpp"

#endif

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
// end copy includes from facets.cpp
//#include "camal_interface.hpp"

// camal edge mesher: need to mesh with even numbers
//#include "CMLEdgeMesher.hpp"
#include "SmoothCurveEval.hpp"
#include "SmoothFaceEval.hpp"
#include "SmoothVertex.hpp"

// these are needed by CAMAL evaluator; 
#include "CAMALGeomEval.hpp"
#include "CAMALSizeEval.hpp"

// the paver is here: 
#include "CMLPaver.hpp"

//#include "iMesh.h"  // already included

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <set>
//#include <map>
#include <assert.h>
#include <string.h>

const bool debugLocal = false;
extern bool debug_surf_eval;

CamalPaveDriver::CamalPaveDriver( MBInterface * mb, MBEntityHandle set,
        MBInterface * out, double angle) {
	_mb = mb;
	_mbo = out;
	_set = set;
	_angle = angle;
	_smthFace= NULL;
	_smthCurve=NULL;
	_smthVertex=NULL;
}

CamalPaveDriver::~CamalPaveDriver()
{
	delete [] _smthFace;
	delete [] _smthCurve;
	delete [] _smthVertex;
	_smthFace= NULL;
	_smthCurve=NULL;
	_smthVertex=NULL;
}
bool CamalPaveDriver::remesh(double mesh_size, int mesh_intervals,
		const bool force_intervals) {

	_mesh_size = mesh_size;
	_mesh_intervals = mesh_intervals;
	_force_intervals = force_intervals;

	clock_t start_time = clock();
	initializeSmoothing();
	clock_t init_time = clock();

#ifdef USE_CGM
	// now, redo the part about cgm, using moab, not iMesh
	prepareCGMEvaluator();// this might be commented out, actually
#endif
	mesh_vertices();

	MBErrorCode rval = find_face_loops();
	// here, we need to establish the mesh count on each curve,
	// the constraints are related to even count on faces
	establish_mesh_curve_count();

	mesh_curves();
	clock_t meshCurves_time = clock();

	mesh_surfaces();

	clock_t mesh_surfaces_time = clock();

	std::cout << "    Initialization time: "
	       << (double) (init_time - start_time)/CLOCKS_PER_SEC
	      << " s\n    Mesh Curves time: "
	       << (double) (meshCurves_time - init_time)/CLOCKS_PER_SEC
	       << "s\n    Mesh Surfaces time: "
	       << (double) (mesh_surfaces_time - meshCurves_time)/CLOCKS_PER_SEC  << " s\n ";

	return true;
}

bool CamalPaveDriver::initializeSmoothing()
{
	//
	// first of all, we need to retrieve all the surfaces from the (root) set
	// in icesheet_test we use iGeom, but maybe that is a stretch
	// get directly the sets with geom dim 2, and from there create the SmoothFaceEval
	MBTag geom_tag, gid_tag;
	MBErrorCode rval = _mb->tag_get_handle(GEOM_DIMENSION_TAG_NAME, geom_tag);
	rval = _mb->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);

		// traverse the model, dimension 2, 1, 0
	MBRange csets, fsets, vsets;
	//std::vector<moab::EntityHandle> sense_ents;
	//std::vector<int> senses, pgids;
	int dim=2;
	void *dim_ptr = &dim;
	//bool sense;

	moab::GeomTopoTool * my_geomTopoTool = new moab::GeomTopoTool(_mb);

	rval = my_geomTopoTool->construct_obb_trees();
	assert(MB_SUCCESS==rval);


	fsets.clear();
	rval = _mb->get_entities_by_type_and_tag(0, MBENTITYSET,
												&geom_tag, &dim_ptr, 1,
												fsets, 1, false);
	int numSurfaces = fsets.size();
	//SmoothFaceEval ** smthFace = new SmoothFaceEval *[numSurfaces];
	_smthFace = new SmoothFaceEval *[numSurfaces];
	// there should also be a map from surfaces to evaluators
	//std::map<MBEntityHandle, SmoothFaceEval*> mapSurfaces;

	int i=0;
	MBRange::iterator it;
	for ( it = fsets.begin(); it!=fsets.end(); it++, i++)
	{
		MBEntityHandle face= *it;
		_smthFace[i]= new SmoothFaceEval(_mb, face, _mbo, my_geomTopoTool);// geom topo tool will be used for searching,
		// among other things; also for senses in edge sets...
		_mapSurfaces[face] = _smthFace[i];
	}

	csets.clear();
	dim = 1; // get now the curves
	rval = _mb->get_entities_by_type_and_tag(0, MBENTITYSET,
												&geom_tag, &dim_ptr, 1,
												csets, 1, false);
	int numCurves = csets.size();
	//SmoothCurveEval ** smthCurve = new SmoothCurveEval *[numCurves];
	_smthCurve = new SmoothCurveEval *[numCurves];
	// there should also be a map from surfaces to evaluators
	//std::map<MBEntityHandle, SmoothCurveEval*> mapCurves;

	i=0;
	for ( it = csets.begin(); it!=csets.end(); it++, i++)
	{
		MBEntityHandle curve= *it;
		_smthCurve[i]= new SmoothCurveEval(_mb, curve, _mbo);
		_mapCurves[curve] = _smthCurve[i];
	}

	// create another mapping for vertices (sets of dimension 0)
	vsets.clear();
	dim = 0; // get now the vertice sets (dimension 0)
	rval = _mb->get_entities_by_type_and_tag(0, MBENTITYSET,
												&geom_tag, &dim_ptr, 1,
												vsets, 1, false);
	int numGeoVertices = vsets.size();
	//SmoothVertex ** smthVertex = new SmoothVertex *[numGeoVertices];
	_smthVertex = new SmoothVertex *[numGeoVertices];
	// there should also be a map from original vertex sets to new vertex sets (or just new vertex??)
	//std::map<MBEntityHandle, SmoothVertex*> mapVertices;

	i=0;
	for ( it = vsets.begin(); it!=vsets.end(); it++, i++)
	{
		MBEntityHandle vertex= *it;
		_smthVertex[i]= new SmoothVertex(_mb, vertex, _mbo);
		_mapVertices[vertex] = _smthVertex[i];
	}

	// _mb, mapSurfaces, mapCurves, mapVertices are characterizing the geometry/topology of the initial mesh

	//SmoothFaceEval geom_eval(_mb, _set);
	// initialize the smooth mesh evaluator
	//
	//geom_eval.Initialize(): it is decomposed in initializing first the gradients
	for (i=0; i<numSurfaces; i++)
	{
		_smthFace[i]->init_gradient();// this will also retrieve the triangles in each surface
		_smthFace[i]->compute_tangents_for_each_edge();// this one will consider all edges internal, so the
		// tangents are all in the direction of the edge; a little bit of waste, as we store
		// one tangent for each edge node , even though they are equal here...
		// no loops are considered
	}

	// this will be used to mark boundary edges, so for them the control points are computed earlier
	unsigned char value = 0; // default value is "not used"=0 for the tag
	// unsigned char def_data_bit = 1;// valid by default
	// rval = mb->tag_create("valid", 1, MB_TAG_BIT, validTag, &def_data_bit);
	MBTag markTag;
	rval = _mb->tag_create("MARKER", 1, MB_TAG_BIT, markTag, &value); // default value : 0 = not computed yet
	// each feature edge will need to have a way to retrieve at every moment the surfaces it belongs to
	// from edge sets, using the sense tag, we can get faces, and from each face, using the map, we can get
	// the SmoothFaceEval (surface evaluator), that has everything, including the normals!!!
	assert(rval==MB_SUCCESS);

	// create the tag also for control points on the edges
	double defCtrlPoints[9] = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };
	MBTag edgeCtrlTag;
	rval = _mb->tag_create("CONTROLEDGE", 9 * sizeof(double),
			MB_TAG_DENSE, edgeCtrlTag, &defCtrlPoints);
	if (MB_SUCCESS != rval)
		return MB_FAILURE;

	MBTag facetCtrlTag;
	double defControls[18] = { 0. };
	rval = _mb->tag_create("CONTROLFACE", 18 * sizeof(double),
			MB_TAG_DENSE, facetCtrlTag, &defControls);
	assert(rval == MB_SUCCESS);

	MBTag facetEdgeCtrlTag;
	double defControls2[27] = { 0. }; // corresponding to 9 control points on edges, in order from edge 0, 1, 2 ( 1-2, 2-0, 0-1 )
	rval = _mb->tag_create("CONTROLEDGEFACE", 27 * sizeof(double),
			MB_TAG_DENSE, facetEdgeCtrlTag, &defControls2);
	assert(rval == MB_SUCCESS);
	// if the
	double min_dot= -1.0; // depends on _angle, but now we ignore it, for the time being
	for (i=0; i<numCurves; i++)
	{
		_smthCurve[i]->compute_tangents_for_each_edge();// do we need surfaces now? or just the chains?
		// the computed edges will be marked accordingly; later one, only internal edges to surfaces are left
		_smthCurve[i]->compute_control_points_on_boundary_edges( min_dot,  _mapSurfaces, edgeCtrlTag, markTag);
	}

	// when done with boundary edges, compute the control points on all edges in the surfaces

	for (i=0;i<numSurfaces; i++)
	{
		// also pass the tags for
		_smthFace[i]->compute_control_points_on_edges(min_dot, edgeCtrlTag, markTag);
	}

	// now we should be able to compute the control points for the facets

	for (i=0;i<numSurfaces; i++)
	{
		// also pass the tags for edge and facet control points
		_smthFace[i]->compute_internal_control_points_on_facets(min_dot, facetCtrlTag, facetEdgeCtrlTag);
	}
	// we will need to compute the tangents for all edges in the model
	// they will be needed for control points for each edge
	// the boundary edges and the feature edges are more complicated
	// the boundary edges need to consider their loops, but feature edges need to consider loops and the normals
	// on each connected surface

	// some control points
	if (debug_surf_eval)
		for (i=0;i<numSurfaces; i++)
				_smthFace[i]->DumpModelControlPoints();

	return true;
}
#if 0
bool CamalPaveDriver::prepareCGMEvaluator()
{
	//
	// we will have to mesh every geo edge separately, and we have to ensure that the number of mesh edges
	// for a face is even.
	// pretty tough to do. Initially, we have to decide loops, number of edges on each face, etc
	// first build
	//int err;
	// get the triangles and the vertices from moab set

	/*iBase_EntityHandle *triangles = NULL;
	int triangles_alloc = 0;
	iBase_EntityHandle *vert_adj = NULL;
	int vert_adj_alloc = 0, vert_adj_size;
	int numTriangles;
	int * offsets = NULL, offsets_alloc = 0, indices_size;
	int * indices = NULL, indices_alloc = 0, offsets_size;
	iMesh_getAdjEntIndices(_meshIface, _set, iBase_FACE, iMesh_TRIANGLE,
			iBase_VERTEX, &triangles, &triangles_alloc, &numTriangles,
			&vert_adj, &vert_adj_alloc, &vert_adj_size, &indices,
			&indices_alloc, &indices_size, &offsets, &offsets_alloc,
			&offsets_size, &err);
	ERRORR("Couldn't get connectivity for triangles.", 1);*/

	MBRange triangles;
	MBErrorCode rval = _mb->get_entities_by_type( 0 /* root set, as above, we know */,
											   MBTRI, triangles);
	// get all the nodes
	MBRange vertices;
	rval = _mb->get_adjacencies(triangles, 0, false, vertices, MBInterface::UNION);

	// first, create CubitPointData list, from the coordinates in one array
	/* get the coordinates in one array */

	/*int vert_coords_alloc = 0, vertex_coord_size;
	double * xyz = NULL;
	iMesh_getVtxArrCoords(_meshIface, vert_adj, vert_adj_size,
			iBase_INTERLEAVED, &xyz, &vert_coords_alloc, &vertex_coord_size,
			&err);
	ERRORR("Couldn't get coordinates for vertices.", 1);*/

	// here, we use Cholla from CGM
	// we need to replace it with something equivalent, but simpler
	// the first try would be some tags in MOAB
	// create the cubit point data
	// initialize CGM
	AppUtil::instance()->startup(0, NULL);
	CGMApp::instance()->startup(0, NULL);

	// Initialize the GeometryTool
	GeometryQueryTool *gqt = GeometryQueryTool::instance();
	FacetModifyEngine *fme = FacetModifyEngine::instance();

	int vert_adj_size =  vertices.size();
	int numTriangles = triangles.size();
	DLIList<CubitFacet*> f_list(numTriangles);
	DLIList<CubitPoint*> p_list(vert_adj_size);
	double * xyz = new double [3*vert_adj_size];
	rval = _mb->  get_coords(vertices, xyz);
	//std::map<MBEntityHandle, CubitPoint *> mapPoints;
	//MBRange::iterator it = vertices.begin();
	for (int i = 0; i < vert_adj_size; i++/*, it++*/) {
		double * pCoord = &xyz[3 * i];
		CubitPointData * newPoint = new CubitPointData(pCoord[0], pCoord[1],
				pCoord[2]);
		p_list.append(newPoint);
		//mapPoints[*it] = newPoint;// or maybe we should use finding the index in MBRange??
	}

	// yes
	// define all the triangles, to see what we have
	for (MBRange::iterator it = triangles.begin(); it!=triangles.end(); it++) {
		MBEntityHandle tri = *it;
		int nnodes;
		const MBEntityHandle * conn3;//
		_mb->get_connectivity(tri, conn3, nnodes);
		assert(nnodes == 3);
		int vtri[3];// indices for triangles
		int ii = 0;
		for (ii = 0; ii < 3; ii++)
			vtri[ii] = vertices.index(conn3[ii]); // vtri[ii] = indices[offsets[j] + ii];
		CubitFacetData * triangle = new CubitFacetData(p_list[vtri[0]],
				p_list[vtri[1]], p_list[vtri[2]]);
		f_list.append(triangle);
	}

	DLIList<LoopSM*> my_loops;

	DLIList<Surface*> surf_list;
	CubitStatus result;
	//double angle = 0.01;// very small, negligible; is this radians or degrees?
	result = fme->build_facet_surface(NULL, f_list, p_list, _angle, 4, true,
			false, surf_list);

	if (surf_list.size() == 0 || result != CUBIT_SUCCESS) {
		PRINT_ERROR("Problems building mesh based surfaces.\n");
		return result;
	} else
		PRINT_INFO("Constructed %d surfaces.\n", surf_list.size());

	//Now build the shell.  If we had it set up right this would be
	//in a loop.  We need to store list of DLBlockSurfaceLists on each
	//blockvolumemesh to store the shell information.  But that will
	//be saved for later.
	ShellSM *shell_ptr;
	result = fme->make_facet_shell(surf_list, shell_ptr);

	if (shell_ptr == NULL || result != CUBIT_SUCCESS) {
		PRINT_ERROR("Problems building mesh based shell entity.\n");
		return result;
	}
#if 1
	DLIList<ShellSM*> shell_list;
	shell_list.append(shell_ptr);
	Lump *lump_ptr;
	result = fme->make_facet_lump(shell_list, lump_ptr);

	if (lump_ptr == NULL || result != CUBIT_SUCCESS) {
		PRINT_ERROR("Problems building mesh based lump entity.\n");
		return result;
	}
	DLIList<Lump*> lump_list;
	lump_list.append(lump_ptr);

	BodySM *bodysm_ptr;
	Body *body_ptr;
	result = fme->make_facet_body(lump_list, bodysm_ptr);

	body_ptr = GeometryQueryTool::instance()->make_Body(bodysm_ptr);

	if (body_ptr == NULL || result != CUBIT_SUCCESS) {
		PRINT_ERROR("Problems building mesh based body entity.\n");
		return result;
	}

	if (!body_ptr) {
		exit(1);
	}

	PRINT_INFO("Body successfully created.\n");
#endif
	PRINT_INFO("Number of vertices = %d\n", gqt->num_ref_vertices());
	PRINT_INFO("Number of edges = %d\n", gqt->num_ref_edges());
	PRINT_INFO("Number of faces = %d\n", gqt->num_ref_faces());

	// print vertex positions
	DLIList<RefVertex*> verts;
	gqt->ref_vertices(verts);

	int i;
	for (i = 0; i < verts.size(); i++) {
		RefVertex * vert = verts[i];
		CubitVector coords = vert->coordinates();
		PRINT_INFO("Vertex %d: %4.2f, %4.2f, %4.2f.\n", vert->id(), coords.x(),
				coords.y(), coords.z());
	}
	// print edges and faces

	DLIList<RefEdge*> refEdges;
	gqt->ref_edges(refEdges);

	for (i = 0; i < refEdges.size(); i++) {
		RefEdge * edg = refEdges[i];
		PRINT_INFO("Edge %d: %d %d\n", edg->id(), edg->start_vertex()->id(), edg->end_vertex ()->id() );
	}

	DLIList<RefFace*> refFaces;
	gqt->ref_faces(refFaces);

	for (i = 0; i < refFaces.size(); i++) {
		RefFace * face = refFaces[i];
		DLIList<  Loop * >  loop_list  ;
		face->ordered_loops (loop_list ) ;
		DLIList<  RefEdge  * > ordered_edge_list;
		loop_list[0]->ordered_ref_edges  (ordered_edge_list);

		//DLIList<  RefVertex* >  *listV = ref_vert_loop_list[0];
		PRINT_INFO("face %d: loop 0 size %d\n", face->id(),  ordered_edge_list.size() );
		for (int j=0; j<ordered_edge_list.size(); j++)
		{
			PRINT_INFO("  %d", ordered_edge_list[j]->id() );
		}
		PRINT_INFO("\n");
	}
	return true;
}
#endif

void CamalPaveDriver::mesh_vertices()
{
	//
	int numVertices = _mapVertices.size();
	for (int i=0; i<numVertices; i++)
	{
		/*if (_smthVertex[i]->isMeshed())
			continue;*/
		_smthVertex[i]-> create_mesh_vertex();
	}
	return;
} //CamalPaveDriver::mesh_vertices()

void CamalPaveDriver::mesh_curves()
{
	// hmmm
	// at this point, mesh counts are established;
	// we need to call meshers for each edge
	int numCurves = _mapCurves.size();
	int num_points;
	for (int i=0; i<numCurves; i++)
		_smthCurve[i]->create_mesh_edges(_mapVertices);
	return;
}// CamalPaveDriver::mesh_curves()

void CamalPaveDriver::mesh_surfaces()
{
	int numSurfaces = _mapSurfaces.size();
	int i=0;
	for (i=0; i<numSurfaces; i++)
	{
		_smthFace[i]->mesh( _mesh_size, _mapCurves);
	}
	for (i=0; i<numSurfaces; i++)
   {
      std::cout << " number of evaluations for surface index " << i <<
            ": "<<_smthFace[i]->eval_counter() << "\n";
   }
	return ;
}// CamalPaveDriver::mesh_surfaces()

void CamalPaveDriver::establish_mesh_curve_count()
{
	// hmmm
	// where do we start?
	// probably get the Edge Mesher from Camal
	// determine edge loops in the faces

	// first, estimate mesh count on every edge
	int numCurves = _mapCurves.size();
	int num_points;
	for (int i=0; i<numCurves; i++)
		_smthCurve[i]->estimate_mesh_count(this->_mesh_size, num_points);

	// check if the total mesh count is even for each face
	int numSurfaces = _mapSurfaces.size();
	int problem = 0;
	for (int i=0; i<numSurfaces; i++)
	{
		int mesh_count = 0;

		MBErrorCode rval = _smthFace[i]->mesh_count_total(_mapCurves, mesh_count);
		std::cout<<"face " << i << " boundary mesh count " << mesh_count<< std::endl;
		if (mesh_count%2 !=0)
		{
			// problem, redistribute mesh count on

			problem = 1;
			std::cout<<"        problem !!! odd boundary"<< std::endl;
		}
	}
	if (problem)
	{
		// need to do something, maybe even set everything to even
		std::cout<<"Set all to even \n";
		for (int i=0; i<numSurfaces; i++)
			_smthFace[i]->evenify(_mapCurves);
		return;
	}
	return;

}//establish_mesh_curve_count();

MBErrorCode CamalPaveDriver::find_face_loops()
{
	//
	int numSurfaces = _mapSurfaces.size();
	for (int i=0; i<numSurfaces; i++)
	{
		MBErrorCode rval = this->_smthFace[i]->find_loops();
		if (MB_SUCCESS!=rval)
			return rval;
	}
	return MB_SUCCESS;
}
