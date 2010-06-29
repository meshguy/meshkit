/**
 this file will prepare the paver call; basically, create the
 smooth surface, and pass it to the Pave routine call
 SmoothMeshEval will be derived from CAMAL eval, and will call the smooth faceter that I will 
 extract from CGM. The arithmetic is explained in Face-Based Surfaces
   for 3D Mesh generation
 // input
 */
//
// mesh something with paver from CAMAL
//
#include "CamalPaveDriver.hpp"

#include "SmoothMeshEval.hpp"

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

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
// end copy includes from facets.cpp
//#include "camal_interface.hpp"

// camal edge mesher: need to mesh with even numbers
#include "CMLEdgeMesher.hpp"
#include "SmoothCurveEval.hpp"

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
#include <assert.h>
#include <string.h>

const bool debugLocal = false;

CamalPaveDriver::CamalPaveDriver(iMesh_Instance mesh,
		iBase_EntitySetHandle set, iMesh_Instance omesh) {
	_meshIface = mesh;
	_meshOutface = omesh;
	_set = set;
}

bool CamalPaveDriver::remesh(double mesh_size, int mesh_intervals,
		const bool force_intervals) {

	SmoothMeshEval geom_eval(_meshIface, _set);
	// initialize the smooth mesh evaluator
	// compute the loops first, then construct the data structure needed for evaluation
	geom_eval.Initialize();
	std::vector<int> loops, loop_sizes;
	std::vector<iBase_EntityHandle> bdy_verts;
	std::vector<double> bdy_coords;
	std::vector<int> connect;
	bool success = true;
	int new_points;
	iMesh_EntityTopology etop = iMesh_ALL_TOPOLOGIES;
	// create the facetted shell

	// start copy from facets
	// first build
	int err;
	// get the triangles and the vertices in one shot
	iBase_EntityHandle *triangles = NULL;
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
	ERRORR("Couldn't get connectivity for triangles.", 1);

	// first, create CubitPointData list, from the coordinates in one array
	/* get the coordinates in one array */

	int vert_coords_alloc = 0, vertex_coord_size;
	double * xyz = NULL;
	iMesh_getVtxArrCoords(_meshIface, vert_adj, vert_adj_size,
			iBase_INTERLEAVED, &xyz, &vert_coords_alloc, &vertex_coord_size,
			&err);
	ERRORR("Couldn't get coordinates for vertices.", 1);

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

	DLIList<CubitFacet*> f_list(numTriangles);
	DLIList<CubitPoint*> p_list(vert_adj_size);
	for (int i = 0; i < vert_adj_size; i++) {
		double * pCoord = &xyz[3 * i];
		CubitPointData * newPoint = new CubitPointData(pCoord[0], pCoord[1],
				pCoord[2]);
		p_list.append(newPoint);
	}

	// define all the triangles, to see what we have
	for (int j = 0; j < numTriangles; j++) {
		int vtri[3];// indices for triangles
		int ii = 0;
		for (ii = 0; ii < 3; ii++)
			vtri[ii] = indices[offsets[j] + ii];
		CubitFacetData * triangle = new CubitFacetData(p_list[vtri[0]],
				p_list[vtri[1]], p_list[vtri[2]]);
		f_list.append(triangle);
	}

	DLIList<LoopSM*> my_loops;

	DLIList<Surface*> surf_list;
	CubitStatus result;
	double angle = 0.01;// very small, negligible
	result = fme->build_facet_surface(NULL, f_list, p_list, angle, 4, true,
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
	for (i = verts.size(); i > 0; i--) {
		CubitVector coords = verts.get_and_step()->coordinates();
		PRINT_INFO("Vertex %d: %4.2f, %4.2f, %4.2f.\n", 8 - i, coords.x(),
				coords.y(), coords.z());
	}

	RefFace *face = gqt->get_first_ref_face();

	RefEdge *edge = gqt->get_first_ref_edge();

	DLIList<RefEdge*> ref_edges;
	gqt->ref_edges(ref_edges);

#if CAMAL_VERSION > 500
	CAMALSizeEval size_eval(mesh_size);
#endif
	// we should have only one edge and one vertex
	// mesh it "uniformly" along parametric direction
	for (i = ref_edges.size(); i > 0; i--) {
		RefEdge * edgeC = ref_edges.get_and_step();
		if (!edgeC)
			continue;
		double mass = edgeC->measure();
		int numMeshEdges = (int) mass / mesh_size;
		numMeshEdges = 2 * (numMeshEdges / 2 + 1); // to make it odd// at least 2
		PRINT_INFO("Edge %d: length: %f .\n", i + 1, mass);
		// how to divide the length
		// now do a meshing on the curve, just to see how we are doing.
		// replicate mesh_curve in camel..
		// we will have to set it to run; we need the boundary coordinates,
		double  start_param, end_param;
		CubitBoolean res = edgeC->get_param_range(  start_param, end_param );
		// at some point we will eliminate ref edge and cgm
		SmoothCurveEval smthCrvEval(edgeC, &geom_eval, i-1);// loop index is 0 in general
		CMLEdgeMesher cmlEdgeMesher(&smthCrvEval, CML::STANDARD);
		cmlEdgeMesher.set_sizing_function(CML::LINEAR_SIZING);
		int num_points_out = 0;
		//! \note The number of points computed will be one more than the number
		//! of edge segments in the mesh.  For closed curves, the first and last
		//! points will be identical.
		double curve_mesh_size = mesh_size;

		CAMALSizeEval size_eval2(curve_mesh_size);
		cmlEdgeMesher.generate_mesh_sizing(&size_eval2, num_points_out);
		if (2 * (num_points_out / 2) == num_points_out) // we want it even nb seg, odd num_points_out
		{

			// now, the number of intervals will be set to num_points_out (even)
			int num_intervals = num_points_out;
			cmlEdgeMesher.generate_mesh_uniform(num_intervals, num_points_out);
			//curve_mesh_size = curve_mesh_size*((3*num_points_out-1.0)/3.0/num_points_out);
		}

		// we are now settled on an odd number of edge segments;
		// get them, and put them in the bdy_loops...
		bdy_coords.resize(num_points_out * 3);// we allocate more, but we will use one les
		bool stat = cmlEdgeMesher.get_mesh(num_points_out, &bdy_coords[0]);
		// create the loops and id
		loops.resize(num_points_out - 1);
		loop_sizes.resize(1);
		loop_sizes[0] = num_points_out - 1;
		loops.resize(num_points_out - 1);
		for (int j = 0; j < num_points_out - 1; j++)
			loops[j] = j;// simple connectivity for the loop
		// do something else;  use position from u
		/*
		double delta_u = (end_param-start_param)/(num_points_out - 1);
		CubitVector output_position;
		for (int i=0; i<num_points_out-1; i++)
		{
			double u_value = start_param + i*delta_u;
			edgeC->position_from_u ( u_value,
					output_position);
			output_position.get_xyz( &bdy_coords[i*3]);
		}*/

	}

	// end copy from facets
	// so at this moment we have the edge, and the face;
	// we need to mesh the ref edges first
#if 0
	success = CAMAL_bdy_loops_coords(cmel, gentity, loops, loop_sizes,
			bdy_verts, bdy_coords);
#endif 

	if (!success) {
		std::cerr << "Couldn't get bounding mesh for surface." << std::endl;
		return success;
	}

	if (mesh_size < 0.) {
		double tot_length = 0.;
		unsigned int start_current_loop = 0;
		for (unsigned int k = 0; k < loop_sizes.size(); k++) {
			// for each loop, compute the edge lengths individually
			int current_loop_size = loop_sizes[k];
			for (unsigned int i = 0; i < current_loop_size; i++) {
				unsigned int i1 = loops[start_current_loop + i];
				unsigned int i2 = loops[start_current_loop + (i + 1)
						% current_loop_size];
				tot_length += sqrt((bdy_coords[3 * i1] - bdy_coords[3 * i2])
						* (bdy_coords[3 * i1] - bdy_coords[3 * i2])
						+ (bdy_coords[3 * i1 + 1] - bdy_coords[3 * i2 + 1])
								* (bdy_coords[3 * i1 + 1] - bdy_coords[3 * i2
										+ 1]) + (bdy_coords[3 * i1 + 2]
						- bdy_coords[3 * i2 + 2]) * (bdy_coords[3 * i1 + 2]
						- bdy_coords[3 * i2 + 2]));
			}
			start_current_loop += current_loop_size;
		}
		mesh_size = tot_length / start_current_loop;
		geom_eval.set_mesh_size(mesh_size);
	}

	if (debugLocal) {
		std::ofstream oFile("points.Point3D", std::ios::out);
		if (!oFile)
			return false;
		oFile.precision(12);
		std::cout << "Surface " << " " // cmel->get_gentity_id(gentity)
				<< ", mesh_size = " << mesh_size << ", boundary mesh: "
				<< std::endl;
		oFile<< "# x y z" << std::endl;
		std::cout << bdy_coords.size() / 3 - 1 << "  " << loop_sizes.size()
				<< std::endl;
		for (unsigned int i = 0; i < bdy_coords.size() / 3 - 1; i++)
		{
			std::cout << bdy_coords[3 * i] << "  " << bdy_coords[3 * i + 1]
					<< "  " << bdy_coords[3 * i + 2] << std::endl;
			oFile << bdy_coords[3 * i] << "  " << bdy_coords[3 * i + 1]
								<< "  " << bdy_coords[3 * i + 2] << std::endl;
		}
		//oFile
		for (std::vector<int>::iterator vit = loop_sizes.begin(); vit
				!= loop_sizes.end(); vit++)
			std::cout << *vit << "  ";

		std::cout << std::endl;
		for (std::vector<int>::iterator vit = loops.begin(); vit != loops.end(); vit++)
			std::cout << *vit << std::endl;
	}

	// pass to CAMAL
	geom_eval.set_ref_face(face);
	CMLPaver pave_mesher(&geom_eval, &size_eval);

	// set only num_points_out -1 , because the last one is repeated
	success = pave_mesher.set_boundary_mesh(bdy_coords.size() / 3 - 1,
			&bdy_coords[0], (int) loop_sizes.size(), &loop_sizes[0], &loops[0]);
	if (!success) {
		std::cerr << "Failed setting boundary mesh" << std::endl;
		return success;
	}

	pave_mesher.set_sizing_function(CML::LINEAR_SIZING);

	// generate the mesh
	int num_quads;
	success = pave_mesher.generate_mesh(new_points, num_quads);
	if (!success) {
		std::cerr << "Failed generating mesh" << std::endl;
		//return success;
	}

	std::cout << "Meshed surface with " << new_points << " new vertices and "
			<< num_quads << " quadrilaterals." << std::endl;

	// get the generated mesh
	//bdy_coords.resize(3*(bdy_verts.size() + new_points));
	std::vector<double> new_coords;
	new_coords.resize(3 * new_points);
	connect.resize(4 * num_quads);
	success = pave_mesher.get_mesh(new_points, &new_coords[0], num_quads,
			&connect[0]);
	if (!success) {
		std::cerr << "Failed get generated mesh" << std::endl;
		//return success;
	}

	// we assume that the first points from pave mesher are the boundary vertices
	// (we should verify that, actually)
	// maybe for the time being, create the whole set again?
	// should we worry about the boundary when we do not know if it fails or not?

	//etop = iMesh_QUADRILATERAL;
	// add to a new set the elements created, anyway, just to see them in visit?
	std::vector<iBase_EntityHandle> qverts;
	qverts.resize (new_points);
	iBase_EntityHandle *verts_ptr = &qverts[0];
	int new_verts_size = new_points, new_verts_alloc =
			new_verts_size;

	// we recreate all the vertices, just because; the original ones will not be used, but this is
	// life; maybe later we will just create the truly new vertices only, and assume that they are
	// exactly at the beginning of the
	iMesh_createVtxArr(_meshOutface, new_verts_size, iBase_INTERLEAVED,
			&new_coords[0], 3 * new_verts_alloc, &verts_ptr,
			&new_verts_alloc, &new_verts_size, &err);
	ERRORR("Couldn't create new vertices.",false);
	iBase_EntitySetHandle outset2;
	iMesh_createEntSet(_meshOutface, false, &outset2, &err);
	iMesh_addEntArrToSet(_meshOutface, &qverts[0], qverts.size(), outset2,
			&err);
	// create the quads
	// start copy
	  // create them
	 // assemble connectivity of new elements
	  std::vector<iBase_EntityHandle> new_connect(connect.size());
	  //int old_verts = bdy_verts.size();
	  for (unsigned int i = 0; i < connect.size(); i++)
	  //{
	    //if (connect[i] >= old_verts)
	      new_connect[i] = qverts[connect[i]];
	    //else
	    //  new_connect[i] = bdy_verts[connect[i]];
	  //}
	  int *status = NULL;
	  int status_size = 0, status_alloc = 0;
	  iBase_EntityHandle *new_ments = NULL;
	  int new_ments_size = 0, new_ments_alloc = 0;
	  iMesh_createEntArr(_meshOutface, iMesh_QUADRILATERAL, &new_connect[0], connect.size(),
	                     &new_ments, &new_ments_alloc, &new_ments_size,
	                     &status, &status_alloc, &status_size,
	                     &err);
	  ERRORR("Couldn't create new elements.", false);

	    // put them into the new_entities vector
	  //new_entities.resize(new_ments_size);
	  //std::copy(new_ments, new_ments+new_ments_size, &new_entities[0]);

	    // add elements to the geometric owner's set
	  iMesh_addEntArrToSet(_meshOutface, new_ments, new_ments_size,
			  outset2, &err);
	// end copy

	return true;
}

