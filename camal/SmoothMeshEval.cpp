#include <iostream>
#include "SmoothMeshEval.hpp"
#include "RefFace.hpp"
//#include "moab/GeomTopoTool.hpp"

#include <algorithm>

// included in the header now
// #include "MBRange.hpp"
// #include "MBCartVect.hpp"


bool within_tolerance(MBCartVect & p1, MBCartVect & p2,
		const double & tolerance) {
	if ((fabs(p1[0] - p2[0]) < tolerance) && (fabs(p1[1] - p2[1]) < tolerance)
			&& (fabs(p1[2] - p2[2]) < tolerance))
		return true;
	return false;
}
int numAdjTriInSet(MBInterface * mb, MBEntityHandle startEdge,
		MBEntityHandle set) {
	std::vector<MBEntityHandle> adjTri;
	mb->get_adjacencies(&startEdge, 1, 2, false, adjTri, MBInterface::UNION);
	// count how many are in the set
	int nbInSet = 0;
	for (int i = 0; i < adjTri.size(); i++) {
		MBEntityHandle tri = adjTri[i];
		if (mb->contains_entities(set, &tri, 1))
			nbInSet++;
	}
	return nbInSet;
}
extern bool debug_surf_eval;
extern bool use_cgm;

void SmoothMeshEval::set_mesh_size(double tmp_size) {
	if (debug_surf_eval) {
		std::cout << "set_mesh_size called." << std::endl;
	}
	meshSize = tmp_size;
}

double SmoothMeshEval::get_mesh_size() {
	if (debug_surf_eval) {
		std::cout << "get_mesh_size called." << std::endl;
	}
	return meshSize;
}

void SmoothMeshEval::set_ref_face(RefFace * refFace) {
	_smooth_face = refFace;
}
SmoothMeshEval::SmoothMeshEval(iMesh_Instance meshIface,
		iBase_EntitySetHandle surface_set) :
	_meshIface(meshIface), _surf_set(surface_set) {
	_smooth_face = NULL;
}

SmoothMeshEval::~SmoothMeshEval() {
}

double SmoothMeshEval::area() {
	// find the area of this entity
	assert(_smooth_face);
	double area1 = _smooth_face->area();
	double totArea = 0.;
	for (MBRange::iterator it = _triangles.begin(); it != _triangles.end(); it++) {
		MBEntityHandle tria = *it;
		const MBEntityHandle * conn3;
		int nnodes;
		_mb->get_connectivity(tria, conn3, nnodes);
		//
		//double coords[9]; // store the coordinates for the nodes
		//_mb->get_coords(conn3, 3, coords);
		MBCartVect p[3];
		_mb->get_coords(conn3, 3, (double*) &p[0]);
		// need to compute the angles
		// compute angles and the normal
		//MBCartVect p1(&coords[0]), p2(&coords[3]), p3(&coords[6]);

		MBCartVect AB(p[1] - p[0]);//(p2 - p1);
		MBCartVect BC(p[2] - p[1]);//(p3 - p2);
		MBCartVect normal = AB * BC;
		totArea += normal.length() / 2;
	}
	return totArea;
}

void SmoothMeshEval::bounding_box(double box_min[3], double box_max[3]) {
	assert(_smooth_face);
	CubitBox box = _smooth_face->bounding_box();
	CubitVector Minim = box.minimum();
	CubitVector Maxim = box.maximum();
	Minim.get_xyz(box_min);
	Maxim.get_xyz(box_max);
	for (int i = 0; i < 3; i++) {
		box_min[i] = _minim[i];
		box_max[i] = _maxim[i];
	}
	// _minim, _maxim
}

void SmoothMeshEval::move_to_surface(double& x, double& y, double& z) {
	assert(_smooth_face);
	CubitVector location(x, y, z);
	MBCartVect loc2(x, y, z);
	bool trim = false;// is it needed?
	bool outside = true;
	MBCartVect closestPoint;
        if (use_cgm)
	{
           _smooth_face->move_to_surface(location);
	    location.get_xyz(x, y, z);
        }
        else
        {
	// try our luck
	// not interested in normal
	      MBErrorCode rval = project_to_facets_main(loc2, trim, outside, &closestPoint, NULL);
               x = closestPoint[0];
               y = closestPoint[1];
               z = closestPoint[2];
         }
#if 0

	int result;
	double ori[3] = {x, y, z};
	iGeom_getEntClosestPt(geomIface, myEnt, x, y, z, &x, &y, &z, &result);
	if (iBase_SUCCESS != result) {
		std::cerr << "Trouble getting closest point on gentity." << std::endl;
	}
	if (debug_surf_eval) {
		std::cout << "myEnt=" << myEnt << " is moved from "
		<< ori[0] << "," << ori[1] << "," << ori[2] << " to "
		<< x << "," << y << "," << z << std::endl;
	}
#endif
}

void SmoothMeshEval::move_to_surface(double& x, double& y, double& z,
		double& u_guess, double& v_guess) {
	if (debug_surf_eval) {
		std::cout << "move_to_surface called." << std::endl;
	}
}

bool SmoothMeshEval::normal_at(double x, double y, double z, double& nx,
		double& ny, double& nz) {
	// normal at
	assert(_smooth_face);
	CubitVector location(x, y, z);
	MBCartVect loc2(x, y, z);
        if (use_cgm)
        { 
	   CubitVector v = _smooth_face->normal_at(location);
	   v.get_xyz(nx, ny, nz);
        }
        else
        {			// try our luck
	   bool trim = false;// is it needed?
	   bool outside = true;
	//MBCartVect closestPoint;// not needed
	// not interested in normal
	   MBCartVect normal;
	   MBErrorCode rval = project_to_facets_main(loc2, trim, outside, NULL, &normal);
		nx = normal[0];
		ny = normal[1];
		nz = normal[2];
	}
	// compare now the normals
	return true;
}

bool SmoothMeshEval::normal_at(double x, double y, double z, double& u_guess,
		double& v_guess, double& nx, double& ny, double& nz) {
	if (debug_surf_eval) {
		std::cout << "normal_at called." << std::endl;
	}
	return true;
}

bool SmoothMeshEval::is_planar() {
	if (debug_surf_eval) {
		std::cout << "is_planar called." << std::endl;
	}
	return false;
}

bool SmoothMeshEval::is_parametric() {
	if (debug_surf_eval) {
		std::cout << "is_parametric called." << std::endl;
	}
	return false;
}

bool SmoothMeshEval::is_periodic_in_u(double& u_period) {
	if (debug_surf_eval) {
		std::cout << "is_periodic_in_u called." << std::endl;
	}
	return false;
}

bool SmoothMeshEval::is_periodic_in_v(double& v_period) {
	if (debug_surf_eval) {
		std::cout << "is_periodic_in_v called." << std::endl;
	}
	return false;
}

void SmoothMeshEval::get_param_range_u(double& u_low, double& u_high) {
	if (debug_surf_eval) {
		std::cout << "is_param_range_u called." << std::endl;
	}
}

void SmoothMeshEval::get_param_range_v(double& v_low, double& v_high) {
}

bool SmoothMeshEval::uv_from_position(double x, double y, double z, double& u,
		double& v)

{
	if (debug_surf_eval) {
		std::cout << "uv_from_position." << std::endl;
	}
	return true;
}

bool SmoothMeshEval::uv_from_position(double x, double y, double z, double& u,
		double& v, double& cx, double& cy, double& cz) {
	if (debug_surf_eval) {
		std::cout << "uv_from_position." << std::endl;
	}
	return true;
}

void SmoothMeshEval::position_from_uv(double u, double v, double& x, double& y,
		double& z) {
	if (debug_surf_eval) {
		std::cout << "position_from_uv." << std::endl;
	}
}

void SmoothMeshEval::distortion_at_uv(double u, double v, double du[3],
		double dv[3]) {
	if (debug_surf_eval) {
		std::cout << "distortion_at_uv." << std::endl;
	}
}

int SmoothMeshEval::Initialize() {
	// first get the MOAB set, create the loops, etc.
	// first, we need to delete all the edges in MOAB, and create new ones only for
	// the set we are interested in; in this way, the orientation of the edges will be fine,
	// counterclockwise; also, check for
	_mb = (MBInterface *) _meshIface;
	_mbset = (MBEntityHandle) this->_surf_set;
	// get all the edges from this subset
	if (NULL == _mb)
		return 1; //fail
	_triangles.clear();
	MBErrorCode rval = _mb->get_entities_by_type(_mbset, MBTRI, _triangles);
	if (MB_SUCCESS != rval)
		return 1;
	// get a new range of edges, and decide the loops from here
	_edges.clear();
	rval = _mb-> get_adjacencies(_triangles, 1, true, _edges,
			MBInterface::UNION);

	rval = _mb->get_adjacencies(_triangles, 0, false, _nodes,
			MBInterface::UNION);

	// mark all edges as unused yet, for the loop creation
	// tag_create( "my_tag", sizeof(double), MB_TAG_DENSE, tag_handle, &value );
	unsigned char value = 0; // default value is "not used"=0 for the tag
	// unsigned char def_data_bit = 1;// valid by default
	// rval = mb->tag_create("valid", 1, MB_TAG_BIT, validTag, &def_data_bit);
	rval = _mb->tag_create("MARKER", 1, MB_TAG_BIT, _markTag, &value);
	// now loop over the edges and create boundary loops (maybe one only)
	// how do we make sure that the loops are ccw?
	// we don't; we just assume that the edges on the boundary will be created in ccw manner
	// by MOAB

	int currentLoopIndex = -1;
	//_loopEnds.push_back(0);// the first loop (current), starts at 0
	for (MBRange::iterator it = _edges.begin(); it != _edges.end(); it++) {
		//
		MBEntityHandle currentEdge = *it;
		unsigned char tagVal = 0;
		_mb->tag_get_data(_markTag, &currentEdge, 1, &tagVal);
		if (tagVal)
			continue; // already used
		// see if it was not marked already
		// determine number of triangles that are part of the set
		// if there is only one triangle, it is a boundary edge
		int nbAdj = numAdjTriInSet(_mb, currentEdge, _mbset);
		if (nbAdj == 1) {
			// start a loop at this edge because it was not marked already
			unsigned char used = 1;
			_mb->tag_set_data(_markTag, &currentEdge, 1, &used);// so it will not be part of
			// other loops in this face
			currentLoopIndex++; // the first time it will be 0;
			_loopEnds.push_back(currentLoopIndex);
			_loopLengths.push_back(0);
			_loopEnds[currentLoopIndex]++;
			_loops.push_back(currentEdge);
			_senses.push_back(1); // forward
			const MBEntityHandle * conn2;
			int nnodes;
			_mb->get_connectivity(currentEdge, conn2, nnodes);
			MBEntityHandle startNode = conn2[0];
			MBEntityHandle currentNode = conn2[1];
			MBCartVect vN[2];
			_mb->get_coords(conn2, 2, (double*) &vN[0]);
			MBCartVect P01(vN[1] - vN[0]);
			double len1 = P01.length();
			_loopLengths[currentLoopIndex] = len1;
			_fractions.push_back(len1);
			int i = _loopEnds[currentLoopIndex] - 1;// the index in the loop

			while (currentNode != startNode) {
				// find the next edge in the loop, with adjacency 1
				std::vector<MBEntityHandle> edgesAdj;
				rval = _mb->get_adjacencies(&currentNode, 1, 1, false,
						edgesAdj, MBInterface::UNION);
				for (int j = 0; j < edgesAdj.size(); j++) {
					MBEntityHandle edg = edgesAdj[j];
					_mb->tag_get_data(_markTag, &edg, 1, &tagVal);
					if (tagVal)
						continue; // the edge was already marked, although
					// it could be only the current edge
					if (edg == currentEdge)
						continue; // we should get rid of the tag maybe

					int nb = numAdjTriInSet(_mb, edg, _mbset);
					if (nb > 1)
						continue;
					// if we get here, this is what we want
					_mb->get_connectivity(edg, conn2, nnodes);
					char senseC = 0;
					if (conn2[0] == currentNode) {
						currentNode = conn2[1]; // next node
						senseC = 1; // forward
					} else {
						currentNode = conn2[0]; // next node
						senseC = 2; // reverse
					}
					_loops.push_back(edg);
					_loopEnds[currentLoopIndex]++;
					_mb->get_coords(conn2, 2, (double*) &vN[0]);
					MBCartVect P01(vN[1] - vN[0]);
					len1 = P01.length();
					// increase the total length
					_loopLengths[currentLoopIndex] += len1;
					double lastLength = _fractions[i];
					_fractions.push_back(lastLength + len1);
					i++; // increment the index in the loop for next accumulated length
					_senses.push_back(senseC);
					currentEdge = edg;
					_mb->tag_set_data(_markTag, &currentEdge, 1, &used);// so it will not be part of
					// other loops in this face
					break; // we should advance from the for{} to another point
				}
				// connected to current vertex, unmarked
			}
			// now we need to modify _fractions, to divide by the length of the loop
			int startCurrentLoop = (currentLoopIndex == 0) ? 0
					: _loopEnds[currentLoopIndex - 1];
			int endCurrentLoop = _loopEnds[currentLoopIndex];
			assert(_loopLengths[currentLoopIndex] > 0);
			for (int j = startCurrentLoop; j < endCurrentLoop; j++) {
				_fractions[j] /= _loopLengths[currentLoopIndex];
			}
		}

	}
	// print the loops size and some other stuff
	if (debug_surf_eval) {
		std::cout << " _loops.size() " << _loops.size() << std::endl;
		int startLoop = 0;
		for (int i = 0; i < _loopEnds.size(); i++) {
			int endLoop = _loopEnds[i];
			std::cout << "loop no " << i << "  start loop " << startLoop
					<< " end loop " << endLoop << std::endl;
			for (int j = startLoop; j < endLoop; j++) {
				MBEntityHandle edgeH = _loops[j];
				const MBEntityHandle * conn2;
				int nnodes;
				_mb->get_connectivity(edgeH, conn2, nnodes);
				std::cout << " index " << j << " Edge id "
						<< _mb->id_from_handle(edgeH) << " sense "
						<< (int) _senses[j] << " nodes: "
						<< _mb->id_from_handle(conn2[0]) << " "
						<< _mb->id_from_handle(conn2[1]) << std::endl;
			}
			startLoop = endLoop;
		}
	}
	// so we have now the boundary edges in order.
	// need to set the smooth facets; start with boundary edges
	// start with initialization of the gradient (normal field)

	init_gradient();

	// all the tags on the boundary have now _markTag 1
	// we will use that information when we compute the tangents for each edge

	// init edge control points
	// some are on the boundary, some are interior, with 2 faces adjacent
	// setup first the boundary edges, and establish their control points

	// we have the loops, and we know that there is only one boundary

	// at some point we should put this back into the game
	double min_dot = -1.; // very small for the time being , allow every angle
	compute_tangents_for_each_edge();

	compute_control_points_on_edges(min_dot);

	// initialize bounding box limits
	MBCartVect vert1;
	MBEntityHandle v1 = *(_nodes.begin());// first vertex
	_mb->get_coords(&v1, 1, (double*) &vert1);
	_minim = vert1;
	_maxim = vert1;
	// compute the internal control points on each facet, from
	// edge control points

	// also adjust the bounding box minim and maxim
	compute_internal_control_points_on_facets(min_dot);

	// initialize geom topo tool, and construct obb tree
	initialize_geom_topo_tool();

	return 0; // success
}
MBErrorCode SmoothMeshEval::initialize_geom_topo_tool() {
	_my_geomTopoTool = new moab::GeomTopoTool(_mb);
	// construct the obb tree!!
	MBErrorCode rval = _my_geomTopoTool->find_geomsets(_gsets);
	assert(rval == MB_SUCCESS);
	rval = _my_geomTopoTool->construct_obb_trees();
	// we know that we have only one surf set; pretty lame
	assert(rval == MB_SUCCESS);
	rval = _my_geomTopoTool->get_root(_gsets[2][0], _obb_root);

	assert(rval == MB_SUCCESS);

	return MB_SUCCESS;
}

MBErrorCode SmoothMeshEval::compute_control_points_on_edges(double min_dot) {
	double defCtrlPoints[9] = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };
	MBErrorCode rval = _mb->tag_create("CONTROLEDGE", 9 * sizeof(double),
			MB_TAG_DENSE, _edgeCtrlTag, &defCtrlPoints);
	if (MB_SUCCESS != rval)
		return MB_FAILURE;
	// now, compute control points for all edges
	for (MBRange::iterator it = _edges.begin(); it != _edges.end(); it++) {
		MBEntityHandle edg = *it;
		//double min_dot;
		init_bezier_edge(edg, min_dot);

	}
	return MB_SUCCESS;
}

int SmoothMeshEval::init_gradient() {
	// first, create a Tag for gradient (or normal)
	// loop over all triangles in set, and modify the normals according to the angle as weight
	double defNormal[] = { 0., 0., 0. };
	MBErrorCode rval = _mb->tag_create("GRADIENT", 3 * sizeof(double),
			MB_TAG_DENSE, _gradientTag, &defNormal);

	double defPlane[4] = { 0., 0., 1., 0. };
	rval = _mb->tag_create("PLANE", 4 * sizeof(double), MB_TAG_DENSE,
			_planeTag, &defPlane);
	// the fourth double is for weight, accumulated at each vertex so far
	// maybe not needed in the end
	for (MBRange::iterator it = _triangles.begin(); it != _triangles.end(); it++) {
		MBEntityHandle tria = *it;
		const MBEntityHandle * conn3;
		int nnodes;
		_mb->get_connectivity(tria, conn3, nnodes);
		if (nnodes != 3)
			return 1; // error
		//double coords[9]; // store the coordinates for the nodes
		//_mb->get_coords(conn3, 3, coords);
		MBCartVect p[3];
		_mb->get_coords(conn3, 3, (double*) &p[0]);
		// need to compute the angles
		// compute angles and the normal
		//MBCartVect p1(&coords[0]), p2(&coords[3]), p3(&coords[6]);

		MBCartVect AB(p[1] - p[0]);//(p2 - p1);
		MBCartVect BC(p[2] - p[1]);//(p3 - p2);
		MBCartVect CA(p[0] - p[2]);//(p1 - p3);
		double a[3];
		a[1] = angle(AB, -BC); // angle at B (p2), etc.
		a[2] = angle(BC, -CA);
		a[0] = angle(CA, -AB);
		MBCartVect normal = -AB * CA;
		normal.normalize();
		double plane[4];

		const double * coordNormal = normal.array();

		plane[0] = coordNormal[0];
		plane[1] = coordNormal[1];
		plane[2] = coordNormal[2];
		plane[3] = -normal % p[0]; // dot product
		//set the plane
		rval = _mb->tag_set_data(_planeTag, &tria, 1, plane);
		assert(rval == MB_SUCCESS);

		// add to each vertex the tag value the normal multiplied by the angle
		double values[9];

		_mb->tag_get_data(_gradientTag, conn3, 3, values);
		for (int i = 0; i < 3; i++) {
			//values[4*i]+=a[i]; // first is the weight, which we do not really need
			values[3 * i + 0] += a[i] * coordNormal[0];
			values[3 * i + 1] += a[i] * coordNormal[1];
			values[3 * i + 2] += a[i] * coordNormal[2];
		}

		// reset those values
		_mb->tag_set_data(_gradientTag, conn3, 3, values);

	}
	// normalize the gradients at each node; maybe not needed here?
	// no we will do it, it is important
	int numNodes = _nodes.size();
	double * normalVal = new double[3 * numNodes];
	_mb->tag_get_data(_gradientTag, _nodes, normalVal);
	for (int i = 0; i < numNodes; i++) {
		MBCartVect p1(&normalVal[3 * i]);
		p1.normalize();
		p1.get(&normalVal[3 * i]);
	}

	_mb->tag_set_data(_gradientTag, _nodes, normalVal);
	// print the loops size and some other stuff
	if (debug_surf_eval) {
		std::cout << " normals at  " << numNodes << " nodes" << std::endl;
		int i = 0;
		for (MBRange::iterator it = _nodes.begin(); it != _nodes.end(); it++, i++) {
			MBEntityHandle node = *it;
			std::cout << " Node id " << _mb->id_from_handle(node) << "  "
					<< normalVal[3 * i] << " " << normalVal[3 * i + 1] << " "
					<< normalVal[3 * i + 2] << std::endl;
		}
		delete[] normalVal; // free this memory
	}
	return 0;
}

// init bezier edges
// start copy
//===========================================================================
//Function Name: init_bezier_edge
//
//Member Type:  PRIVATE
//Description:  compute the control points for an edge
//===========================================================================
MBErrorCode SmoothMeshEval::init_bezier_edge(MBEntityHandle edge,
		double min_dot) {
	// min dot was used for angle here
	//int stat = 0; // CUBIT_SUCCESS;
	// all boundaries will be simple, initially
	// we may complicate them afterwards
#if 0
	TDFacetBoundaryEdge *td_bfe =
	TDFacetBoundaryEdge::get_facet_boundary_edge( edge );
	if (td_bfe)
	{
		stat = td_bfe->init_control_points( min_dot );
		if (stat != CUBIT_SUCCESS)
		return stat;
		stat = td_bfe->merge_control_points();
	}
	else
	{
#endif

	MBCartVect ctrl_pts[3];
	int nnodes;
	const MBEntityHandle * conn2;
	MBErrorCode rval = _mb->get_connectivity(edge, conn2, nnodes);
	assert(rval == MB_SUCCESS);
	assert(2 == nnodes);
	//double coords[6]; // store the coordinates for the nodes
	MBCartVect P[2];
	//MBErrorCode rval = _mb->get_coords(conn2, 2, coords);
	rval = _mb->get_coords(conn2, 2, (double*) &P[0]);
	assert(rval == MB_SUCCESS);

	//MBCartVect P0(&coords[0]);
	//MBCartVect P3(&coords[3]);

	//double normalVec[6];
	MBCartVect N[2];
	//_mb->tag_get_data(_gradientTag, conn2, 2, normalVec);
	rval = _mb->tag_get_data(_gradientTag, conn2, 2, (double*) &N[0]);
	assert(rval == MB_SUCCESS);

	//MBCartVect N0(&normalVec[0]);
	//MBCartVect N3(&normalVec[3]);
	MBCartVect T[2]; // T0, T3
	//if (edge->num_adj_facets() <= 1) {
	//stat = compute_curve_tangent(edge, min_dot, T0, T3);
	//  if (stat != CUBIT_SUCCESS)
	//	return stat;
	//} else {
	//}
	rval = _mb->tag_get_data(_tangentsTag, &edge, 1, &T[0]);
	assert(rval == MB_SUCCESS);

	rval = init_edge_control_points(P[0], P[1], N[0], N[1], T[0], T[1],
			ctrl_pts);
	assert(rval == MB_SUCCESS);

	// set  a tag with control points on the edge
	//edge->control_points(ctrl_pts, 4);
	rval = _mb->tag_set_data(_edgeCtrlTag, &edge, 1, &ctrl_pts[0]);
	assert(rval == MB_SUCCESS);

	if (debug_surf_eval) {
		std::cout << "edge: " << _mb-> id_from_handle(edge) << " tangents: "
				<< T[0] << T[1] << std::endl;
		std::cout << "  points: " << P[0] << " " << P[1] << std::endl;
		std::cout << "  normals: " << N[0] << " " << N[1] << std::endl;
		std::cout << "  Control points  " << ctrl_pts[0] << " " << ctrl_pts[1]
				<< " " << ctrl_pts[2] << std::endl;
	}

	return MB_SUCCESS;
}
// end copy

MBErrorCode SmoothMeshEval::compute_tangents_for_each_edge()
// they will be used for control points
{
	double defTangents[6] = { 0., 0., 0., 0., 0., 0. };
	MBErrorCode rval = _mb->tag_create("TANGENTS", 6 * sizeof(double),
			MB_TAG_DENSE, _tangentsTag, &defTangents);
	if (MB_SUCCESS != rval)
		return MB_FAILURE;
	// for each edge that is marked
	// first, edges in a loop, they are marked also
	int startLoop = 0;
	for (int i = 0; i < _loopEnds.size(); i++) {
		int endLoop = _loopEnds[i];
		std::cout << "loop no " << i << "  start loop " << startLoop
				<< " end loop " << endLoop << std::endl;
		// we will look at next 2 edges (past the end too)
		// we also assume that the senses are the same in the edges, along the boundary
		int sizeLoop = (i == 0) ? endLoop : (_loopEnds[i] - _loopEnds[i - 1]);
		for (int j = startLoop; j < endLoop; j++) {
			int index1 = (j + 1) % sizeLoop;
			int index2 = (j + 2) % sizeLoop;
			MBEntityHandle prevEdge = _loops[j];
			MBEntityHandle currentEdge = _loops[index1];
			MBEntityHandle nextEdge = _loops[index2];
			assert(_senses[j] == _senses[index1]);
			assert(_senses[j] == _senses[index2]);

			const MBEntityHandle * conn2;// we are truly interested in 4 points position, for each edge
			int nnodes;
			_mb->get_connectivity(prevEdge, conn2, nnodes);
			assert(nnodes == 2);
			const MBEntityHandle * conn4;//
			_mb->get_connectivity(nextEdge, conn4, nnodes);
			assert(nnodes == 2);
			// this
			/*std::cout << " index " << j << " Edge id "
			 << _mb->id_from_handle(edgeH) << " sense "
			 << (int) _senses[j] << " nodes: "
			 << _mb->id_from_handle(conn2[0]) << " "
			 << _mb->id_from_handle(conn2[1]) << std::endl;*/
			//start copy

			MBCartVect PP[4]; // store the coordinates for the nodes
			rval = _mb->get_coords(conn2, 2, (double*) &PP[0]); // populate with vertex 0 and 1
			assert(MB_SUCCESS == rval);
			rval = _mb->get_coords(conn4, 2, (double*) &PP[2]);// the nodes for next edge
			assert(MB_SUCCESS == rval);

			MBCartVect T[2]; // we will set the tangents tag for it
			T[0] = PP[2] - PP[0];
			T[0].normalize();
			T[1] = PP[3] - PP[1];
			T[1].normalize();
			// if the senses are reversed, the tangents are negative
			if (2 == (int) _senses[j]) {
				T[0] = -T[0];
				T[1] = -T[1];
			}
			// set those values
			_mb->tag_set_data(_tangentsTag, &currentEdge, 1, &T[0]);
			/*CubitPoint *p0 = edge->point( 0 );
			 CubitPoint *p1 = edge->point( 1 );
			 CubitFacetEdge *prev_edge = next_boundary_edge( edge, p0 );

			 CubitPoint *p2 = prev_edge->other_point( p0 );
			 CubitVector e0 = p0->coordinates() - p2->coordinates();
			 CubitVector e1 = p1->coordinates() - p0->coordinates();
			 e0.normalize();
			 e1.normalize();
			 if (e0 % e1 >= min_dot)
			 {
			 T0 = (p0->coordinates() - p2->coordinates()) +
			 (p1->coordinates() - p0->coordinates());
			 T0.normalize();
			 }
			 else
			 {
			 T0 = e1;
			 }*/
			// end copy
		}
		startLoop = endLoop;
	}
	// now, compute Tangents for all edges that are not on boundary, so they are not marked
	for (MBRange::iterator it = _edges.begin(); it != _edges.end(); it++) {
		MBEntityHandle edg = *it;
		unsigned char tagVal = 0;
		_mb->tag_get_data(_markTag, &edg, 1, &tagVal);
		if (tagVal)
			continue; // it must be on the boundary, the tangents must have been already computed

		int nnodes;
		const MBEntityHandle * conn2;//
		_mb->get_connectivity(edg, conn2, nnodes);
		assert(nnodes == 2);
		double coords[6]; // store the coordinates for the nodes
		MBErrorCode rval = _mb->get_coords(conn2, 2, coords);
		MBCartVect P0(coords);
		MBCartVect P1(&coords[3]);
		MBCartVect T[2];
		T[0] = P1 - P0;
		T[0].normalize();
		T[1] = T[0]; //
		_mb->tag_set_data(_tangentsTag, &edg, 1, &T[0]);// set the tangents computed at every edge

	}
	return MB_SUCCESS;
}
// start copy
//===========================================================================
//Function Name: init_edge_control_points
//
//Member Type:  PRIVATE
//Description:  compute the control points for an edge
//===========================================================================
MBErrorCode SmoothMeshEval::init_edge_control_points(MBCartVect &P0,
		MBCartVect &P3, MBCartVect &N0, MBCartVect &N3, MBCartVect &T0,
		MBCartVect &T3, MBCartVect * Pi) {
	MBCartVect Vi[4];
	Vi[0] = P0;
	Vi[3] = P3;
	MBCartVect P03(P3 - P0);
	double di = P03.length();
	double ai = N0 % N3;// this is the dot operator, the same as in cgm for CubitVector
	double ai0 = N0 % T0;
	double ai3 = N3 % T3;
	double denom = 4 - ai * ai;
	if (fabs(denom) < 1e-20) {
		return MB_FAILURE; // CUBIT_FAILURE;
	}
	double row = 6.0e0 * (2.0e0 * ai0 + ai * ai3) / denom;
	double omega = 6.0e0 * (2.0e0 * ai3 + ai * ai0) / denom;
	Vi[1] = Vi[0] + (di * (((6.0e0 * T0) - ((2.0e0 * row) * N0) + (omega * N3))
			/ 18.0e0));
	Vi[2] = Vi[3] - (di * (((6.0e0 * T3) + (row * N0) - ((2.0e0 * omega) * N3))
			/ 18.0e0));
	MBCartVect Wi[3];
	Wi[0] = Vi[1] - Vi[0];
	Wi[1] = Vi[2] - Vi[1];
	Wi[2] = Vi[3] - Vi[2];

	Pi[0] = 0.25 * Vi[0] + 0.75 * Vi[1];
	Pi[1] = 0.50 * Vi[1] + 0.50 * Vi[2];
	Pi[2] = 0.75 * Vi[2] + 0.25 * Vi[3];

	return MB_SUCCESS;
}
// end copy
#if 0
//

for (i = 0; i < myPointList.size(); i++)
{
	CubitPoint* point = myPointList.get_and_step();

	DLIList<CubitFacet*> adj_facet_list;
	point->facets(adj_facet_list);
	if (adj_facet_list.size() > 0) {
		CubitVector avg_normal(0.0e0, 0.0e0, 0.0e0);
		double totangle = 0.0e0;

		// weight the normal by the spanning angle at the point

		for (j = 0; j < adj_facet_list.size(); j++)
		{
			CubitFacet* facet = adj_facet_list.get_and_step();
			double angle = facet->angle( point );
			facet->weight( angle );
			totangle += angle;
		}
		for (j = 0; j < adj_facet_list.size(); j++)
		{
			CubitFacet* facet = adj_facet_list.get_and_step();
			CubitVector normal = facet->normal();
			normal.normalize();
			avg_normal += (facet->weight() / totangle) * normal;
		}
		avg_normal.normalize();
		point->normal(avg_normal);
		double coefd = -(point->coordinates()%avg_normal);
		point->d_coef( coefd );
	}
}
#endif
#if 0
CubitStatus FacetEvalTool::init_gradient()
{
	int i,j;

	// retrieve all faces attached to the points in point_list

	for (i = 0; i < myPointList.size(); i++)
	{
		CubitPoint* point = myPointList.get_and_step();

		DLIList<CubitFacet*> adj_facet_list;
		point->facets(adj_facet_list);
		if (adj_facet_list.size() > 0) {
			CubitVector avg_normal(0.0e0, 0.0e0, 0.0e0);
			double totangle = 0.0e0;

			// weight the normal by the spanning angle at the point

			for (j = 0; j < adj_facet_list.size(); j++)
			{
				CubitFacet* facet = adj_facet_list.get_and_step();
				double angle = facet->angle( point );
				facet->weight( angle );
				totangle += angle;
			}
			for (j = 0; j < adj_facet_list.size(); j++)
			{
				CubitFacet* facet = adj_facet_list.get_and_step();
				CubitVector normal = facet->normal();
				normal.normalize();
				avg_normal += (facet->weight() / totangle) * normal;
			}
			avg_normal.normalize();
			point->normal(avg_normal);
			double coefd = -(point->coordinates()%avg_normal);
			point->d_coef( coefd );
		}
	}
	return CUBIT_SUCCESS;
}
#endif
MBErrorCode SmoothMeshEval::find_edges_orientations(MBEntityHandle edges[3],
		const MBEntityHandle * conn3, int orient[3])// maybe we will set it?
{
	// find the edge that is adjacent to 2 vertices at a time
	for (int i = 0; i < 3; i++) {
		// edge 0 is 1-2, 1 is 3-1, 2 is 0-1
		MBEntityHandle v[2];
		v[0] = conn3[(i + 1) % 3];
		v[1] = conn3[(i + 2) % 3];
		std::vector<MBEntityHandle> adjacencies;
		// generate all edges for these two hexes
		MBErrorCode rval = _mb->get_adjacencies(v, 2, 1, false, adjacencies,
				MBInterface::INTERSECT);
		// find the edge connected to both vertices, and then see its orientation
		assert(adjacencies.size() == 1);
		const MBEntityHandle * conn2;
		int nnodes;
		rval = _mb->get_connectivity(adjacencies[0], conn2, nnodes);
		assert(rval == MB_SUCCESS);
		assert(2 == nnodes);
		edges[i] = adjacencies[0];
		// what is the story morning glory?
		if (conn2[0] == v[0] && conn2[1] == v[1])
			orient[i] = 1;
		else if (conn2[0] == v[1] && conn2[1] == v[0])
			orient[i] = -1;
		else
			return MB_FAILURE;
	}
	return MB_SUCCESS;
}
MBErrorCode SmoothMeshEval::compute_internal_control_points_on_facets(
		double min_dot) {
	// collect from each triangle the control points in order
	//
	double defControls[18] = { 0. };
	MBErrorCode rval = _mb->tag_create("CONTROLFACE", 18 * sizeof(double),
			MB_TAG_DENSE, _facetCtrlTag, &defControls);
	assert(rval == MB_SUCCESS);

	double defControls2[27] = { 0. }; // corresponding to 9 control points on edges, in order from edge 0, 1, 2 ( 1-2, 2-0, 0-1 )
	rval = _mb->tag_create("CONTROLEDGEFACE", 27 * sizeof(double),
			MB_TAG_DENSE, _facetEdgeCtrlTag, &defControls2);
	assert(rval == MB_SUCCESS);

	for (MBRange::iterator it = _triangles.begin(); it != _triangles.end(); it++) {
		MBEntityHandle tri = *it;
		// first get connectivity, and the edges
		// we need a fast method to retrieve the adjacent edges to each triangle
		const MBEntityHandle * conn3;
		int nnodes;
		rval = _mb->get_connectivity(tri, conn3, nnodes);
		assert(MB_SUCCESS == rval);
		assert(3 == nnodes);
		/*double coords[9]; // get the coordinates for the nodes
		 _mb->get_coords(conn3, 3, coords);

		 // position of nodes
		 MBCartVect vN[3];
		 int i=0;
		 for (i=0; i<3; i++)
		 vN[i] = MBCartVect(&coords[3*i]);*/
		// would it be easier to do
		MBCartVect vNode[3];// position at nodes
		rval = _mb->get_coords(conn3, 3, (double*) &vNode[0]);
		assert(MB_SUCCESS == rval);

		// get gradients (normal) at each node of triangle
		MBCartVect NN[3];
		rval = _mb->tag_get_data(_gradientTag, conn3, 3, &NN[0]);
		assert(MB_SUCCESS == rval);

		MBEntityHandle edges[3];
		int orient[3]; // + 1 or -1, if the edge is positive or negative within the face
		rval = find_edges_orientations(edges, conn3, orient);// maybe we will set it?
		assert(MB_SUCCESS == rval);
		// maybe we will store some tags with edges and their orientation with respect to
		// a triangle;
		MBCartVect P[3][5];
		MBCartVect N[6], G[6];
		// create the linear array for control points on edges, for storage (expensive !!!)
		MBCartVect CP[9];
		int index = 0;
		// maybe store a tag / entity handle for edges?
		for (int i = 0; i < 3; i++) {
			// populate P and N with the right vectors
			int i1 = (i + 1) % 3; // the first node of the edge
			int i2 = (i + 2) % 3; // the second node of the edge
			N[2 * i] = NN[i1];
			N[2 * i + 1] = NN[i2];
			P[i][0] = vNode[i1];
			rval = _mb->tag_get_data(_edgeCtrlTag, &edges[i], 1, &(P[i][1]));
			// if sense is -1, swap 1 and 3 control points
			if (orient[i] == -1) {
				MBCartVect tmp;
				tmp = P[i][1];
				P[i][1] = P[i][3];
				P[i][3] = tmp;
			}
			P[i][4] = vNode[i2];
			for (int j = 1; j < 4; j++)
				CP[index++] = P[i][j];

			// the first edge control points
		}

		//  stat = facet->get_edge_control_points( P );
		init_facet_control_points(N, P, G);
		// what do we need to store in the tag control points?
		rval = _mb->tag_set_data(_facetCtrlTag, &tri, 1, &G[0]);
		assert(MB_SUCCESS == rval);

		// store here again the 9 control points on the edges
		rval = _mb->tag_set_data(_facetEdgeCtrlTag, &tri, 1, &CP[0]);
		assert(MB_SUCCESS == rval);
		// look at what we retrieve later

		// adjust the bounding box
		int j = 0;
		for (j = 0; j < 3; j++)
			adjust_bounding_box(vNode[j]);
		// edge control points
		for (j = 0; j < 9; j++)
			adjust_bounding_box(CP[j]);
		// internal facet control points
		for (j = 0; j < 6; j++)
			adjust_bounding_box(G[j]);

	}
	return MB_SUCCESS;
}
void SmoothMeshEval::adjust_bounding_box(MBCartVect & vect) {
	// _minim, _maxim
	for (int j = 0; j < 3; j++) {
		if (_minim[j] > vect[j])
			_minim[j] = vect[j];
		if (_maxim[j] < vect[j])
			_maxim[j] = vect[j];
	}
}
//===============================================================
////Function Name: init_facet_control_points
////
////Member Type:  PRIVATE
////Descriptoin:  compute the control points for a facet
////===============================================================
MBErrorCode SmoothMeshEval::init_facet_control_points(MBCartVect N[6], // vertex normals (per edge)
		MBCartVect P[3][5], // edge control points
		MBCartVect G[6]) // return internal control points
{
	MBCartVect Di[4], Ai[3], N0, N3, Vi[4], Wi[3];
	double denom;
	double lambda[2], mu[2];

	MBErrorCode rval = MB_SUCCESS;

	for (int i = 0; i < 3; i++) {
		N0 = N[i * 2];
		N3 = N[i * 2 + 1];
		Vi[0] = P[i][0];
		Vi[1] = (P[i][1] - 0.25 * P[i][0]) / 0.75;
		Vi[2] = (P[i][3] - 0.25 * P[i][4]) / 0.75;
		Vi[3] = P[i][4];
		Wi[0] = Vi[1] - Vi[0];
		Wi[1] = Vi[2] - Vi[1];
		Wi[2] = Vi[3] - Vi[2];
		Di[0] = P[(i + 2) % 3][3] - 0.5 * (P[i][1] + P[i][0]);
		Di[3] = P[(i + 1) % 3][1] - 0.5 * (P[i][4] + P[i][3]);
		Ai[0] = (N0 * Wi[0]) / Wi[0].length();
		Ai[2] = (N3 * Wi[2]) / Wi[2].length();
		Ai[1] = Ai[0] + Ai[2];
		denom = Ai[1].length();
		Ai[1] /= denom;
		lambda[0] = (Di[0] % Wi[0]) / (Wi[0] % Wi[0]);
		lambda[1] = (Di[3] % Wi[2]) / (Wi[2] % Wi[2]);
		mu[0] = (Di[0] % Ai[0]);
		mu[1] = (Di[3] % Ai[2]);
		G[i * 2] = 0.5 * (P[i][1] + P[i][2]) + 0.66666666666666 * lambda[0]
				* Wi[1] + 0.33333333333333 * lambda[1] * Wi[0]
				+ 0.66666666666666 * mu[0] * Ai[1] + 0.33333333333333 * mu[1]
				* Ai[0];
		G[i * 2 + 1] = 0.5 * (P[i][2] + P[i][3]) + 0.33333333333333 * lambda[0]
				* Wi[2] + 0.66666666666666 * lambda[1] * Wi[1]
				+ 0.33333333333333 * mu[0] * Ai[2] + 0.66666666666666 * mu[1]
				* Ai[1];
	}
	return rval;
}

//===========================================================================
//Function Name: evaluate_single
//
//Member Type:  PUBLIC
//Description:  evaluate edge not associated with a facet (this is used
// by camal edge mesher!!!)
//Note:         t is a value from 0 to 1, for us
//===========================================================================
MBErrorCode SmoothMeshEval::evaluate_smooth_edge(MBEntityHandle eh, double &tt,
		MBCartVect & outv) {
	MBCartVect P[2]; // P0 and P1
	MBCartVect controlPoints[3]; // edge control points
	double t4, t3, t2, one_minus_t, one_minus_t2, one_minus_t3, one_minus_t4;

	// project the position to the linear edge
	// t is from 0 to 1 only!!
	//double tt = (t + 1) * 0.5;
	if (tt <= 0.0)
		tt = 0.0;
	if (tt >= 1.0)
		tt = 1.0;

	int nnodes;
	const MBEntityHandle * conn2;
	MBErrorCode rval = _mb->get_connectivity(eh, conn2, nnodes);
	assert(rval == MB_SUCCESS);

	rval = _mb->get_coords(conn2, 2, (double*) &P[0]);
	assert(rval == MB_SUCCESS);

	rval = _mb->tag_get_data(_edgeCtrlTag, &eh, 1, (double*) &controlPoints[0]);
	assert(rval == MB_SUCCESS);

	t2 = tt * tt;
	t3 = t2 * tt;
	t4 = t3 * tt;
	one_minus_t = 1. - tt;
	one_minus_t2 = one_minus_t * one_minus_t;
	one_minus_t3 = one_minus_t2 * one_minus_t;
	one_minus_t4 = one_minus_t3 * one_minus_t;

	outv = one_minus_t4 * P[0] + 4. * one_minus_t3 * tt * controlPoints[0] + 6.
			* one_minus_t2 * t2 * controlPoints[1] + 4. * one_minus_t * t3
			* controlPoints[2] + t4 * P[1];

	return MB_SUCCESS;
}

//this is the
double SmoothMeshEval::get_length_loop(int loopIndex) {
	assert(0 <= loopIndex && loopIndex < _loopLengths.size());
	return _loopLengths[loopIndex];
}

MBErrorCode SmoothMeshEval::evaluate_loop_at_u(int loopIndex, double u,
		MBCartVect & position) {
	// first
	assert(loopIndex < _loopLengths.size());
	assert(loopIndex < _loopLengths.size());
	int startCurrentLoop = (loopIndex == 0) ? 0 : _loopEnds[loopIndex - 1];
	int endCurrentLoop = _loopEnds[loopIndex] - 1;
	if (u > 1.)
		u = 1.;//
	if (u < 0.)
		u = 0.;
	// _fractions are increasing, so find the
	double * ptr = std::lower_bound(&_fractions[startCurrentLoop],
			&_fractions[endCurrentLoop], u);
	int index = ptr - &_fractions[startCurrentLoop];
	double nextFraction = _fractions[startCurrentLoop + index];
	double prevFraction = 0;
	if (index > 0) {
		prevFraction = _fractions[startCurrentLoop + index - 1];
	}
	double t = (u - prevFraction) / (nextFraction - prevFraction);

	MBEntityHandle edge = _loops[startCurrentLoop + index];
	if (2 == this->_senses[startCurrentLoop + index])
		t = 1 - t; // reverse

	MBErrorCode rval = evaluate_smooth_edge(edge, t, position);
	assert(MB_SUCCESS == rval);

	return MB_SUCCESS;

}

MBErrorCode SmoothMeshEval::eval_bezier_patch(MBEntityHandle tri,
		MBCartVect &areacoord, MBCartVect &pt) {
	//
	// interpolate internal control points

	MBCartVect gctrl_pts[6];
	// get the control points  facet->get_control_points( gctrl_pts );
	//init_facet_control_points( N, P, G) ;
	// what do we need to store in the tag control points?
	MBErrorCode rval = _mb->tag_get_data(_facetCtrlTag, &tri, 1, &gctrl_pts[0]);// get all 6 control points
	assert(MB_SUCCESS == rval);
	const MBEntityHandle * conn3;
	int nnodes;
	rval = _mb->get_connectivity(tri, conn3, nnodes);
	assert(MB_SUCCESS == rval);

	MBCartVect vN[3];
	_mb->get_coords(conn3, 3, (double*) &vN[0]); // fill the coordinates of the vertices

	if (fabs(areacoord[1] + areacoord[2]) < 1.0e-6) {
		pt = vN[0];
		return MB_SUCCESS;
	}
	if (fabs(areacoord[0] + areacoord[2]) < 1.0e-6) {
		pt = vN[0];
		return MB_SUCCESS;
	}
	if (fabs(areacoord[0] + areacoord[1]) < 1.0e-6) {
		pt = vN[0];
		return MB_SUCCESS;
	}

	MBCartVect P_facet[3];
	//2,1,1
	P_facet[0] = (1.0e0 / (areacoord[1] + areacoord[2])) * (areacoord[1]
			* gctrl_pts[3] + areacoord[2] * gctrl_pts[4]);
	//1,2,1
	P_facet[1] = (1.0e0 / (areacoord[0] + areacoord[2])) * (areacoord[0]
			* gctrl_pts[0] + areacoord[2] * gctrl_pts[5]);
	//1,1,2
	P_facet[2] = (1.0e0 / (areacoord[0] + areacoord[1])) * (areacoord[0]
			* gctrl_pts[1] + areacoord[1] * gctrl_pts[2]);

	// sum the contribution from each of the control points

	pt = MBCartVect(0.);// set all to 0, we start adding / accumulating different parts
	// first edge is from node 0 to 1, index 2 in

	// retrieve the points, in order, and the control points on edges

	// store here again the 9 control points on the edges
	MBCartVect CP[9];
	rval = _mb->tag_get_data(_facetEdgeCtrlTag, &tri, 1, &CP[0]);
	assert(MB_SUCCESS == rval);

	//CubitFacetEdge *edge;
	//edge = facet->edge(2);! start with edge 2, from 0-1
	int k = 0;
	MBCartVect ctrl_pts[5];
	//edge->control_points(facet, ctrl_pts);
	ctrl_pts[0] = vN[0]; //
	for (k = 1; k < 4; k++)
		ctrl_pts[k] = CP[k + 5]; // for edge index 2
	ctrl_pts[4] = vN[1]; //

	//i=4; j=0; k=0;
	double B = quart(areacoord[0]);
	pt += B * ctrl_pts[0];

	//i=3; j=1; k=0;
	B = 4.0 * cube(areacoord[0]) * areacoord[1];
	pt += B * ctrl_pts[1];

	//i=2; j=2; k=0;
	B = 6.0 * sqr(areacoord[0]) * sqr(areacoord[1]);
	pt += B * ctrl_pts[2];

	//i=1; j=3; k=0;
	B = 4.0 * areacoord[0] * cube(areacoord[1]);
	pt += B * ctrl_pts[3];

	//edge = facet->edge(0);
	//edge->control_points(facet, ctrl_pts);
	// edge index 0, from 1 to 2
	ctrl_pts[0] = vN[1]; //
	for (k = 1; k < 4; k++)
		ctrl_pts[k] = CP[k - 1]; // for edge index 0
	ctrl_pts[4] = vN[2]; //

	//i=0; j=4; k=0;
	B = quart(areacoord[1]);
	pt += B * ctrl_pts[0];

	//i=0; j=3; k=1;
	B = 4.0 * cube(areacoord[1]) * areacoord[2];
	pt += B * ctrl_pts[1];

	//i=0; j=2; k=2;
	B = 6.0 * sqr(areacoord[1]) * sqr(areacoord[2]);
	pt += B * ctrl_pts[2];

	//i=0; j=1; k=3;
	B = 4.0 * areacoord[1] * cube(areacoord[2]);
	pt += B * ctrl_pts[3];

	//edge = facet->edge(1);
	//edge->control_points(facet, ctrl_pts);
	// edge index 1, from 2 to 0
	ctrl_pts[0] = vN[2]; //
	for (k = 1; k < 4; k++)
		ctrl_pts[k] = CP[k + 2]; // for edge index 0
	ctrl_pts[4] = vN[0]; //

	//i=0; j=0; k=4;
	B = quart(areacoord[2]);
	pt += B * ctrl_pts[0];

	//i=1; j=0; k=3;
	B = 4.0 * areacoord[0] * cube(areacoord[2]);
	pt += B * ctrl_pts[1];

	//i=2; j=0; k=2;
	B = 6.0 * sqr(areacoord[0]) * sqr(areacoord[2]);
	pt += B * ctrl_pts[2];

	//i=3; j=0; k=1;
	B = 4.0 * cube(areacoord[0]) * areacoord[2];
	pt += B * ctrl_pts[3];

	//i=2; j=1; k=1;
	B = 12.0 * sqr(areacoord[0]) * areacoord[1] * areacoord[2];
	pt += B * P_facet[0];

	//i=1; j=2; k=1;
	B = 12.0 * areacoord[0] * sqr(areacoord[1]) * areacoord[2];
	pt += B * P_facet[1];

	//i=1; j=1; k=2;
	B = 12.0 * areacoord[0] * areacoord[1] * sqr(areacoord[2]);
	pt += B * P_facet[2];

	return MB_SUCCESS;

}

//===========================================================================
//Function Name: project_to_facet_plane
//
//Member Type:  PUBLIC
//Descriptoin:  Project a point to the plane of a facet
//===========================================================================
void SmoothMeshEval::project_to_facet_plane(MBEntityHandle tri, MBCartVect &pt,
		MBCartVect &point_on_plane, double &dist_to_plane) {
	double plane[4];

	MBErrorCode rval = _mb->tag_get_data(_planeTag, &tri, 1, plane);
	assert(rval == MB_SUCCESS);
	// _planeTag
	MBCartVect normal(&plane[0]);// just first 3 components are used
	// normal.normalize();
	/* CubitPoint *facPoint = facet->point(0);
	 double coefd = -(normal.x() * facPoint->x() +
	 normal.y() * facPoint->y() +
	 normal.z() * facPoint->z());
	 double dist = normal.x() * pt.x() +
	 normal.y() * pt.y() +
	 normal.z() * pt.z() + coefd;*/
	double dist = normal % pt + plane[3]; // coeff d is saved!!!
	dist_to_plane = fabs(dist);
	point_on_plane = pt - dist * normal;
	/* point_on_plane.set( pt.x() - normal.x() * dist,
	 pt.y() - normal.y() * dist,
	 pt.z() - normal.z() * dist );*/
	return;
}

//===========================================================================
//Function Name: facet_area_coordinate
//
//Member Type:  PUBLIC
//Descriptoin:  Determine the area coordinates of a point on the plane
//              of a facet
//===========================================================================
void SmoothMeshEval::facet_area_coordinate(MBEntityHandle facet,
		MBCartVect & pt_on_plane, MBCartVect & areacoord) {

	const MBEntityHandle * conn3;
	int nnodes;
	MBErrorCode rval = _mb->get_connectivity(facet, conn3, nnodes);
	assert(MB_SUCCESS == rval);
	//double coords[9]; // store the coordinates for the nodes
	//_mb->get_coords(conn3, 3, coords);
	MBCartVect p[3];
	rval = _mb->get_coords(conn3, 3, (double*) &p[0]);
	assert(MB_SUCCESS == rval);
	double plane[4];

	rval = _mb->tag_get_data(_planeTag, &facet, 1, plane);
	assert(rval == MB_SUCCESS);
	MBCartVect normal(&plane[0]);// just first 3 components are used

	double area2;
	/*CubitVector normal = facet->normal();
	 CubitPoint *pt0 = facet->point(0);
	 CubitPoint *pt1 = facet->point(1);
	 CubitPoint *pt2 = facet->point(2);*/
	double tol = GEOMETRY_RESABS * 1.e-5;// 1.e-11;
	/*CubitVector v1( pt1->x() - pt0->x(),
	 pt1->y() - pt0->y(),
	 pt1->z() - pt0->z());//(*p1-*p0);
	 CubitVector v2( pt2->x() - pt0->x(),
	 pt2->y() - pt0->y(),
	 pt2->z() - pt0->z());// = (*p2-*p0);*/
	MBCartVect v1(p[1] - p[0]);
	MBCartVect v2(p[2] - p[0]);

	area2 = (v1 * v2).length_squared();// the same for MBCartVect
	if (area2 < 100 * tol) {
		tol = .01 * area2;
	}
	MBCartVect absnorm(fabs(normal[0]), fabs(normal[1]), fabs(normal[2]));

	// project to the closest coordinate plane so we only have to do this in 2D

	if (absnorm[0] >= absnorm[1] && absnorm[0] >= absnorm[2]) {
		/*area2 = determ3(pt0->y(), pt0->z(),
		 pt1->y(), pt1->z(),
		 pt2->y(), pt2->z());*/
		area2 = determ3( p[0][1], p[0][2],
				p[1][1], p[1][2],
				p[2][1], p[2][2]);
		if (fabs(area2) < tol) {
			areacoord = MBCartVect(-CUBIT_DBL_MAX);// .set( -CUBIT_DBL_MAX, -CUBIT_DBL_MAX, -CUBIT_DBL_MAX );
		} else if (within_tolerance(p[0], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(1., 0., 0.);//.set( 1.0, 0.0, 0.0 );
		} else if (within_tolerance(p[1], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(0., 1., 0.);//.set( 0.0, 1.0, 0.0 );
		} else if (within_tolerance(p[2], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(0., 0., 1.);//.set( 0.0, 0.0, 1.0 );
		} else {
			/*areacoord.x(
			 determ3( pt_on_plane.y(), pt_on_plane.z(),
			 pt1->y(), pt1->z(), pt2->y(), pt2->z() ) / area2 );*/
			areacoord[0] = determ3(pt_on_plane[1], pt_on_plane[2],
					p[1][1], p[1][2], p[2][1], p[2][2] ) / area2;
			/*areacoord.y(
			 determ3( pt0->y(), pt0->z(),
			 pt_on_plane.y(), pt_on_plane.z(),
			 pt2->y(), pt2->z() ) / area2 );*/
			areacoord[1] = determ3( p[0][1], p[0][2],
					pt_on_plane[1], pt_on_plane[2],
					p[2][1], p[2][2]) / area2;

			/* areacoord.z(
			 determ3( pt0->y(), pt0->z(), pt1->y(), pt1->z(),
			 pt_on_plane.y(), pt_on_plane.z() ) / area2 );*/
			areacoord[2] = determ3( p[0][1], p[0][2], p[1][1], p[1][2],
					pt_on_plane[1], pt_on_plane[2]) / area2;
		}
	} else if (absnorm[1] >= absnorm[0] && absnorm[1] >= absnorm[2]) {
		/*area2 = determ3(pt0->x(), pt0->z(),
		 pt1->x(), pt1->z(),
		 pt2->x(), pt2->z());*/
		area2 = determ3(p[0][0], p[0][2],
				p[1][0], p[1][2],
				p[2][0], p[2][2]);
		if (fabs(area2) < tol) {
			areacoord = MBCartVect(-CUBIT_DBL_MAX);//.set( -CUBIT_DBL_MAX, -CUBIT_DBL_MAX, -CUBIT_DBL_MAX );
		} else if (within_tolerance(p[0], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(1., 0., 0.);//.set( 1.0, 0.0, 0.0 );
		} else if (within_tolerance(p[1], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(0., 1., 0.);//.set( 0.0, 1.0, 0.0 );
		} else if (within_tolerance(p[2], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(0., 0., 1.);//.set( 0.0, 0.0, 1.0 );
		} else {
			/*areacoord.x(
			 determ3( pt_on_plane.x(), pt_on_plane.z(),
			 pt1->x(), pt1->z(), pt2->x(), pt2->z() ) / area2 );
			 areacoord.y(
			 determ3( pt0->x(), pt0->z(),
			 pt_on_plane.x(), pt_on_plane.z(),
			 pt2->x(), pt2->z() ) / area2 );
			 areacoord.z(
			 determ3( pt0->x(), pt0->z(), pt1->x(), pt1->z(),
			 pt_on_plane.x(), pt_on_plane.z() ) / area2 );*/
			areacoord[0] = determ3(pt_on_plane[0], pt_on_plane[2],
					p[1][0], p[1][2], p[2][0], p[2][2] ) / area2;

			areacoord[1] = determ3( p[0][0], p[0][2],
					pt_on_plane[0], pt_on_plane[2],
					p[2][0], p[2][2]) / area2;

			areacoord[2] = determ3( p[0][0], p[0][2], p[1][0], p[1][2],
					pt_on_plane[0], pt_on_plane[2]) / area2;

		}
	} else {
		/*area2 = determ3(pt0->x(), pt0->y(),
		 pt1->x(), pt1->y(),
		 pt2->x(), pt2->y());*/
		area2 = determ3( p[0][0], p[0][1],
				p[1][0], p[1][1],
				p[2][0], p[2][1]);
		if (fabs(area2) < tol) {
			areacoord = MBCartVect(-CUBIT_DBL_MAX);//.set( -CUBIT_DBL_MAX, -CUBIT_DBL_MAX, -CUBIT_DBL_MAX );
		} else if (within_tolerance(p[0], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(1., 0., 0.);//.set( 1.0, 0.0, 0.0 );
		} else if (within_tolerance(p[1], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(0., 1., 0.);//.set( 0.0, 1.0, 0.0 );
		} else if (within_tolerance(p[2], pt_on_plane, GEOMETRY_RESABS)) {
			areacoord = MBCartVect(0., 0., 1.);//.set( 0.0, 0.0, 1.0 );
		} else {
			/*areacoord.x(
			 determ3( pt_on_plane.x(), pt_on_plane.y(),
			 pt1->x(), pt1->y(), pt2->x(), pt2->y() ) / area2 );
			 areacoord.y(
			 determ3( pt0->x(), pt0->y(),
			 pt_on_plane.x(), pt_on_plane.y(),
			 pt2->x(), pt2->y() ) / area2 );
			 areacoord.z(
			 determ3( pt0->x(), pt0->y(), pt1->x(), pt1->y(),
			 pt_on_plane.x(), pt_on_plane.y() ) / area2 );*/
			areacoord[0] = determ3(pt_on_plane[0], pt_on_plane[1],
					p[1][0], p[1][1], p[2][0], p[2][1] ) / area2;

			areacoord[1] = determ3( p[0][0], p[0][1],
					pt_on_plane[0], pt_on_plane[1],
					p[2][0], p[2][1]) / area2;

			areacoord[2] = determ3( p[0][0], p[0][1], p[1][0], p[1][1],
					pt_on_plane[0], pt_on_plane[1]) / area2;
		}
	}
}

MBErrorCode SmoothMeshEval::project_to_facets_main(MBCartVect &this_point,
		bool trim, bool & outside, MBCartVect * closest_point_ptr,
		MBCartVect * normal_ptr) {

	// if there are a lot of facets on this surface - use the OBB search first
	// to narrow the selection

	double tolerance = 1.e-5;
	std::vector<MBEntityHandle> facets_out;
	// we will start with a list of facets anyway, the best among them wins
	MBErrorCode rval = _my_geomTopoTool->obb_tree()->closest_to_location(
			(double*) &this_point, _obb_root, tolerance, facets_out);

	int interpOrder = 4;
	double compareTol = 1.e-5;
	MBEntityHandle lastFacet = facets_out.front();
	rval = project_to_facets(facets_out, lastFacet, interpOrder, compareTol,
			this_point, trim, outside, closest_point_ptr, normal_ptr);

	return rval;
}
MBErrorCode SmoothMeshEval::project_to_facets(
		std::vector<MBEntityHandle> & facet_list, MBEntityHandle & lastFacet,
		int interpOrder, double compareTol, MBCartVect &this_point, bool trim,
		bool & outside, MBCartVect *closest_point_ptr, MBCartVect * normal_ptr) {

	int ncheck, ii, nincr = 0;
	static int calls = 0;
	static int nncheck = 0;
	static int ntol = 0;
	static int mydebug = 0;
	bool outside_facet, best_outside_facet = true;
	MBCartVect close_point, best_point, best_areacoord;
	MBEntityHandle best_facet= 0L;// no best facet found yet
	MBEntityHandle facet;
	assert(facet_list.size() > 0);

	double big_dist = compareTol * 1.0e3;

	int nevald = 0;
  double tol = compareTol * 10;
  const double atol = 0.001;
  double mindist = 1.e20;
  bool done = false; // maybe not use this

	// from the list of close facets, determine the closest point
	for (int i=0; i<facet_list.size(); i++)
	{
		facet = facet_list[i];
		MBCartVect pt_on_plane;
		double dist_to_plane;
		project_to_facet_plane(facet, this_point, pt_on_plane,
				dist_to_plane);

		MBCartVect areacoord;
		//MBCartVect close_point;
		facet_area_coordinate(facet, pt_on_plane, areacoord);
		if (interpOrder != 0) {

			// modify the areacoord - project to the bezier patch- snaps to the
			// edge of the patch if necessary



			if (project_to_facet(facet, this_point, areacoord, close_point,
					outside_facet, compareTol) != MB_SUCCESS) {
				return MB_FAILURE;
			}
			//if (closest_point_ptr)
				//*closest_point_ptr = close_point;
		}
		// keep track of the minimum distance

		double dist = (close_point-this_point).length();//close_point.distance_between(this_point);
		if ( (best_outside_facet == outside_facet && dist < mindist)
				|| (best_outside_facet && !outside_facet && (dist
						< big_dist || best_facet==0L/*!best_facet*/))) {
			mindist = dist;
			best_point = close_point;
			best_facet = facet;
			best_areacoord = areacoord;
			best_outside_facet = outside_facet;

			if (dist < compareTol) {
				done = true;
				break;
			}
			big_dist = 10.0 * mindist;
		}
		//facet->marked(1);
		//used_facet_list.append(facet);

	}

	if (normal_ptr) {
		MBCartVect normal;
		if (eval_bezier_patch_normal(best_facet, best_areacoord, normal)
				!= MB_SUCCESS) {
			return MB_FAILURE;
		}
		*normal_ptr = normal;
	}

	if (closest_point_ptr) {
		*closest_point_ptr = best_point;
	}

	outside = best_outside_facet;
	lastFacet = best_facet;

#if 0
	// big loop here;
	// start copy
	int ncheck, ii, nincr = 0;
	static int calls = 0;
	static int nncheck = 0;
	static int ntol = 0;
	static int mydebug = 0;
	CubitBoolean outside_facet, best_outside_facet;
	CubitVector close_point, best_point, best_areacoord;
	CubitFacet *best_facet = NULL;
	CubitFacet *facet;
	assert(facet_list.size() > 0);
	double big_dist = compare_tol * 1.0e3;

	// set the first facet to be checked as the last one located

	if (last_facet) {
		if (!facet_list.move_to(last_facet)) {
			facet_list.reset();
		}
	} else {
		facet_list.reset();
	}

	// so we don't evaluate a facet more than once - mark the facets
	// as we evaluate them.  Put the evaluated facets on a used_facet_list
	// so we clear the marks off when we are done.  Note: this assumes
	// theat marks are initially cleared.

	DLIList<CubitFacet *> used_facet_list;
	for (ii = 0; ii < facet_list.size(); ii++)
		facet_list.get_and_step()->marked(0);

	int nfacets = facet_list.size();
	int nevald = 0;
	double tol = compare_tol * 10;
	const double atol = 0.001;
	double mindist = CUBIT_DBL_MAX;
	CubitBoolean eval_all = CUBIT_FALSE;
	CubitBoolean done = CUBIT_FALSE;
	best_outside_facet = CUBIT_TRUE;

	while (!done) {

		// define a bounding box around the point

		double ptmin_x = this_point.x() - tol;
		double ptmin_y = this_point.y() - tol;
		double ptmin_z = this_point.z() - tol;
		double ptmax_x = this_point.x() + tol;
		double ptmax_y = this_point.y() + tol;
		double ptmax_z = this_point.z() + tol;

		ncheck = 0;
		for (ii = facet_list.size(); ii > 0 && !done; ii--) {
			facet = facet_list.get_and_step();
			if (facet->marked())
				continue;

			// Try to trivially reject this facet with a bounding box test
			// (Does the bounding box of the facet intersect with the
			// bounding box of the point)

			if (!eval_all) {
				const CubitBox &bbox = facet->bounding_box();
				if (ptmax_x < bbox.min_x() || ptmin_x > bbox.max_x()) {
					continue;
				}
				if (ptmax_y < bbox.min_y() || ptmin_y > bbox.max_y()) {
					continue;
				}
				if (ptmax_z < bbox.min_z() || ptmin_z > bbox.max_z()) {
					continue;
				}
			}

			// Only facets that pass the bounding box test will get past here!

			// Project point to plane of the facet and determine its area coordinates

			ncheck++;
			nevald++;
			CubitVector pt_on_plane;
			double dist_to_plane;
			project_to_facet_plane(facet, this_point, pt_on_plane,
					dist_to_plane);

			CubitVector areacoord;
			facet_area_coordinate(facet, pt_on_plane, areacoord);
			if (interp_order != 0) {

				// modify the areacoord - project to the bezier patch- snaps to the
				// edge of the patch if necessary

				if (project_to_facet(facet, this_point, areacoord, close_point,
						outside_facet, compare_tol) != CUBIT_SUCCESS) {
					return CUBIT_FAILURE;
				}
			} else {

				// If sign of areacoords are all positive then its inside the triangle
				// and we are done - go interpolate the point. (use an absolute
				// tolerance since the coordinates arenormalized)
				if (areacoord.x() > -atol && areacoord.y() > -atol
						&& areacoord.z() > -atol) {
					if (dist_to_plane < compare_tol) {
						close_point = this_point;
					} else {
						close_point = pt_on_plane;
					}
					outside_facet = CUBIT_FALSE;
				}

				// otherwise find the closest vertex or edge to the projected point

				else if (areacoord.x() < atol) {
					outside_facet = CUBIT_TRUE;
					if (areacoord.y() < atol) {
						if (eval_point(facet, 2, close_point) != CUBIT_SUCCESS) {
							return CUBIT_FAILURE;
						}
					} else if (areacoord.z() < atol) {
						if (eval_point(facet, 1, close_point) != CUBIT_SUCCESS) {
							return CUBIT_FAILURE;
						}
					} else {
						if (project_to_facetedge(facet, 1, 2, this_point,
								pt_on_plane, close_point, outside_facet, trim)
								!= CUBIT_SUCCESS) {
							return CUBIT_FAILURE;
						}
					}
				} else if (areacoord.y() < atol) {
					outside_facet = CUBIT_TRUE;
					if (areacoord.z() < atol) {
						if (eval_point(facet, 0, close_point) != CUBIT_SUCCESS) {
							return CUBIT_FAILURE;
						}
					} else {
						if (project_to_facetedge(facet, 2, 0, this_point,
								pt_on_plane, close_point, outside_facet, trim)
								!= CUBIT_SUCCESS) {
							return CUBIT_FAILURE;
						}
					}
				} else {
					outside_facet = CUBIT_TRUE;
					if (project_to_facetedge(facet, 0, 1, this_point,
							pt_on_plane, close_point, outside_facet, trim)
							!= CUBIT_SUCCESS) {
						return CUBIT_FAILURE;
					}
				}
			}

			// keep track of the minimum distance

			double dist = close_point.distance_between(this_point);
			if ((best_outside_facet == outside_facet && dist < mindist)
					|| (best_outside_facet && !outside_facet && (dist
							< big_dist || !best_facet))) {
				mindist = dist;
				best_point = close_point;
				best_facet = facet;
				best_areacoord = areacoord;
				best_outside_facet = outside_facet;

				if (dist < compare_tol) {
					done = CUBIT_TRUE;
				}
				big_dist = 10.0 * mindist;
			}
			facet->marked(1);
			used_facet_list.append(facet);
		}

		// We are done if we found at least one triangle.  Otherwise
		// increase the tolerance and try again

		if (!done) {
			if (nevald == nfacets) {
				done = CUBIT_TRUE;
			} else {
				nincr++;
				if (ncheck > 0) {
					if (best_outside_facet) {
						if (nincr < 10) {
							tol *= 2.0;
							ntol++;
						} else
						// getting here means that the compare_tol probably is too small
						// just try all the remaining facets
						{
							eval_all = CUBIT_TRUE;
						}
					} else {
						done = CUBIT_TRUE;
					}
				} else {
					if (nincr >= 10) {
						eval_all = CUBIT_TRUE;
					} else {
						tol *= 2.0e0;
						ntol++;
					}
				}
			}
		}
	} // while(!done)
	if (best_facet == NULL) {
		PRINT_ERROR("Unable to determine facet correctly.\n");
		return CUBIT_FAILURE;
	}
	// make sure we actually got something
	assert(best_facet != NULL);

	// if the closest point is outside of a facet, then evaluate the point
	// on the facet using its area coordinates (otherwise it would be
	// trimmed to an edge or point)


	if (!trim && best_outside_facet && interp_order != 4) {
		if (project_to_facet(best_facet, this_point, best_areacoord,
				best_point, best_outside_facet, compare_tol) != CUBIT_SUCCESS) {
			return CUBIT_FAILURE;
		}

		// see if its really outside (it could just be on an edge where the
		// curvature is convex)

		best_outside_facet = is_outside(best_facet, best_areacoord);
	}

	// evaluate the normal if required

	if (normal_ptr) {
		CubitVector normal;
		if (eval_facet_normal(best_facet, best_areacoord, normal)
				!= CUBIT_SUCCESS) {
			return CUBIT_FAILURE;
		}
		*normal_ptr = normal;
	}

	if (closest_point_ptr) {
		*closest_point_ptr = best_point;
	}

	*outside = best_outside_facet;
	last_facet = best_facet;

	// clear the marks from the used facets

	for (ii = 0; ii < used_facet_list.size(); ii++) {
		facet = used_facet_list.get_and_step();
		facet->marked(0);
	}

	if (mydebug) {
		nncheck += ncheck;
		calls++;
		if (calls % 100 == 0) {
			char message[100];
			sprintf(message, "calls = %d, ckecks = %d, ntol = %d\n", calls,
					nncheck, ntol);
			PRINT_INFO(message);
		}
	}

	return CUBIT_SUCCESS;
#endif
	return MB_SUCCESS;
	//end copy
}


//===========================================================================
//Function Name: project_to_patch
//
//Member Type:  PUBLIC
//Descriptoin:  Project a point to a bezier patch. Pass in the areacoord
//              of the point projected to the linear facet.  Function
//              assumes that the point is contained within the patch -
//              if not, it will project to one of its edges.
//===========================================================================
MBErrorCode SmoothMeshEval::project_to_patch(
  MBEntityHandle facet,     // (IN) the facet where the patch is defined
  MBCartVect &ac,       // (IN) area coordinate initial guess (from linear facet)
  MBCartVect &pt,       // (IN) point we are projecting to patch
  MBCartVect &eval_pt,  // (OUT) The projected point
  MBCartVect *eval_norm, // (OUT) normal at evaluated point
  bool &outside, // (OUT) the closest point on patch to pt is on an edge
  double compare_tol,    // (IN) comparison tolerance
  int edge_id )          // (IN) only used if this is to be projected to one
                         //      of the edges.  Otherwise, should be -1
{
  MBErrorCode status = MB_SUCCESS;

  // see if we are at a vertex

#define INCR 0.01
  const double tol = compare_tol;

  if (is_at_vertex( facet, pt, ac, compare_tol, eval_pt, eval_norm ))
  {
    outside = false;
    return MB_SUCCESS;
  }

  // check if the start ac is inside the patch -if not, then move it there

  int nout = 0;
  const double atol = 0.001;
  if(move_ac_inside( ac, atol ))
    nout++;

  int diverge = 0;
  int iter = 0;
  MBCartVect  newpt;
  eval_bezier_patch(facet, ac, newpt);
  MBCartVect move = pt - newpt;
  double lastdist = move.length();
  double bestdist = lastdist;
  MBCartVect bestac = ac;
  MBCartVect bestpt = newpt;
  MBCartVect bestnorm;

  // If we are already close enough, then return now

  if (lastdist <= tol && !eval_norm && nout == 0) {
    eval_pt = pt;
    outside = false;
    return status;
  }

  double ratio, mag, umove, vmove, det, distnew, movedist;
  MBCartVect lastpt = newpt;
  MBCartVect lastac = ac;
  MBCartVect norm;
  MBCartVect xpt, ypt, zpt, xac, yac, zac, xvec, yvec, zvec;
  MBCartVect du, dv, newac;
  bool done = false;
  while (!done) {

    // We will be locating the projected point within the u,v,w coordinate
    // system of the triangular bezier patch.  Since u+v+w=1, the components
    // are not linearly independent.  We will choose only two of the
    // coordinates to use and compute the third.

    int system;
    if (lastac[0] >= lastac[1] && lastac[0] >= lastac[2]) {
      system = 0;
    }
    else if (lastac[1] >= lastac[2]) {
      system = 1;
    }
    else {
      system = 2;
    }

    // compute the surface derivatives with respect to each
    // of the barycentric coordinates


    if (system == 1 || system == 2) {
      xac[0] = lastac[0]+INCR; // xac.x( lastac.x() + INCR );
      if (lastac[1] + lastac[2] == 0.0)
        return MB_FAILURE;
      ratio = lastac[2] / (lastac[1] + lastac[2]);//ratio = lastac.z() / (lastac.y() + lastac.z());
      xac[1] =  (1.0 - xac[0]) * (1.0 - ratio) ;//xac.y( (1.0 - xac.x()) * (1.0 - ratio) );
      xac[2] =  1.0 - xac[0] - xac[1]; // xac.z( 1.0 - xac.x() - xac.y() );
      eval_bezier_patch(facet, xac, xpt);
      xvec = xpt - lastpt;
      xvec /= INCR;
    }
    if (system == 0 || system == 2) {
       yac[1] = ( lastac[1] + INCR );//yac.y( lastac.y() + INCR );
       if (lastac[0] + lastac[2] == 0.0)//if (lastac.x() + lastac.z() == 0.0)
        return MB_FAILURE;
       ratio = lastac[2] / (lastac[0] + lastac[2]);//ratio = lastac.z() / (lastac.x() + lastac.z());
       yac[0] = ( (1.0 - yac[1]) * (1.0 - ratio) );//yac.x( (1.0 - yac.y()) * (1.0 - ratio) );
       yac[2] = ( 1.0 - yac[0] - yac[1] );//yac.z( 1.0 - yac.x() - yac.y() );
       eval_bezier_patch(facet, yac, ypt);
       yvec = ypt - lastpt;
       yvec /= INCR;
    }
    if (system == 0 || system == 1) {
    	zac[2] = ( lastac[2] + INCR );//zac.z( lastac.z() + INCR );
    	if (lastac[0] + lastac[1] == 0.0)//if (lastac.x() + lastac.y() == 0.0)
          return MB_FAILURE;
    	ratio = lastac[1] / (lastac[0] + lastac[1]);//ratio = lastac.y() / (lastac.x() + lastac.y());
    	zac[0] = ( (1.0 - zac[2]) * (1.0 - ratio) );//zac.x( (1.0 - zac.z()) * (1.0 - ratio) );
    	zac[1] = ( 1.0 - zac[0] - zac[2] );//zac.y( 1.0 - zac.x() - zac.z() );
      eval_bezier_patch(facet, zac, zpt);
      zvec = zpt - lastpt;
      zvec /= INCR;
    }

    // compute the surface normal

     switch (system) {
     case 0:
       du = yvec;
       dv = zvec;
       break;
     case 1:
       du = zvec;
       dv = xvec;
       break;
     case 2:
       du = xvec;
       dv = yvec;
       break;
     }
     norm = du * dv;
     mag = norm.length();
     if (mag < DBL_EPSILON) {
       return MB_FAILURE;
       // do something else here (it is likely a flat triangle -
       // so try evaluating just an edge of the bezier patch)
     }
     norm /= mag;
     if (iter == 0)
       bestnorm = norm;

     // project the move vector to the tangent plane

     move = (norm * move) * norm;

     // compute an equivalent u-v-w vector

     MBCartVect absnorm( fabs(norm[0]), fabs(norm[1]), fabs(norm[2]) );
     if (absnorm[2] >= absnorm[1] && absnorm[2] >= absnorm[0]) {
       det = du[0] * dv[1] - dv[0] * du[1];
       if (fabs(det) <= DBL_EPSILON) {
         return MB_FAILURE;  // do something else here
       }
       umove = (move[0] * dv[1] - dv[0] * move[1]) / det;
       vmove = (du[0] * move[1] - move[0] * du[1]) / det;
     }
     else if (absnorm[1] >= absnorm[2] && absnorm[1] >= absnorm[0]) {
       det = du[0] * dv[2] - dv[0] * du[2];
       if (fabs(det) <= DBL_EPSILON) {
         return MB_FAILURE;
       }
       umove = (move[0] * dv[2] - dv[0] * move[2]) / det;
       vmove = (du[0] * move[2] - move[0] * du[2]) / det;
     }
     else {
       det = du[1] * dv[2] - dv[1] * du[2];
       if (fabs(det) <= DBL_EPSILON) {
         return MB_FAILURE;
       }
       umove = (move[1] * dv[2] - dv[1] * move[2]) / det;
       vmove = (du[1] * move[2] - move[1] * du[2]) / det;
     }

    /* === compute the new u-v coords and evaluate surface at new location */

    switch (system) {
    case 0:
    	newac[1] = ( lastac[1] + umove );//newac.y( lastac.y() + umove );
    	newac[2] = ( lastac[2] + vmove );//newac.z( lastac.z() + vmove );
    	newac[0] = ( 1.0 - newac[1] - newac[2] );//newac.x( 1.0 - newac.y() - newac.z() );
      break;
    case 1:
    	newac[2] = ( lastac[2] + umove );//newac.z( lastac.z() + umove );
    	newac[0] = ( lastac[0] + vmove );//newac.x( lastac.x() + vmove );
    	newac[1] =( 1.0 - newac[2] - newac[0] );//newac.y( 1.0 - newac.z() - newac.x() );
      break;
    case 2:
    	newac[0] =( lastac[0] + umove );//newac.x( lastac.x() + umove );
    	newac[1] = ( lastac[1] + vmove );//newac.y( lastac.y() + vmove );
    	newac[2] = ( 1.0 - newac[0] - newac[1] );//newac.z( 1.0 - newac.x() - newac.y() );
      break;
    }

    // Keep it inside the patch

    if ( newac[0] >= -atol &&
         newac[1] >= -atol &&
         newac[2] >= -atol) {
      nout = 0;
    }
    else {
      if (move_ac_inside( newac, atol ) )
        nout++;
    }

    // Evaluate at the new location

    if (edge_id != -1)
      ac_at_edge( newac, newac, edge_id );  // move to edge first
    eval_bezier_patch(facet, newac, newpt);

    // Check for convergence

    distnew = (pt-newpt).length();//pt.distance_between(newpt);
    move = newpt - lastpt;
    movedist = move.length();
    if (movedist < tol ||
        distnew < tol ) {
      done = true;
      if (distnew < bestdist)
      {
        bestdist = distnew;
        bestac = newac;
        bestpt = newpt;
        bestnorm = norm;
      }
    }
    else {

      // don't allow more than 30 iterations

      iter++;
      if (iter > 30) {
        //if (movedist > tol * 100.0) nout=1;
        done = true;
      }

      // Check for divergence - don't allow more than 5 divergent
      // iterations

      if (distnew > lastdist) {
        diverge++;
        if (diverge > 10) {
          done = true;
          //if (movedist > tol * 100.0) nout=1;
        }
      }

      // Check if we are continuing to project outside the facet.
      // If so, then stop now

      if (nout > 3) {
        done = true;
      }

      // set up for next iteration

      if (!done) {
        if (distnew < bestdist)
        {
          bestdist = distnew;
          bestac = newac;
          bestpt = newpt;
          bestnorm = norm;
        }
        lastdist = distnew;
        lastpt = newpt;
        lastac = newac;
        move = pt - lastpt;
      }
    }
  }

  eval_pt = bestpt;
  if (eval_norm) {
    *eval_norm = bestnorm;
  }
  outside = (nout > 0) ? true : false;
  ac = bestac;

  return status;
}


//===========================================================================
//Function Name: ac_at_edge
//
//Member Type:  PRIVATE
//Description:  determine the area coordinate of the facet at the edge
//===========================================================================
void SmoothMeshEval::ac_at_edge( MBCartVect &fac,  // facet area coordinate
								MBCartVect &eac,   // edge area coordinate
                                int edge_id )      // id of edge
{
  double u, v, w;
  switch (edge_id)
  {
  case 0:
    u = 0.0;
    v = fac[1] / (fac[1] + fac[2]);//v = fac.y() / (fac.y() + fac.z());
    w = 1.0 - v;
    break;
  case 1:
	u = fac[0] / (fac[0] + fac[2]);// u = fac.x() / (fac.x() + fac.z());
    v = 0.0;
    w = 1.0 - u;
    break;
  case 2:
	u = fac[0] / (fac[0] + fac[1]);//u = fac.x() / (fac.x() + fac.y());
    v = 1.0 - u;
    w = 0.0;
    break;
  default:
    assert(0);
    break;
  }
  eac[0] = u;
  eac[1] = v;
  eac[2] = w; //= MBCartVect(u, v, w);
}


//===========================================================================
//Function Name: project_to_facet
//
//Member Type:  PUBLIC
//Description:  project to a single facet.  Uses the input areacoord as
//              a starting guess.
//===========================================================================
MBErrorCode SmoothMeshEval::project_to_facet( MBEntityHandle facet,
                                       MBCartVect &pt,
                                       MBCartVect &areacoord,
                                       MBCartVect &close_point,
                                       bool &outside_facet,
                                       double compare_tol)
{

  MBErrorCode stat = MB_SUCCESS;
  const MBEntityHandle * conn3;
  	int nnodes;
  	_mb->get_connectivity(facet, conn3, nnodes);
  	//
  	//double coords[9]; // store the coordinates for the nodes
  	//_mb->get_coords(conn3, 3, coords);
  	MBCartVect p[3];
  	_mb->get_coords(conn3, 3, (double*) &p[0]);
  /*CubitPoint *pt0 = facet->point(0);
  CubitPoint *pt1 = facet->point(1);
  CubitPoint *pt2 = facet->point(2);

  int interp_order = facet->eval_order();
  if (facet->is_flat())
  {
    interp_order = 0;
  }

  switch(interp_order) {
  case 0:
    close_point.x( areacoord.x() * pt0->x() +
                   areacoord.y() * pt1->x() +
                   areacoord.z() * pt2->x() );
    close_point.y( areacoord.x() * pt0->y() +
                   areacoord.y() * pt1->y() +
                   areacoord.z() * pt2->y() );
    close_point.z( areacoord.x() * pt0->z() +
                   areacoord.y() * pt1->z() +
                   areacoord.z() * pt2->z() );
    outside_facet = CUBIT_FALSE;
    break;
  case 1:
  case 2:
  case 3:
    assert(0);  // not available from this function
    break;
  case 4:
    {*/
      int edge_id = -1;
      stat = project_to_patch(facet, areacoord, pt, close_point, NULL,
                              outside_facet, compare_tol, edge_id );
   /* }
    break;
  }*/

  return stat;
}


//===========================================================================
//Function Name: is_at_vertex
//
//Member Type:  PRIVATE
//Description:  determine if the point is at one of the facet's vertices
//===========================================================================
bool SmoothMeshEval::is_at_vertex(
  MBEntityHandle facet,  // (IN) facet we are evaluating
  MBCartVect &pt,    // (IN) the point
  MBCartVect &ac,    // (IN) the ac of the point on the facet plane
  double compare_tol, // (IN) return TRUE of closer than this
  MBCartVect &eval_pt, // (OUT) location at vertex if TRUE
  MBCartVect *eval_norm_ptr ) // (OUT) normal at vertex if TRUE
{
  double dist;
  MBCartVect vert_loc;
  const double actol = 0.1;
  // get coordinates get_coords
  const MBEntityHandle * conn3;
	int nnodes;
	_mb->get_connectivity(facet, conn3, nnodes);
	//
	//double coords[9]; // store the coordinates for the nodes
	//_mb->get_coords(conn3, 3, coords);
	MBCartVect p[3];
	_mb->get_coords(conn3, 3, (double*) &p[0]);
	// also get the normals at nodes
	MBCartVect NN[3];
	_mb->tag_get_data(_gradientTag, conn3, 3, (double*)&NN[0]);
  if (fabs(ac[0]) < actol && fabs(ac[1]) < actol) {
    vert_loc = p[2];
    dist =(pt-vert_loc).length(); //pt.distance_between( vert_loc );
    if (dist <= compare_tol)
    {
      eval_pt = vert_loc;
      if (eval_norm_ptr)
      {
        *eval_norm_ptr = NN[2];
      }
      return true;
    }
  }

  if (fabs(ac[0]) < actol && fabs(ac[2]) < actol) {
    vert_loc = p[1];
    dist = (pt-vert_loc).length();//pt.distance_between( vert_loc );
    if (dist <= compare_tol)
    {
      eval_pt = vert_loc;
      if (eval_norm_ptr)
      {
        *eval_norm_ptr = NN[1];//facet->point(1)->normal( facet );
      }
      return true;
    }
  }

  if (fabs(ac[1]) < actol && fabs(ac[2]) < actol) {
    vert_loc = p[0];
    dist = (pt-vert_loc).length();//pt.distance_between( vert_loc );
    if (dist <= compare_tol)
    {
      eval_pt = vert_loc;
      if (eval_norm_ptr)
      {
        *eval_norm_ptr = NN[0];
      }
      return true;
    }
  }

  return false;
}


//===========================================================================
//Function Name: move_ac_inside
//
//Member Type:  PRIVATE
//Description:  find the closest area coordinate to the boundary of the
//              patch if any of its components are < 0
//              Return if the ac was modified.
//===========================================================================
bool SmoothMeshEval::move_ac_inside( MBCartVect &ac, double tol )
{
  int nout = 0;
  if (ac[0] < -tol) {
    ac[0]= 0.0 ;
    ac[1] =  ac[1] / (ac[1] + ac[2]) ; //( ac.y() / (ac.y() + ac.z()) ;
    ac[2] = 1.-ac[1]; //ac.z( 1.0 - ac.y() );
    nout++;
  }
  if (ac[1] < -tol) {
    ac[1] = 0.;//ac.y( 0.0 );
    ac[0] = ac[0]/(ac[0]+ac[2]);//ac.x( ac.x() / (ac.x() + ac.z()) );
    ac[2] = 1.-ac[0];//ac.z( 1.0 - ac.x() );
    nout++;
  }
  if (ac[2] < -tol) {
    ac[2] = 0.;// ac.z( 0.0 );
    ac[0] = ac[0]/(ac[0]+ac[1]);//ac.x( ac.x() / (ac.x() + ac.y()) );
    ac[1] = 1.-ac[0]; // ac.y( 1.0 - ac.x() );
    nout++;
  }
  return (nout > 0) ? true : false;
}


//===========================================================================
//Function Name: hodograph
//
//Member Type:  PUBLIC
//Description:  get the hodograph control points for the facet
//Note:  This is a triangle cubic patch that is defined by the
//       normals of quartic facet control point lattice.  Returned coordinates
//       in Nijk are defined by the following diagram
//
//
//                         *9               index  polar
//                        / \                 0     300    point(0)
//                       /   \                1     210
//                     7*-----*8              2     120
//                     / \   / \              3     030    point(1)
//                    /   \ /   \             4     201
//                  4*----5*-----*6           5     111
//                  / \   / \   / \           6     021
//                 /   \ /   \ /   \          7     102
//                *-----*-----*-----*         8     012
//                0     1     2     3         9     003    point(2)
//

//===========================================================================
//Function Name: eval_bezier_patch_normal
//
//Member Type:  PRIVATE
//Descriptoin:  evaluate the bezier patch defined at a facet
//===========================================================================
MBErrorCode SmoothMeshEval::eval_bezier_patch_normal( MBEntityHandle facet,
		MBCartVect &areacoord,
		MBCartVect &normal )
{
  // interpolate internal control points

  MBCartVect gctrl_pts[6];
  //facet->get_control_points( gctrl_pts );
  MBErrorCode rval = _mb->tag_get_data(_facetCtrlTag, &facet, 1, &gctrl_pts[0]);
  assert(rval == MB_SUCCESS);
  // _gradientTag
  // get normals at points
  const MBEntityHandle * conn3;
  		int nnodes;
  rval = _mb->get_connectivity(facet, conn3, nnodes);

  MBCartVect NN[3];
  rval =  _mb->tag_get_data(_gradientTag, conn3, 3, &NN[0]);

  assert(rval == MB_SUCCESS);

  if (fabs(areacoord[1] + areacoord[2]) < 1.0e-6) {
    normal = NN[0];
    return MB_SUCCESS;
  }
  if (fabs(areacoord[0] + areacoord[2]) < 1.0e-6) {
    normal = NN[1];//facet->point(1)->normal(facet);
    return MB_SUCCESS;
  }
  if (fabs(areacoord[0] + areacoord[1]) < 1.0e-6) {
    normal = NN[2]; //facet->point(2)->normal(facet);
    return MB_SUCCESS;
  }

  // compute the hodograph of the quartic Gregory patch

  MBCartVect Nijk[10];
  //hodograph(facet,areacoord,Nijk);
  // start copy from hodograph
  //CubitVector gctrl_pts[6];
  // facet->get_control_points( gctrl_pts );
   MBCartVect P_facet[3];

   //2,1,1
   /*P_facet[0] = (1.0e0 / (areacoord.y() + areacoord.z())) *
                (areacoord.y() * gctrl_pts[3] +
                 areacoord.z() * gctrl_pts[4]);*/
   P_facet[0] = (1.0e0 / (areacoord[1] + areacoord[2])) *
                   (areacoord[1] * gctrl_pts[3] +
                    areacoord[2] * gctrl_pts[4]);
   //1,2,1
   /*P_facet[1] = (1.0e0 / (areacoord.x() + areacoord.z())) *
                (areacoord.x() * gctrl_pts[0] +
                 areacoord.z() * gctrl_pts[5]);*/
   P_facet[1] = (1.0e0 / (areacoord[0] + areacoord[2])) *
                   (areacoord[0] * gctrl_pts[0] +
                    areacoord[2] * gctrl_pts[5]);
   //1,1,2
   /*P_facet[2] = (1.0e0 / (areacoord.x() + areacoord.y())) *
                (areacoord.x() * gctrl_pts[1] +
                 areacoord.y() * gctrl_pts[2]);*/
   P_facet[2] = (1.0e0 / (areacoord[0] + areacoord[1])) *
                   (areacoord[0] * gctrl_pts[1] +
                    areacoord[1] * gctrl_pts[2]);

   // corner control points are just the normals at the points

   //3, 0, 0
   Nijk[0] = NN[0];
   //0, 3, 0
   Nijk[3] = NN[1];
   //0, 0, 3
   Nijk[9] = NN[2];//facet->point(2)->normal(facet);

   // fill in the boundary control points.  Define as the normal to the local
   // triangle formed by the quartic control point lattice

   // store here again the 9 control points on the edges
   MBCartVect CP[9]; // 9 control points on the edges,
   rval = _mb->tag_get_data(_facetEdgeCtrlTag, &facet, 1, &CP[0]);
   // there are 3 CP for each edge, 0, 1, 2; first edge is 1-2
   //CubitFacetEdge *edge;
   //edge = facet->edge( 2 );
   //CubitVector ctrl_pts[5];
   //edge->control_points(facet, ctrl_pts);

   //2, 1, 0
   //Nijk[1] = (ctrl_pts[2] - ctrl_pts[1]) * (P_facet[0] - ctrl_pts[1]);
   Nijk[1] = (CP[7] - CP[6]) * (P_facet[0] - CP[6]);
   Nijk[1].normalize();

   //1, 2, 0
   //Nijk[2] = (ctrl_pts[3] - ctrl_pts[2]) * (P_facet[1] - ctrl_pts[2]);
   Nijk[2] = (CP[8] - CP[7]) * (P_facet[1] - CP[7]);
   Nijk[2].normalize();

   //edge = facet->edge( 0 );
   //edge->control_points(facet, ctrl_pts);

   //0, 2, 1
   //Nijk[6] = (ctrl_pts[1] - P_facet[1]) * (ctrl_pts[2] - P_facet[1]);
   Nijk[6] = (CP[0] - P_facet[1]) * (CP[1] - P_facet[1]);
   Nijk[6].normalize();

   //0, 1, 2
   //Nijk[8] = (ctrl_pts[2] - P_facet[2]) * (ctrl_pts[3] - P_facet[2]);
   Nijk[8] = (CP[1] - P_facet[2]) * (CP[2] - P_facet[2]);
   Nijk[8].normalize();

   //edge = facet->edge( 1 );
   //edge->control_points(facet, ctrl_pts);

   //1, 0, 2
   //Nijk[7] = (P_facet[2] - ctrl_pts[2]) * (ctrl_pts[1] - ctrl_pts[2]);
   Nijk[7] = (P_facet[2] - CP[4]) * (CP[3] - CP[4]);
   Nijk[7].normalize();

   //2, 0, 1
   //Nijk[4] = (P_facet[0] - ctrl_pts[3]) * (ctrl_pts[2] - ctrl_pts[3]);
   Nijk[4] = (P_facet[0] - CP[5]) * (CP[4] - CP[5]);
   Nijk[4].normalize();

   //1, 1, 1
   Nijk[5] = (P_facet[1] - P_facet[0]) * (P_facet[2] - P_facet[0]);
   Nijk[5].normalize();
  // end copy from hodograph

  // sum the contribution from each of the control points

  normal = MBCartVect(0.0e0, 0.0e0, 0.0e0);

  //i=3; j=0; k=0;
  double Bsum = 0.0;
  double B = cube(areacoord[0]);
  Bsum += B;
  normal += B * Nijk[0];

  //i=2; j=1; k=0;
  B = 3.0 * sqr(areacoord[0]) * areacoord[1];
  Bsum += B;
  normal += B * Nijk[1];

  //i=1; j=2; k=0;
  B = 3.0 * areacoord[0] * sqr(areacoord[1]);
  Bsum += B;
  normal += B * Nijk[2];

  //i=0; j=3; k=0;
  B = cube(areacoord[1]);
  Bsum += B;
  normal += B * Nijk[3];

  //i=2; j=0; k=1;
  B = 3.0 * sqr(areacoord[0]) * areacoord[2];
  Bsum += B;
  normal += B * Nijk[4];

    //i=1; j=1; k=1;
  B = 6.0 * areacoord[0] * areacoord[1] * areacoord[2];
  Bsum += B;
  normal += B * Nijk[5];

  //i=0; j=2; k=1;
  B = 3.0 * sqr(areacoord[1]) * areacoord[2];
  Bsum += B;
  normal += B * Nijk[6];

  //i=1; j=0; k=2;
  B = 3.0 * areacoord[0] * sqr(areacoord[2]);
  Bsum += B;
  normal += B * Nijk[7];

  //i=0; j=1; k=2;
  B = 3.0 * areacoord[1] * sqr(areacoord[2]);
  Bsum += B;
  normal += B * Nijk[8];

  //i=0; j=0; k=3;
  B = cube(areacoord[2]);
  Bsum += B;
  normal += B * Nijk[9];

  //assert(fabs(Bsum - 1.0) < 1e-9);

  normal.normalize();

  return MB_SUCCESS;
}
