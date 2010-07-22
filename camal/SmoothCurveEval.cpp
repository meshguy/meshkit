/*
 * SmoothCurveEval.cpp
 *
 *  Created on: Jun 9, 2010
 *      Author: iulian
 */

#include "SmoothCurveEval.hpp"
#include "SmoothVertex.hpp"
#include "CubitVector.hpp"
#include "SmoothFaceEval.hpp"
#include "CAMALSizeEval.hpp"

// the edge mesher from Camal; we do not really need it
#include "CMLEdgeMesher.hpp"

#include "moab/GeomTopoTool.hpp"


SmoothCurveEval::SmoothCurveEval(MBInterface * mb, MBEntityHandle curve, MBInterface * mbo) :
		SmoothBase(mb, curve, mbo)
{
	_mbOut->create_meshset(MESHSET_ORDERED, _oSet);
	_cmlEdgeMesher = new CMLEdgeMesher (this, CML::STANDARD);
	_cmlEdgeMesher->set_sizing_function(CML::LINEAR_SIZING);
	_leng = 0; // not initialized
	_edgeTag = 0; // not initialized
}
SmoothCurveEval::~SmoothCurveEval() {
	// TODO Auto-generated destructor stub
}

double SmoothCurveEval::arc_length()
{

	return _leng;
}

    //! \brief Get the parametric status of the curve.
    //!
    //! \return \a true if curve is parametric, \a false otherwise.
bool SmoothCurveEval::is_parametric()
{
	return true;
}

    //! \brief Get the periodic status of the curve.
    //!
    //! \param period The period of the curve if periodic.
    //!
    //! \return \a true if curve is periodic, \a false otherwise.
bool SmoothCurveEval::is_periodic(double& period)
{
	//assert(_ref_edge);
	//return _ref_edge->is_periodic(   period);
	MBRange vsets;
	MBErrorCode rval = _mb->get_child_meshsets(_set, vsets);// num_hops =1
	if (vsets.size() == 1)
	{
		period = _leng;
		return true;//  true , especially for ice sheet data
	}
	return false;
}

    //! \brief Get the parameter range of the curve.
    //!
    //! \param u_start The beginning curve parameter
    //! \param u_end The ending curve parameter
    //!
    //! \note The numerical value of \a u_start may be greater
    //! than the numerical value of \a u_end.
void SmoothCurveEval::get_param_range(double& u_start, double& u_end)
{
    //assert(_ref_edge);
    u_start = 0;
    u_end = 1.;

    return;
}

    //! Compute the parameter value at a specified distance along the curve.
    //!
    //! \param u_root The start parameter from which to compute the distance
    //! along the curve.
    //! \param arc_length The distance to move along the curve.
    //!
    //! \note For positive values of \a arc_length the distance will be
    //! computed in the direction of increasing parameter value along the
    //! curve.  For negative values of \a arc_length the distance will be
    //! computed in the direction of decreasing parameter value along the
    //! curve.
    //!
    //! \return The parametric coordinate u along the curve
double SmoothCurveEval::u_from_arc_length(double u_root, double arc_leng)
{

	if (_leng<=0)
		return 0;
	return u_root + arc_leng/_leng;
}


    //! \brief Evaluate the curve at a specified parameter value.
    //!
    //! \param u The parameter at which to evaluate the curve
    //! \param x The x coordinate of the evaluated point
    //! \param y The y coordinate of the evaluated point
    //! \param z The z coordinate of the evaluated point
bool SmoothCurveEval::position_from_u(double u,
                               double& x, double& y, double& z )
{

		// _fractions are increasing, so find the
		double * ptr = std::lower_bound(&_fractions[0],
				&_fractions[_fractions.size()], u);
		int index = ptr - &_fractions[0];
		double nextFraction = _fractions[index];
		double prevFraction = 0;
		if (index > 0) {
			prevFraction = _fractions[index - 1];
		}
		double t = (u - prevFraction) / (nextFraction - prevFraction);

		MBEntityHandle edge = _entities[index];

		MBCartVect position;
		MBErrorCode rval = evaluate_smooth_edge(edge, t, position);
		x = position[0];
		y = position[1];
		z = position[2];

		return true;



}

    //! \brief Move a point near the curve to the closest point on the curve.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
void SmoothCurveEval::move_to_curve(double& x, double& y, double& z)
{

	return;
}

    //! Get the u parameter value on the curve closest to x,y,z
    //! and the point on the curve.
    //!
    //! \param x The x coordinate of the point
    //! \param y The y coordinate of the point
    //! \param z The z coordinate of the point
    //!
    //! \return The parametric coordinate u on the curve
double SmoothCurveEval::u_from_position(double x, double y, double z)
{
	return 0;// not implemented
}

    //! \brief Get the starting point of the curve.
    //!
    //! \param x The x coordinate of the start point
    //! \param y The y coordinate of the start point
    //! \param z The z coordinate of the start point
void SmoothCurveEval::start_coordinates(double& x, double& y, double& z)
{


	int nnodes;
	const MBEntityHandle * conn2;
	_mb->get_connectivity(_entities[0], conn2, nnodes);
	double c[3];
	_mb->get_coords(conn2, 1, c);

	x = c[0];
	y = c[1];
	z = c[2];

	return;
}

    //! \brief Get the ending point of the curve.
    //!
    //! \param x The x coordinate of the start point
    //! \param y The y coordinate of the start point
    //! \param z The z coordinate of the start point
void SmoothCurveEval::end_coordinates(double& x, double& y, double& z)
{

	int nnodes;
	const MBEntityHandle * conn2;
	_mb->get_connectivity(_entities[_entities.size()-1], conn2, nnodes);
	double c[3];
	// careful, the second node here
	_mb->get_coords(&conn2[1], 1, c);

	x = c[0];
	y = c[1];
	z = c[2];
	return;
}

// this will recompute the 2 tangents for each edge, considering the geo edge they are into
void SmoothCurveEval::compute_tangents_for_each_edge()
{
	// will retrieve the edges in each set; they are retrieved in order they were put into, because
	// these sets are "MESHSET_ORDERED"
	// retrieve the tag handle for the tangents; it should have been created already
	// this tangents are computed for the chain of edges that form a geometric edge
	// some checks should be performed on the vertices, but we trust the correctness of the model completely
	// (like the vertices should match in the chain...)
	MBTag tangentsTag;
	MBErrorCode rval = _mb->tag_get_handle("TANGENTS", tangentsTag);
	if (rval != MB_SUCCESS)
		return; // some error should be thrown
	std::vector<MBEntityHandle> entities;
	_mb->get_entities_by_type(_set, MBEDGE, entities); // no recursion!!
	// basically, each tangent at a node will depend on previous tangent
	int nbEdges = entities.size();
	// now, we can advance in the loop
	// the only special problem is if the first node coincides with the last node, then we should
	// consider the closed loop; or maybe we should look at angles in that case too?
	// also, should we look at the 2 semi-circles case? How to decide if we need to continue the "tangents"
	// maybe we can do that later, and we can alter the tangents at the feature nodes, in the directions of the loops
	// again, do we need to decide the "closed" loop or not? Not yet...
	MBEntityHandle previousEdge = entities[0]; // this is the first edge in the chain
	MBCartVect TP[2]; // tangents for the previous edge
	rval = _mb->tag_get_data(tangentsTag, &previousEdge, 1, &TP[0]);// tangents for previous edge
	if (rval != MB_SUCCESS)
		return; // some error should be thrown
	MBCartVect TC[2]; // tangents for the current edge
	MBEntityHandle currentEdge;
	for (int i=1; i<nbEdges; i++)
	{
		// current edge will start after first one
		currentEdge = entities[i];
		rval = _mb->tag_get_data(tangentsTag, &currentEdge, 1, &TC[0]);//
		// now compute the new tangent at common vertex; reset tangents for previous edge and current edge
		// a little bit of CPU and memory waste, but this is life
		MBCartVect T = 0.5*TC[0]+0.5*TP[1]; //
		T.normalize();
		TP[1] = T;
		rval = _mb->tag_set_data(tangentsTag, &previousEdge, 1, &TP[0]);//
		TC[0] = T;
		rval = _mb->tag_set_data(tangentsTag, &currentEdge, 1, &TC[0]);//
		// now set the next edge
		previousEdge = currentEdge;
		TP[0] = TC[0];
		TP[1] = TC[1];
	}
	return;
}

void SmoothCurveEval::compute_control_points_on_boundary_edges(double min_dot,
			  std::map<MBEntityHandle, SmoothFaceEval*> & mapSurfaces, MBTag controlPointsTag, MBTag markTag)
{
	// these points really need the surfaces they belong to, because the control points on edges
	// depend on the normals on surfaces
	// the control points are averaged from different surfaces, by simple mean average
	// the surfaces have
	// do we really need min_dot here?
	// first of all, find out the SmoothFaceEval for each surface set that is adjacent here
	moab::GeomTopoTool  gTopoTool(_mb);
	std::vector<MBEntityHandle> faces;
	std::vector<int> senses;
	MBErrorCode rval = gTopoTool.get_senses(_set , faces, senses);
	if (MB_SUCCESS!=rval)
		return;

	// need to find the smooth face attached
	int numSurfacesAdjacent = faces.size();
	// get the edges, and then get the
	//std::vector<MBEntityHandle> entities;
	_mb->get_entities_by_type(_set, MBEDGE, _entities); // no recursion!!
	// each edge has the tangent computed already
	MBTag tangentsTag;
	rval = _mb->tag_get_handle("TANGENTS", tangentsTag);
	if (rval != MB_SUCCESS)
		return; // some error should be thrown

	// we do not want to search every time
	std::vector<SmoothFaceEval*> smoothFaceArray;
	int i=0;
	for (i=0; i<numSurfacesAdjacent; i++)
	{
		SmoothFaceEval * sms = mapSurfaces[faces[i]];
		smoothFaceArray.push_back(sms);
	}

	int e=0;
	for (e=0; e<_entities.size(); e++)
	{
		MBCartVect zero(0.);
		MBCartVect  ctrlP[3]= {zero, zero, zero};// null positions initially
		// the control points are averaged from connected faces
		MBEntityHandle edge = _entities[e]; // the edge in the chain

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

		MBCartVect chord = P[1]-P[0];
		_leng+= chord.length();
		_fractions.push_back(_leng);
		MBCartVect N[2];


		//MBCartVect N0(&normalVec[0]);
		//MBCartVect N3(&normalVec[3]);
		MBCartVect T[2]; // T0, T3
		//if (edge->num_adj_facets() <= 1) {
		//stat = compute_curve_tangent(edge, min_dot, T0, T3);
		//  if (stat != CUBIT_SUCCESS)
		//	return stat;
		//} else {
		//}
		rval = _mb->tag_get_data(tangentsTag, &edge, 1, &T[0]);
		assert(rval == MB_SUCCESS);

		for (i=0; i<numSurfacesAdjacent; i++)
		{
			MBCartVect controlForEdge[3];
			rval = smoothFaceArray[i]->get_normals_for_vertices( conn2,  N );
			assert(rval == MB_SUCCESS);
			rval = smoothFaceArray[i]->init_edge_control_points(P[0], P[1], N[0], N[1], T[0], T[1],
					controlForEdge);
			// accumulate those over faces!!!
			for (int j=0; j<3; j++)
			{
				ctrlP[j]+=controlForEdge[j];
			}
		}
		// now divide them for the average position!
		for (int j=0; j<3; j++)
		{
			ctrlP[j]/=numSurfacesAdjacent;
		}
		// we are done, set the control points now!
		//edge->control_points(ctrl_pts, 4);
		rval = _mb->tag_set_data(controlPointsTag, &edge, 1, &ctrlP[0]);
		assert(rval == MB_SUCCESS);
		this->_edgeTag = controlPointsTag;// this is a tag that will be stored with the edge
		// is that a waste of memory or not...
		// also mark the edge for later on
		unsigned char used = 1;
		_mb->tag_set_data(markTag, &edge, 1, &used);
	}
	// now divide fractions, to make them vary from 0 to 1
	assert(_leng>0.);
	for (e=0; e<_entities.size(); e++)
		_fractions[e]/=_leng;


}

void SmoothCurveEval::estimate_mesh_count(double curve_mesh_size, int & num_points_out) // this is based on initial settings

{
	CAMALSizeEval size_eval2(curve_mesh_size);
	_cmlEdgeMesher->generate_mesh_sizing(&size_eval2, num_points_out);
	mesh_count = num_points_out-1;
	return;
}

MBErrorCode SmoothCurveEval::evaluate_smooth_edge(MBEntityHandle eh, double &tt,
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

	if (0==_edgeTag)
	{
		rval = _mb->tag_get_handle("CONTROLEDGE", _edgeTag);
		if (rval != MB_SUCCESS)
			return rval;
	}
	rval = _mb->tag_get_data(_edgeTag, &eh, 1, (double*) &controlPoints[0]);
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

void SmoothCurveEval::create_mesh_edges(std::map<MBEntityHandle, SmoothVertex*>  & mapVertices)
{
	// first, we need to find the nodes already created in _mbo, at mesh vertices
	// also, get the nodes at start and end; compare with the nodes from children of the
	// edge set
	// mesh_set
	MBRange vsets;
	MBErrorCode rval = _mb->get_child_meshsets(_set, vsets);// num_hops =1
	int numVertexSets = vsets.size();
	if (numVertexSets!=2 && numVertexSets!=1)
		return;// error
	// there should be exactly 2 children sets to each edge
	MBEntityHandle startVert, endVert; // those are node handles for _mb
	const MBEntityHandle * conn2;
	int nnodes ;

	rval = _mb->get_connectivity(_entities[0], conn2, nnodes);
	if (nnodes!=2 )
		return ;// failure
	startVert = conn2[0]; // first node of first edge
	rval = _mb->get_connectivity(_entities[_entities.size()-1], conn2, nnodes);
	if (nnodes!=2 )
		return ;// failure
	endVert = conn2[1]; // second node of last edge

	// get the first entity of each vsets; it should be either startVert or endVert;
	// basically we need the vertex from SmoothVertex that was already created in _mbo

	std::vector<MBEntityHandle>  nodes;
	rval = _mb->get_entities_by_handle(vsets[0], nodes);

	if (nodes.size()<1)
		return; // failure
	//MBEntityHandle newNodeStart, newNodeEnd;

	SmoothVertex * smV0 = mapVertices[vsets[0]];
	SmoothVertex * smV1 = NULL;
	if (numVertexSets==1)
		smV1 = smV0; // periodic
	else
		smV1 = mapVertices[vsets[1]];

	MBEntityHandle newNode1, newNode2;// node handles on new mesh instance
	if (nodes[0]==startVert)
	{
		newNode1 = smV0->get_new_node();
		newNode2 = smV1->get_new_node();
	}
	else
	{
		newNode1 = smV1->get_new_node(); // sort of reversed
		newNode2 = smV0->get_new_node();
	}

	// so now we know the start and end nodes; we just need to create some
	// more nodes in _mbOut, and new edges, that will be added in order to _oSet
	// invoke the _cmlEdgeMesher

	// start copy
	// now, the number of intervals will be set to num_points_out (even)
	int num_intervals = mesh_count;
	int num_points_out;
	_cmlEdgeMesher->generate_mesh_uniform(num_intervals, num_points_out);


	// we are now settled on an odd number of edge segments;
	// get them, and put them in the bdy_loops...
	std::vector<double> coords;
	coords.resize(num_points_out * 3);// we allocate more, but we will use one less
	bool stat = _cmlEdgeMesher->get_mesh(num_points_out, &coords[0]);
	// now we can create new nodes in _mbOut, for the new internal nodes
	int numInteriorNodes = num_points_out-2; // exclude the extreme nodes

	MBRange newNodes;
	newNodes.clear();
	if (numInteriorNodes > 0)
	{
		// start at the second node
		_mbOut->create_vertices(&coords[3], numInteriorNodes, newNodes);
	}
	MBEntityHandle conn2New[2];
	conn2New[0] = newNode1;
	MBEntityHandle newEdge;
	for (MBRange::iterator it = newNodes.begin(); it!=newNodes.end(); it++)
	{
		conn2New[1] = *it; // the second new node in an edge
		// create the new edge
		_mbOut->create_element(MBEDGE, conn2New, 2, newEdge);
		_mbOut-> add_entities(_oSet, &newEdge, 1);
		conn2New[0] = *it;
	}
	conn2New [1]=newNode2;
	_mbOut->create_element(MBEDGE, conn2New, 2, newEdge);
	_mbOut-> add_entities(_oSet, &newEdge, 1);
	// the last edge created will have the last node
	// end copy
	return ;
}
