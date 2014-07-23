/*
 * SmoothVertex.cpp
 *
 *  Created on: Jul 19, 2010
 *      Author: iulian
 */

#include "SmoothVertex.hpp"
#include <vector>
#include "assert.h"

SmoothVertex::SmoothVertex(MBInterface * mb, MBEntityHandle vset, MBInterface * mbo):
		SmoothBase(mb, vset, mbo)
{
	_mbOut->create_meshset(MESHSET_SET, _oSet);
}

SmoothVertex::~SmoothVertex() {
	// TODO Auto-generated destructor stub
}

MBErrorCode SmoothVertex::create_mesh_vertex()
{
	// just get the coordinates from original vertex, then create the new vertex in the new set

	std::vector<MBEntityHandle>  nodes;
	MBErrorCode rval = _mb->get_entities_by_handle(_set, nodes);
	// there should be only one node, but we do not care to check
	assert(nodes.size()>=1);

	double coords[3];
	rval = _mb->get_coords(&nodes[0], 1, coords);

	MBEntityHandle newNode;
	rval = _mbOut -> create_vertex(coords, newNode );
	_mbOut->add_entities(_oSet, &newNode, 1);

	return rval;
}
// this will be used to generate new mesh edges on the geo edges
// this belongs to _mbOut instance
MBEntityHandle SmoothVertex::get_new_node()
{
	std::vector<MBEntityHandle>  nodes;
	MBErrorCode rval = _mbOut->get_entities_by_handle(_oSet, nodes);
	return nodes[0];// first node
}
