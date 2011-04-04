/**
 * \file Geometrize.cpp
 *
 * \brief add geometry sets to a bland surface mesh file, so it can be 
 *  imported easily into cubit; the preffered out is exo
 *
 *  typical scenario: geometrize an smf file, add boundary edges, and geometry sets
 *   ( surface, one periodic boundary and one vertex)
 *
 */

//#include "iMesh.h"

#include "MBRange.hpp"
#include "MBSkinner.hpp"
#include "MBInterface.hpp"
#include "MBCore.hpp"
#include "MBTagConventions.hpp"

#include "moab/GeomTopoTool.hpp"

#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <assert.h>
#include <string.h>

bool debug = false;

#define ERRORR(a) {if (iBase_SUCCESS != err) {printf(a); return 1;}}
#define MBERRORR(a) {if (MB_SUCCESS != rval) {std::cerr << a << std::endl; return 1;}}

int main(int argc, char *argv[]) {
	// Check command line arg
	std::string input_filename, output_filename;

	if (argc < 3) {
		std::cout << "Usage: " << argv[0]
				<< " <input_mesh_filename> <output_mesh_filename> "
				<< std::endl;
		return 0;
	} else {
		input_filename = argv[1];
		output_filename = argv[2];
	}

	MBCore moab;
	MBInterface * MBI = &moab;

	MBErrorCode rval = MBI->load_mesh(input_filename.c_str(), NULL, 0);
	if (rval != MB_SUCCESS)
		return 1;
	MBRange surface_ents, edge_ents, loop_range;

	rval = MBI->get_entities_by_type(0, MBTRI, surface_ents);

	// mb
	// get surface sets
	MBTag geom_tag;
	rval = MBI->tag_create(GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE,
			MB_TYPE_INTEGER, geom_tag, 0, true);
	assert(MB_SUCCESS==rval);

	MBTag gid_tag;
	rval = MBI->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
			MB_TYPE_INTEGER, gid_tag, 0, true);
	assert(MB_SUCCESS==rval);

	MBEntityHandle face ;
	rval = MBI->create_meshset(MESHSET_SET, face);
	assert(MB_SUCCESS==rval);
	// set some tags on face
	int dimension = 2;
	rval = MBI->tag_set_data(geom_tag, &face, 1, &dimension);
	assert(MB_SUCCESS==rval);

	int gid = 1;
	rval = MBI->tag_set_data(gid_tag, &face, 1, &gid);

	// create an edge set
	MBEntityHandle edge;
	rval = MBI->create_meshset(MESHSET_ORDERED, edge);
	assert(MB_SUCCESS==rval);
	dimension =1;
	rval = MBI->tag_set_data(geom_tag, &edge, 1, &dimension);
	assert(MB_SUCCESS==rval);
	gid = 1;
	rval = MBI->tag_set_data(gid_tag, &edge, 1, &gid);
	assert(MB_SUCCESS==rval);

	// create a vertex set
	MBEntityHandle vertex;
	rval = MBI->create_meshset(MESHSET_SET, vertex);
	assert(MB_SUCCESS==rval);
	dimension =0;
	rval = MBI->tag_set_data(geom_tag, &vertex, 1, &dimension);
	assert(MB_SUCCESS==rval);
	gid = 1;
	rval = MBI->tag_set_data(gid_tag, &edge, 1, &gid);
	assert(MB_SUCCESS==rval);

	// add triangles to the face

	rval = MBI->add_entities(face, surface_ents);
	assert(MB_SUCCESS==rval);

	MBSkinner tool(MBI);
	rval = tool.find_skin(surface_ents, 1, edge_ents);
	MBERRORR("Skinner failure.");

	// arrange edges according to loops; one edge per loop
	//
	rval = MBI ->add_parent_child( face, edge);
	assert(MB_SUCCESS==rval);
	rval = MBI ->add_parent_child( edge, vertex);
	assert(MB_SUCCESS==rval);

	// need to arrange the edges; first node in the edge list is vertex1
	std::vector<MBEntityHandle> firstV;
	std::vector<MBEntityHandle> endV;
	std::vector<MBEntityHandle> loops; // edge loops, edges need to be chained in loops

	std::vector<int> marks;
	for (MBRange::iterator it =edge_ents.begin(); it!=edge_ents.end(); it++ )
	{
		MBEntityHandle e = *it;
		int nnodes ;
		const MBEntityHandle * conn2;
		rval = MBI->get_connectivity(e, conn2, nnodes);
		assert (nnodes==2);
		firstV.push_back(conn2[0]);
		endV.push_back(conn2[1]);
		marks.push_back(0);
	}
	// add an edge, one at a time, until all exhausted
	// first node is added to the vertex set
	MBI->add_entities(vertex, &firstV[0], 1);
	loops.push_back(edge_ents[0]);// first edge
	marks[0] = 1; // marked already
	MBEntityHandle curVert = endV[0];
	int numEdges = firstV.size();
	while (curVert!=firstV[0])
	{
		for (int i=1; i<numEdges; i++)
		{
			if (marks[i])
				continue;
			if (curVert == firstV[i])
			{
				loops.push_back(edge_ents[i]);
				curVert = endV[i];
				marks[i] = 1;// to not add it again
				break;// from for loop
			}
		}
	}

	// add loops edges to the edge set
	MBI->add_entities(edge, &loops[0], loops.size());// only one edge
	// add the sense tags to the edges on the boundary
	// get the GeomTopoTool on the mesh

	moab::GeomTopoTool geomTool(MBI);
	int sense = 1;//
	std::vector<int> senses;
	senses.push_back(sense);
	std::vector<MBEntityHandle> sets1;
	//MBEntityHandle setFace((MBEntityHandle) sets[2]);
	sets1.push_back(face);
	//MBEntityHandle setCurve((MBEntityHandle) sets[1]);
	// GeomTopoTool::set_senses(EntityHandle edge,  std::vector<EntityHandle> &faces,  std::vector<int> &senses)
	geomTool.set_senses(edge, sets1, senses);// only one set of dimension 1
	// write the mesh as an h5m file, because it can save the sets too, and their relationship
	// then look at it with mbsize -ll

	rval = MBI-> write_mesh(output_filename.c_str());
	MBERRORR("Couldn't write geometrized mesh file\n");

	return 0;
}
