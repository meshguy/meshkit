/*
 * CreateShell.cpp
 *
 *  Created on: Jul 12, 2010
 *      Author: iulian
 */
#include "iMesh.h"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"

#include "moab/GeomTopoTool.hpp"

#include <assert.h>

int main (int argc, char * argv[])
{
	// initialize mesh interface instances
	int err;

	iMesh_Instance mesh;
	iMesh_newMesh(0, &mesh, &err, 0);
	assert(iBase_SUCCESS == err);
	MBInterface * mb = (MBInterface *) mesh;
		//_mbset = (MBEntityHandle) this->_surf_set;
	MBEntityHandle root = mb->get_root_set();
	// create some vertices
	double coords [] = { 0, 0, 0,
			             1, 0, 0.1,
			             2, 0, 0,
			             3, 0, -0.1,
			             0, 1, 0,
			             1, 1, 0,
			             2, 1, 0,
			             3, 1, -0.1,
			             0, 2, 0,
			             1, 2, -0.1,
			             2, 2, -0.1,
			             3, 2, -0.2,
			             0, 0, 1,
			             1, 0, 0.9,
			             2, 0.1, 0.85,
			             3, 0.2, 0.8,
			             0, 0.1, 2,
			             1, 0.1, 2,
			             2.1, 0.2, 2.1,
			             3.1, 0.2, 2.1 };

	int nvert = 20;
	MBRange verts;
	MBErrorCode rval = mb->create_vertices(coords, nvert, verts);
	assert(MB_SUCCESS==rval);

	MBEntityHandle connec [] = { 1, 2, 5,
	                  5, 2, 6,
	                  2, 3, 6,
	                  6, 3, 7,
	                  3, 4, 7,
	                  7, 4, 8,
	                  5, 6, 9,
	                  9, 6, 10,
	                  6, 7, 10,
	                  10, 7, 11,
	                  7, 8, 11,
	                  11, 8, 12,  // first face, 1-12
	                  13, 14, 1,
	                  1, 14, 2,
	                  14, 15, 2,
	                  2, 15, 3,
	                  15, 16, 3,
	                  3, 16, 4,
	                  17, 18, 13,
	                  13, 18, 14,
	                  18, 19, 14,
	                  14, 19, 15,
	                  19, 20, 15,
	                  15, 20, 16 // second face, 13-24
	                 };
	MBEntityHandle elem;
	int nbTri = sizeof(connec)/3/sizeof(MBEntityHandle);
	int i = 0;
	std::vector<MBEntityHandle> tris;
	for (i=0; i<nbTri; i++)
	{
		mb->create_element(MBTRI, &connec[3*i], 3, elem);
		tris.push_back(elem);
	}


	// create some edges too
	MBEntityHandle edges [] = { 1, 2,
			                    2, 3,
			                    3, 4,  // geo edge 1  1:3
			                    4, 8,
			                    8, 12,  // geo 2         4:5
			                    12, 11,
			                    11, 10,
			                    10, 9, // geo 3 6:8
			                    9, 5,
			                    5, 1,  // geo 4  9:10
								1, 13,
	                            13, 17, // geo 5  11:12
	                            17, 18,
	                            18, 19,
	                            19, 20, // geo 6  13:15
	                            20, 16,
	                            16, 4   // geo 7  16:17
	                             };
	int nbEdges = sizeof(edges)/2/sizeof(MBEntityHandle);
	std::vector<MBEntityHandle> edgs;
	for (i=0; i<nbEdges; i++)
	{
		mb->create_element(MBEDGE, &edges[2*i], 2, elem);
		edgs.push_back(elem);
	}
	// create some sets, and create some ordered sets for edges
	MBEntityHandle face1, face2;
	rval = mb->create_meshset(MESHSET_SET, face1);
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(face1, &tris[0], 12);
	assert(MB_SUCCESS==rval);

	rval = mb->create_meshset(MESHSET_SET, face2);
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(face2, &tris[12], 12); // next 12 triangles
	assert(MB_SUCCESS==rval);

	// mb
	// get surface sets
	MBTag geom_tag;
	rval = mb->tag_create( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE,
	                         MB_TYPE_INTEGER, geom_tag, 0, true );
	assert(MB_SUCCESS==rval);

	MBTag gid_tag;
	rval = mb->tag_create( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
							 MB_TYPE_INTEGER, gid_tag, 0, true );
	assert(MB_SUCCESS==rval);

	int dimension = 2;
	rval = mb->tag_set_data(geom_tag, &face1, 1, &dimension);
	assert(MB_SUCCESS==rval);

	rval = mb->tag_set_data(geom_tag, &face2, 1, &dimension);
	assert(MB_SUCCESS==rval);

	int gid = 1;
	rval = mb->tag_set_data(gid_tag, &face1, 1, &gid);
	assert(MB_SUCCESS==rval);
	gid = 2;
	rval = mb->tag_set_data(gid_tag, &face2, 1, &gid);
	assert(MB_SUCCESS==rval);

	// create some edges
	MBEntityHandle edge[7]; //edge[0] has EH 1...
	dimension = 1;
	for (i=0; i<7; i++)
	{
		rval = mb->create_meshset(MESHSET_ORDERED, edge[i]);
		assert(MB_SUCCESS==rval);
		rval = mb->tag_set_data(geom_tag, &edge[i], 1, &dimension);
		assert(MB_SUCCESS==rval);
		gid = i+1;
		rval = mb->tag_set_data(gid_tag, &edge[i], 1, &gid);
		assert(MB_SUCCESS==rval);

	}


	rval = mb->add_entities(edge[0], &edgs[0], 3); // first 3 mesh edges...
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(edge[1], &edgs[3], 2); //
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(edge[2], &edgs[5], 3); //
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(edge[3], &edgs[8], 2); //
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(edge[4], &edgs[10], 2); //
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(edge[5], &edgs[12], 3); //
	assert(MB_SUCCESS==rval);
	rval = mb->add_entities(edge[6], &edgs[15], 2);
	assert(MB_SUCCESS==rval);

	// create some sets for vertices; also need to create some for parent/child relationships
	MBEntityHandle vertSets[6];// start from 0
	dimension =0;
	for (i=0; i<6; i++)
	{
		rval = mb->create_meshset(MESHSET_SET, vertSets[i]);
		assert(MB_SUCCESS==rval);
		rval = mb->tag_set_data(geom_tag, &vertSets[i], 1, &dimension);
		assert(MB_SUCCESS==rval);
		gid = i+1;
		rval = mb->tag_set_data(gid_tag, &vertSets[i], 1, &gid);
		assert(MB_SUCCESS==rval);
	}



	MBEntityHandle v(1); // first vertex;
	rval = mb->add_entities(vertSets[0], &v, 1);
	assert(MB_SUCCESS==rval);
	v = MBEntityHandle (4);
	rval = mb->add_entities(vertSets[1], &v, 1);
	assert(MB_SUCCESS==rval);
	v = MBEntityHandle (9);
	rval = mb->add_entities(vertSets[2], &v, 1);
	assert(MB_SUCCESS==rval);
	v = MBEntityHandle (12);
	rval = mb->add_entities(vertSets[3], &v, 1);
	assert(MB_SUCCESS==rval);
	v = MBEntityHandle (17);
	rval = mb->add_entities(vertSets[4], &v, 1);
	assert(MB_SUCCESS==rval);
	v = MBEntityHandle (20);
	rval = mb->add_entities(vertSets[5], &v, 1);
	assert(MB_SUCCESS==rval);

	// need to add parent-child relations between sets
	// edge 1 : 1-2
	rval = mb ->add_parent_child( edge[0], vertSets[0]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( edge[0], vertSets[1]);
	assert(MB_SUCCESS==rval);
	// edge 2 : 2-4
	rval = mb ->add_parent_child( edge[1], vertSets[1]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( edge[1], vertSets[3]);
	assert(MB_SUCCESS==rval);

	// edge 3 : 4-3
	rval = mb ->add_parent_child( edge[2], vertSets[3]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( edge[2], vertSets[2]);
	assert(MB_SUCCESS==rval);

	// edge 4 : 4-1
	rval = mb ->add_parent_child( edge[3], vertSets[2]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( edge[3], vertSets[0]);
	assert(MB_SUCCESS==rval);

	// edge 5 : 1-5
	rval = mb ->add_parent_child( edge[4], vertSets[0]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( edge[4], vertSets[4]);
	assert(MB_SUCCESS==rval);

	// edge 6 : 5-6
	rval = mb ->add_parent_child( edge[5], vertSets[4]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( edge[5], vertSets[5]);
	assert(MB_SUCCESS==rval);

	// edge 7 : 6-2
	rval = mb ->add_parent_child( edge[6], vertSets[5]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( edge[6], vertSets[1]);
	assert(MB_SUCCESS==rval);

	// face 1: edges 1, 2, 3, 4
	rval = mb ->add_parent_child( face1, edge[0]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( face1, edge[1]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( face1, edge[2]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( face1, edge[3]);
	assert(MB_SUCCESS==rval);

	// face 2: edges 1, 5, 6, 7
	rval = mb ->add_parent_child( face2, edge[0]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( face2, edge[4]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( face2, edge[5]);
	assert(MB_SUCCESS==rval);
	rval = mb ->add_parent_child( face2, edge[6]);
	assert(MB_SUCCESS==rval);

	// the orientation and senses need to be set for face edges

	moab::GeomTopoTool gTopoTool(mb, false);

	// set senses !!
	std::vector<MBEntityHandle> faces;
	faces.push_back(face1); // the face1 has all edges oriented positively
	std::vector<int> senses;
	senses.push_back(0); //

	//faces.push_back(face1);
	//faces.push_back(face2);
	gTopoTool.set_senses(edge[1], faces, senses);
	gTopoTool.set_senses(edge[2], faces, senses);
	gTopoTool.set_senses(edge[3], faces, senses);

	faces[0]=face2; // set senses for edges for face2
	gTopoTool.set_senses(edge[4], faces, senses);
	gTopoTool.set_senses(edge[5], faces, senses);
	gTopoTool.set_senses(edge[6], faces, senses);

	// the only complication is edge1 (edge[0]), that has face1 forward and face 2 reverse
	faces[0]=face1;
	faces.push_back(face2);
	senses.push_back(1); // 1 is reverse; face 2 is reverse for edge1 (0)
	// forward == 0, reverse ==1
	gTopoTool.set_senses(edge[0], faces, senses);


	mb->write_mesh("shell.h5m");



	mb-> write_mesh("shell.vtk");

	return 0;

}
