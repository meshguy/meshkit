/*
 * primitives.cpp
 *
 *  Created on: Mar 22, 2010
 *      Author: iulian
 */

#include "primitives.h"

void filterValid(MBInterface * mb, std::vector<MBEntityHandle> & io) {
	int next = 0, size = io.size();
	std::vector<unsigned char> tags(size);
	MBErrorCode rval = mb->tag_get_data(validTag, &io[0], size, &tags[0] );
	if (MB_SUCCESS != rval)
			return;
	for (int i = 0; i < size; i++) {
		//MBEntityHandle eh = io[i];
		if (tags[i]) {
			io[next++] = io[i];
		}
	}
	if (next<size)
		io.resize(next);
	return;
}

MBErrorCode contractionRegion(MBInterface * mb, MBEntityHandle v1,
		MBEntityHandle v2, std::vector<MBEntityHandle>& changed) {
	MBEntityHandle vlist[2];
	vlist[0] = v1;
	vlist[1] = v2;
	// it makes more sense to use vector, than range
	//  the entities could be very disjoint at some point
	// MBRange adj_ents;
	MBErrorCode rval = mb->get_adjacencies(vlist, 2, 2, false, changed,
			MBInterface::UNION);
	if (opts.useDelayedDeletion)
		filterValid(mb, changed);
	return rval;
}

//
int classifyVertex(MBInterface * mb, MBEntityHandle v1) {
	// return 0, 1, or 2 if the vertex is interior, on the border, or
	// border only (corner)
	// to do
	std::vector<MBEntityHandle> adjEdges;
	MBErrorCode rval = mb->get_adjacencies(&v1, 1, 1, false, adjEdges,
			MBInterface::UNION);
	if (MB_SUCCESS != rval)
		return 0; // interior??
	if (opts.useDelayedDeletion)
		filterValid(mb, adjEdges);
	int nBorder = 0;
	for (int i = 0; i < adjEdges.size(); i++) {
		MBEntityHandle edg = adjEdges[i];
		std::vector<MBEntityHandle> adjFaces;
		rval = mb->get_adjacencies(&edg, 1, 2, false, adjFaces,
				MBInterface::UNION);
		if (opts.useDelayedDeletion)
			filterValid(mb, adjFaces);
		if (adjFaces.size() == 1)
			nBorder++;
	}
	if (nBorder == 0)
		return 0;// everything interior
	else if (nBorder == adjEdges.size())
		return 2;
	else
		return 1; // some edges are interior
}

Vec3 getVec3FromMBVertex(MBInterface * mbi, MBEntityHandle v) {
	double c[3];
	mbi->get_coords(&v, 1, c);
	return Vec3(c[0], c[1], c[2]);
}
// every time we are getting the normal, we compute a new plane
// maybe we should store it??
//  No debate needed
// it will be much cheaper to store it, for "-m" option
// there, we will need it a lot
Plane trianglePlane(MBInterface * mb, MBEntityHandle tri) {
	// get connectivity of triangle
	const MBEntityHandle * conn;
	int num_nodes;
	MBErrorCode rval = mb->get_connectivity(tri, conn, num_nodes);
	assert(3==num_nodes && rval == MB_SUCCESS);
	Vec3 ve1 = getVec3FromMBVertex(mb, conn[0]);
	Vec3 ve2 = getVec3FromMBVertex(mb, conn[1]);
	Vec3 ve3 = getVec3FromMBVertex(mb, conn[2]);
	return Plane(ve1, ve2, ve3);

}

MBErrorCode contract(MBInterface * mb, MBEntityHandle v0, MBEntityHandle v1,
		Vec3 & vnew, std::vector<MBEntityHandle> & changed) {

	//
	//// Collect all the faces that are going to be changed
	////
	std::vector<MBEntityHandle> adj_entities;
	contractionRegion(mb, v0, v1, adj_entities);
	// those are all triangles that are affected
	// find also all edges that are affect
	MBEntityHandle vlist[2];
	vlist[0] = v0;
	vlist[1] = v1;
	// it makes more sense to use vector, than range
	//  the entities could be very disjoint at some point
	// MBRange adj_ents;
	std::vector<MBEntityHandle> edges;
	MBErrorCode rval = mb->get_adjacencies(vlist, 2, 1, false, edges,
			MBInterface::UNION);
	if (opts.useDelayedDeletion)
		filterValid(mb, edges);
	if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
		*opts.logfile << "Edges Adjacent:" << edges.size();
		for (int i = 0; i < edges.size(); i++)
			*opts.logfile << " " << mb->id_from_handle(edges[i]);
		*opts.logfile << std::endl;
	}
	// we have the edges and the triangles that are affected
	// 2 situations
	//   1) edge v0 v1 is existing
	//      we will delete edge (v0, v1), and triangles formed
	//       with edge (v0, v1)
	//   2) edge v0 v1 is not existing, but due to proximity
	//      only edges v2 v1 , v2, v0 will need to merge
	// more important is case 1)

	// first, find edge v0, v1
	MBEntityHandle ev0v1;
	int foundEdge = 0;
	for (int i = 0; i < edges.size(); i++) {
		MBEntityHandle e = edges[i];
		int nnodes;
		const MBEntityHandle * conn2;
		mb->get_connectivity(e, conn2, nnodes);
		if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
			*opts.logfile << "Edge: " << mb->id_from_handle(e) << " nnodes:"
					<< nnodes << " vertices:" << mb->id_from_handle(conn2[0])
					<< " " << mb->id_from_handle(conn2[1]) << std::endl;
		if ((conn2[0] == v0 && conn2[1] == v1) || (conn2[0] == v1 && conn2[1]
				== v0)) {
			foundEdge = 1;
			ev0v1 = e; // could be ev1v0, but who cares?
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
				*opts.logfile << "Edge found " << mb->id_from_handle(e)
						<< std::endl;
			break;
		}

	}
	// set the position of new vertices in vnew
	double newCoords[3];
	newCoords[0] = vnew[0];
	newCoords[1] = vnew[1];
	newCoords[2] = vnew[2];
	mb->set_coords(&v0, 1, newCoords);
	mb->set_coords(&v1, 1, newCoords);
	// although, vertex v1 will be deleted in the end; do we really need to set it?
	// yes, for merging purposes
	//

	if (opts.useDelayedDeletion) {
		// big copy from version 3512
		// although vertex v1 will be merged!!
		std::vector<MBEntityHandle> edgePairs; // the one that has v0 will
		// be kept
		std::vector<MBEntityHandle> edgesWithV1;
		if (foundEdge) {
			// this is case 1, the most complicated
			// get triangles connected to edge ev0v1
			std::vector<MBEntityHandle> tris;
			rval = mb->get_adjacencies(&ev0v1, 1, 2, false, tris,
					MBInterface::UNION);
			filterValid(mb, tris);
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
				*opts.logfile << "Triangles adjacent to found edge:"
						<< tris.size() << ":";
				for (int i = 0; i < tris.size(); i++)
					*opts.logfile << " " << mb->id_from_handle(tris[i]);
				*opts.logfile << std::endl;
			}
			for (int i = 0; i < tris.size(); i++) {
				MBEntityHandle triangleThatCollapses = tris[i];
				std::vector<MBEntityHandle> localEdges;
				rval = mb->get_adjacencies(&triangleThatCollapses, 1, 1, false,
						localEdges, MBInterface::UNION);
				filterValid(mb, localEdges);
				if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
					*opts.logfile << "Triangle " << mb->id_from_handle(
							triangleThatCollapses) << " Edges: "
							<< localEdges.size();

					for (int i = 0; i < localEdges.size(); i++)
						*opts.logfile << " " << mb->id_from_handle(
								localEdges[i]);
					*opts.logfile << std::endl;
				}
				// find the edges that contains v0
				MBEntityHandle e[2];// the 2 edges e0, e1, that are not ev0v1;
				if (localEdges.size() != 3)
					return MB_FAILURE; // failure
				int index = 0;
				for (int k = 0; k < 3; k++)
					if (localEdges[k] != ev0v1)
						e[index++] = localEdges[k];
				// among those 2 edges, find out which one has v0, and which one v1
				if (index != 2)
					return MB_FAILURE; // failure
				for (int j = 0; j < 2; j++) {
					int nn;
					const MBEntityHandle * conn2;
					mb->get_connectivity(e[j], conn2, nn);
					if (conn2[0] == v0 || conn2[1] == v0) {
						// this is the edge that will be kept, the other one collapsed
						edgePairs.push_back(e[j]);
						j = (j + 1) % 2;// the other one
						edgePairs.push_back(e[j]);
						edgesWithV1.push_back(e[j]);
						break; // no need to check the other one. it
						// will contain v1
					}
				}
			}
			// look at all triangles that are adjacent
			// invalidate first tris
			unsigned char invalid = 0;
			rval
					= mb->tag_set_data(validTag, &(tris[0]), tris.size(),
							&invalid);
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
				*opts.logfile << "Triangles invalidated: " << tris.size()
						<< ":";
				for (int i = 0; i < tris.size(); i++)
					*opts.logfile << " " << mb->id_from_handle(tris[i]);
				*opts.logfile << std::endl;
			}
			// then look at all triangles that are not in tris (adj_entities), and then
			// replace connectivity of v1 with v0
			std::vector<MBEntityHandle> trisChanged;
			for (int i = 0; i < adj_entities.size(); i++) {
				MBEntityHandle tr = adj_entities[i];
				if (!ehIsValid(tr))
					continue;
				// see the connectivity of tr; if there is a v1, replace it with v0
				// will that propagate to edges? or do we have to set it separately?
				int nn;
				const MBEntityHandle * conn3;
				mb->get_connectivity(tr, conn3, nn);

				assert(3==nn);
				for (int j = 0; j < 3; j++) {
					if (conn3[j] == v1) {
						// replace it with v0, and reset it
						MBEntityHandle connNew[3];
						connNew[0] = conn3[0];
						connNew[1] = conn3[1];
						connNew[2] = conn3[2];
						connNew[j] = v0;
						if (opts.logfile && opts.selected_output
								& OUTPUT_CONTRACTIONS) {
							std::vector<MBEntityHandle> localEdges;
							rval = mb->get_adjacencies(&tr, 1, 1, false,
									localEdges, MBInterface::UNION);
							filterValid(mb, localEdges);
							*opts.logfile << "Triangle t"
									<< mb->id_from_handle(tr)
									<< "  filtered : " << localEdges.size();
							for (int j = 0; j < localEdges.size(); j++)
								*opts.logfile << " e" << mb->id_from_handle(
										localEdges[j]);
							*opts.logfile << std::endl;
						}
						if (opts.logfile && opts.selected_output
								& OUTPUT_CONTRACTIONS) {
							*opts.logfile << "replace connectivity t"
									<< mb->id_from_handle(tr) << " v"
									<< mb->id_from_handle(conn3[0]) << " v"
									<< mb->id_from_handle(conn3[1]) << " v"
									<< mb->id_from_handle(conn3[2])
									<< "  to: v";
						}
						rval = mb->set_connectivity(tr, connNew, 3);
						if (opts.logfile && opts.selected_output
								& OUTPUT_CONTRACTIONS) {
							*opts.logfile << mb->id_from_handle(connNew[0])
									<< " v" << mb->id_from_handle(connNew[1])
									<< " v" << mb->id_from_handle(connNew[2])
									<< std::endl;
						}
						trisChanged.push_back(tr);
					}
				}

			}
			validFaceCount -= tris.size();
			rval = mb->tag_set_data(validTag, &ev0v1, 1, &invalid);
			// invalidate the edges connected for sure to v1
			rval = mb->tag_set_data(validTag, &edgesWithV1[0],
					edgesWithV1.size(), &invalid);
			// reset the connectivity of some edges (from v1 to v0)
			for (int i = 0; i < edges.size(); i++) {
				MBEntityHandle e1 = edges[i];
				if (!ehIsValid(e1)) // it could be the ones invalidated
					continue;

				int nn;
				const MBEntityHandle * conn;
				mb->get_connectivity(e1, conn, nn);

				assert(2==nn);
				for (int j = 0; j < 2; j++) {
					if (conn[j] == v1) {
						// replace it with v0, and reset it
						MBEntityHandle connNew[2];
						connNew[0] = conn[0];
						connNew[1] = conn[1];
						connNew[j] = v0;
						if (opts.logfile && opts.selected_output
								& OUTPUT_CONTRACTIONS) {
							*opts.logfile << "replace connectivity edge: "
									<< mb->id_from_handle(e1) << " "
									<< mb->id_from_handle(conn[0]) << " "
									<< mb->id_from_handle(conn[1]) << "  to: ";
						}
						rval = mb->set_connectivity(e1, connNew, 2);
						if (opts.logfile && opts.selected_output
								& OUTPUT_CONTRACTIONS) {
							*opts.logfile << mb->id_from_handle(connNew[0])
									<< " " << mb->id_from_handle(connNew[1])
									<< std::endl;
						}
					}
				}
			}
			// the question: is the adjacency between triangles and edges restored?
			// yes, it is; check set_connectivity logic
			// we need to remove adjacencies between triangles and edges that are not valid
			// and add some more adjacencies to the edges that are now part of the triangle
			for (int i = 0; i < trisChanged.size(); i++) {
				// get adjacencies now, and bail out early; maybe we should create if missing
				std::vector<MBEntityHandle> localEdges;
				MBEntityHandle tr = trisChanged[i];
				rval = mb->get_adjacencies(&tr, 1, 1, false, localEdges,
						MBInterface::UNION);
				filterValid(mb, localEdges);
				if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
					*opts.logfile << "Triangle t" << mb->id_from_handle(tr)
							<< "  filtered : " << localEdges.size();
					for (int j = 0; j < localEdges.size(); j++)
						*opts.logfile << " e" << mb->id_from_handle(
								localEdges[j]);
					*opts.logfile << std::endl;
				}
				assert(localEdges.size()==3);
			}

			filterValid(mb, adj_entities);
			changed = adj_entities; // deep copy
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
				*opts.logfile << "Triangles changed:" << changed.size();
				for (int i = 0; i < changed.size(); i++)
					*opts.logfile << " t" << mb->id_from_handle(changed[i]);
				*opts.logfile << std::endl;
			}
		} else // this should appear only when proximity limit > 0
		{
			// in this case, we only need to worry if vertex 0 and 1 are connected to the same
			// vertex 2; then we need to merge edges v0-v2 and v1 - v2
			// no triangles get deleted, only those edges;
			// the crack v0-v2-v1 gets seamed
			// get edges connected to vertex v0
			std::vector<MBEntityHandle> edges0;
			rval = mb->get_adjacencies(&v0, 1, 1, false, edges0,
					MBInterface::UNION);
			filterValid(mb, edges0);

			// get edges connected to vertex v1
			std::vector<MBEntityHandle> edges1;
			rval = mb->get_adjacencies(&v1, 1, 1, false, edges1,
					MBInterface::UNION);
			filterValid(mb, edges1);
			// find all edges that will be merged, of type v0-v2 v1-v2 (so that they have a
			// common vertex v2
			// in that case, we will have to merge them as before, and delete the
			// one that contains v1 before
			// for sure, there is no edge from v0 to v1 !!!
			// keep all the vertices from edges1, different from v1
			std::vector<MBEntityHandle> otherVertices1;
			int i = 0;
			for (i = 0; i < edges1.size(); i++) {
				MBEntityHandle edgeThatGoes = edges1[1];
				// find the other vertex, not v1
				int nn;
				const MBEntityHandle * conn2;
				mb->get_connectivity(edgeThatGoes, conn2, nn);
				int other = 0;
				if (conn2[0] == v1)
					other = 1;
				else if (conn2[1] == v1)// is this really necessary in a correct model?
					other = 0;
				otherVertices1.push_back(conn2[other]);
			}
			for (i = 0; i < edges0.size(); i++) {
				MBEntityHandle edgeThatStays = edges0[i];
				// find the other vertex, not v0
				int nn;
				const MBEntityHandle * conn2;
				mb->get_connectivity(edgeThatStays, conn2, nn);
				int other = 0;
				if (conn2[0] == v1)
					other = 1;
				else if (conn2[1] == v1)// is this really necessary in a correct model?
					other = 0;
				MBEntityHandle v2 = conn2[other];
				// let's see now if v2 is among vertices from otherVertices1
				for (int i1 = 0; i1 < otherVertices1.size(); i1++) {
					if (v2 == otherVertices1[i1]) {
						// we have a match, some work to do
						// invalidate the edge edges1[i1]
						unsigned char invalid = 0;
						mb->tag_set_data(validTag, &(edges1[i1]), 1, &invalid);
						break; // we stop looking for a match, only one possible
					}
				}
			}
			// triangles that need reconnected
			std::vector<MBEntityHandle> tri1;
			rval = mb->get_adjacencies(&v1, 1, 2, false, tri1,
					MBInterface::UNION);
			// start copy tri reconnect
			for (int i = 0; i < tri1.size(); i++) {
				MBEntityHandle tr = tri1[i];
				if (!ehIsValid(tr))
					continue;
				// see the connectivity of tr; if there is a v1, replace it with v0
				// will that propagate to edges? or do we have to set it separately?
				int nn;
				const MBEntityHandle * conn3;
				mb->get_connectivity(tr, conn3, nn);

				assert(3==nn);
				for (int j = 0; j < 3; j++) {
					if (conn3[j] == v1) {
						// replace it with v0, and reset it
						MBEntityHandle connNew[3];
						connNew[0] = conn3[0];
						connNew[1] = conn3[1];
						connNew[2] = conn3[2];
						connNew[j] = v0;
						rval = mb->set_connectivity(tr, connNew, 3);
						if (opts.logfile && opts.selected_output
								& OUTPUT_CONTRACTIONS) {
							*opts.logfile << mb->id_from_handle(connNew[0])
									<< " v" << mb->id_from_handle(connNew[1])
									<< " v" << mb->id_from_handle(connNew[2])
									<< std::endl;
						}
					}
				}

			}
			// now reconnect edges1 that are still valid
			for (int i = 0; i < edges1.size(); i++) {
				MBEntityHandle e1 = edges1[i];
				if (!ehIsValid(e1))
					continue;
				// see the connectivity of e1; if there is a v1, replace it with v0
				// will that propagate to edges? or do we have to set it separately?
				int nn;
				const MBEntityHandle * conn3;
				mb->get_connectivity(e1, conn3, nn);

				assert(2==nn);
				for (int j = 0; j < 2; j++) {
					if (conn3[j] == v1) {
						// replace it with v0, and reset it
						MBEntityHandle connNew[2];
						connNew[0] = conn3[0];
						connNew[1] = conn3[1];
						connNew[j] = v0;
						rval = mb->set_connectivity(e1, connNew, 2);

					}
				}

			}
			// all the triangles connected to v0 are now changed, need recomputation
			//
			rval = mb->get_adjacencies(&v0, 1, 2, false, changed,
					MBInterface::UNION);
			filterValid(mb, changed);

		}
		// end big copy from version 3512
	} else {
		// vertex v1 will be merged!!
		std::vector<MBEntityHandle> edgePairsToMerge; // the one that has v0 will
		// be kept
		if (foundEdge) {
			// this is case 1, the most complicated
			// get triangles connected to edge ev0v1
			std::vector<MBEntityHandle> tris;
			rval = mb->get_adjacencies(&ev0v1, 1, 2, false, tris,
					MBInterface::UNION);
			// find all edges that will be merged ( xv0, xv1, etc)
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
				*opts.logfile << "Triangles adjacent to found edge:"
						<< tris.size();
				for (int i = 0; i < tris.size(); i++)
					*opts.logfile << " " << mb->id_from_handle(tris[i]);
				*opts.logfile << std::endl;
			}
			for (int i = 0; i < tris.size(); i++) {
				MBEntityHandle triangleThatCollapses = tris[i];
				std::vector<MBEntityHandle> localEdges;
				rval = mb->get_adjacencies(&triangleThatCollapses, 1, 1, false,
						localEdges, MBInterface::UNION);
				if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
					*opts.logfile << "Triangle " << mb->id_from_handle(
							triangleThatCollapses) << " Edges: "
							<< localEdges.size();

					for (int i = 0; i < localEdges.size(); i++)
						*opts.logfile << " " << mb->id_from_handle(
								localEdges[i]);
					*opts.logfile << std::endl;
				}
				// find the edges that contains v0
				MBEntityHandle e[2];// the 2 edges e0, e1, that are not ev0v1;
				if (localEdges.size() != 3)
					return MB_FAILURE; // failure
				int index = 0;
				for (int k = 0; k < 3; k++)
					if (localEdges[k] != ev0v1)
						e[index++] = localEdges[k];
				// among those 2 edges, find out which one has v0, and which one v1
				if (index != 2)
					return MB_FAILURE; // failure
				for (int j = 0; j < 2; j++) {
					int nn;
					const MBEntityHandle * conn2;
					mb->get_connectivity(e[j], conn2, nn);
					if (conn2[0] == v0 || conn2[1] == v0) {
						// this is the edge that will be kept, the other one collapsed
						edgePairsToMerge.push_back(e[j]);
						j = (j + 1) % 2;// the other one
						edgePairsToMerge.push_back(e[j]);
						break; // no need to check the other one. it
						// will contain v1
					}
				}
			}
			// first merge vertices v0 and v1 : will also NOT delete v1 (yet)
			// the tag on v1 will be deleted too, and we do not want that, at least until
			// after the merging of edges, and deleting the pair
			rval = mb->merge_entities(v0, v1, false, false);
			// merge edgePairsToMerge // now, v0 and v1 should be collapsed!
			for (int j = 0; j < edgePairsToMerge.size(); j += 2) {
				// will also delete edges that contained v1 before
				if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
					*opts.logfile << "Edges merged:" << mb->id_from_handle(
							edgePairsToMerge[j]) << " " << mb->id_from_handle(
							edgePairsToMerge[j + 1]) << std::endl;
				mb->merge_entities(edgePairsToMerge[j],
						edgePairsToMerge[j + 1], false, true);

			}
			// the only things that need deleted are triangles
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
				*opts.logfile << "Triangles invalidated:" << tris.size();
				for (int i = 0; i < tris.size(); i++)
					*opts.logfile << " " << mb->id_from_handle(tris[i]);
				*opts.logfile << std::endl;
			}
			rval = mb->delete_entities(&(tris[0]), tris.size());
			validFaceCount -= tris.size();
			// hopefully, all adjacencies are preserved
			// delete now the edge ev0v1
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS)
				*opts.logfile << "Edge invalidated"
						<< mb->id_from_handle(ev0v1) << std::endl;
			rval = mb->delete_entities(&ev0v1, 1);
			// among adj_entities, remove the tris triangles, and return them
			for (int j = 0; j < adj_entities.size(); j++) {
				MBEntityHandle F = adj_entities[j];
				int inTris = 0;
				for (int k = 0; k < tris.size(); k++)
					if (F == tris[k]) {
						inTris = 1;
						break;
					}
				if (!inTris)
					changed.push_back(F);
			}
			if (opts.logfile && opts.selected_output & OUTPUT_CONTRACTIONS) {
				*opts.logfile << "Triangles changed:" << changed.size();
				for (int i = 0; i < changed.size(); i++)
					*opts.logfile << " " << mb->id_from_handle(changed[i]);
				*opts.logfile << std::endl;
			}

		} else // this should appear only when proximity limit > 0
		{
			// in this case, we only need to worry if vertex 0 and 1 are connected to the same
			// vertex 2; then we need to merge edges v0-v2 and v1 - v2
			// no triangles get deleted, only those edges
			// get edges connected to vertex v0
			std::vector<MBEntityHandle> edges0;
			rval = mb->get_adjacencies(&v0, 1, 1, false, edges0,
					MBInterface::UNION);

			// get edges connected to vertex v1
			std::vector<MBEntityHandle> edges1;
			rval = mb->get_adjacencies(&v1, 1, 1, false, edges1,
					MBInterface::UNION);
			// find all edges that will be merged, of type v0-v2 v1-v2 (so that they have a
			// common vertex v2
			// in that case, we will have to merge them as before, and delete the
			// one that contains v1 before
			// for sure, there is no edge from v0 to v1 !!!
			// keep all the vertices from edges1, different from v1
			std::vector<MBEntityHandle> otherVertices1;
			for (int i1 = 0; i1 < edges1.size(); i1++) {
				MBEntityHandle edgeThatGoes = edges1[i1];
				// find the other vertex, not v1
				int nn;
				const MBEntityHandle * conn2;
				mb->get_connectivity(edgeThatGoes, conn2, nn);
				int other = 0;
				if (conn2[0] == v1)
					other = 1;
				else if (conn2[1] == v1)// is this really necessary in a correct model?
					other = 0;
				otherVertices1.push_back(conn2[other]);
			}
			for (int i = 0; i < edges0.size(); i++) {
				MBEntityHandle edgeThatStays = edges0[i];
				// find the other vertex, not v0
				int nn;
				const MBEntityHandle * conn2;
				mb->get_connectivity(edgeThatStays, conn2, nn);
				int other = 0;
				if (conn2[0] == v1)
					other = 1;
				else if (conn2[1] == v1)// is this really necessary in a correct model?
					other = 0;
				MBEntityHandle v2 = conn2[other];
				// let's see now if v2 is among vertices from otherVertices1
				// if yes, then we have a match, edges that need collapsed
				for (int i1 = 0; i1 < otherVertices1.size(); i1++) {
					if (v2 == otherVertices1[i1]) {
						// we have a match, some work to do
						edgePairsToMerge.push_back(edgeThatStays);
						edgePairsToMerge.push_back(edges1[i1]);
						break; // we stopp looking for a match, only one possible
					}
				}
			}

			// first merge vertices v0 and v1 : will also NOT delete v1 (yet)
			// the tag on v1 will be deleted too, and we do not want that, at least until
			// after the merging of edges, and deleting the pair
			rval = mb->merge_entities(v0, v1, false, false);
			// merge edgePairsToMerge // now, v0 and v1 should be collapsed!
			for (int j = 0; j < edgePairsToMerge.size(); j += 2) {
				// will also delete edges that contained v1 before
				mb->merge_entities(edgePairsToMerge[j],
						edgePairsToMerge[j + 1], false, true);
			}
			// all the triangles connected to v0 are now changed, need recomputation
			//
			rval = mb->get_adjacencies(&v0, 1, 2, false, changed,
					MBInterface::UNION);

		}
	}
	return MB_SUCCESS;

}

