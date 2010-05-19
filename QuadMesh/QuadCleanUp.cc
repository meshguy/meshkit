#include "QuadCleanUp.h"
#include <sstream>
#include <assert.h>

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
void set_diamond_tag(Mesh *mesh, vector<FaceClose> &diamonds)
{
   size_t numfaces = mesh->getSize(2);
   for( size_t i = 0; i <  numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        face->setTag( 0  );
   }

   set<Face*>::const_iterator it;
   for( size_t i = 0; i < diamonds.size(); i++) {
        for( it = diamonds[i].shield.begin(); it != diamonds[i].shield.end(); ++it){
	     Face *face = *it;
	     face->setTag(1);
        }
   }

   for( size_t i = 0; i < diamonds.size(); i++) {
        Face *face = diamonds[i].face;
        face->setTag(2);
   }
}

////////////////////////////////////////////////////////////////////

void set_doublet_tag(Mesh *mesh, vector<Doublet> &doublet)
{
/*
   size_t numfaces = mesh->getSize(2);
   for( size_t i = 0; i <  numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        face->setTag( 0  );
   }

   set<Face*>::const_iterator it;
   for( size_t i = 0; i < diamonds.size(); i++) {
        for( it = diamonds[i].shield.begin(); it != diamonds[i].shield.end(); ++it){
	     Face *face = *it;
	     face->setTag(1);
        }
   }

   for( size_t i = 0; i < diamonds.size(); i++) {
        Face *face = diamonds[i].face;
        face->setTag(2);
   }
*/
}


////////////////////////////////////////////////////////////////////

bool FaceClose :: isSafe() const 
{
   if( vertex0->isBoundary() ) return 0;
   if( vertex1->isBoundary() ) return 0;

   assert( face->hasNode( vertex0 ) );
   assert( face->hasNode( vertex1 ) );

   Vertex *v0 = face->getNodeAt(0);
   Vertex *v1 = face->getNodeAt(1);
   Vertex *v2 = face->getNodeAt(2);
   Vertex *v3 = face->getNodeAt(3);

   vector<Face*> neighs;

   neighs = v0->getRelations2();
   for (size_t j = 0; j < neighs.size(); j++)
       if( neighs[j]->isVisited() ) return 0;

   // Every neighbour must share only one diamond vertex.
   for (size_t j = 0; j < neighs.size(); j++) {
        if( neighs[j] != face ) {
	    int val0 = neighs[j]->hasNode(vertex0);
	    int val1 = neighs[j]->hasNode(vertex1);
	    if( val0 + val1 == 2 ) return 0;
        }
   }

   neighs = v1->getRelations2();
   for (size_t j = 0; j < neighs.size(); j++)
       if( neighs[j]->isVisited() ) return 0;

   // Every neighbour must share only one diamond vertex.
   for (size_t j = 0; j < neighs.size(); j++) {
        if( neighs[j] != face ) {
	    int val0 = neighs[j]->hasNode(vertex0);
	    int val1 = neighs[j]->hasNode(vertex1);
	    if( val0 + val1 == 2 ) return 0;
        }
   }

   neighs = v2->getRelations2();
   for (size_t j = 0; j < neighs.size(); j++)
       if( neighs[j]->isVisited() ) return 0;

   // Every neighbour must share only one diamond vertex.
   for (size_t j = 0; j < neighs.size(); j++) {
        if( neighs[j] != face ) {
	    int val0 = neighs[j]->hasNode(vertex0);
	    int val1 = neighs[j]->hasNode(vertex1);
	    if( val0 + val1 == 2 ) return 0;
        }
   }

   neighs = v3->getRelations2();
   for (size_t j = 0; j < neighs.size(); j++)
       if( neighs[j]->isVisited() ) return 0;

   // Every neighbour must share only one diamond vertex.
   for (size_t j = 0; j < neighs.size(); j++) {
        if( neighs[j] != face ) {
	    int val0 = neighs[j]->hasNode(vertex0);
	    int val1 = neighs[j]->hasNode(vertex1);
	    if( val0 + val1 == 2 ) return 0;
        }
   }

   return 1;
}
////////////////////////////////////////////////////////////////////
void FaceClose :: makeShield() 
{
   Vertex *v0 = face->getNodeAt(0);
   Vertex *v1 = face->getNodeAt(1);
   Vertex *v2 = face->getNodeAt(2);
   Vertex *v3 = face->getNodeAt(3);

   vector<Face*> neighs;
   shield.clear();

   neighs = vertex0->getRelations2(); assert( !neighs.empty() );
   for (size_t j = 0; j < neighs.size(); j++) {
       neighs[j]->setVisitMark(1);
       shield.insert( neighs[j] );
   }

   neighs = vertex1->getRelations2(); assert( !neighs.empty() );
   for (size_t j = 0; j < neighs.size(); j++) {
       neighs[j]->setVisitMark(1);
       shield.insert( neighs[j] );
   }
}

////////////////////////////////////////////////////////////////////

bool Doublet::isSafe() const 
{
  assert( vertex  );

  if( vertex->isRemoved() ) return 0;

  vector<Face*> apexfaces = vertex->getRelations2();
  for (size_t i = 0; i < apexfaces.size(); i++) {
        Face *face = apexfaces[i];
	if( face->isRemoved() ) return 0;
	if( face->isVisited() ) return 0;
  }
  return 1;
}

////////////////////////////////////////////////////////////////////

void Doublet::makeShield() 
{
  vector<Face*> apexfaces = vertex->getRelations2();
  for (size_t i = 0; i < apexfaces.size(); i++) {
        Face *face = apexfaces[i];
	face->setVisitMark(1);
  }
}

////////////////////////////////////////////////////////////////////

void
YRing::create_boundary() {
    vector<Vertex*> inner_nodes;
    inner_nodes.reserve(7);

    vector<Face*> apexfaces = apex->getRelations2();
    for (size_t i = 0; i < apexfaces.size(); i++) {
        Face *face = apexfaces[i];
        for (int j = 0; j < 4; j++) {
            Vertex *v = face->getNodeAt(j);
            if (std::find(inner_nodes.begin(), inner_nodes.end(), v) == inner_nodes.end())
                inner_nodes.push_back(v);
        }
    }

    assert(inner_nodes.size() == 7);

    sort(inner_nodes.begin(), inner_nodes.end());
    boundary.clear();

    boundary.reserve(14);
    for (size_t i = 0; i < faces.size(); i++) {
        Face *face = faces[i];
        for (int j = 0; j < 4; j++) {
            Vertex *v0 = face->getNodeAt((j + 0) % 4);
            Vertex *v1 = face->getNodeAt((j + 1) % 4);
            if (std::binary_search(inner_nodes.begin(), inner_nodes.end(), v0)) continue;
            if (std::binary_search(inner_nodes.begin(), inner_nodes.end(), v1)) continue;
            Edge newedge(v0, v1);
            boundary.push_back(newedge);
        }
    }
}

////////////////////////////////////////////////////////////////////

bool
YRing::isValid() {
    if (faces.size() != 13) return 0;

    std::map<Vertex*, set<Face*> > relations02;
    for (size_t iface = 0; iface < faces.size(); iface++) {
        Face *face = faces[iface];
        if (face->getSize(0) != 4) return 0;
        for (int j = 0; j < face->getSize(0); j++) {
            Vertex *vertex = face->getNodeAt(j);
            relations02[vertex].insert(face);
        }
    }

    // We need exact four corners in the patch.
    std::set<Vertex*> corners;
    std::map<Vertex*, set<Face*> >::const_iterator it;
    for (it = relations02.begin(); it != relations02.end(); ++it) {
        set<Face*> neigh = it->second;
        if (neigh.size() == 1) corners.insert(it->first);
    }
    if (corners.size() != 4) return 0;

    // Create the boundary of the patch...
    create_boundary();

    // Sequence the chain and start from one of the corner...
    int err;
    err = Mesh::make_chain(boundary);
    if (err) return 0;

    err = Mesh::rotate_chain(boundary, *corners.begin());
    if (err) return 0;

    err = Mesh::is_closed_chain(boundary);
    if (err == 0) return 0;

    int nSize = boundary.size();
    boundnodes.resize(nSize);

    for (int i = 0; i < nSize; i++)
        boundnodes[i] = boundary[i].getNodeAt(0);

    int gridSize[] = {0, 0, 0, 0};
    int lastpos = 0;
    int currside = 0;
    for (int j = 1; j < nSize; j++) {
        if (corners.find(boundnodes[j]) != corners.end()) {
            gridSize[currside++] = j - lastpos + 1;
            lastpos = j;
        }
    }
    gridSize[3] = nSize - lastpos + 1;

    for (int i = 0; i < 4; i++) {
        if (gridSize[i] < 4 || gridSize[i] > 5) {
            return 0;
        }
    }
    cout << " Valid Pattern "
            << gridSize[0] << " "
            << gridSize[1] << " "
            << gridSize[2] << " "
            << gridSize[3] << endl;

    return 1;
}

////////////////////////////////////////////////////////////////////

void
QuadCleanUp::initialize_wavefront() {
    size_t numnodes = mesh->getSize(0);

    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
            vertex->setConstrainedMark(1);
        else
            vertex->setConstrainedMark(0);
    }
}

////////////////////////////////////////////////////////////////////

vector<Vertex*>
QuadCleanUp::next_front_nodes() const {
    set<Vertex*> vset;
    vector<Vertex*> neighs;

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        if (!vertex->isConstrained()) {
            neighs = vertex->getRelations0();
            for (size_t j = 0; j < neighs.size(); i++) {
                if (neighs[j]->isConstrained()) {
                    vset.insert(vertex);
                    break;
                }
            }
        }
    }
    neighs.clear();
    if (vset.empty()) return neighs;

    neighs.resize(vset.size());
    set<Vertex*>::const_iterator it;

    int index = 0;
    for (it = vset.begin(); it != vset.end(); ++it)
        neighs[index++] = *it;
    return neighs;

}

////////////////////////////////////////////////////////////////////

vector<int>
QuadCleanUp::getVertexFaceDegrees() {
    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    size_t numnodes = mesh->getSize(0);

    vector<int> degree(numnodes);
    vector<Face*> neighs;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNodeAt(i);
        neighs = v->getRelations2();
        degree[i] = neighs.size();
    }

    int mindegree = *min_element(degree.begin(), degree.end());
    int maxdegree = *max_element(degree.begin(), degree.end());

    cout << "Topological Information :" << endl;
    cout << "Min Vertex-Face Degree : " << mindegree << endl;
    cout << "Max Vertex-Face Degree : " << maxdegree << endl;
    cout << "Node-Face Degree  |  "  <<  " Count " << endl;
    for (int i = mindegree; i <= maxdegree; i++) {
        int ncount = 0;
        for (size_t j = 0; j < degree.size(); j++)
            if (degree[j] == i)
                ncount++;
        cout <<  i <<  "        " << ncount << endl;
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    return degree;
}

////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::get_VertexOf_FaceDegree(int n) {
    size_t numnodes = mesh->getSize(0);
    vector<Face*> neighs;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNodeAt(i);
        neighs = v->getRelations2();
        if (neighs.size() == size_t(n))
            return v;
    }
    return NULL;
}

////////////////////////////////////////////////////////////////////

vector<FaceClose>
QuadCleanUp::search_diamonds(bool both_sides, bool allow_boundary_faces) {
    vDiamonds.clear();

    size_t numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    if (!allow_boundary_faces)
        mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    for (size_t iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
        face->setVisitMark(0);
    }

    vector<Face*>  neighs;
    Vertex *vertex0, *vertex1;

    for (size_t iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
            Vertex *v0 = face->getNodeAt(0);
            Vertex *v1 = face->getNodeAt(1);
            Vertex *v2 = face->getNodeAt(2);
            Vertex *v3 = face->getNodeAt(3);

            neighs = v0->getRelations2();
            size_t d0 = neighs.size();

            neighs = v1->getRelations2();
            size_t d1 = neighs.size();

            neighs = v2->getRelations2();
            size_t d2 = neighs.size();

            neighs = v3->getRelations2();
            size_t d3 = neighs.size();

            bool enlist = 0;
            int boundCount = 0;
            if (allow_boundary_faces) {
                boundCount += v0->isBoundary();
                boundCount += v1->isBoundary();
                boundCount += v2->isBoundary();
                boundCount += v3->isBoundary();
            }

            if (boundCount == 1) {
                enlist = 1;
                if (v0->isBoundary()) {
                    vertex0 = v1;
                    vertex1 = v3;
                }
                if (v1->isBoundary()) {
                    vertex0 = v2;
                    vertex1 = v0;
                }
                if (v2->isBoundary()) {
                    vertex0 = v1;
                    vertex1 = v3;
                }
                if (v3->isBoundary()) {
                    vertex0 = v0;
                    vertex1 = v2;
                }
            }

            if (boundCount == 0) {
                if (both_sides) {
                    if ((d0 == 3 && d2 == 3)) {
                        enlist = 1;
                        vertex0 = v0;
                        vertex1 = v2;
                    } else if (d1 == 3 && d3 == 3) {
                        enlist = 1;
                        vertex0 = v1;
                        vertex1 = v3;
                    }
                } else {
                    if ((d0 == 3 || d2 == 3) || (d1 == 3 || d3 == 3)) {
                        if (d0 + d2 < d1 + d3) {
                            enlist = 1;
                            vertex0 = v0;
                            vertex1 = v2;
                        } else {
                            enlist = 1;
                            vertex0 = v1;
                            vertex1 = v3;
                        }
                    }
                }
            }

            if (enlist) {
	        FaceClose diamond(mesh, face, vertex0, vertex1);
		if( diamond.isSafe() ) {
		    diamond.makeShield();
                    vDiamonds.push_back(diamond);
                }
            }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    return vDiamonds;
}

/////////////////////////////////////////////////////////////////////////////

vector<Bridge>
QuadCleanUp::search_bridges(bool allow_boundary_nodes) {
    //
    // This module identifies independent set of bridges. A bridge is defined
    // as an edge whoose end points have three neighbouring cells.
    // So if an edge is selected as a bridge, then all its neighbours are locked,
    // therefore, they cann't be part of another bridge. Therefore, operations
    // on bridge can executed independently..
    //
    // A Bridge is surrounded by four quads and enclosed by six edges and six
    // nodes.
    //
    vBridges.clear();

    size_t numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    for (size_t iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
        face->setVisitMark(0);
        face->setTag(0);
    }

    Bridge bridge;
    vector<Face*> diamonds, neighs;
    for (size_t iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
        if (face->getSize(0) == 4) {
            Vertex *v0 = face->getNodeAt(0);
            Vertex *v1 = face->getNodeAt(1);
            Vertex *v2 = face->getNodeAt(2);
            Vertex *v3 = face->getNodeAt(3);

            if (v0->isBoundary())
                continue;
            if (v1->isBoundary())
                continue;
            if (v2->isBoundary())
                continue;
            if (v3->isBoundary())
                continue;

            int visitCount = 0;
            int hitBoundary = 0;

            neighs = v0->getRelations2();
            for (size_t j = 0; j < neighs.size(); j++)
                if (neighs[j]->hasBoundaryNode()) hitBoundary = 1;
            if (hitBoundary) continue;

            for (size_t j = 0; j < neighs.size(); j++)
                visitCount += neighs[j]->isVisited();
            if (visitCount) continue;

            int d0 = neighs.size();

            neighs = v1->getRelations2();
            for (size_t j = 0; j < neighs.size(); j++)
                if (neighs[j]->hasBoundaryNode()) hitBoundary = 1;
            if (hitBoundary) continue;

            for (size_t j = 0; j < neighs.size(); j++)
                visitCount += neighs[j]->isVisited();
            if (visitCount) continue;
            int d1 = neighs.size();

            neighs = v2->getRelations2();
            for (size_t j = 0; j < neighs.size(); j++)
                if (neighs[j]->hasBoundaryNode()) hitBoundary = 1;
            if (hitBoundary) continue;

            for (size_t j = 0; j < neighs.size(); j++)
                visitCount += neighs[j]->isVisited();
            if (visitCount)
                continue;
            int d2 = neighs.size();

            neighs = v3->getRelations2();
            for (size_t j = 0; j < neighs.size(); j++)
                if (neighs[j]->hasBoundaryNode()) hitBoundary = 1;
            if (hitBoundary) continue;

            for (size_t j = 0; j < neighs.size(); j++)
                visitCount += neighs[j]->isVisited();
            if (visitCount)
                continue;
            int d3 = neighs.size();

            bool enlist = 0;

            if (enlist == 0 && (d0 == 3 && d1 == 3)) {
                enlist = 1;
                bridge.vertex0 = v0;
                bridge.vertex1 = v1;
                vBridges.push_back(bridge);
            }

            if (enlist == 0 && (d1 == 3 && d2 == 3)) {
                enlist = 1;
                bridge.vertex0 = v1;
                bridge.vertex1 = v2;
                vBridges.push_back(bridge);
            }

            if (enlist == 0 && (d2 == 3 && d3 == 3)) {
                enlist = 1;
                bridge.vertex0 = v2;
                bridge.vertex1 = v3;
                vBridges.push_back(bridge);
            }

            if (enlist == 0 && (d3 == 3 && d0 == 3)) {
                enlist = 1;
                bridge.vertex0 = v3;
                bridge.vertex1 = v0;
                vBridges.push_back(bridge);
            }

            if (enlist) {
                neighs = bridge.vertex0->getRelations2();
                for (size_t j = 0; j < neighs.size(); j++) {
                    neighs[j]->setVisitMark(1);
                    neighs[j]->setTag(1);
                }
                neighs = bridge.vertex1->getRelations2();
                for (size_t j = 0; j < neighs.size(); j++) {
                    neighs[j]->setVisitMark(1);
                    neighs[j]->setTag(1);
                }

            }

        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    cout << "Info: Number of Bridges in the mesh " << vBridges.size() << endl;
    return vBridges;
}

/////////////////////////////////////////////////////////////////////////////

vector<YRing>
QuadCleanUp::search_yrings() {
    //
    // This module identifies independent set of bridges. A bridge is defined
    // as an edge whoose end points have three neighbouring cells.
    // So if an edge is selected as a bridge, then all its neighbours are locked,
    // therefore, they cann't be part of another bridge. Therefore, operations
    // on bridge can executed independently..
    //
    // A Bridge is surrounded by four quads and enclosed by six edges and six
    // nodes.
    //
    vyRings.clear();


    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    for (size_t inode = 0; inode < numnodes; inode++) {
        Vertex *vertex = mesh->getNodeAt(inode);
        vertex->setVisitMark(0);
    }

    for (size_t iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
        face->setVisitMark(0);
        face->setTag(0);
    }

    vector<Face*> neighs1, neighs2;
    set<Face*> ringneighs;

    set<Face*>::const_iterator iter;

    YRing newring;

    int index = 0;
    for (size_t inode = 0; inode < numnodes; inode++) {
        Vertex *vertex = mesh->getNodeAt(inode);

        if (vertex->isBoundary()) continue;

        neighs1 = vertex->getRelations2();

        if (neighs1.size() != 3) continue;

        ringneighs.clear();
        for (size_t j = 0; j < neighs1.size(); j++) {
            neighs2 = neighs1[j]->getRelations202();
            for (size_t k = 0; k < neighs2.size(); k++)
                ringneighs.insert(neighs2[k]);
        }

        if (ringneighs.size() != 13) continue;

        int enlist = 1;
        for (iter = ringneighs.begin(); iter != ringneighs.end(); ++iter) {
            Face *face = *iter;
            if (face->hasBoundaryNode() || face->isVisited()) enlist = 0;
            for (int j = 0; j < face->getSize(0); j++) {
                Vertex *vertex = face->getNodeAt(j);
                if (vertex->isVisited()) {
                    enlist = 0;
                    break;
                }
            }
            if (enlist == 0) break;
        }

        if (enlist) {
            newring.apex = vertex;
            newring.faces.clear();
            for (iter = ringneighs.begin(); iter != ringneighs.end(); ++iter) {
                Face *face = *iter;
                newring.faces.push_back(face);
            }

            if (newring.isValid()) {
                index++;
                for (iter = ringneighs.begin(); iter != ringneighs.end(); ++iter) {
                    Face *face = *iter;
                    face->setVisitMark(1);
                    face->setTag(index);
                    for (int j = 0; j < face->getSize(0); j++) {
                        Vertex *vertex = face->getNodeAt(j);
                        vertex->setVisitMark(1);
                    }
                }
                vyRings.push_back(newring);
            }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    cout << "Info: Number of YRings in the mesh " << vyRings.size() << endl;
    return vyRings;
}

/////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_yring(const YRing &aring) {
    //
    // The original idea of Y Ring came from Guy Bunin's paper " Non-Local
    // Topological Clean-Up. But I have modified the algorithm slightly.
    // Input Mesh:
    //           #Faces :  13    #Nodes : 7
    // Output original algorithm:
    //           #Faces :  10    #Nodes : 4
    // This implementation
    //           #Faces :  12    #Nodes : 6
    //
    //Therefore, the complexity of this implementation is close to the
    //input mesh. Whether it is a good idea or not, can be verified by the
    //experiements and evaluation.
    //
    // We will try to reuse objects as much as possible i.e  only one face
    // and one node will be marked for remove.
    //
    //

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::remove_yrings(bool recursive) {
    search_yrings();
}

/////////////////////////////////////////////////////////////////////////////

vector<Face*>
QuadCleanUp::search_flat_quads() {
    //  Public Function ...
    size_t numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    vector<Face*> flatQ, neighs;

    int edgefaces[4];
    for (size_t iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
        assert(face);
        if (face->getSize(0) == 4) {
            int boundary = 1;
            for (int j = 0; j < 4; j++) {
                Vertex *vb = face->getNodeAt(j);
                if (!vb->isBoundary()) {
                    boundary = 0;
                    break;
                }
            }

            if (boundary) {
                int bound_edges = 0;
                for (int j = 0; j < 4; j++) {
                    Vertex *v0 = face->getNodeAt((j + 0) % 4);
                    Vertex *v1 = face->getNodeAt((j + 1) % 4);
                    neighs = Mesh::getRelations112(v0, v1);
                    if (neighs.size() == 1)
                        bound_edges++;
                    edgefaces[j] = neighs.size();
                }

                if (bound_edges == 3) {
                    /*
                     Point3D v1 = make_vector( neighs[0], node );
                     Point3D v2 = make_vector( neighs[1], node );
                     double  angle = getVectorAngle(v1,v2);
                     if( angle > cutOffAngle )
                     degree2nodes.push_back(node);
                     flatQ.push_back(face);
                     */
                }

            }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    cout << "Number of flat Quads " << flatQ.size() << endl;
    return flatQ;
}

///////////////////////////////////////////////////////////////////////////////

bool QuadCleanUp :: isDoublet(const Vertex *v)
{
   if (!v->isBoundary() && (v->getRelations2().size() == 2)) return 1;
   return 0;
}
///////////////////////////////////////////////////////////////////////////////

vector<Doublet>
QuadCleanUp::search_interior_doublets() {
    //
    // An interior doublet is a vertex, which is shared by two
    // face neighbours. They are undesirables in the quadmesh.
    // as it would mean the angle is 180 between some adjacent
    // edges...

    // This module finds doublets in the interiour mesh...

    size_t numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    size_t numfaces = mesh->getSize(2);

    for( int i = 0; i < numfaces; i++) {
         Face *face = mesh->getFaceAt(i);
	 face->setVisitMark(0);
    }

    vDoublets.clear();

    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        if( isDoublet( vertex ) ) {
            Doublet newdoublet(mesh, vertex);
	    if( newdoublet.isSafe() ) {
	        newdoublet.makeShield();
	        vDoublets.push_back(newdoublet);
	    }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    cout << "Number of interior doublets Detected : " << vDoublets.size() << endl;
    return vDoublets;
}

//////////////////////////////////////////////////////////////////////////

vector<Vertex*>
QuadCleanUp::search_boundary_singlets() {
    //
    // A boundary singlet is a vertex which is shared by only one face.
    // They are undesirables in the quad mesh as that would mean large
    // angle on some of the edges..
    //
    // For the flat singlet ( angle closer to 180 degree ). it is easy
    // to remove the neighbouring quad from the mesh.
    //

    size_t numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    vector<Vertex*> singlets;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNodeAt(i);
        if (v->isBoundary() && (v->getRelations2().size() == 1))
            singlets.push_back(v);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    cout << "Number of boundary singlets Detected : " << singlets.size() << endl;
    return singlets;
}

////////////////////////////////////////////////////////////////////

int
FaceClose::apply()
{
    assert(face->getSize(0) == 4);

    Vertex *v0 = vertex0;
    Vertex *v2 = vertex1;

    //
    // We shouldn't modify the boundary, therefore skip this operation. It is
    // advisable to close the face, if the diagonal distance between (v0,v2)
    // is much less than the diagonal distance between (v1,v3).
    //
    if (v0->isBoundary() || v2->isBoundary())
        return 1;

    int vpos = -1;
    for (int i = 0; i < 4; i++) {
        if (face->getNodeAt(i) == v0) {
            vpos = i;
            break;
        }
    }
    assert(vpos >= 0);
    if (face->getNodeAt((vpos + 2) % 4) != v2) {
        cout << "Warning: Face-open requires opposite vertices " << endl;
        cout << "Debug  : Face is : " << face->getNodeAt(0)->getID() << " "
                << face->getNodeAt(1)->getID() << " "
                << face->getNodeAt(2)->getID() << " "
                << face->getNodeAt(3)->getID() << endl;
        cout << "Opposite ends are " << v0->getID() << " " << v2->getID() << endl;
        return 1;
    }

    Vertex *v1 = face->getNodeAt((vpos + 1) % 4);
    Vertex *v3 = face->getNodeAt((vpos + 3) % 4);

    // Add a new vertex in the mesh...
    Point3D p3d = Vertex::mid_point(v0, v2);
    replacedNode  = Vertex::newObject();
    replacedNode->setXYZCoords(p3d);
    mesh->addNode( replacedNode );

    // This face will go away, remove this from vertex relationship.
    v0->removeRelation2(face);
    v1->removeRelation2(face);
    v2->removeRelation2(face);
    v3->removeRelation2(face);

    vector<Face*> vr0 = v0->getRelations2();
    vector<Face*> vr1 = v1->getRelations2();
    vector<Face*> vr2 = v2->getRelations2();
    vector<Face*> vr3 = v3->getRelations2();

    for (size_t i = 0; i < vr0.size(); i++) {
         Face *face = vr0[i];
	 assert( face->isValid() );
    }

    for (size_t i = 0; i < vr1.size(); i++) {
         Face *face = vr1[i];
	 assert( face->isValid() );
    }

    for (size_t i = 0; i < vr2.size(); i++) {
         Face *face = vr2[i];
	 assert( face->isValid() );
    }

    for (size_t i = 0; i < vr3.size(); i++) {
         Face *face = vr3[i];
	 assert( face->isValid() );
    }

    // Affected neighbouring faces must replace the v0/v2 the new
    // vertex and the new vertex has new neighbours.
    for (size_t i = 0; i < vr0.size(); i++) {
        assert( !vr0[i]->isRemoved() );
        vr0[i]->replaceNode(v0, replacedNode );
        replacedNode->addRelation2(vr0[i]);
    }

    for (size_t i = 0; i < vr2.size(); i++) {
        assert( !vr2[i]->isRemoved() );
	assert(  vr2[i]->hasNode(v2 ) );
        vr2[i]->replaceNode(v2, replacedNode);
        replacedNode->addRelation2(vr2[i]);
	assert( !vr2[i]->hasNode(v2 ) );
    }

    if( mesh->getAdjTable(0,0) ) {
        vector<Vertex*> rel0;
        rel0 = vertex0->getRelations0();
        for( size_t i = 0; i < rel0.size(); i++) {
             replacedNode->addRelation0( rel0[i] );
             rel0[i]->addRelation0( replacedNode );
	     rel0[i]->removeRelation0( vertex0 );
        }

       rel0 = vertex1->getRelations0();
       for( size_t i = 0; i < rel0.size(); i++) {
            replacedNode->addRelation0( rel0[i] );
            rel0[i]->addRelation0( replacedNode );
	    rel0[i]->removeRelation0( vertex1 );
       }
    }

    // Two nodes and face go away from the mesh..
    v0->setRemoveMark(1);
    v2->setRemoveMark(1);
    face->setRemoveMark(1);

    for (size_t i = 0; i < vr0.size(); i++) {
         Face *face = vr0[i];
	 assert( face->isValid() );
    }

    for (size_t i = 0; i < vr1.size(); i++) {
         Face *face = vr1[i];
	 assert( face->isValid() );
    }

    for (size_t i = 0; i < vr2.size(); i++) {
         Face *face = vr2[i];
	 assert( face->isValid() );
    }

    for (size_t i = 0; i < vr3.size(); i++) {
         Face *face = vr3[i];
	 assert( face->isValid() );
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::set_regular_node_tag() {
    size_t numnodes = mesh->getSize(0);
    int relexist = mesh->build_relations(0, 2);

    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->getRelations2().size() == 4)
            vertex->setTag(0);
        else
            vertex->setTag(1);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////////////////
int
QuadCleanUp::clean_layer_once(int layernum ) 
{
    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    for( size_t i = 0; i < numfaces; i++) {
         Face *face = mesh->getFaceAt(i);
	 face->setVisitMark(0);
	 assert( !face->isRemoved() );
    }

    // Here we are using a narrow domain search for doublets for efficiency
    // purpose and therefore don't call the member functions which searches
    // doublets in the entire mesh...
    //
    while(1) {

    vDoublets.clear();
    for( size_t i = 0; i < numnodes; i++) {
         Vertex *vertex = mesh->getNodeAt(i);
	 int lid = vertex->getLayerID();
	 if( lid == layernum-1 || lid == layernum || lid == layernum+1) {
             if( isDoublet(vertex) ) {
	         Doublet newdoublet(mesh, vertex);
		 if( newdoublet.isSafe() ) {
		     newdoublet.makeShield();
		     vDoublets.push_back(newdoublet);
                 }
             }
	 }
    }

    if( vDoublets.empty() )  break;

    remove_doublets_once();

    int layerid[4];
    for( size_t i = 0; i < vDoublets.size(); i++) {
         Face *face = vDoublets[i].replacedFace;
	 if( face ) {
	    for( int j = 0; j < 4; j++) 
	         layerid[j] = face->getNodeAt(j)->getLayerID();
            int minid = *min_element( layerid, layerid + 4 );
            face->setLayerID(minid);
         }
    }
    }

    // At this stage, the mesh must be free of any doublets.

    vDiamonds.clear();

    vector<Face*> neighs;
    Vertex *currnode, *oppnode;

    for (size_t i = 0; i < numnodes; i++) {
        Vertex *currnode = mesh->getNodeAt(i);
        int id = currnode->getLayerID();
        if (id == layernum) {
            neighs = currnode->getRelations2();
            for (size_t j = 0; j < neighs.size(); j++) {
                Face *face = neighs[j];
                if (face->getLayerID() == layernum + 1) {
                    int k = face->queryNodeAt(currnode);
                        Vertex *vside0 = face->getNodeAt((k + 1) % 4);
                        Vertex *vside1 = face->getNodeAt((k + 3) % 4);
                        if ( (vside0->getLayerID() == layernum + 1) && 
			     (vside1->getLayerID() == layernum + 1) &&
			     (!vside0->isBoundary()  && !vside1->isBoundary() ) ) 
		        {
                             FaceClose diamond(mesh, face,vside0, vside1);
			     if( diamond.isSafe() ) {
				 diamond.makeShield();
                                 vDiamonds.push_back(diamond);
                                 break;
                             }
                        }
		}
             }
         }
    }

    static int niter = 0;
    if( niter == 0) {
        Jaal::set_layer_tag( mesh );
        mesh->saveAs("debug0");
    }

    int ncount = remove_diamonds_once();

    if( niter == 0) {
        Jaal::set_layer_tag( mesh );
        mesh->saveAs("debug1");
	exit(0);
    }
    niter++;

    for( size_t i = 0; i < vDiamonds.size(); i++) {
         Vertex *vertex = vDiamonds[i].replacedNode;
	 if( vertex ) 
	     vertex->setLayerID( layernum + 1);
    }

    return ncount;
}

////////////////////////////////////////////////////////////////////////////////
int 
QuadCleanUp::clean_layer(int layernum ) 
{
   if( mesh->isHomogeneous() != 4 ) {
        cout << "Warning: layer cleaning only for quad mesh " << endl;
	return 1;
   }

   int relexist2 = mesh->build_relations(0, 2);
   int ncount = 0;

   assert( !mesh->hasDuplicates(2) );

   mesh->setWavefront(0);
   mesh->setWavefront(2);

   assert( mesh->getAdjTable(0,2) );

   static int id = 0;
   while(1) {
       int n = clean_layer_once( layernum );
       if( mesh->hasDuplicates(2) ) {
           cout << "Error: Mesh has Duplicates " << endl;
	   exit(0);
       }
       if( n == 0) break;
       ncount += n;
       id++;
  }

  cout << " CSV " << endl;
  Jaal::set_layer_tag(mesh);
  mesh->saveAs("debug");

  exit(0);

  if (!relexist2) mesh->clear_relations(0, 2);
  return ncount;
}
////////////////////////////////////////////////////////////////////////////////
void QuadCleanUp :: advancing_front_cleanup()
{
   search_and_remove_doublets( 1 );

   int numlayers = mesh->setWavefront(2);

   for( int i = 0; i < numlayers; i++) 
        clean_layer(i);
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face, Vertex *v0, Vertex *v2) {

    //  Public Function ...
    Point3D p3d = Vertex::mid_point(v0, v2);

    Vertex *doublet = Vertex::newObject();
    doublet->setXYZCoords(p3d);

    Vertex *o1 = NULL, *o2 = NULL;

    vector<Vertex*> connect = face->getNodes();
    for (size_t i = 0; i < connect.size(); i++) {
        if (connect[i] == v0) {
            o1 = connect[(i + 1) % 4];
            o2 = connect[(i + 3) % 4];
            break;
        }
    }

    mesh->addNode(doublet);

    connect[0] = doublet;
    connect[1] = v0;
    connect[2] = o1;
    connect[3] = v2;

    Face *newquad1 = new Face;
    newquad1->setNodes(connect);
    mesh->addFace(newquad1);

    face->setRemoveMark(1);

    connect[0] = doublet;
    connect[1] = v2;
    connect[2] = o2;
    connect[3] = v0;

    Face *newquad2 = new Face;
    newquad2->setNodes(connect);
    mesh->addFace(newquad2);

    if (mesh->getAdjTable(0, 2)) {
        v0->removeRelation2(face);
        o1->removeRelation2(face);
        v2->removeRelation2(face);
        o2->removeRelation2(face);

        v0->addRelation2(newquad1);
        v0->addRelation2(newquad2);

        v2->addRelation2(newquad1);
        v2->addRelation2(newquad2);

        doublet->addRelation2(newquad1);
        doublet->addRelation2(newquad2);

        o1->addRelation2(newquad1);
        o2->addRelation2(newquad2);
    }

    return doublet;
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face) {

    //  Public Function ...
    if (face->getSize(0) != 4)
        return NULL;
    Vertex *v0 = face->getNodeAt(0);
    Vertex *v2 = face->getNodeAt(2);
    return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_boundary_doublet(Face *face) {

    //  Public Function ...
    if (!face->isBoundary())
        return NULL;

    if (face->getSize(0) != 4)
        return NULL;

    int ncount = 0;
    for (int i = 0; i < 4; i++) {
        Vertex *v = face->getNodeAt(i);
        ncount += v->isBoundary();
    }

    if (ncount != 3)
        return NULL;

    Vertex *v0 = NULL, *v2 = NULL;
    for (int i = 0; i < 4; i++) {
        Vertex *v = face->getNodeAt(i);
        if (!v->isBoundary()) {
            v0 = face->getNodeAt((i + 0) % 4);
            v2 = face->getNodeAt((i + 2) % 4);
        }
    }

    return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::diamond_collapse(FaceClose &closeface) {

    //  Private function ...

    Face *face = closeface.face;
    if (face->getSize(0) != 4)
        return 1;

    assert(mesh->getAdjTable(0, 2));

    Vertex *v0 = face->getNodeAt(0);
    Vertex *v1 = face->getNodeAt(1);
    Vertex *v2 = face->getNodeAt(2);
    Vertex *v3 = face->getNodeAt(3);

    vector<Face*> vr0 = v0->getRelations2();
    for (size_t i = 0; i < vr0.size(); i++)
        if (vr0[i]->isRemoved())
            return 1;

    vector<Face*> vr1 = v1->getRelations2();
    for (size_t i = 0; i < vr1.size(); i++)
        if (vr1[i]->isRemoved())
            return 1;

    vector<Face*> vr2 = v2->getRelations2();
    for (size_t i = 0; i < vr2.size(); i++)
        if (vr2[i]->isRemoved())
            return 1;

    vector<Face*> vr3 = v3->getRelations2();
    for (size_t i = 0; i < vr3.size(); i++)
        if (vr3[i]->isRemoved())
            return 1;

   closeface.apply();

   return 0;
}

////////////////////////////////////////////////////////////////////

int
Doublet ::remove() {

    if (vertex->isBoundary())
        return 1;

    vector<Face*> neighs = vertex->getRelations2();

    if (neighs.size() != 2)
        return 1;

    assert( neighs[0] != neighs[1] );

    Vertex *d1 = NULL, *d2 = NULL, *o1 = NULL, *o2 = NULL;

    vector<Vertex*> connect = neighs[0]->getNodes();
    if (connect.size() != 4)
        return 1;

    if (neighs[0]->isRemoved())
        return 2;
    if (neighs[1]->isRemoved())
        return 2;

    for (size_t i = 0; i < connect.size(); i++) {
        if (connect[i] == vertex) {
            d1 = connect[(i + 1) % 4];
            o1 = connect[(i + 2) % 4];
            d2 = connect[(i + 3) % 4];
            break;
        }
    }

    connect = neighs[1]->getNodes();
    if (connect.size() != 4)
        return 1;

    for (size_t i = 0; i < connect.size(); i++) {
        if (connect[i] == vertex) {
            o2 = connect[(i + 2) % 4];
            break;
        }
    }

    assert(d1);
    assert(d2);
    assert(o1);
    assert(o2);
    assert( o1 != o2 );

    d1->removeRelation2(neighs[0]);
    d1->removeRelation2(neighs[1]);

    d2->removeRelation2(neighs[0]);
    d2->removeRelation2(neighs[1]);

    o1->removeRelation2(neighs[0]);
    o2->removeRelation2(neighs[1]);

    connect[0] = d1;
    connect[1] = o1;
    connect[2] = d2;
    connect[3] = o2;

    Face *newquad = new Face;
    newquad->setNodes(connect);
    mesh->addFace(newquad);

    d1->addRelation2(newquad);
    d2->addRelation2(newquad);

    o1->addRelation2(newquad);
    o2->addRelation2(newquad);

    neighs[0]->setRemoveMark(1);
    neighs[1]->setRemoveMark(1);
    vertex->setRemoveMark(1);

    return 0;
}

///////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_boundary_singlet(Vertex *vertex) {
    if (!vertex->isBoundary())
        return 1;

    vector<Face*> neighs = vertex->getRelations2();

    if (neighs.size() != 1)
        return 1;

    Vertex *d1 = NULL, *d2 = NULL, *o1 = NULL;

    vector<Vertex*> connect = neighs[0]->getNodes();
    if (connect.size() != 4)
        return 1;

    if (neighs[0]->isRemoved())
        return 2;

    for (size_t i = 0; i < connect.size(); i++) {
        if (connect[i] == vertex) {
            d1 = connect[(i + 1) % 4];
            o1 = connect[(i + 2) % 4];
            d2 = connect[(i + 3) % 4];
            break;
        }
    }

    if (o1->isBoundary())
        return 3;

    /*
     Point3D p3d = vertex->getXYZCoords();
     o1->setXYZCoords(p3d);
     o1->setBoundaryMark( vertex->getBoundaryMark() );

     d1->removeRelation2(neighs[0]);
     d2->removeRelation2(neighs[0]);
     o1->removeRelation2(neighs[0]);

     neighs[0]->setRemoveMark(1);
     vertex->setRemoveMark(1);
     */

    return 0;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridge(const Bridge &bridge) {
    //
    // Presently this code does not check the convexity of the hexagon
    // which is split into two quadrilateral. So in some sense, this
    // is not a robust code and probably smoothing might solve the
    // problem, but we can not guarantee.
    //
    // One workaround will be to ensure all the enclosing vertices
    // are on a circle.
    //
    Vertex *v0 = bridge.vertex0;
    Vertex *v1 = bridge.vertex1;

    if (v0->isRemoved() || v1->isRemoved()) return 1;

    vector<Face*> neighs = Mesh::getRelations102(v0, v1);

    for (size_t i = 0; i < neighs.size(); i++)
        if (neighs[i]->isRemoved()) return 1;

    assert(neighs.size() == 4);

    // Create a closed chain of bounding edges ...

    list<Edge> boundedges;
    for (size_t i = 0; i < neighs.size(); i++) {
        Face *face = neighs[i];
        for (int j = 0; j < 4; j++) {
            Vertex *ev0 = face->getNodeAt((j + 0) % 4);
            Vertex *ev1 = face->getNodeAt((j + 1) % 4);
            if (ev0 == v0 || ev0 == v1) continue;
            if (ev1 == v0 || ev1 == v1) continue;
            Edge edge(ev0, ev1);
            boundedges.push_back(edge);
        }
    }
    assert(boundedges.size() == 6);

    // Create a closed chain of bounding nodess ...

    vector<Vertex*> chain_nodes(6);

    Edge edge = boundedges.front();
    boundedges.pop_front();
    chain_nodes[0] = edge.getNodeAt(0);
    chain_nodes[1] = edge.getNodeAt(1);

    Vertex *curr_vertex = chain_nodes[1];

    int index = 2;
    list<Edge>::iterator it;
    for (int i = 0; i < 4; i++) {
        for (it = boundedges.begin(); it != boundedges.end(); ++it) {
            edge = *it;
            if (edge.getNodeAt(0) == curr_vertex) {
                curr_vertex = edge.getNodeAt(1);
                chain_nodes[index++] = curr_vertex;
                break;
            }
            if (edge.getNodeAt(1) == curr_vertex) {
                curr_vertex = edge.getNodeAt(0);
                chain_nodes[index++] = curr_vertex;
                break;
            }
        }

        if (it != boundedges.end()) boundedges.erase(it);
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++)
            chain_nodes[i]->removeRelation2(neighs[j]);
    }

    //
    // Only for the convex polygons, we can be sure about joining two diagonal.
    // May be there are some other elegant ways, but for now, just try to adjust
    // coordinates of all the six vertices so that they lie on the circle/sphere.
    //
    Point3D pCenter = Vertex::mid_point(v0, v1);
    v0->setConstrainedMark(1);
    v1->setConstrainedMark(1);

    double radialdist[6];
    for (int i = 0; i < 6; i++)
        radialdist[i] = Math::length2(pCenter, chain_nodes[i]->getXYZCoords());
    sort(radialdist, radialdist + 6);
    double radius2 = radialdist[4]; // Bias towards large length.
    double radius = sqrt(radius2);

    // Keep these nodes fixed and smooth the mesh.
    for (int j = 0; j < 6; j++)
        chain_nodes[j]->setConstrainedMark(1);

    Point3D oldPos, newPos, posVec;
    double alpha;
    // Try bringing the nodes closer the circle/sphere boundary

    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 6; j++) {
            oldPos = chain_nodes[j]->getXYZCoords();
            double d = Math::length2(pCenter, oldPos);
            if (d / radius2 > 1.10) {
                // Attract the node towards center.
                alpha = 0.90;
                newPos[0] = (1.0 - alpha) * pCenter[0] + alpha * oldPos[0];
                newPos[1] = (1.0 - alpha) * pCenter[1] + alpha * oldPos[1];
                newPos[2] = (1.0 - alpha) * pCenter[2] + alpha * oldPos[2];
            } else if (d / radius2 < 0.90) {
                // Repeal the node from the center.
                alpha = 1.10;
                newPos[0] = (1.0 - alpha) * pCenter[0] + alpha * oldPos[0];
                newPos[1] = (1.0 - alpha) * pCenter[1] + alpha * oldPos[1];
                newPos[2] = (1.0 - alpha) * pCenter[2] + alpha * oldPos[2];
            } else {
                //  Force the node to come on the boundary...
                posVec[0] = oldPos[0] - pCenter[0];
                posVec[1] = oldPos[1] - pCenter[1];
                posVec[2] = oldPos[2] - pCenter[2];
                alpha = radius / Math::magnitude(posVec);
                newPos[0] = pCenter[0] + alpha * posVec[0];
                newPos[1] = pCenter[1] + alpha * posVec[1];
                newPos[2] = pCenter[2] + alpha * posVec[2];
            }
            chain_nodes[j]->setXYZCoords(newPos);
        }
        // Carry out laplacian smoothing 5 times (arbitrary choice).
        laplacian_smoothing(mesh, chain_nodes, 5);
    }

    // Unmark the nodes..
    for (int j = 0; j < 6; j++)
        chain_nodes[j]->setConstrainedMark(0);

    // Hopefully by now the enclosing polygon is convex. Now add two
    // quadrlaterals along two diagonally min-degree vertices.

    int degree[3];

    for (int i = 0; i < 3; i++) {
        int d0 = chain_nodes[(i + 0) % 6]->getRelations2().size();
        int d1 = chain_nodes[(i + 3) % 6]->getRelations2().size();
        degree[i] = d0 + d1;
    }

    int minat = 0;
    int mindegree = degree[0];

    for (int i = 0; i < 3; i++) {
        if (degree[i] < mindegree) {
            minat = i;
            mindegree = degree[i];
        }
    }

    // Finally update the data strucutures...
    minat = 0;
    vector<Vertex*> qConnect(4);
    for (int i = 0; i < 4; i++)
        qConnect[i] = chain_nodes[(minat + i) % 6];
    Face *quad1 = new Face;
    quad1->setNodes(qConnect);
    for (int i = 0; i < 4; i++)
        qConnect[i]->addRelation2(quad1);
    mesh->addFace(quad1);

    for (int i = 0; i < 4; i++)
        qConnect[i] = chain_nodes[(minat + 3 + i) % 6];
    Face *quad2 = new Face;
    quad2->setNodes(qConnect);
    for (int i = 0; i < 4; i++)
        qConnect[i]->addRelation2(quad2);
    mesh->addFace(quad2);

    for (size_t i = 0; i < neighs.size(); i++)
        neighs[i]->setRemoveMark(1);
    v0->setRemoveMark(1);
    v1->setRemoveMark(1);

    return 0;
}
////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridges_once(bool allow_boundary_nodes) 
{
    // These relationship are used in "remove_bridge" function call
    // and it is optimization step to call them here.
    //
    int rel0exist = mesh->build_relations(0, 0);
    int rel2exist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    search_bridges(allow_boundary_nodes);

    int ncount = 0;
    for (size_t i = 0; i < vBridges.size(); i++) {
        int err = remove_bridge(vBridges[i]);
        if (!err) ncount++;
    }

    if (ncount) {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);

    cout << "Info: number of bridges removed " << ncount << endl;

    return ncount;
}

////////////////////////////////////////////////////////////////////

void
QuadCleanUp::remove_bridges(bool recursive, bool allow_boundary_nodes) {
    int ncount = remove_bridges_once(allow_boundary_nodes);

    if (recursive) {
        while (1) {
            ncount = remove_bridges_once(allow_boundary_nodes);
            if (ncount == 0) break;
        }
    }
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds_once() 
{
    int rel0exist = mesh->build_relations(0, 0);
    int rel2exist = mesh->build_relations(0, 2);

    vector<Vertex*> vertexQ;
    for (size_t i = 0; i < vDiamonds.size(); i++) {
        Vertex *v0 = vDiamonds[i].vertex0;
        Vertex *v1 = vDiamonds[i].vertex1;
        v0->setConstrainedMark(1);
        v1->setConstrainedMark(1);
        vertexQ.push_back(v0);
        vertexQ.push_back(v1);
    }

    Point3D p3d0, p3d1;
    for (size_t j = 0; j < 4; j++) {
        for (size_t i = 0; i < vDiamonds.size(); i++) {
            Vertex *v0 = vDiamonds[i].vertex0;
            Vertex *v1 = vDiamonds[i].vertex1;
            p3d0 = Vertex::mid_point(v0, v1, 0.25);
            p3d1 = Vertex::mid_point(v1, v0, 0.25);
            v0->setXYZCoords(p3d0);
            v1->setXYZCoords(p3d1);
        }
        laplacian_smoothing(mesh, vertexQ, 5);
    }

    for (size_t i = 0; i < vDiamonds.size(); i++) {
        Vertex *v0 = vDiamonds[i].vertex0;
        Vertex *v1 = vDiamonds[i].vertex1;
        v0->setConstrainedMark(0);
        v1->setConstrainedMark(0);
    }

    int err, ncount = 0;
    for (size_t i = 0; i < vDiamonds.size() ; i++) {
        err = diamond_collapse(vDiamonds[i]);
        if (!err) ncount++;
    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);

    if (ncount) {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
    }

    return ncount;
}

////////////////////////////////////////////////////////////////////

void
QuadCleanUp::remove_diamonds(bool recursive, bool both_sides,
        bool allow_boundary_faces) {

    int rel0exist = mesh->build_relations(0, 0);
    int rel2exist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    search_diamonds(both_sides, allow_boundary_faces);
    int ncount = remove_diamonds_once();

    if (recursive) {
        while (1) {
            search_diamonds(both_sides, allow_boundary_faces);
            ncount = remove_diamonds_once();
            if (ncount == 0) break;
        }
    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_interior_doublets_once() 
{
    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    vector<Doublet> doublets = search_interior_doublets();

    int ncount = 0;
    for (size_t i = 0; i < doublets.size(); i++) {
        int err = doublets[i].remove();
        if (!err)
            ncount++;
    }

    if (ncount) {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    return ncount;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_doublets_once() 
{
    int ncount = 0;
    for (size_t i = 0; i < vDoublets.size(); i++) {
        int err = vDoublets[i].remove();
        if (!err)
            ncount++;
    }

    if (ncount) {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
    }

    return ncount;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::search_and_remove_doublets_once() 
{
    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();
    vector<Doublet> doublets = search_interior_doublets();

    int ncount = remove_doublets_once();

    if (!relexist)
        mesh->clear_relations(0, 2);

    return ncount;
}

////////////////////////////////////////////////////////////////////

void
QuadCleanUp::search_and_remove_doublets(bool recursive) 
{
    int ncount = search_and_remove_doublets_once();

    if (recursive) {
        while (1) {
            ncount = search_and_remove_doublets_once();
            if (ncount == 0) break;
        }
    }
}

////////////////////////////////////////////////////////////////////

void
QuadCleanUp::cleanup_internal_boundary_face() {
    int relexist = mesh->build_relations(0, 2);

    size_t numfaces = mesh->getSize(2);

    vector<Face*> boundfaces;

    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        if (face->hasBoundaryNode()) {
            int nfnodes = face->getSize(0);
            Vertex *v0 = NULL;
            Vertex *v1 = NULL;
            for (int j = 0; j < nfnodes; j++) {
                Vertex *v = face->getNodeAt(j);
                if (v->isBoundary()) {
                    v0 = face->getNodeAt((j + 1) % nfnodes);
                    v1 = face->getNodeAt((j + 3) % nfnodes);
                    if (!v0->isBoundary() && !v1->isBoundary()) {
                        FaceClose closeface(mesh, face,v0,v1);
                        break;
                    }
                }
            }
        }
    }

    mesh->prune();

    if (!relexist)
        mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////

void
QuadCleanUp::cleanup_boundary(double cutOffAngle) {
    mesh->search_boundary();

    ///////////////////////////////////////////////////////////////
    // First try to handle flat node...
    ///////////////////////////////////////////////////////////////
    //  vector<Vertex*> flatnodes = search_flat_nodes();
    //  cleanup_internal_boundary_face() ;

    int relexist = mesh->build_relations(0, 0);

    vector<Vertex*> degree2nodes;

    size_t numnodes = mesh->getSize(0);

    Vertex* node, *onode;
    vector<Vertex*> neighs;
    for (size_t i = 0; i < numnodes; i++) {
        node = mesh->getNodeAt(i);
        if (node->isBoundary()) {
            neighs = node->getRelations0();
            if (neighs.size() == 2) {
                if (neighs[0]->isBoundary() && neighs[1]->isBoundary()) {
                    Point3D v1 = make_vector(neighs[0], node);
                    Point3D v2 = make_vector(neighs[1], node);
                    double angle = Math::getVectorAngle(v1, v2);
                    if (angle > cutOffAngle)
                        degree2nodes.push_back(node);
                }
            }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 0);

    relexist = mesh->build_relations(0, 2);

    Face *boundface;
    vector<Face*> faceneighs;
    for (size_t i = 0; i < degree2nodes.size(); i++) {
        node = degree2nodes[i];
        faceneighs = node->getRelations2();
        if (faceneighs.size() == 1) {
            boundface = faceneighs[0];
            if (boundface->getSize(0) == 4) {
                int j = boundface->queryNodeAt(node);
                onode = boundface->getNodeAt((j + 2) % 4);
                if (!onode->isBoundary())
                    insert_doublet(boundface, node, onode);
            }
        }
    }

    mesh->prune();
    mesh->enumerate(0);
    mesh->enumerate(2);

    if (!relexist)
        mesh->clear_relations(0, 2);

    mesh->setBoundaryStatus(0);
}

///////////////////////////////////////////////////////////////////////
