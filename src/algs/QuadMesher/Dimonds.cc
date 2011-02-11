#include "QuadCleanUp.h"
#include <sstream>
#include <assert.h>

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

/*
bool
SwapEdge::is_topologically_valid_swap(int d1, int d2, int d3, int d4)
{
    if (d1 < 4 || d2 < 4) return 0;
    if ((d1 > 4) && (d2 > 4) && (d3 == 3) && (d4 == 3)) return 1;
    if ((d1 == 5) && (d2 == 5) && (d3 == 3) && (d4 == 4)) return 1;
    if ((d1 == 5) && (d2 == 5) && (d3 == 4) && (d4 == 3)) return 1;
    if (max(d1, d2) > max(d3 + 1, d4 + 1)) return 1;
    return 0;
}
 */
///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_diamond_tag(Mesh *mesh)
{
    QuadCleanUp qClean(mesh);
    vector<Diamond> diamonds = qClean.search_diamonds();

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(1);
    }

    set<Face*>::const_iterator it;
    for (size_t i = 0; i < diamonds.size(); i++)
    {
        Face *face = diamonds[i].face;
        face->setTag(0);
    }
}

////////////////////////////////////////////////////////////////////

void
Jaal::set_doublet_tag(Mesh *mesh)
{
    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (QuadCleanUp::isDoublet(vertex))
        {
            vertex->setTag(0);
        }
        else
            vertex->setTag(1);
    }

    if (relexist)
        mesh->clear_relations(0, 2);
}

//////////////////////////////////////////////////////////////////////

void
Jaal::set_bridge_tag(Mesh *mesh)
{
    int relexist = mesh->build_relations(0, 2);

    QuadCleanUp qClean(mesh);
    vector<Edge33> bridges = qClean.search_edges33();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        v->setTag(1);
    }

    for (size_t i = 0; i < bridges.size(); i++)
    {
        Vertex *v0 = bridges[i].connect[0];
        Vertex *v1 = bridges[i].connect[1];
        if ((v0 != NULL) && (v1 != NULL))
        {
            v0->setTag(0);
            v1->setTag(0);
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);
}

//////////////////////////////////////////////////////////////////////

void
Jaal::set_singlet_tag(Mesh *mesh)
{
    size_t numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 2);

    if (!mesh->isBoundaryKnown());
    mesh->search_boundary();

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (QuadCleanUp::isSinglet(vertex))
        {
            vertex->setTag(0);
        }
        else
            vertex->setTag(1);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////

bool
Doublet::isSafe() const
{
    assert(vertex);

    if (vertex->isRemoved()) return 0;

    FaceSequence apexfaces = vertex->getRelations2();
    for (size_t i = 0; i < apexfaces.size(); i++)
    {
        Face *face = apexfaces[i];
        if (face->isRemoved()) return 0;
        if (face->isVisited()) return 0;
    }
    return 1;
}

////////////////////////////////////////////////////////////////////

void
Doublet::makeShield()
{
    FaceSequence apexfaces = vertex->getRelations2();
    for (size_t i = 0; i < apexfaces.size(); i++)
    {
        Face *face = apexfaces[i];
        face->setVisitMark(1);
    }
}

///////////////////////////////////////////////////////////////////////////////

int QuadCleanUp::has_interior_nodes_degree_345()
{
    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    size_t numnodes = mesh->getSize(0);

    vector<int> degree(numnodes);
    vector<Face*> neighs;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        neighs = v->getRelations2();
        if (neighs.size() < 3 || neighs.size() > 5) return 0;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

bool
QuadCleanUp::isDiamond(const Face *face, int &pos, int type)
{
    pos = -1;
    const Vertex *v0 = face->getNodeAt(0);
    const Vertex *v1 = face->getNodeAt(1);
    const Vertex *v2 = face->getNodeAt(2);
    const Vertex *v3 = face->getNodeAt(3);

    FaceSequence neighs;

    neighs = v0->getRelations2();
    size_t d0 = neighs.size();

    neighs = v1->getRelations2();
    size_t d1 = neighs.size();

    neighs = v2->getRelations2();
    size_t d2 = neighs.size();

    neighs = v3->getRelations2();
    size_t d3 = neighs.size();

    // Boundary Cases ...
    if (v0->isBoundary() || v2->isBoundary())
    {
        if (!v1->isBoundary() && !v3->isBoundary())
        {
            if (d1 == 3 && d3 == 3)
            {
                pos = 1;
                return 1;
            }
        }
    }

    if (v1->isBoundary() || v3->isBoundary())
    {
        if (!v0->isBoundary() && !v2->isBoundary())
        {
            if ((d0 == 3 && d2 == 3))
            {
                pos = 0;
                return 1;
            }
        }
    }

    if (v0->isBoundary()) return 0;
    if (v1->isBoundary()) return 0;
    if (v2->isBoundary()) return 0;
    if (v3->isBoundary()) return 0;

    if ((d0 == 3 && d2 == 3))
    {
        pos = 0;
        return 1;
    }


    if (d1 == 3 && d3 == 3)
    {
        pos = 1;
        return 1;
    }

    if (d1 == 5 && d3 == 5)
    {
        if ((d0 == 3 && d2 == 4) || (d0 == 4 && d2 == 3))
        {
            pos = 0;
            return 1;
        }
    }

    if (d0 == 5 && d2 == 5)
    {
        if ((d1 == 3 && d3 == 4) || (d1 == 4 && d3 == 3))
        {
            pos = 1;
            return 1;
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Diamond::remove()
{
    cout << " Hello CSV " << endl;
    abort();
    assert(face);
    assert(mesh);
    assert(vertex0);
    assert(vertex2);

    if (face->isRemoved()) return 1;
    if (vertex0->isRemoved()) return 1;
    if (vertex2->isRemoved()) return 1;

    vector<Face*> neighs;
    neighs = vertex0->getRelations2();
    if (neighs.size() != 3) return 2;

    neighs = vertex2->getRelations2();
    if (neighs.size() != 3) return 2;

    FaceClose op(mesh, face, vertex0, vertex2);
    int err = op.remove();

    return err;
}


////////////////////////////////////////////////////////////////////

vector<Diamond>
QuadCleanUp::search_diamonds(int type)
{
    if (!mesh->isPruned()) mesh->prune();

    size_t numfaces = mesh->getSize(2);

    mesh->search_boundary();

    int relexist = mesh->build_relations(0, 2);

    int pos;

    assert(mesh->isBoundaryKnown());

    vDiamonds.clear();
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        if (isDiamond(face, pos))
        {
            Diamond diamond(mesh, face, pos);
            vDiamonds.push_back(diamond);
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    if (vDiamonds.size())
        cout << "Info: number of diamonds detected " << vDiamonds.size() << endl;

    return vDiamonds;
}

/////////////////////////////////////////////////////////////////////////////

vector<Edge33>
QuadCleanUp::search_edges33_in_layer(int layerid)
{
    int relexist0 = mesh->build_relations(0, 0);

    vector<Edge*> medges = mesh->getEdges(0);

    mesh->search_boundary();

    vEdges33.clear();
    for (size_t iedge = 0; iedge < medges.size(); iedge++)
    {
        Edge *edge = medges[iedge];
        Vertex *v0 = edge->getNodeAt(0);
        Vertex *v1 = edge->getNodeAt(1);
        int l0 = v0->getLayerID();
        int l1 = v1->getLayerID();
        if (l0 < layerid - 2 || l0 > layerid + 2) continue;
        if (l1 < layerid - 2 || l1 > layerid + 2) continue;
        if (isEdge33(edge))
        {
            Edge33 bridge(mesh, v0, v1);
            vEdges33.push_back(bridge);
        }
    }

    // We don't need edge pointers anymore
    for (size_t iedge = 0; iedge < medges.size(); iedge++)
        delete medges[iedge];

    if (!relexist0)
        mesh->clear_relations(0, 0);

    return vEdges33;
}

/////////////////////////////////////////////////////////////////////////////

vector<Edge33>
QuadCleanUp::search_edges33()
{
    vEdges33.clear();

    int relexist0 = mesh->build_relations(0, 0);

    vector<Edge*> medges = mesh->getEdges(0);

    mesh->search_boundary();

    for (size_t iedge = 0; iedge < medges.size(); iedge++)
    {
        Edge *edge = medges[iedge];
        if (isEdge33(edge))
        {
            Vertex *v0 = edge->getNodeAt(0);
            Vertex *v1 = edge->getNodeAt(1);
            Edge33 bridge(mesh, v0, v1);
            vEdges33.push_back(bridge);
        }
        delete edge; // No more needed..
    }

    if (!relexist0)
        mesh->clear_relations(0, 0);

    return vEdges33;
}

/////////////////////////////////////////////////////////////////////////////

int
Edge33::build_boundary()
{
    assert(!connect[0]->isBoundary());
    assert(!connect[1]->isBoundary());

    adjFaces = Mesh::getRelations102(connect[0], connect[1]);
    assert(adjFaces.size() == 4);

    // Create a closed chain of bounding edges ...
    vector<Edge> boundedges;
    for (size_t i = 0; i < adjFaces.size(); i++)
    {
        Face *face = adjFaces[i];
        for (int j = 0; j < 4; j++)
        {
            Vertex *ev0 = face->getNodeAt((j + 0) % 4);
            Vertex *ev1 = face->getNodeAt((j + 1) % 4);
            if (ev0 == connect[0] || ev0 == connect[1]) continue;
            if (ev1 == connect[0] || ev1 == connect[1]) continue;
            Edge edge(ev0, ev1);
            boundedges.push_back(edge);
        }
    }

    assert(boundedges.size() == 6);
    Mesh::make_chain(boundedges);

    // Create a closed chain of bounding nodess ...
    bound_nodes = Mesh::chain_nodes(boundedges);
    assert(bound_nodes.size() == 6);

    // It is important that the first node of is attached to one of
    // the two nodes of the bridge.
    assert(mesh->getAdjTable(0, 0));
    vector<Vertex*> vneighs = connect[0]->getRelations0();
    assert(vneighs.size() == 3);
    vector<Vertex*> rotated(6);
    for (int i = 0; i < 3; i++)
    {
        if (vneighs[i] != connect[1])
        {
            int pos = -1;
            Vertex *start_vertex = vneighs[i];
            for (int j = 0; j < 6; j++)
            {
                if (bound_nodes[j] == start_vertex)
                {
                    pos = j;
                    break;
                }
            }
            assert(pos >= 0);
            for (int j = 0; j < 6; j++)
                rotated[j] = bound_nodes[ (pos + j) % 6 ];
            bound_nodes = rotated;
        }
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

int
Edge33::build()
{
    Vertex *v0 = connect[0];
    Vertex *v1 = connect[1];

    if (v0->isVisited() || v1->isVisited()) return 1;

    if (v0->getRelations2().size() != 3) return 1;
    if (v1->getRelations2().size() != 3) return 1;

    assert(mesh->getAdjTable(0, 2));

    // Our edges are assumed to be simple i.e. shared by at the most
    // two faces.
    adjFaces = Mesh::getRelations112(v0, v1);
    if (adjFaces.size() != 2) return 1;

    // Check if the bridge touches boundary.
    boundary = 0;
    if (adjFaces[0]->hasBoundaryEdge())
        boundary = 1;

    if (adjFaces[1]->hasBoundaryEdge())
        boundary = 1;

    int err;
    err = build_boundary();
    if (err) return 2;

    // If boundary nothing is done here.
    if (boundary) return 0;

    // Build the outer shell of the bridge. Must contain closed
    // contour with six line segments. There are many ways to
    // convert a hexagon into quads ( mostly 2, 3 or 4 quads).
    // For general case, it is a hard problem.
    //

    int ncount[6];
    for (int i = 0; i < 6; i++)
    {
        vector<Face*> vfaces = bound_nodes[i]->getRelations2();
        if (vfaces.size() == 3) return 3;
        ncount[i] = vfaces.size();
    }

    if (ncount[0] == 4 && ncount[3] == 4)
    {
        err = Face::hexagon_2_quads(bound_nodes, newFaces, 0);
        if (err) return 4;
    }
    else if (ncount[1] == 4 && ncount[4] == 4)
    {
        err = Face::hexagon_2_quads(bound_nodes, newFaces, 1);
        if (err) return 4;
    }

    //
    // Check for "area invariance". Before and after the decomposition
    // area must be same. When there is concave faces, such violations
    // may occur. This method just avoid overlapping cases and produce
    // simple polygons.
    //

    double area0 = 0.0;
    for (size_t i = 0; i < adjFaces.size(); i++)
        area0 += adjFaces[i]->getArea();

    // We have decomposed a hexagon into 2 quads...
    double area1 = 0.0;
    for (size_t i = 0; i < 2; i++)
        area1 += newFaces[i]->getArea();

    if (fabs(area1 - area0) < 1.0E-06) return 0;

    return 1;
}

////////////////////////////////////////////////////////////////////

int
Edge33::remove_internal_one()
{
    if (newFaces.empty()) return 1;

    Vertex *v0 = connect[0];
    Vertex *v1 = connect[1];

    if (v0->isVisited() || v1->isVisited()) return 1;
    if (v0->isRemoved() || v1->isRemoved()) return 1;

    // Check for double removal..
    assert(bound_nodes.size() == 6);
    for (int i = 0; i < 6; i++)
        if (bound_nodes[i]->isVisited()) return 1;

    assert(adjFaces.size() == 4);

    // In this case, no new vertex is created, only for faces are
    // deleted and two new inserted.
    assert(mesh);
    mesh->addFace(newFaces[0]);
    mesh->addFace(newFaces[1]);

    // The bridge edge goes away along with the neighboring faces..
    for (int i = 0; i < 4; i++)
        mesh->remove(adjFaces[i]);

    mesh->remove(connect[0]);
    mesh->remove(connect[1]);

    Point3D backupCoords[6];
    for (int i = 0; i < 6; i++)
        backupCoords[i] = bound_nodes[i]->getXYZCoords();

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(bound_nodes);

    set<Face*> faces_to_check;
    for (int i = 0; i < 6; i++)
    {
        vector<Face*> vfaces = bound_nodes[i]->getRelations2();
        for (size_t j = 0; j < vfaces.size(); j++)
            faces_to_check.insert(vfaces[j]);
    }

    bool pass = 1;
    set<Face*>::const_iterator siter;
    for (siter = faces_to_check.begin(); siter != faces_to_check.end(); ++siter)
    {
        Face *f = *siter;
        if (f->invertedAt() >= 0)
        {
            pass = 0;
            break;
        }
    }

    if (!pass)
    {
        mesh->remove(newFaces[0]);
        mesh->remove(newFaces[1]);
        for (int i = 0; i < 4; i++)
            mesh->addFace(adjFaces[i]);
        mesh->addNode(connect[0]);
        mesh->addNode(connect[1]);
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setXYZCoords(backupCoords[i]);
    }
    else
    {
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setVisitMark(1);
        v0->setVisitMark(1);
        v1->setVisitMark(1);
    }

    // Destructor must not deallocate newFaces since they are not kept by
    // the mesh object.
    newFaces[0] = NULL;
    newFaces[1] = NULL;

    return 0;
}

///////////////////////////////////////////////////////////////////////////

int
Edge33::remove_boundary_one()
{
    Vertex *v0 = connect[0];
    Vertex *v1 = connect[1];

    if (v0->isVisited() || v1->isVisited()) return 1;
    if (v0->isRemoved() || v1->isRemoved()) return 1;

    // Our assumption is that all the edges are simple.
    adjFaces = Mesh::getRelations112(v0, v1);
    if (adjFaces.size() != 2) return 1;

    // Check for the internal face.
    Face *internal_face = NULL;

    int ncount_boundfaces = 0;
    if (internal_face == NULL)
    {
        if (!adjFaces[0]->hasBoundaryEdge())
            internal_face = adjFaces[0];
        else
            ncount_boundfaces++;
    }
    if (internal_face == NULL)
    {
        if (!adjFaces[1]->hasBoundaryEdge())
            internal_face = adjFaces[1];
        else
            ncount_boundfaces++;
    }

    // A valid boundary bridge must have one boundary face and one internal
    // face.
    if (internal_face == NULL) return 1;
    if (ncount_boundfaces == 2) return 2;

    // Swap the opposite edge:
    Vertex *v2, *v3;
    Face::opposite_nodes(internal_face, v0, v1, v2, v3);

    // Try swapping at the first node.
    int err;
    SwapEdge edge1(mesh, v2, v3, internal_face);
    err = edge1.apply_deficient_rule(connect[0]);
    if (!err)
    {
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setVisitMark(1);
        v0->setVisitMark(1);
        v1->setVisitMark(1);
        return 0;
    }

    // Try swapping at the second node.
    SwapEdge edge2(mesh, v2, v3, internal_face);
    err = edge2.apply_deficient_rule(connect[1]);
    if (!err)
    {
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setVisitMark(1);
        v0->setVisitMark(1);
        v1->setVisitMark(1);
        return 0;
    }

    return 3;
}
///////////////////////////////////////////////////////////////////////////////

int
Edge33::commit()
{
    int err = 1;

    if (boundary)
        err = remove_boundary_one();
    else
        err = remove_internal_one();

    return err;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridges_in_layer(int layerid)
{
    cout << "CSV : " << endl;
    abort();
    int rel2exist = mesh->build_relations(0, 2);
    mesh->search_boundary();

    if (vEdges33.empty())
        search_edges33_in_layer(layerid);

    return remove_bridges_once();
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridges_once()
{
    int rel0exist = mesh->build_relations(0, 0);
    int rel2exist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    int ncount = 0;

    if (vEdges33.empty())
        search_edges33();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
    }

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    for (size_t i = 0; i < vEdges33.size(); i++)
    {
        int err = vEdges33[i].remove();
        if (!err) ncount++;
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: number of bridges removed " << ncount << endl;
        mesh->collect_garbage();
    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);

    vEdges33.clear();

    return ncount;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridges()
{
    if (!mesh->isPruned()) mesh->prune();

    while (1)
    {
        int ncount = remove_bridges_once();
        if (ncount == 0) break;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

/*
int
SwapEdge::getPosOf(const Vertex *vertex) const
{
    for (int i = 0; i < 6; i++)
        if (bound_nodes[i] == vertex) return i;

    return -1;
}
/////////////////////////////////////////////////////////////////////////////

int
SwapEdge::build_boundary()
{
    assert(mesh);
    assert(connect[0]);
    assert(connect[1]);
    adjFaces.resize(2);
    adjFaces[0] = NULL;
    adjFaces[1] = NULL;

    assert(mesh->getAdjTable(0, 2));

    vector<Face*> nghs;
    // Since the degree of each node of an existing edge decreases by
    // one, by restricting the swapping to degree greater than 3 will
    // ensure that no doublet is created.
    //
    nghs = connect[0]->getRelations2();
    int nsize1 = nghs.size();
    if (connect[0]->isBoundary())
    {
        if (nsize1 < 3) return 1;
    }
    else
    {
        if (nsize1 < 4) return 1;
    }

    nghs = connect[1]->getRelations2();
    int nsize2 = nghs.size();

    if (connect[1]->isBoundary())
    {
        if (nsize2 < 3) return 1;
    }
    else
    {
        if (nsize2 < 4) return 1;
    }

    nghs = Mesh::getRelations112(connect[0], connect[1]);

    if (nghs.size() != 2) return 1;

    if (firstFace == NULL)
    {
        adjFaces[0] = nghs[0];
        adjFaces[1] = nghs[1];
    }
    if (nghs[0] == firstFace)
    {
        adjFaces[0] = nghs[0];
        adjFaces[1] = nghs[1];
    }
    if (nghs[1] == firstFace)
    {
        adjFaces[0] = nghs[1];
        adjFaces[1] = nghs[0];
    }

    int err = get_boundary_nodes_chain();

    if (err) return 2;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
SwapEdge::get_boundary_nodes_chain()
{
    vector<Edge> bndedges;
    bndedges.reserve(6);

    Edge sharededge(connect[0], connect[1]);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Vertex *vf0 = adjFaces[i]->getNodeAt((j + 0) % 4);
            Vertex *vf1 = adjFaces[i]->getNodeAt((j + 1) % 4);
            Edge edge(vf0, vf1);
            if (!edge.isSameAs(sharededge))
                bndedges.push_back(edge);
        }
    }

    Mesh::make_chain(bndedges);
    bound_nodes = Mesh::chain_nodes(bndedges);
    assert(bound_nodes.size() == 6);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void SwapEdge::update_front()
{
    int layerid = connect[0]->getLayerID();
    connect[1]->setLayerID(layerid + 1);
    adjFaces[0]->setLayerID(layerid);
    adjFaces[1]->setLayerID(layerid);
}

////////////////////////////////////////////////////////////////////////////////

void
SwapEdge::backup()
{
    bkp_data.diagonalConnect[0] = connect[0];
    bkp_data.diagonalConnect[1] = connect[1];
    bkp_data.face1Connect = adjFaces[0]->getNodes();
    bkp_data.face2Connect = adjFaces[1]->getNodes();
    for (int i = 0; i < 6; i++)
    {
        Vertex *v = bound_nodes[i];
        bkp_data.pCoords[v] = v->getXYZCoords();
    }
}
////////////////////////////////////////////////////////////////////////////////

void
SwapEdge::rollback()
{
    connect[0]->removeRelation0(connect[1]);
    connect[1]->removeRelation0(connect[0]);

    // Remove the vertex-face relationship.
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->removeRelation2(adjFaces[0]);
    }

    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->removeRelation2(adjFaces[1]);
    }

    // Change the connectivity of the two quads.
    adjFaces[0]->setNodes(bkp_data.face1Connect);
    adjFaces[1]->setNodes(bkp_data.face2Connect);
    connect[0] = bkp_data.diagonalConnect[0];
    connect[1] = bkp_data.diagonalConnect[1];

    // Update the vertex-vertex relationship.
    connect[0]->addRelation0(connect[1]);
    connect[1]->addRelation0(connect[0]);

    // Update the vertex-face relationships.
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->addRelation2(adjFaces[0]);
    }

    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->addRelation2(adjFaces[1]);
    }

    for (int i = 0; i < 6; i++)
    {
        Vertex *v = bound_nodes[i];
        v->setXYZCoords(bkp_data.pCoords[v]);
    }

    update_front();
}
////////////////////////////////////////////////////////////////////////////////

int
SwapEdge::make_new_diagonal_at(int pos, bool bound_check)
{
    //
    // Given closednodes ( exactly six nodes) in an order ( clockwise or
    // counter-clockwise) and create two quads having common diagonal between
    // (pos, pos+3).
    //
    assert(pos >= 0 && pos < 6);

    if (bound_check)
    {
        if (bound_nodes[(pos + 0) % 6]->isBoundary()) return 1;
        if (bound_nodes[(pos + 3) % 6]->isBoundary()) return 1;
    }

    // First perform full backup of the data structures.
    backup();

    // Remove the vertex-vertex relationship.
    connect[0]->removeRelation0(connect[1]);
    connect[1]->removeRelation0(connect[0]);

    // Remove the vertex-face relationship.
    assert(adjFaces[0]);
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->removeRelation2(adjFaces[0]);
    }

    assert(adjFaces[1]);
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->removeRelation2(adjFaces[1]);
    }

    // Change the connectivity of the two quads.
    vector<Vertex*> qConnect(4);
    qConnect[0] = bound_nodes[(pos + 0) % 6];
    qConnect[1] = bound_nodes[(pos + 1) % 6];
    qConnect[2] = bound_nodes[(pos + 2) % 6];
    qConnect[3] = bound_nodes[(pos + 3) % 6];
    connect[0] = qConnect[0];
    adjFaces[0]->setNodes(qConnect);

    qConnect[0] = bound_nodes[(pos + 3) % 6];
    qConnect[1] = bound_nodes[(pos + 4) % 6];
    qConnect[2] = bound_nodes[(pos + 5) % 6];
    qConnect[3] = bound_nodes[(pos + 6) % 6];
    connect[1] = qConnect[0];
    adjFaces[1]->setNodes(qConnect);

    // Update the vertex-vertex relationship.
    connect[0]->addRelation0(connect[1]);
    connect[1]->addRelation0(connect[0]);

    // Update the vertex-face relationships.
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->addRelation2(adjFaces[0]);
    }

    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->addRelation2(adjFaces[1]);
    }

    // The change may produce concave or tangled mesh. Do some local
    // Laplacian smoothing at the six nodes only. Hopefully, the elements
    // will be acceptable. If not, we need to rollback this operation.

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(bound_nodes);

    // Look at the neighbors, if some elements is inverted. Here we create
    // a set of neighbours, to avoid duplicate checking.
    set<Face*> neighSet;
    for (int i = 0; i < 6; i++)
    {
        Vertex *v = bound_nodes[i];
        if (!v->isBoundary())
        {
            vector<Face*> vfaces = v->getRelations2();
            for (int j = 0; j < vfaces.size(); j++)
                neighSet.insert(vfaces[j]);
        }
    }
    assert(!neighSet.empty());

    set<Face*>::const_iterator it;
    for (it = neighSet.begin(); it != neighSet.end(); ++it)
    {
        Face *face = *it;
        if (face->invertedAt() >= 0)
        {
            rollback();
            return 1;
        }
    }
    update_front();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
QuadEdge::commit()
{
    if (newFaces.empty()) return 1;

    for (int i = 0; i < newFaces.size(); i++)
        if (newFaces[i] == NULL) return 2;

    if (adjFaces[0]->isRemoved() || adjFaces[1]->isRemoved())
    {
        for (int i = 0; i < newNodes.size(); i++)
            delete newNodes[i];
        newNodes.clear();

        for (int i = 0; i < newFaces.size(); i++)
            delete newFaces[i];
        newFaces.clear();
        return 3;
    }

    Vertex *vertex;
    for (int i = 0; i < 4; i++)
    {
        vertex = adjFaces[0]->getNodeAt(i);
        if (vertex->isVisited()) return 4;
        vertex = adjFaces[1]->getNodeAt(i);
        if (vertex->isVisited()) return 4;
    }

    adjFaces[0]->setRemoveMark(1);
    adjFaces[1]->setRemoveMark(1);

    assert(mesh);

    for (int i = 0; i < newNodes.size(); i++)
        mesh->addNode(newNodes[i]);

    for (int i = 0; i < newFaces.size(); i++)
        mesh->addFace(newFaces[i]);

    for (int i = 0; i < 4; i++)
    {
        vertex = adjFaces[0]->getNodeAt(i);
        vertex->setVisitMark(1);
        vertex = adjFaces[1]->getNodeAt(i);
        vertex->setVisitMark(1);
    }

    return 0;
}
 */

/////////////////////////////////////////////////////////////////////

/*
int
SwapEdge::apply_reduce_degree_rule()
{
    if (build_boundary() != 0) return 1;

    vector<Vertex*> neighs;

    neighs = connect[0]->getRelations0();
    int d1 = neighs.size();
    if (d1 < 4) return 1;

    neighs = connect[1]->getRelations0();
    int d2 = neighs.size();
    if (d2 < 4) return 1;

    int start_pos = this->getPosOf(connect[0]);
    assert(start_pos >= 0);

    for (int i = 0; i < 2; i++)
    {
        int pos = (start_pos + i + 1) % 6;
        Vertex *v0 = bound_nodes[ pos ];
        if (v0->isBoundary()) continue;
        neighs = v0->getRelations0();
        int d3 = neighs.size();

        Vertex *v3 = bound_nodes[(pos + 3) % 6];
        if (v3->isBoundary()) continue;
        neighs = v3->getRelations0();
        int d4 = neighs.size();

        if (SwapEdge::is_topologically_valid_swap(d1, d2, d3, d4))
        {
            int err = make_new_diagonal_at(pos);
            if (!err) return 0;
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_concave_rule()
{
    if (build_boundary() != 0) return 1;

    int start_pos = getPosOf(connect[0]);
    for (int i = 0; i < 2; i++)
    {
        int err = make_new_diagonal_at((start_pos + 1 + i) % 6);
        if (!err) return 0;
    }
    return 1;
}

/////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_bound_rule()
{
    if (!connect[0]->isBoundary()) return 1;
    if (connect[1]->isBoundary()) return 1;

    if (build_boundary() != 0) return 2;

    for (int i = 0; i < 3; i++)
    {
        Vertex *v0 = bound_nodes[(i + 0) % 6];
        Vertex *v3 = bound_nodes[(i + 3) % 6];
        if ((!v0->isBoundary()) && (!v3->isBoundary()))
        {
            int err = make_new_diagonal_at(i, 0);
            if (!err) return 0;
        }
    }
    return 3;
}

/////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_advance_front_rule()
{
    if (build_boundary() != 0) return 1;

    int layerid = connect[0]->getLayerID();
    for (int i = 0; i < 3; i++)
    {
        Vertex *v0 = bound_nodes[(i + 0) % 6];
        Vertex *v3 = bound_nodes[(i + 3) % 6];
        if ((v0->getLayerID() > layerid) && (v3->getLayerID() > layerid))
        {
            int err = make_new_diagonal_at(i);
            if (!err) return 0;
        }
    }
    return 3;
}

/////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_singlet_rule(Vertex *singlet)
{
    assert(connect[0]->isBoundary());

    if (connect[1]->isBoundary()) return 1;

    if (build_boundary() != 0) return 1;
    ////////////////////////////////////////////////////////////////////////////
    // Objective :: Swap quads common diagonal such that the resulting diagonal
    //              contains the Singlet node.
    ////////////////////////////////////////////////////////////////////////////
    assert(QuadCleanUp::isSinglet(singlet));

    // First Vertex must be boundary

    if (adjFaces[0] == NULL || adjFaces[1] == NULL) return 2;

    assert(QuadCleanUp::hasSinglet(adjFaces[0]));

    // Make sure that other face doesn't have a singlet
    if (QuadCleanUp::hasSinglet(adjFaces[1])) return 3;

    //Find the position of singlet in the contour ...
    int pos = this->getPosOf(singlet);
    assert(pos >= 0);

    // Create new quads whose common diagonal must contain the singlet.
    return make_new_diagonal_at(pos, 0);
}

/////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_deficient_rule(Vertex *vertex)
{
    assert(firstFace);
    if (vertex->isBoundary()) return 1;

    int d1 = connect[0]->getRelations2().size();
    int d2 = connect[1]->getRelations2().size();

    // Note 1: Just near the deficient node, there must be atleast 5 nodes in
    // order to swap the edge,otherwise, things may not converge, on the
    // same level, the deficient node will keep changing the location.
    //
    // Note 2: The node connect[1] is opposite to the deficient node and
    // most probably on the next level. If we keep moving the deficient
    // node inside the domain, then it is fine to have its degree to 4, so
    // that the deficient node will be generated there.
    vector<Face*> faces;
    if (d1 == 4 && d2 > 3)
    {
        int layerid = connect[0]->getLayerID();
        faces = vertex->getRelations2();
        if (faces.size() == 3)
        {
            Face *f0 = firstFace;
            faces = Mesh::getRelations112(connect[0], connect[1]);
            assert(faces.size() == 2);
            Face *f1 = NULL;
            if (faces[0] == f0) f1 = faces[1];
            if (faces[1] == f0) f1 = faces[0];
            assert(f1);
            int pos = f1->queryPosOf(connect[0]);
            Vertex *vopp = f1->getNodeAt((pos + 2) % 4);
            if (connect[1]->getLayerID() == layerid + 1 && vopp->getLayerID() == layerid + 1)
            {
                SwapEdge edge(mesh, connect[1], vopp, f1);
                int err = edge.apply_deficient_rule(connect[0]);
                if (err) return 1;
            }
        }
    }

    d1 = connect[0]->getRelations2().size();
    d2 = connect[1]->getRelations2().size();

    if (d1 < 5 && d2 < 3) return 2;

    // Having these conditions set, now build the structure.
    if (build_boundary() != 0)
    {
        return 3;
    }

    //Find the position of triplet in the contour ...
    int pos = this->getPosOf(vertex);
    assert(pos >= 0);

    // Create new quads whose common diagonal must contain the singlet.
    return make_new_diagonal_at(pos);
}
 */

/////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::boundary_vertex_degree_reduction_once()
{
    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();
    mesh->setFeatureAngles();

    vector<Doublet> doublets = search_interior_doublets();
    assert(doublets.empty());

    size_t numnodes = mesh->getSize(0);

    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
        {
            if (vertex->getFeatureAngle() >= mesh->getMaxFeatureAngle())
            {
                vector<Vertex*> vneighs = vertex->getRelations0();
                if (vneighs.size() > 3) // Two lie on the boundary and only one must be internal
                {
                    for (int k = 0; k < vneighs.size(); k++)
                    {
                        vector<Vertex*> wneighs = vneighs[k]->getRelations0();
                        if (!vneighs[k]->isBoundary() && wneighs.size() > 3)
                        {
                            SwapEdge edge(mesh, vertex, vneighs[k]);
                            int err = edge.apply_bound_rule();
                            if (!err)
                            {
                                ncount++;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    if (ncount)
        cout << "Info: number of boundary edges swapped " << ncount << endl;

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    doublets = search_interior_doublets();
    assert(doublets.empty());


    return ncount;
}
///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::internal_vertex_degree_reduction_once()
{
    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();

    vector<Doublet> doublets = search_interior_doublets();
    assert(doublets.empty());

    size_t numnodes = mesh->getSize(0);

    vector<Vertex*> nodes = mesh->getNodes();
    sort(nodes.begin(), nodes.end(), HighVertexDegreeCompare());

    size_t ncount = 0;

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = nodes[i];
        if (!vertex->isBoundary())
        {
            vector<Vertex*> vneighs = vertex->getRelations0();
            if (vneighs.size() > 4)
            {
                sort(vneighs.begin(), vneighs.end(), HighVertexDegreeCompare());
                for (int k = 0; k < vneighs.size(); k++)
                {
                    vector<Vertex*> wneighs = vneighs[k]->getRelations0();
                    if (wneighs.size() > 3)
                    {
                        SwapEdge edge(mesh, vertex, vneighs[k]);
                        int err = edge.apply_reduce_degree_rule();
                        if (!err)
                        {
                            ncount++;
                            break;
                        }
                    }
                }
            }

        }
    }

    if (ncount)
        cout << "Info: number of internal edges swapped " << ncount << endl;

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    doublets = search_interior_doublets();
    assert(doublets.empty());

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

/*
int
QuadCleanUp::swap_concave_faces()
{
    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
    }

    int ncount = 0;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        int pos = face->invertedAt();
        if (pos >= 0)
        {
            for (int i = 0; i < 2; i++)
            {
                Vertex *v0 = face->getNodeAt((pos + 1 + i) % 4);
                Vertex *v1 = face->getNodeAt((pos + 2 + i) % 4);
                SwapEdge edge(mesh, v0, v1);
                int err = edge.apply_concave_rule();
                if (!err)
                {
                    ncount++;
                    break;
                }
            }
        }
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: # of swapped edges " << ncount << endl;
    }

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    return ncount;
}
 */

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::vertex_degree_reduction()
{
    int ncount1 = 0;
    int ncount2 = 0;
    while (1)
    {
        ncount1 = boundary_vertex_degree_reduction_once();
        ncount2 = internal_vertex_degree_reduction_once();
        int err = lapsmooth->execute();
        if (err || (ncount1 + ncount2) == 0) break;
    }

}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_singlet_rule(Vertex *singlet)
{
    int err;
    if (singlet->getLayerID() != 0) return 1;
    if (!QuadCleanUp::isSinglet(singlet)) return 1;

    vector<Face*> vfaces = singlet->getRelations2();
    assert(vfaces.size() == 1);

    Face *face = vfaces[0];
    int pos = face->queryPosOf(singlet);

    Vertex *v0 = face->getNodeAt((pos + 0) % 4);
    Vertex *v1 = face->getNodeAt((pos + 1) % 4);
    Vertex *v2 = face->getNodeAt((pos + 2) % 4);
    Vertex *v3 = face->getNodeAt((pos + 3) % 4);

    SwapEdge edge1(mesh, v1, v2);
    err = edge1.apply_singlet_rule(singlet);
    if (!err) return 0;

    SwapEdge edge2(mesh, v3, v2);
    err = edge2.apply_singlet_rule(singlet);
    if (!err) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_triplet_rule(Vertex *vertex)
{
    //****************************************************************************
    // Implementation ideas:
    //
    // This is the case of 5-434 transition, Where 5 degree node is opposite
    // to 3 degree node and with first swap, it becomes 4-335 or 4-534 and with
    // the second swap, it becomes 3-444 and therefore, the 3 degree node moves
    // inside the domain.
    //
    // Having degree 5 is a somewhat strict condition, because degree four is
    // sufficient, but it will create doublet, and as per our design goals, we
    // do not want to introduce any new doublet in the mesh. In future, we may
    // change the conditions, if that is useful.
    //
    //****************************************************************************
    if (vertex->isBoundary()) return 1;

    vector<Vertex*> vnodes = vertex->getRelations0();
    if (vnodes.size() != 3) return 1;

    int layerid = vertex->getLayerID();
    vector<Face*> vfaces = vertex->getRelations2();
    for (int j = 0; j < 3; j++)
    {
        Face *face = vfaces[j];
        if (face->getLayerID() == layerid)
        {
            int pos = face->queryPosOf(vertex);
            Vertex *v1 = face->getNodeAt((pos + 1) % 4);
            Vertex *v2 = face->getNodeAt((pos + 2) % 4);
            Vertex *v3 = face->getNodeAt((pos + 3) % 4);
            // If v1 and v3 are at the lower level, they are already processed.
            if (v1->getLayerID() < layerid || v3->getLayerID() < layerid) return 1;
            int d1 = v1->getRelations2().size();
            int d2 = v2->getRelations2().size();
            int d3 = v3->getRelations2().size();
            if (d3 > d1)
                std::swap(v1, v3);

            //  The opposite node must have degree at least five.
            if (max(d1, d3) <= 4 && d2 < 5) return 1;

            //  The opposite node must be at higher level than at (v0, v1, v3)
            if (v2->getLayerID() > layerid)
            {
                SwapEdge edge1(mesh, v1, v2, face);
                int err1 = edge1.apply_deficient_rule(vertex);
                if (!err1)
                {
                    return 0;
                }

                SwapEdge edge2(mesh, v3, v2, face);
                int err2 = edge2.apply_deficient_rule(vertex);
                if (!err2 == 0)
                {
                    return 0;
                }
            }
        }
    }

    return 1;

}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_bridge_rule(Vertex *v0, Vertex *v1)
{
    int layerid = v0->getLayerID();
    if (v1->getLayerID() != layerid) return 1;

    // Although it is checked once again, but it is essential.
    vector<Face*> adjFaces;
    adjFaces = v0->getRelations2();
    if (adjFaces.size() != 3) return 2;

    adjFaces = v1->getRelations2();
    if (adjFaces.size() != 3) return 3;

    // Our assumption is that all the edges are simple.
    adjFaces = Mesh::getRelations112(v0, v1);
    if (adjFaces.size() != 2) return 4;

    if (adjFaces[0]->getLayerID() == adjFaces[1]->getLayerID()) return 5;

    // Check for the internal face.
    Face *internal_face = NULL;
    internal_face = adjFaces[0]->getLayerID() > adjFaces[1]->getLayerID() ? adjFaces[0] : adjFaces[1];

    // Swap the opposite edge:
    Vertex *v2, *v3;
    Face::opposite_nodes(internal_face, v0, v1, v2, v3);

    // Try swapping at the first node.
    int err;
    SwapEdge edge1(mesh, v2, v3, internal_face);
    err = edge1.apply_deficient_rule(v0);
    if (!err) return 0;

    // Try swapping at the second node.
    SwapEdge edge2(mesh, v2, v3, internal_face);
    err = edge2.apply_deficient_rule(v1);
    if (!err) return 0;

    return 1;

}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_excess_rule(Vertex *vertex)
{
    if (vertex->isBoundary()) return 1;

    int layerid = vertex->getLayerID();
    vector<Vertex*> vneighs = vertex->getRelations0();
    int degree = vneighs.size();

    if (degree < 5) return 1;

    int ncount = 0;
    for (int k = 0; k < degree; k++)
        if (vneighs[k]->getLayerID() > layerid) ncount++;

    if (ncount < 2) return 1;

    for (int k = 0; k < degree; k++)
    {
        if (vneighs[k]->getLayerID() > layerid)
        {
            SwapEdge edge(mesh, vertex, vneighs[k]);
            int err = edge.apply_advance_front_rule();
            if (!err) return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::advance_front_edges_swap_once(int layerid)
{
    int err;
    if (mesh->getAdjTable(0, 0))
        mesh->clear_relations(0, 0);

    if (mesh->getAdjTable(0, 2))
        mesh->clear_relations(0, 2);

    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    //
    // The input mesh should be doublet free and the output mesh will
    // be doublet free.
    //
    vector<Doublet> doublets = search_interior_doublets();
    assert(doublets.empty());

    // If the boundary is unknown ...
    if (layerid == 0)
        mesh->search_boundary();

    // We need atleast two layer to work with.
    size_t numnodes = mesh->getSize(0);

    // Check how many iiregular nodes in the present layer...
    size_t num_irregular_nodes = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
        if (vertex->getLayerID() == layerid)
        {
            if (vertex->getRelations2().size() != 4)
                num_irregular_nodes++;
        }
    }

    size_t ncount = 0;

    //
    // If the layer is boundary, then highest priority must be given to remove
    // the singlets. Many singlets can be removed by swapping the edges. There are
    // some cases, where "Swapping" may not remove singlets. There are two such
    // scenarios, which are handled by calling "remove_boundary_singlet" method.
    //

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
        {
            err = apply_advance_front_singlet_rule(vertex);
            if (!err) ncount++;
        }
    }

    // Second priority is given to bridges in the mesh. Removal of bridges
    // create a three degree node, which is the next case.

    if (ncount == 0)
    {
        /*
           vector<Edge33> bridges = search_edges33();
           for (size_t i = 0; i < bridges.size(); i++){
                Vertex *v0 = bridges[i].getNodeAt(0);
                Vertex *v1 = bridges[i].getNodeAt(1);
                if( v0->getLayerID() == layerid  && v1->getLayerID() == layerid ) {
                    err = apply_advance_front_bridge_rule(v0, v1);
                    if( !err ) ncount++;
                    v0->setTag(3);
                    v1->setTag(3);
                    mesh->saveAs("dbg.dat");
                    exit(0);
                }
           }
         */
    }

    // Third Priority is given to degree three nodes. if the adjacent nodes are
    // regular, then two swaps are required. In the first step, the degree of the
    // adjacent node is increased by swapping from outer layer, and then swapping
    // is done to make 3 degree node to regular node.

    if (ncount == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);

            if (vertex->getLayerID() == layerid)
            {
                int err = apply_advance_front_triplet_rule(vertex);
                if (!err) ncount++;
            }
        }
    }

    // Final case, this is the most intuitive and general; swapping is done to
    // reduce the vertex degree.
    if (ncount == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);
            int lid = vertex->getLayerID();
            if (lid >= layerid)
            {
                while (1)
                {
                    int err = apply_advance_front_excess_rule(vertex);
                    if (!err)
                        ncount++;
                    else
                        break;
                }
            }
        }
    }

#ifdef DEBUG
    if (ncount)
        cout << "Info: Layer " << layerid << " number of layer edges swapped " << ncount << endl;
#endif

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    mesh->collect_garbage();

    doublets = search_interior_doublets();
    assert(doublets.empty());

    return ncount;
}

////////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::advancing_front_edges_swap()
{
    //****************************************************************************
    // Implementation ideas:
    // In general, irregular nodes are scattered around in the domain ( even after
    // initial cleanup operations. With the advance front swapping, we try to
    // (1) clean the mesh starting from the boundary (2) cluster irregular nodes
    // deep inside the mesh ( far from the boundary ) with the hope that clustering
    // will provide more opportunities for cleanup operations ( for example, enable
    // diamonds, bridges etc).
    //
    // There are two side-effects of this procedure.
    // 1)  It may create very high valance nodes in the last layers, which we may
    //     not be able to clean. ( this is not a stopper, as we can call vertex
    //     reduction modules, which will again scatter the irregular nodes.
    // 2)  The position of irregular nodes may depends on the density of mesh.
    //
    //
    //***************************************************************************

    size_t ncount = 0;
    int nfronts;

    nfronts = mesh->setWavefront(0);
    nfronts = mesh->setWavefront(2);

    int istart = 1;
    int iend = nfronts;

    cout << "Temporay Changed " << endl;
    iend = 1;

    for (int ilayer = istart; ilayer < iend; ilayer++)
    {
        // All nodes below a current level never move or modified, therefore
        // make them inactive for smooting operations.
        //
        size_t numnodes = mesh->getSize(0);
        for (size_t j = 0; j < numnodes; j++)
        {
            Vertex *v = mesh->getNodeAt(j);
            if (v->getLayerID() < ilayer)
                v->setConstrainedMark(1);
            else
                v->setConstrainedMark(0);
        }

        int ncount1 = advance_front_edges_swap_once(ilayer);
        //          lapsmooth->execute ();
        ncount += ncount1;
    }

}
////////////////////////////////////////////////////////////////////

vector<Face*>
QuadCleanUp::search_flat_quads()
{
    //  Public Function ...
    size_t numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    vector<Face*> flatQ, neighs;

    int edgefaces[4];
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        assert(face);
        if (face->getSize(0) == 4)
        {
            int boundary = 1;
            for (int j = 0; j < 4; j++)
            {
                Vertex *vb = face->getNodeAt(j);
                if (!vb->isBoundary())
                {
                    boundary = 0;
                    break;
                }
            }

            if (boundary)
            {
                int bound_edges = 0;
                for (int j = 0; j < 4; j++)
                {
                    Vertex *v0 = face->getNodeAt((j + 0) % 4);
                    Vertex *v1 = face->getNodeAt((j + 1) % 4);
                    neighs = Mesh::getRelations112(v0, v1);
                    if (neighs.size() == 1)
                        bound_edges++;
                    edgefaces[j] = neighs.size();
                }

                if (bound_edges == 3)
                {
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

vector<Vertex*>
QuadCleanUp::search_restricted_nodes()
{
    size_t numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 0);

    mesh->search_boundary();

    vector<Vertex*> restricted_nodes, vneighs;

    int ncount;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        ncount = 0;
        if (!vertex->isBoundary())
        {
            vneighs = vertex->getRelations0();
            for (size_t j = 0; j < vneighs.size(); j++)
                if (vneighs[j]->isBoundary()) ncount++;
        }

        vertex->setTag(1);
        if (ncount > 1)
        {
            restricted_nodes.push_back(vertex);
            vertex->setTag(0);
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 0);

    if (restricted_nodes.size())
        cout << "Info: Number of restricted nodes detected : " << restricted_nodes.size() << endl;

    return restricted_nodes;
}

///////////////////////////////////////////////////////////////////////////////

vector<Face*>
QuadCleanUp::search_restricted_faces()
{
    size_t numfaces = mesh->getSize(2);

    mesh->search_boundary();

    vector<Face*> restricted_faces;

    int ncount;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(0);
        Vertex *v0 = face->getNodeAt(0);
        Vertex *v1 = face->getNodeAt(1);
        Vertex *v2 = face->getNodeAt(2);
        Vertex *v3 = face->getNodeAt(3);
        if ((v0->isBoundary() && v2->isBoundary()) || (v1->isBoundary() && v3->isBoundary()))
        {
            restricted_faces.push_back(face);
            face->setTag(1);
        }
    }

    cout << "Info: Number of restricted faces detected : " << restricted_faces.size() << endl;

    return restricted_faces;
}

////////////////////////////////////////////////////////////////////////////////

vector<Doublet>
QuadCleanUp::search_interior_doublets()
{
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    vDoublets.clear();

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (isDoublet(vertex))
        {
            Doublet newdoublet(mesh, vertex);
            if (newdoublet.isSafe())
            {
                newdoublet.makeShield();
                vDoublets.push_back(newdoublet);
            }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    if (vDoublets.size())
    {
        cout << "Info: Number of interior doublets Detected : "
                << vDoublets.size() << endl;
    }

    return vDoublets;
}

//////////////////////////////////////////////////////////////////////////

vector<Singlet>
QuadCleanUp::search_boundary_singlets()
{
    //
    // A boundary singlet is a vertex which is shared by only one face.
    // They are undesirables in the quad mesh as that would mean large
    // angle on some of the edges..
    //
    // For the flat singlet ( angle closer to 180 degree ). it is easy
    // to remove the neighbouring quad from the mesh.
    //
    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();
    mesh->setFeatureAngles();

    assert(mesh->getAdjTable(0, 2));

    for (size_t i = 0; i < mesh->getSize(2); i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }


    size_t numnodes = mesh->getSize(0);
    vSinglets.clear();
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        if (isSinglet(v))
        {
            if (v->getFeatureAngle() >= mesh->getMaxFeatureAngle())
            {
                Singlet newsinglet(mesh, v);
                vSinglets.push_back(newsinglet);
            }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    if (vSinglets.size())
        cout << "Info: Number of Singlets detected " << vSinglets.size() << endl;

    return vSinglets;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_boundary_singlets_once()
{
    int relexist = mesh->build_relations(0, 2);

    if (vSinglets.empty())
        search_boundary_singlets();

    // Simplest case: just remove it by refinement;
    int ncount = 0;
    for (int i = 0; i < vSinglets.size(); i++)
    {
        int err = vSinglets[i].remove();
        if (!err) ncount++;
    }

    /* Reconsider the following code later.

    // Necessary for type2 changes.
    size_t numfaces = mesh->getSize (2);
    for (size_t i = 0; i < numfaces; i++)
      {
        Face *face = mesh->getFaceAt (i);
        face->setVisitMark (0);
      }

    for (int i = 0; i < vSinglets.size (); i++)
      vSinglets[i].update_type1 ();

    for (int i = 0; i < vSinglets.size (); i++)
      vSinglets[i].update_type2 ();

      for (int i = 0; i < vSinglets.size (); i++)
        vSinglets[i].update_type3 ();

      size_t numnodes = mesh->getSize (0);
      for (size_t i = 0; i < numnodes; i++)
        {
          Vertex *vertex = mesh->getNodeAt (i);
          vertex->setVisitMark (0);
        }

      for (int i = 0; i < vSinglets.size (); i++)
        {
          int err = vSinglets[i].commit ();
          if (!err) ncount++;
        }
     */

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        mesh->collect_garbage();
        mesh->setBoundaryKnown(0);

    }

    vSinglets.clear();

    if (!relexist)
        mesh->clear_relations(0, 2);


    return ncount;
}
//////////////////////////////////////////////////////////////////////////

int Singlet::remove()
{
    vector<Face*> vfaces = vertex->getRelations2();
    if (vfaces.size() > 1) return 1;

    return mesh->refine_quad15(vfaces[0]);
}
//////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_boundary_singlets()
{
    mesh->setFeatureAngles();

    int ncount = 0;
    while (1)
    {
        size_t nremoved = remove_boundary_singlets_once();
        if (nremoved == 0) break;
        if (!nremoved) return ncount;
        ncount += nremoved;
    }

    if (!mesh->isConsistentlyOriented())
        mesh->makeConsistentlyOriented();

    return 0;
}

////////////////////////////////////////////////////////////////////

int FaceClose::remove()
{
    return 1;
}

////////////////////////////////////////////////////////////////////

bool
FaceClose::isSafe() const
{
    if (face->isRemoved()) return 0;

    if (vertex0->isRemoved()) return 0;
    if (vertex2->isRemoved()) return 0;

    if (vertex0->isBoundary()) return 0;
    if (vertex2->isBoundary()) return 0;

    int vpos = face->queryPosOf(vertex0);
    assert(vpos >= 0);

    if (face->getNodeAt((vpos + 2) % 4) != vertex2)
    {
        cout << "Warning: Face-open requires opposite vertices " << endl;
        cout << "Debug  : Face is : " << face->getNodeAt(0)->getID() << " "
                << face->getNodeAt(1)->getID() << " "
                << face->getNodeAt(2)->getID() << " "
                << face->getNodeAt(3)->getID() << endl;
        cout << "Opposite ends are " << vertex0->getID() << " " << vertex2->getID() << endl;
        return 1;
    }

    assert(face->hasNode(vertex0));
    assert(face->hasNode(vertex2));

    vector<Face*> neighs;
    neighs = vertex0->getRelations2();
    for (size_t j = 0; j < neighs.size(); j++)
    {
        if (neighs[j] != face)
        {
            int val0 = neighs[j]->hasNode(vertex0);
            int val1 = neighs[j]->hasNode(vertex2);
            if (val0 + val1 == 2) return 0;
        }
    }

    neighs = vertex2->getRelations2();
    for (size_t j = 0; j < neighs.size(); j++)
    {
        if (neighs[j] != face)
        {
            int val0 = neighs[j]->hasNode(vertex0);
            int val1 = neighs[j]->hasNode(vertex2);
            if (val0 + val1 == 2) return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
FaceClose::build()
{
    if (!isSafe()) return 1;

    assert(face->getSize(0) == 4);

    Vertex *v0 = vertex0;
    Vertex *v2 = vertex2;

    FaceSequence neighs = face->getRelations202();
    NodeSet nodeset;

    for (size_t i = 0; i < neighs.size(); i++)
    {
        Face *f = neighs[i];
        if (f->isRemoved()) return 1;
        for (int j = 0; j < 4; j++)
        {
            Vertex *v = neighs[i]->getNodeAt(j);
            if (v->isRemoved()) return 1;
            nodeset.insert(v);
        }
    }

    NodeSequence localnodes;
    localnodes.reserve(nodeset.size());

    FaceSet faces_to_check;
    map<Vertex*, Point3D> backupCoords;
    set<Vertex*> ::const_iterator it;
    for (it = nodeset.begin(); it != nodeset.end(); ++it)
    {
        Vertex *v = *it;
        backupCoords[v] = v->getXYZCoords();
        localnodes.push_back(v);
        FaceSequence vneighs = v->getRelations2();
        for (int i = 0; i < vneighs.size(); i++)
        {
            Face *f = vneighs[i];
            if (f->isRemoved()) return 1;
            faces_to_check.insert(f);
        }
    }
    faces_to_check.erase(face);

    // Add a new vertex in the mesh...
    Point3D p3d = Vertex::mid_point(v0, v2);

    // Both nodes will merge at the same point...
    vertex0->setXYZCoords(p3d);
    vertex2->setXYZCoords(p3d);

    // Don't move vertex0 and vertex2, because they will be merged.
    int fixnodes[] = {0, 0};
    if (!vertex0->isConstrained())
    {
        vertex0->setConstrainedMark(1);
        fixnodes[0] = 1;
    }

    if (!vertex2->isConstrained())
    {
        vertex2->setConstrainedMark(1);
        fixnodes[1] = 1;
    }

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(localnodes);

    bool pass = 1;
    set<Face*>::const_iterator siter;
    for (siter = faces_to_check.begin(); siter != faces_to_check.end(); ++siter)
    {
        Face *f = *siter;
        if (f->invertedAt() >= 0)
        {
            pass = 0;
            break;
        }
    }

    // Roll back the positions
    if (!pass)
    {
        if (fixnodes[0]) vertex0->setConstrainedMark(0);
        if (fixnodes[1]) vertex2->setConstrainedMark(0);
        map<Vertex*, Point3D>::const_iterator miter;
        for (miter = backupCoords.begin(); miter != backupCoords.end(); ++miter)
        {
            Vertex *v = miter->first;
            v->setXYZCoords(miter->second);
        }
        return 1;
    }

    replacedNode = Vertex::newObject();
    replacedNode->setXYZCoords(p3d);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
FaceClose::commit()
{
    if (replacedNode == NULL) return 1;

    vector<Face*> vr0 = vertex0->getRelations2();
    vector<Face*> vr2 = vertex2->getRelations2();

    for (size_t i = 0; i < vr0.size(); i++)
        if (vr0[i]->isVisited()) return 1;

    for (size_t i = 0; i < vr2.size(); i++)
        if (vr2[i]->isVisited()) return 1;

    cout << " Face Close 2 " << endl;

    int err;
    for (size_t i = 0; i < vr0.size(); i++)
    {
        if (vr0[i] != face)
        {
            err = vr0[i]->replaceNode(vertex0, replacedNode);
            if (!err) mesh->reactivate(vr0[i]);
        }
        vr0[i]->setVisitMark(1);
    }

    for (size_t i = 0; i < vr2.size(); i++)
    {
        if (vr2[i] != face)
        {
            err = vr2[i]->replaceNode(vertex2, replacedNode);
            if (!err) mesh->reactivate(vr2[i]);
        }
        vr2[i]->setVisitMark(1);
    }

    assert(replacedNode);
    mesh->addNode(replacedNode);

    // Two nodes and face go away from the mesh..
    mesh->remove(face);
    mesh->remove(vertex0);
    mesh->remove(vertex2);

    replacedNode = NULL; // So that destructor can delete if this is not used.
    cout << "   CSV Hfdsfsdfsd " << endl;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int QuadCleanUp::remove_tunnels()
{
}

////////////////////////////////////////////////////////////////////////////////

void
Jaal::set_regular_node_tag(Mesh *mesh)
{
    int relexist = mesh->build_relations(0, 2);

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        int nsize = vertex->getRelations2().size();
        if (nsize >= 6) vertex->setTag(0);
        if (nsize == 4) vertex->setTag(1);
        if (nsize < 4) vertex->setTag(2);
        if (nsize == 5) vertex->setTag(3);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face, Vertex *v0, Vertex *v2)
{

    //  Public Function ...
    Point3D p3d = Vertex::mid_point(v0, v2);

    Vertex *doublet = Vertex::newObject();
    doublet->setXYZCoords(p3d);

    Vertex *o1 = NULL, *o2 = NULL;

    vector<Vertex*> connect = face->getNodes();
    for (size_t i = 0; i < connect.size(); i++)
    {
        if (connect[i] == v0)
        {
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

    Face *newquad1 = Face::newObject();
    newquad1->setNodes(connect);
    mesh->addFace(newquad1);

    face->setRemoveMark(1);

    connect[0] = doublet;
    connect[1] = v2;
    connect[2] = o2;
    connect[3] = v0;

    Face *newquad2 = Face::newObject();
    newquad2->setNodes(connect);
    mesh->addFace(newquad2);

    if (mesh->getAdjTable(0, 2))
    {
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
QuadCleanUp::insert_doublet(Face *face)
{

    //  Public Function ...
    if (face->getSize(0) != 4)
        return NULL;
    Vertex *v0 = face->getNodeAt(0);
    Vertex *v2 = face->getNodeAt(2);
    return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_boundary_doublet(Face *face)
{

    //  Public Function ...
    if (!face->isBoundary())
        return NULL;

    if (face->getSize(0) != 4)
        return NULL;

    int ncount = 0;
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = face->getNodeAt(i);
        ncount += v->isBoundary();
    }

    if (ncount != 3)
        return NULL;

    Vertex *v0 = NULL, *v2 = NULL;
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = face->getNodeAt(i);
        if (!v->isBoundary())
        {
            v0 = face->getNodeAt((i + 0) % 4);
            v2 = face->getNodeAt((i + 2) % 4);
        }
    }

    return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

int
Doublet::remove()
{
    if (vertex->isBoundary())
        return 1;

    vector<Face*> neighs = vertex->getRelations2();
    assert(neighs.size() > 0);

    if (neighs.size() != 2)
        return 1;

    assert(neighs[0] != neighs[1]);

    Vertex *d1 = NULL, *d2 = NULL, *o1 = NULL, *o2 = NULL;

    vector<Vertex*> connect = neighs[0]->getNodes();
    if (connect.size() != 4)
        return 1;

    if (neighs[0]->isRemoved())
        return 2;
    if (neighs[1]->isRemoved())
        return 2;

    for (size_t i = 0; i < connect.size(); i++)
    {
        if (connect[i] == vertex)
        {
            d1 = connect[(i + 1) % 4];
            o1 = connect[(i + 2) % 4];
            d2 = connect[(i + 3) % 4];
            break;
        }
    }

    connect = neighs[1]->getNodes();
    if (connect.size() != 4)
        return 1;

    for (size_t i = 0; i < connect.size(); i++)
    {
        if (connect[i] == vertex)
        {
            o2 = connect[(i + 2) % 4];
            break;
        }
    }

    assert(d1);
    assert(d2);
    assert(o1);
    assert(o2);
    assert(o1 != o2);

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

    Face *newquad = Face::newObject();
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
Singlet::update_type1()
{
    if (vertex->getFeatureAngle() < mesh->getMaxFeatureAngle()) return 1;

    if (!active) return 1;

    vector<Face*> vfaces = vertex->getRelations2();

    if (vfaces.size() != 1) return 1;

    Face *f0 = vfaces[0];
    int pos = f0->queryPosOf(vertex);

    Vertex *nextvertex0, *nextvertex1;

    nextvertex0 = f0->getNodeAt((pos + 1) % 4);
    nextvertex1 = f0->getNodeAt((pos + 2) % 4);
    assert(nextvertex0->isBoundary());
    vfaces = nextvertex0->getRelations2();
    if (vfaces.size() > 2)
    {
        type = 1;
        SwapEdge edge(mesh, nextvertex0, nextvertex1, f0);
        int err = edge.apply_singlet_rule(vertex);
        if (!err)
        {
            active = 0;
            return 0;
        }
    }

    nextvertex0 = f0->getNodeAt((pos + 3) % 4);
    nextvertex1 = f0->getNodeAt((pos + 2) % 4);
    assert(nextvertex0->isBoundary());
    vfaces = nextvertex0->getRelations2();

    if (vfaces.size() > 2)
    {
        type = 1;
        SwapEdge edge(mesh, nextvertex0, nextvertex1, f0);
        int err = edge.apply_singlet_rule(vertex);
        if (!err)
        {
            active = 0;
            return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
Singlet::update_type2()
{
    if (vertex->getFeatureAngle() < mesh->getMaxFeatureAngle()) return 1;

    if (!active) return 1;

    if (!vertex->isBoundary()) return 1;

    vector<Face*> vfaces = vertex->getRelations2();

    if (vfaces.size() != 1) return 1;

    Face *f0 = vfaces[0];
    Face *f1 = NULL;

    Vertex * v[6];
    v[0] = vertex;
    v[1] = NULL;
    v[2] = NULL;
    v[3] = NULL;
    v[4] = NULL;
    v[5] = NULL;

    int pos = f0->queryPosOf(vertex);
    Vertex *v1 = f0->getNodeAt((pos + 1) % 4);
    Vertex *v2 = f0->getNodeAt((pos + 2) % 4);
    Vertex *v3 = f0->getNodeAt((pos + 3) % 4);

    vector<Face*> neighs;

    bool found = 0;

    if (!found)
    {
        vfaces = v1->getRelations2();
        if (vfaces.size() == 2)
        {
            neighs = Mesh::getRelations112(v1, v2);
            assert(neighs.size() > 0);
            if (neighs.size() == 2)
            {
                if (neighs[0] == f0) f1 = neighs[1];
                if (neighs[1] == f0) f1 = neighs[0];
                if (!QuadCleanUp::hasSinglet(f1))
                {
                    type = 2;
                    pos = f1->queryPosOf(v1);

                    v[1] = v1;
                    v[2] = v2;
                    v[3] = v3;
                    v[5] = f1->getNodeAt((pos + 2) % 4);

                    if (f1->getNodeAt((pos + 1) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 3) % 4);
                        found = 1;
                    }

                    if (f1->getNodeAt((pos + 3) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 1) % 4);
                        found = 1;
                    }
                }
            } // End checking first vertex
        }
    }

    if (!found)
    {
        vfaces = v3->getRelations2();
        if (vfaces.size() == 2)
        {
            neighs = Mesh::getRelations112(v2, v3);
            assert(neighs.size() > 0);
            if (neighs.size() == 2)
            {
                if (neighs[0] == f0) f1 = neighs[1];
                if (neighs[1] == f0) f1 = neighs[0];
                if (!QuadCleanUp::hasSinglet(f1))
                {
                    type = 2;
                    int pos = f1->queryPosOf(v3);
                    v[1] = v3;
                    v[2] = v2;
                    v[3] = v1;
                    v[5] = f1->getNodeAt((pos + 2) % 4);

                    if (f1->getNodeAt((pos + 1) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 3) % 4);
                        found = 1;
                    }

                    if (f1->getNodeAt((pos + 3) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 1) % 4);
                        found = 1;
                    }
                }
            } // End checking second vertex
        }
    }

    if (!found) return 1;

    if (f0->isVisited() || f1->isVisited()) return 2;

    assert(v[0]->isBoundary());
    assert(v[1]->isBoundary());
    assert(v[3]->isBoundary());

    Point3D xyz;

    // Two new nodes will be created ...
    Vertex *doublet1 = Vertex::newObject();
    xyz = Vertex::mid_point(v[0], v[2]);
    doublet1->setXYZCoords(xyz);
    mesh->addNode(doublet1);

    Vertex *doublet2 = Vertex::newObject();
    xyz = Vertex::mid_point(v[1], v[5]);
    doublet2->setXYZCoords(xyz);
    mesh->addNode(doublet2);

    // Two old faces will be updated and two new faces will be created ...

    Face *newface;
    vector<Vertex*> connect(4);

    connect[0] = v[3];
    connect[1] = v[0];
    connect[2] = doublet1;
    connect[3] = v[2];
    f0->setNodes(connect);
    f0->setVisitMark(1);

    connect[0] = v[1];
    connect[1] = v[4];
    connect[2] = v[5];
    connect[3] = doublet2;
    f1->setNodes(connect);
    f1->setVisitMark(1);

    connect[0] = v[0];
    connect[1] = v[1];
    connect[2] = doublet2;
    connect[3] = doublet1;

    newface = Face::newObject();
    newface->setNodes(connect);
    mesh->addFace(newface);

    connect[0] = doublet1;
    connect[1] = doublet2;
    connect[2] = v[5];
    connect[3] = v[2];

    newface = Face::newObject();
    newface->setNodes(connect);
    mesh->addFace(newface);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Singlet::update_type3()
{
    if (vertex->getFeatureAngle() < mesh->getMaxFeatureAngle()) return 1;

    if (!active) return 1;

    if (!vertex->isBoundary()) return 1;

    vector<Face*> vfaces = vertex->getRelations2();

    if (vfaces.size() != 1) return 1;

    Face *f0 = vfaces[0];

    int pos = f0->queryPosOf(vertex);

    Vertex *v0 = f0->getNodeAt((pos + 0) % 4);
    Vertex *v2 = f0->getNodeAt((pos + 2) % 4);

    if (v2->isBoundary()) return 1;

    Point3D p0 = v0->getXYZCoords();
    Point3D p2 = v2->getXYZCoords();

    vfaces = v2->getRelations2();
    double area0 = 0.0;
    for (size_t i = 0; i < vfaces.size(); i++)
        area0 += vfaces[i]->getArea();

    v2->setXYZCoords(p0);

    bool pass = 1;
    double area1 = 0.0;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        if (vfaces[i] != f0)
        {
            if (!vfaces[i]->isConvex())
            {
                pass = 0;
                break;
            }
            area1 += vfaces[i]->getArea();
        }
    }

    if (pass && fabs(area1 - area0) > 1.0E-06) pass = 0;

    if (!pass)
    {
        v2->setXYZCoords(p2);
        return 2;
    }

    Vertex *vnew = Vertex::newObject();
    vnew->setXYZCoords(p0);
    vnew->setBoundaryMark(vertex->getBoundaryMark());
    newNodes.resize(1);
    newNodes[0] = vnew;

    newFaces.reserve(vfaces.size() - 1);
    vector<Vertex*> connect;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        if (vfaces[i] != f0)
        {
            connect = vfaces[i]->getNodes();
            for (int j = 0; j < 4; j++)
            {
                if (connect[j] == v2)
                    connect[j] = vnew;
            }
            Face *newface = Face::newObject();
            newface->setNodes(connect);
            newFaces.push_back(newface);
        }
    }
    oldFaces = vfaces;

    oldNodes.resize(1);
    oldNodes[0] = vertex;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
Singlet::clear()
{

    for (size_t j = 0; j < newNodes.size(); j++)
        delete newNodes[j];
    newNodes.clear();
    for (size_t j = 0; j < newFaces.size(); j++)
        delete newFaces[j];
    newFaces.clear();

}

///////////////////////////////////////////////////////////////////////////////

int
Singlet::commit()
{
    for (size_t i = 0; i < oldFaces.size(); i++)
    {
        if (oldFaces[i]->isRemoved())
        {
            clear();
            return 1;
        }
    }

    if (!active) return 2;

    for (size_t i = 0; i < oldNodes.size(); i++)
        oldNodes[i]->setRemoveMark(1);
    oldNodes.clear();

    for (size_t i = 0; i < oldFaces.size(); i++)
        oldFaces[i]->setRemoveMark(1);
    oldFaces.clear();

    for (size_t i = 0; i < newNodes.size(); i++)
        mesh->addNode(newNodes[i]);
    newNodes.clear();

    for (size_t i = 0; i < newFaces.size(); i++)
        mesh->addFace(newFaces[i]);
    newFaces.clear();

    return 0;

}

///////////////////////////////////////////////////////////////////////////////

int
RestrictedEdge::build()
{
    assert(connect[0] != NULL);
    assert(connect[1] != NULL);

    Vertex *resnode = connect[0];
    Vertex *bndnode = connect[1];

    assert(!resnode->isBoundary());
    assert(bndnode->isBoundary());

    vector<Face*> vneighs = bndnode->getRelations2();
    if (vneighs.size() != 2) return 2;

    adjFaces[0] = vneighs[0];
    adjFaces[1] = vneighs[1];

    Face *face0 = vneighs[0];
    if (face0->isRemoved()) return 3;

    int pos0 = face0->queryPosOf(bndnode);
    assert(pos0 >= 0);
    Vertex *v0 = face0->getNodeAt((pos0 + 2) % 4);

    Face *face1 = vneighs[1];
    if (face1->isRemoved()) return 4;
    int pos1 = face1->queryPosOf(bndnode);
    assert(pos1 >= 0);
    Vertex *v2 = face1->getNodeAt((pos1 + 2) % 4);

    double area_before = face0->getArea() + face1->getArea();

    // One new node will be inserted ::
    Vertex *vnew = Vertex::newObject();
    Point3D xyz = Vertex::mid_point(v0, v2);
    vnew->setXYZCoords(xyz);
    newNodes.resize(1);
    newNodes[0] = vnew;

    // Three new faces will be created ...
    newFaces.resize(3);

    Face *newface;
    vector<Vertex*> connect(4);

    // Face : 1
    connect[0] = resnode;
    connect[1] = v0;
    connect[2] = vnew;
    connect[3] = v2;

    newface = Face::newObject();
    newface->setNodes(connect);
    newFaces[0] = newface;

    // Face : 2
    Vertex *v3 = NULL;
    if (face0->getNodeAt((pos0 + 1) % 4) == resnode)
        v3 = face0->getNodeAt((pos0 + 3) % 4);

    if (face0->getNodeAt((pos0 + 3) % 4) == resnode)
        v3 = face0->getNodeAt((pos0 + 1) % 4);

    assert(v3);
    connect[0] = vnew;
    connect[1] = v0;
    connect[2] = v3;
    connect[3] = bndnode;

    newface = Face::newObject();
    newface->setNodes(connect);
    newFaces[1] = newface;

    // Face : 3
    Vertex *v4 = NULL;
    if (face1->getNodeAt((pos1 + 1) % 4) == resnode)
        v4 = face1->getNodeAt((pos1 + 3) % 4);

    if (face1->getNodeAt((pos1 + 3) % 4) == resnode)
        v4 = face1->getNodeAt((pos1 + 1) % 4);

    assert(v4);
    connect[0] = vnew;
    connect[1] = bndnode;
    connect[2] = v4;
    connect[3] = v2;

    newface = Face::newObject();
    newface->setNodes(connect);
    newFaces[2] = newface;

    double area_after = 0.0;
    for (int i = 0; i < 3; i++)
        area_after += newFaces[i]->getArea();

    if (fabs(area_after - area_before) > 1.0E-10)
    {
        delete newFaces[0];
        delete newFaces[1];
        delete newFaces[2];
        newFaces.clear();
        return 4;
    }

    for (int i = 0; i < 3; i++)
    {
        if (!newFaces[i]->isConvex())
        {
            delete newFaces[0];
            delete newFaces[1];
            delete newFaces[2];
            newFaces.clear();
            return 5;
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::free_restricted_nodes_once()
{
    /*
      vector<Vertex*> vrestrict = search_restricted_nodes ();

      if (vrestrict.size () == 0) return 0;

      int relexist0 = mesh->build_relations (0, 0);
      int relexist2 = mesh->build_relations (0, 2);

      vector<Vertex*> vneighs;
      vector<RestrictedEdge> redges;

      for (size_t i = 0; i < vrestrict.size (); i++)
        {
          Vertex *vertex = vrestrict[i];
          vneighs = vertex->getRelations0 ();
          for (size_t j = 0; j < vneighs.size (); j++)
            {
              if (vneighs[j]->isBoundary ())
                {
                  RestrictedEdge edge (mesh, vertex, vneighs[j]);
                  int err = edge.build ();
                  if (!err)
                    {
                      redges.push_back (edge);
                    }
                }
            }
        }

      if (redges.size ())
        cout << "Info:  # of valid restricted edges : " << redges.size () << endl;

      int ncount = 0;
      for (size_t i = 0; i < redges.size (); i++)
        {
          int err = redges[i].commit ();
          if (!err) ncount++;
        }

      if (ncount)
        {
          mesh->prune ();
          mesh->enumerate (0);
          mesh->enumerate (2);
          cout << "Info: # of restricted edges committed : " << ncount << endl;
        }

      if (!relexist0)
        mesh->clear_relations (0, 0);

      if (!relexist2)
        mesh->clear_relations (0, 2);

      lapsmooth->execute ();

      return ncount;
     */
    return 1;
}

////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::free_restricted_nodes()
{
    int ncount;
    while (1)
    {
        ncount = free_restricted_nodes_once();
        if (ncount == 0) break;
    }

    if (!mesh->isConsistentlyOriented())
        mesh->makeConsistentlyOriented();
}

/////////////////////////////////////////////////////////////////////////

int
Diamond::isSafe()
{
    vector<Face*> vr0 = vertex0->getRelations2();
    for (size_t i = 0; i < vr0.size(); i++)
        if (vr0[i]->isVisited()) return 0;

    vector<Face*> vr2 = vertex2->getRelations2();
    for (size_t i = 0; i < vr2.size(); i++)
        if (vr2[i]->isVisited()) return 0;

    return 1;

}

/////////////////////////////////////////////////////////////////////////

int
Diamond::makeShield()
{
    vector<Face*> vr0 = vertex0->getRelations2();
    for (size_t i = 0; i < vr0.size(); i++)
        vr0[i]->setVisitMark(1);

    vector<Face*> vr2 = vertex2->getRelations2();
    for (size_t i = 0; i < vr2.size(); i++)
        vr2[i]->setVisitMark(1);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
Diamond::build()
{
    faceclose = new FaceClose(mesh, face, vertex0, vertex2);
    return faceclose->build();
}

///////////////////////////////////////////////////////////////////////////////

int
Diamond::commit()
{
    if (faceclose)
    {
        int err = faceclose->commit();
        delete faceclose;
        faceclose = NULL;
        return err;
    }

    return 1;
}
///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds_once()
{
    int rel0exist = mesh->build_relations(0, 0); // Need for laplace smoothing
    int rel2exist = mesh->build_relations(0, 2);

    // Essential step because only the faces are updated only.
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    int err, pos, ncount = 0;

    size_t faceid = 0;
    while (1)
    {
        if (faceid >= mesh->getSize(2)) break;
        Face *face = mesh->getFaceAt(faceid);
        if (isDiamond(face, pos))
        {
            Diamond diamond(mesh, face, pos);
            err = diamond.build();
            if (!err)
            {
                err = diamond.commit();
                if (!err) ncount++;
            }
        }
        faceid++;
    }
    cout << " NCount " << ncount << endl;

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: number of diamonds removed : " << ncount << endl;
        mesh->collect_garbage();
    }

    if (!rel0exist) mesh->clear_relations(0, 0);
    if (!rel2exist) mesh->clear_relations(0, 2);

    vDiamonds.clear();

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds_in_layer(int layerid)
{
    int rel0exist = mesh->build_relations(0, 0); // Need for laplace smoothing
    int rel2exist = mesh->build_relations(0, 2);

    // Essential step because only the faces are updated only.
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    int err, pos, ncount = 0;

    size_t faceid = 0;
    while (1)
    {
        if (faceid >= mesh->getSize(2)) break;
        Face *face = mesh->getFaceAt(faceid);
        if (face->getLayerID() == layerid)
        {
            if (isDiamond(face, pos))
            {
                Diamond diamond(mesh, face, pos);
                err = diamond.build();
                if (!err)
                {
                    err = diamond.commit();
                    if (!err) ncount++;
                }
            }
        }
        faceid++;
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: number of diamonds removed : " << ncount << endl;
        mesh->collect_garbage();
    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);

    vDiamonds.clear();

    return ncount;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds()
{
    if (!mesh->isPruned()) mesh->prune();

    mesh->search_boundary();

    vector<Doublet> doublets;
    doublets = search_interior_doublets();
    size_t num_doublets_start = doublets.size();
    if (num_doublets_start)
    {
        cout << "Info: Mesh has interior doublets. Removing them " << endl;
        remove_interior_doublets();
        doublets = search_interior_doublets();
        num_doublets_start = doublets.size();
    }

    while (1)
    {
        int ncount = remove_diamonds_once();
        if (ncount == 0) break;
    }

    doublets = search_interior_doublets();
    size_t num_doublets_end = doublets.size();
    if (num_doublets_end > num_doublets_start)
    {
        cout << "Warning: diamonds removal created doublet in the mesh: Now removing them " << endl;
        remove_interior_doublets();
    }
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_doublets_once()
{
    int relexist = mesh->build_relations(0, 2);

    search_interior_doublets();

    int ncount = 0;
    for (size_t i = 0; i < vDoublets.size(); i++)
    {
        int err = vDoublets[i].remove();
        if (!err)
            ncount++;
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: Number of doublets committed : " << ncount << endl;
        mesh->collect_garbage();
    }

    mesh->setBoundaryKnown(0);

    if (!relexist)
        mesh->clear_relations(0, 2);

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_interior_doublets()
{
    if (!mesh->isPruned()) mesh->prune();

    mesh->search_boundary();

    // It is possible that removal of doublet may create singlets on the boundary
    // so, it is better to call doublets first and then call to singlet removal next.

    while (1)
    {
        int ncount = remove_doublets_once();
        if (ncount == 0) break;
    }

    if (vDoublets.empty()) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::cleanup_internal_boundary_face()
{
    int relexist = mesh->build_relations(0, 2);

    size_t numfaces = mesh->getSize(2);

    vector<Face*> boundfaces;

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if (face->hasBoundaryNode())
        {
            int nfnodes = face->getSize(0);
            Vertex *v0 = NULL;
            Vertex *v1 = NULL;
            for (int j = 0; j < nfnodes; j++)
            {
                Vertex *v = face->getNodeAt(j);
                if (v->isBoundary())
                {
                    v0 = face->getNodeAt((j + 1) % nfnodes);
                    v1 = face->getNodeAt((j + 3) % nfnodes);
                    if (!v0->isBoundary() && !v1->isBoundary())
                    {
                        FaceClose closeface(mesh, face, v0, v1);
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
QuadCleanUp::cleanup_boundary(double cutOffAngle)
{
    int rel0exist = mesh->build_relations(0, 0);
    int rel2exist = mesh->build_relations(0, 2);

    mesh->search_boundary();
    mesh->setFeatureLength();

    Point3D pf, pb, pn;
    double flen, len, t;

    size_t numnodes = mesh->getSize(0);

    vector<Vertex*> bound_neighs, free_neighs;
    vector<Face*> vfaces;

    vector<Vertex*> updated;

    for (int iter = 0; iter < 10; iter++)
    {
        updated.clear();

        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *bndnode = mesh->getNodeAt(i);
            if (bndnode->isBoundary())
            {
                bound_neighs = bndnode->getRelations0();
                for (size_t j = 0; j < bound_neighs.size(); j++)
                {
                    Vertex *freenode = bound_neighs[j];
                    if (!freenode->isBoundary())
                    {
                        free_neighs = freenode->getRelations0();
                        int ncount = 0;
                        for (size_t k = 0; k < free_neighs.size(); k++)
                            if (free_neighs[k]->isBoundary()) ncount++;
                        if (ncount == 1)
                        {
                            flen = bndnode->getFeatureLength();
                            len = Vertex::length(bndnode, freenode);
                            if (len > 2.0 * flen)
                            {
                                pf = freenode->getXYZCoords();
                                pb = bndnode->getXYZCoords();
                                t = 0.75;
                                pn[0] = t * pf[0] + (1 - t) * pb[0];
                                pn[1] = t * pf[1] + (1 - t) * pb[1];
                                pn[2] = t * pf[2] + (1 - t) * pb[2];
                                freenode->setXYZCoords(pn);
                                vfaces = freenode->getRelations2();
                                bool pass = 1;
                                for (size_t iface = 0; iface < vfaces.size(); iface++)
                                {
                                    if (!vfaces[iface]->isConvex())
                                    {
                                        pass = 0.0;
                                        break;
                                    }
                                }
                                if (!pass)
                                {
                                    freenode->setXYZCoords(pf);
                                }
                                else
                                    updated.push_back(freenode);
                            }
                        }
                    }
                }
            }
        }

        if (updated.empty()) break;

        for (size_t i = 0; i < updated.size(); i++)
            updated[i]->setConstrainedMark(1);

        lapsmooth->execute();

        for (size_t i = 0; i < updated.size(); i++)
            updated[i]->setConstrainedMark(0);

    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);

    return;


#ifdef CSV

    ///////////////////////////////////////////////////////////////
    // First try to handle flat node...
    ///////////////////////////////////////////////////////////////

    vector<Vertex*> degree2nodes;


    Vertex* node, *onode;
    for (size_t i = 0; i < numnodes; i++)
    {
        node = mesh->getNodeAt(i);
        if (node->isBoundary())
        {
            neighs = node->getRelations0();
            if (neighs.size() == 2)
            {
                if (neighs[0]->isBoundary() && neighs[1]->isBoundary())
                {
                    Point3D v1 = make_vector(neighs[0], node);
                    Point3D v2 = make_vector(neighs[1], node);
                    double angle = Math::getVectorAngle(v1, v2);
                    if (angle > cutOffAngle)
                        degree2nodes.push_back(node);
                }
            }
        }
    }


    relexist = mesh->build_relations(0, 2);

    Face *boundface;
    vector<Face*> faceneighs;
    for (size_t i = 0; i < degree2nodes.size(); i++)
    {
        node = degree2nodes[i];
        faceneighs = node->getRelations2();
        if (faceneighs.size() == 1)
        {
            boundface = faceneighs[0];
            if (boundface->getSize(0) == 4)
            {
                int j = boundface->queryPosOf(node);
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

    mesh->setBoundaryKnown(0);
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

int QuadCleanUp::apply_shift_node3_rule(Vertex *vertex)
{
    int layerID = vertex->getLayerID();

    vector<Face*> vfaces = vertex->getRelations2();
    if (vfaces.size() != 3) return 1;

    Face *face = NULL;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        if (vfaces[i]->getLayerID() == layerID)
        {
            face = vfaces[i];
            break;
        }
    }
    if (face == NULL) return 1;

    int pos = face->queryPosOf(vertex);
    assert(pos >= 0);

    Vertex *opp_vertex = face->getNodeAt((pos + 2) % 4);

    vfaces = opp_vertex->getRelations2();
    if (vfaces.size() == 4)
        return refine_3444_pattern(face, pos);

    if (vfaces.size() == 5)
        return refine_3454_pattern(face, pos);

    if (vfaces.size() > 5)
    {
        while (1)
        {
            int err = apply_advance_front_excess_rule(opp_vertex);
            if (err) break;
        }
        vfaces = opp_vertex->getRelations2();
        if (vfaces.size() < 6)
            return apply_shift_node3_rule(vertex);
    }

    return 2;
}

///////////////////////////////////////////////////////////////////////

int QuadCleanUp::refine_3454_pattern(const Face *face, int pos)
{
    Point3D xyz;
    vector<Vertex*> localnodes(13);

    localnodes[0] = face->getNodeAt((pos + 0) % 4);
    localnodes[1] = face->getNodeAt((pos + 1) % 4);
    localnodes[2] = face->getNodeAt((pos + 2) % 4);
    localnodes[3] = face->getNodeAt((pos + 3) % 4);

    int layerid = localnodes[0]->getLayerID();

    if (localnodes[1]->getLayerID() != layerid) return 1;
    if (localnodes[2]->getLayerID() <= layerid) return 1;
    if (localnodes[3]->getLayerID() != layerid) return 1;

    vector<Face*> vfaces;
    vfaces = localnodes[0]->getRelations2();
    if (vfaces.size() != 3) return 1;

    vfaces = localnodes[1]->getRelations2();
    if (vfaces.size() != 4) return 1;

    vfaces = localnodes[2]->getRelations2();
    if (vfaces.size() != 5) return 1;

    vfaces = localnodes[3]->getRelations2();
    if (vfaces.size() != 4) return 1;

    Face *neigh1 = NULL;
    Face *neigh2 = NULL;

    vector<Face*> adjFaces;

    adjFaces = Mesh::getRelations112(localnodes[1], localnodes[2]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[4] = Face::opposite_node(neigh1, localnodes[2]);
    localnodes[5] = Face::opposite_node(neigh1, localnodes[1]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[5]);
    localnodes[11] = Vertex::newObject();
    localnodes[11]->setXYZCoords(xyz);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[5]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == neigh1) neigh2 = adjFaces[1];
    if (adjFaces[1] == neigh1) neigh2 = adjFaces[0];
    localnodes[8] = Face::opposite_node(neigh2, localnodes[2]);
    localnodes[9] = Face::opposite_node(neigh2, localnodes[5]);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[3]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[6] = Face::opposite_node(neigh1, localnodes[3]);
    localnodes[7] = Face::opposite_node(neigh1, localnodes[2]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[6]);
    localnodes[12] = Vertex::newObject();
    localnodes[12]->setXYZCoords(xyz);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[6]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == neigh1) neigh2 = adjFaces[1];
    if (adjFaces[1] == neigh1) neigh2 = adjFaces[0];
    localnodes[10] = Face::opposite_node(neigh2, localnodes[2]);

    Face * newfaces[7];
    vector<Vertex*> connect(4);

    connect[0] = localnodes[0];
    connect[1] = localnodes[1];
    connect[2] = localnodes[11];
    connect[3] = localnodes[2];
    newfaces[0] = Face::newObject();
    newfaces[0]->setNodes(connect);

    connect[0] = localnodes[1];
    connect[1] = localnodes[4];
    connect[2] = localnodes[5];
    connect[3] = localnodes[11];
    newfaces[1] = Face::newObject();
    newfaces[1]->setNodes(connect);

    connect[0] = localnodes[11];
    connect[1] = localnodes[5];
    connect[2] = localnodes[8];
    connect[3] = localnodes[9];
    newfaces[2] = Face::newObject();
    newfaces[2]->setNodes(connect);

    connect[0] = localnodes[0];
    connect[1] = localnodes[2];
    connect[2] = localnodes[12];
    connect[3] = localnodes[3];
    newfaces[3] = Face::newObject();
    newfaces[3]->setNodes(connect);

    connect[0] = localnodes[3];
    connect[1] = localnodes[12];
    connect[2] = localnodes[6];
    connect[3] = localnodes[7];
    newfaces[4] = Face::newObject();
    newfaces[4]->setNodes(connect);

    connect[0] = localnodes[6];
    connect[1] = localnodes[12];
    connect[2] = localnodes[9];
    connect[3] = localnodes[10];
    newfaces[5] = Face::newObject();
    newfaces[5]->setNodes(connect);

    connect[0] = localnodes[2];
    connect[1] = localnodes[11];
    connect[2] = localnodes[9];
    connect[3] = localnodes[12];
    newfaces[6] = Face::newObject();
    newfaces[6]->setNodes(connect);

    // Do backup of Coordinates. Only 11 nodes to be backed up.
    Point3D backupCoords[11];
    for (int i = 0; i < 11; i++)
        backupCoords[i] = localnodes[i]->getXYZCoords();

    Face * backupFaces[5];
    vfaces = localnodes[2]->getRelations2();
    for (int i = 0; i < 5; i++)
    {
        backupFaces[i] = vfaces[i];
        mesh->remove(vfaces[i]); // Goes to garbage but not deallocated.
    }

    // Update new nodes and faces ...
    mesh->addNode(localnodes[11]);
    mesh->addNode(localnodes[12]);
    for (int i = 0; i < 7; i++)
        mesh->addFace(newfaces[i]);

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(localnodes);

    set<Face*> faces_to_check;
    for (int i = 0; i < 13; i++)
    {
        vector<Face*> vfaces = localnodes[i]->getRelations2();
        for (int j = 0; j < vfaces.size(); j++)
            faces_to_check.insert(vfaces[j]);
    }

    bool pass = 1;
    set<Face*>::const_iterator siter;
    for (siter = faces_to_check.begin(); siter != faces_to_check.end(); ++siter)
    {
        Face *f = *siter;
        if (f->invertedAt() >= 0)
        {
            pass = 0;
            break;
        }
    }

    if (!pass)
    {
        for (int i = 0; i < 7; i++)
            mesh->remove(newfaces[i]);

        for (int i = 0; i < 5; i++)
            mesh->addFace(backupFaces[i]); // Reactivated now ...

        for (int i = 0; i < 11; i++)
            localnodes[i]->setXYZCoords(backupCoords[i]);

        mesh->remove(localnodes[11]);
        mesh->remove(localnodes[12]);
        return 1;
    }

    // Update front levels ...
    localnodes[11]->setLayerID(layerid + 1);
    localnodes[12]->setLayerID(layerid + 1);

    newfaces[0]->setLayerID(layerid);
    newfaces[1]->setLayerID(layerid);
    newfaces[2]->setLayerID(layerid + 1);

    newfaces[3]->setLayerID(layerid);
    newfaces[4]->setLayerID(layerid);
    newfaces[5]->setLayerID(layerid + 1);

    newfaces[6]->setLayerID(layerid + 1);

    return 0;
}

///////////////////////////////////////////////////////////////////////

int QuadCleanUp::refine_3444_pattern(const Face *face, int pos)
{
    Point3D xyz;
    vector<Vertex*> localnodes(12);

    localnodes[0] = face->getNodeAt((pos + 0) % 4);
    localnodes[1] = face->getNodeAt((pos + 1) % 4);
    localnodes[2] = face->getNodeAt((pos + 2) % 4);
    localnodes[3] = face->getNodeAt((pos + 3) % 4);

    int layerid = localnodes[0]->getLayerID();

    if (localnodes[1]->getLayerID() != layerid) return 1;
    if (localnodes[2]->getLayerID() <= layerid) return 1;
    if (localnodes[3]->getLayerID() != layerid) return 1;

    vector<Face*> vfaces;
    vfaces = localnodes[0]->getRelations2();
    if (vfaces.size() != 3) return 1;

    vfaces = localnodes[1]->getRelations2();
    if (vfaces.size() != 4) return 1;

    vfaces = localnodes[2]->getRelations2();
    if (vfaces.size() != 4) return 1;

    vfaces = localnodes[3]->getRelations2();
    if (vfaces.size() != 4) return 1;

    Face *neigh1 = NULL;
    Face *neigh2 = NULL;

    FaceSequence adjFaces;

    adjFaces = Mesh::getRelations112(localnodes[1], localnodes[2]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[4] = Face::opposite_node(neigh1, localnodes[2]);
    localnodes[5] = Face::opposite_node(neigh1, localnodes[1]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[5]);
    localnodes[9] = Vertex::newObject();
    localnodes[9]->setXYZCoords(xyz);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[5]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == neigh1) neigh2 = adjFaces[1];
    if (adjFaces[1] == neigh1) neigh2 = adjFaces[0];
    localnodes[8] = Face::opposite_node(neigh2, localnodes[2]);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[3]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[6] = Face::opposite_node(neigh1, localnodes[3]);
    localnodes[7] = Face::opposite_node(neigh1, localnodes[2]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[6]);
    localnodes[11] = Vertex::newObject();
    localnodes[11]->setXYZCoords(xyz);

    xyz = Vertex::mid_point(localnodes[2], localnodes[8]);
    localnodes[10] = Vertex::newObject();
    localnodes[10]->setXYZCoords(xyz);

    Face * newfaces[7];
    vector<Vertex*> connect(4);

    connect[0] = localnodes[0];
    connect[1] = localnodes[1];
    connect[2] = localnodes[9];
    connect[3] = localnodes[2];
    newfaces[0] = Face::newObject();
    newfaces[0]->setNodes(connect);

    connect[0] = localnodes[1];
    connect[1] = localnodes[4];
    connect[2] = localnodes[5];
    connect[3] = localnodes[9];
    newfaces[1] = Face::newObject();
    newfaces[1]->setNodes(connect);

    connect[0] = localnodes[9];
    connect[1] = localnodes[5];
    connect[2] = localnodes[8];
    connect[3] = localnodes[10];
    newfaces[2] = Face::newObject();
    newfaces[2]->setNodes(connect);

    connect[0] = localnodes[0];
    connect[1] = localnodes[2];
    connect[2] = localnodes[11];
    connect[3] = localnodes[3];
    newfaces[3] = Face::newObject();
    newfaces[3]->setNodes(connect);

    connect[0] = localnodes[3];
    connect[1] = localnodes[11];
    connect[2] = localnodes[6];
    connect[3] = localnodes[7];
    newfaces[4] = Face::newObject();
    newfaces[4]->setNodes(connect);

    connect[0] = localnodes[6];
    connect[1] = localnodes[11];
    connect[2] = localnodes[10];
    connect[3] = localnodes[8];
    newfaces[5] = Face::newObject();
    newfaces[5]->setNodes(connect);

    connect[0] = localnodes[2];
    connect[1] = localnodes[9];
    connect[2] = localnodes[10];
    connect[3] = localnodes[11];
    newfaces[6] = Face::newObject();
    newfaces[6]->setNodes(connect);

    Point3D backupCoords[9];
    for (int i = 0; i < 9; i++)
        backupCoords[i] = localnodes[i]->getXYZCoords();

    Face * backupFaces[4];
    vfaces = localnodes[2]->getRelations2();
    for (int i = 0; i < 4; i++)
    {
        backupFaces[i] = vfaces[i];
        mesh->remove(vfaces[i]); // Send to garbage, but don't delete ..
    }

    // Update the data structures ...
    mesh->addNode(localnodes[9]);
    mesh->addNode(localnodes[10]);
    mesh->addNode(localnodes[11]);
    for (int i = 0; i < 7; i++)
        mesh->addFace(newfaces[i]);

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(localnodes);

    FaceSet faces_to_check;
    for (int i = 0; i < 12; i++)
    {
        vector<Face*> vfaces = localnodes[i]->getRelations2();
        for (int j = 0; j < vfaces.size(); j++)
            faces_to_check.insert(vfaces[j]);
    }

    bool pass = 1;
    set<Face*>::const_iterator siter;
    for (siter = faces_to_check.begin(); siter != faces_to_check.end(); ++siter)
    {
        Face *f = *siter;
        if (f->invertedAt() >= 0)
        {
            pass = 0;
            break;
        }
    }

    if (!pass)
    {
        for (int i = 0; i < 7; i++)
            mesh->remove(newfaces[i]);

        for (int i = 0; i < 4; i++)
            mesh->addFace(backupFaces[i]);

        for (int i = 0; i < 9; i++)
            localnodes[i]->setXYZCoords(backupCoords[i]);

        mesh->remove(localnodes[9]);
        mesh->remove(localnodes[10]);
        mesh->remove(localnodes[11]);

        return 1;
    }

    // Update front levels ...
    localnodes[9]->setLayerID(layerid + 1);
    localnodes[10]->setLayerID(layerid + 2);
    localnodes[11]->setLayerID(layerid + 1);

    newfaces[0]->setLayerID(layerid);
    newfaces[1]->setLayerID(layerid);
    newfaces[2]->setLayerID(layerid + 1);

    newfaces[3]->setLayerID(layerid);
    newfaces[4]->setLayerID(layerid);
    newfaces[5]->setLayerID(layerid + 1);

    newfaces[6]->setLayerID(layerid + 1);

    return 0;
}

///////////////////////////////////////////////////////////////////////

int QuadCleanUp::refine_degree3_faces()
{
    int relexist2 = mesh->build_relations(0, 2);

    size_t numnodes = mesh->getSize(0);

    Vertex *vertex;
    vector<Face*> vfaces;
    for (size_t i = 0; i < numnodes; i++)
    {
        vertex = mesh->getNodeAt(i);
        vfaces = vertex->getRelations2();
        if (vfaces.size() == 3)
        {
            for (int j = 0; j < vfaces.size(); j++)
            {
                vertex = vfaces[j]->getNodeAt(0);
                if (vertex->getRelations2().size() > 4) continue;

                vertex = vfaces[j]->getNodeAt(1);
                if (vertex->getRelations2().size() > 4) continue;

                vertex = vfaces[j]->getNodeAt(2);
                if (vertex->getRelations2().size() > 4) continue;

                vertex = vfaces[j]->getNodeAt(3);
                if (vertex->getRelations2().size() > 4) continue;

                mesh->refine_quad15(vfaces[j]);
            }
        }
    }
    mesh->prune();
    mesh->collect_garbage();

    if (!relexist2)
        mesh->clear_relations(0, 2);

    return 0;
}

///////////////////////////////////////////////////////////////////////

void QuadCleanUp::report()
{
    if (mesh == NULL) return;
    cout << "Info: Reporting mesh information " << endl;

}

///////////////////////////////////////////////////////////////////////

void QuadCleanUp::build_irregular_nodes_set()
{
    // Get all the irregular nodes from the mesh ( only the internal ones );
    irregular_nodes_set.clear();
    NodeSequence nset = mesh->get_irregular_nodes(4);
    for (size_t i = 0; i < nset.size(); i++)
        irregular_nodes_set.insert(nset[i]);
}
///////////////////////////////////////////////////////////////////////

int
QuadCleanUp::automatic()
{
    int err, stage = 1;

    MeshOpt mopt;
    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(100);
    //
    // If the mesh is triangular first convert it into Quad mesh ....
    //
    if (mesh->isHomogeneous() == 3)
    {
        Tri2Quads t2quad;
        t2quad.getQuadMesh(mesh, 1);
    }

    // Throughout the cleaning process, euler characteristic should remain same
    int euler1, euler0 = mesh->getEulerCharacteristic(); // Invariant
    // Total area of the domain must be same ...
    double area0 = mesh->getSurfaceArea();

    //  Ensure that there are irregular nodes in the mesh. if not, you are lucky and done.
    NodeSequence irreg_nodes = mesh->get_irregular_nodes(4);

    assert(mesh->isHomogeneous() == 4);
    cout << " Input Mesh :    " << endl;
    cout << "      # Nodes : " << mesh->getSize(0) << endl;
    cout << "      # Faces : " << mesh->getSize(2) << endl;
    cout << "      # Components : " << mesh->getNumOfComponents() << endl;
    cout << "      # of irregular nodes " << irreg_nodes.size() << endl;
    cout << "      # of concave faces   " << mesh->count_concave_faces() << endl;
    cout << "      Euler Characteristics : " << euler0 << endl;
    cout << "      Surface Area : " << area0 << endl;
    cout << "******************************************************************" << endl;

    //  Mesh connectivity must be consistent:  Condition: Strict
    if (!mesh->isConsistentlyOriented())
        mesh->makeConsistentlyOriented();

    //  Initial mesh may have different connectivity. Condition Strict
    size_t ninvert = mesh->count_inverted_faces();
    size_t numfaces = mesh->getSize(2);
    if (ninvert > 0.5 * numfaces) mesh->reverse();

    //  Initial mesh may have  doublets, remove them, they are troublesome: Condition Strict
    err = remove_interior_doublets();
    if (err)
    {
        cout << "Fatal Error: There are interior doublets in the mesh " << endl;
        mesh->saveAs("dbg.dat");
        exit(0);
        return stage;
    }
    cout << "Info:  Stage " << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces() << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert(euler0 == euler1);
    cout << "******************************************************************" << endl;

    //  Check the boundary nodes connectivity, and ensure that all elements adjacent
    //  to them are convex. Singlet must be called after doublets removal ...
    err = remove_boundary_singlets();
    swap_concave_faces();
    lapsmooth.execute();
    mopt.shape_optimize(mesh);
    swap_concave_faces();

    cout << "Info:  Stage " << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces() << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert(euler0 == euler1);
    cout << "******************************************************************" << endl;

    //  Tunnels may put constraints in the movements, so it is better to remove them,
    //  Condition: Soft
    err = remove_tunnels();
    cout << "Info:  Stage " << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces() << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert(euler0 == euler1);
    cout << "******************************************************************" << endl;

    //  Perform some local operations: Condition: Soft.
    err = remove_diamonds();
    cout << "Info:  Stage " << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces() << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert(euler0 == euler1);
    cout << "******************************************************************" << endl;
    mesh->saveAs("dbg.dat");
    exit(0);

    //  Perform  vertex degree reduction with local edge swapping: Soft.
    err = vertex_degree_reduction();
    cout << "Info:  Stage " << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces() << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert(euler0 == euler1);
    cout << "******************************************************************" << endl;

    //  Now it is time to improve geometry;
    if (mesh->count_concave_faces())
        mopt.untangle(mesh);
    mopt.shape_optimize(mesh);

    if (!has_interior_nodes_degree_345())
    {
        cout << "Warning: Some of the nodes has vertex degree outside range" << endl;
        cout << "         threfore, global optimization Skipped " << endl;
        return 1;
    }

    cout << "Info:  Stage " << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces() << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert(euler0 == euler1);
    cout << "******************************************************************" << endl;

    // Perform Global remeshing to elimate 3 and 5 degree nodes ..
    err = remesh_defective_patches();

    return 0;
}


/*
int
QuadCleanUp::clean_layer_once (int layernum)
{
  vector<Doublet> doublets = search_interior_doublets ();
  assert (doublets.empty ());

  mesh->search_boundary ();

  assert (!mesh->getAdjTable (0, 0));
  assert (!mesh->getAdjTable (0, 2));

  mesh->setWavefront (0);
  mesh->setWavefront (2);

  size_t numnodes = mesh->getSize (0);
  size_t numfaces = mesh->getSize (2);

  for (size_t i = 0; i < numfaces; i++)
    {
      Face *face = mesh->getFaceAt (i);
      face->setVisitMark (0);
      assert (!face->isRemoved ());
    }

  int relexist = mesh->build_relations (0, 2);

  vDiamonds.clear ();

  vector<Face*> neighs;
  Vertex *currnode, *oppnode;

  numnodes = mesh->getSize (0);

  int ncount_removable = 0;
  for (size_t i = 0; i < numnodes; i++)
    {
      Vertex *currnode = mesh->getNodeAt (i);

      int id = currnode->getLayerID ();
      if (id == layernum)
        {
          if (currnode->isBoundary ())
            if (currnode->getFeatureAngle () < 120) continue;

          neighs = currnode->getRelations2 ();

          for (size_t j = 0; j < neighs.size (); j++)
            {
              Face *face = neighs[j];
              if (face->getLayerID () >= layernum + 1)
                {
                  int k = face->queryPosOf (currnode);
                  Vertex *vside0 = face->getNodeAt ((k + 1) % 4);
                  Vertex *vside1 = face->getNodeAt ((k + 3) % 4);
                  if ((vside0->getLayerID () == layernum + 1) &&
                      (vside1->getLayerID () == layernum + 1) &&
                      (!vside0->isBoundary () && !vside1->isBoundary ()))
                    {
                      Diamond diamond (mesh, face, k + 1);
                      if (diamond.isSafe ())
                        {
                          diamond.makeShield ();
                          vDiamonds.push_back (diamond);
                        }
                    }
                }
            }
        }
    }

  int ncount = remove_diamonds_once ();

  for (size_t i = 0; i < vDiamonds.size (); i++)
    {
      Vertex *vertex = vDiamonds[i].getNewNode ();
      if (vertex)
        vertex->setLayerID (layernum + 1);
    }

  if (!relexist) mesh->clear_relations (0, 2);

  doublets = search_interior_doublets ();
  assert (doublets.empty ());

  lapsmooth->execute ();

  return ncount;
}

////////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::clean_layer (int layernum)
{
  if (mesh->isHomogeneous () != 4)
    {
      cout << "Warning: layer cleaning only for quad mesh " << endl;
      return 1;
    }

  int ncount = 0;
  while (1)
    {
      int n = clean_layer_once (layernum);
      if (n == 0) break;
      ncount += n;
    }

  return ncount;
}
 */
////////////////////////////////////////////////////////////////////////////////

/*  To discard 
void
QuadCleanUp::advancing_front_cleanup ()
{
  remove_interior_doublets ();

  int numlayers = mesh->setWavefront (2);

  double area0 = mesh->getSurfaceArea ();

  double minlen;
  size_t numnodes;
  Point3D head, tail, uvec;
  map<Vertex*, double> vertex_feature_len;

  numlayers = 2;

  for (int layernum = 0; layernum < numlayers; layernum++)
    {
     clean_layer (layernum);
      mesh->build_relations (0, 0);

      vertex_feature_len.clear ();

      numnodes = mesh->getSize (0);
      for (size_t j = 0; j < numnodes; j++)
        {
          Vertex *vertex = mesh->getNodeAt (j);
          if (vertex->getLayerID () == layernum)
            {
              vector<Vertex*> vneighs = vertex->getRelations0 ();
              minlen = MAXDOUBLE;
              for (size_t k = 0; k < vneighs.size (); k++)
                {
                  Vertex *vneigh = vneighs[k];
                  if (vneighs[k]->getLayerID () == layernum)
                    minlen = std::min (minlen, Vertex::length2 (vertex, vneigh));
                }
              vertex_feature_len[vertex] = sqrt (minlen);
            }
        }

      for (size_t j = 0; j < numnodes; j++)
        {
          Vertex *vertex = mesh->getNodeAt (j);
          if (vertex->getLayerID () == layernum + 1)
            {
              vector<Vertex*> vneighs = vertex->getRelations0 ();

              int ncount = 0;
              for (size_t k = 0; k < vneighs.size (); k++)
                {
                  Vertex *vneigh = vneighs[k];
                  if (vneigh->getLayerID () == layernum) ncount++;
                }

              if (ncount == 1)
                {
                  head = vertex->getXYZCoords ();
                  for (size_t k = 0; k < vneighs.size (); k++)
                    {
                      Vertex *vneigh = vneighs[k];
                      if (vneigh->getLayerID () == layernum)
                        {
                          double len = Vertex::length (vertex, vneigh);
                          double minlen = vertex_feature_len[vneigh];
                          if (len > 0.7 * minlen)
                            {
                              tail = vneigh->getXYZCoords ();
                              uvec = Math::unit_vector (head, tail);
                              head[0] = tail[0] + 0.7 * minlen * uvec[0];
                              head[1] = tail[1] + 0.7 * minlen * uvec[1];
                              head[2] = tail[2] + 0.7 * minlen * uvec[2];
                              vertex->setXYZCoords (head);
                            }
                        }
                    }
                  vertex->setConstrainedMark (1);
                }

            }
        }
      lapsmooth->execute ();
      mesh->clear_relations (0, 0);
    }
}
 */

