#include "QuadCleanUp.hpp"

using namespace Jaal;
///////////////////////////////////////////////////////////////////////////////

int
QuadEdge::commit()
{
    if (newFaces.empty()) return 1;

    int numfaces = newFaces.size();
    int numnodes = newNodes.size();

    for (int i = 0; i < numfaces; i++)
        if (newFaces[i] == NULL) return 2;

    if (adjFaces[0]->isRemoved() || adjFaces[1]->isRemoved())
    {
        for (int i = 0; i < numnodes; i++)
            delete newNodes[i];
        newNodes.clear();

        for (int i = 0; i < numfaces; i++)
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

    adjFaces[0]->setStatus(MeshEntity::REMOVE);
    adjFaces[1]->setStatus(MeshEntity::REMOVE);

    assert(mesh);

    for (int i = 0; i < numnodes; i++)
        mesh->addNode(newNodes[i]);

    for (int i = 0; i < numfaces; i++)
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

///////////////////////////////////////////////////////////////////////////////

bool
SwapQuadEdge::is_topologically_valid_swap(int d1, int d2, int d3, int d4)
{
    if (d1 < 4 || d2 < 4) return 0;
    if ((d1 > 4) && (d2 > 4) && (d3 == 3) && (d4 == 3)) return 1;
    if ((d1 == 5) && (d2 == 5) && (d3 == 3) && (d4 == 4)) return 1;
    if ((d1 == 5) && (d2 == 5) && (d3 == 4) && (d4 == 3)) return 1;
    if (max(d1, d2) > max(d3 + 1, d4 + 1)) return 1;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int
SwapQuadEdge::getPosOf(const Vertex *vertex) const
{
    for (int i = 0; i < 6; i++)
        if (bound_nodes[i] == vertex) return i;

    return -1;
}
/////////////////////////////////////////////////////////////////////////////

int
SwapQuadEdge::build_boundary()
{
    assert(mesh);

    assert(connect[0]);
    assert(connect[1]);

    adjFaces.resize(2);
    adjFaces[0] = NULL;
    adjFaces[1] = NULL;

    assert(mesh->getAdjTable(0, 2));

    FaceSequence nghs;
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

    // Possibility of doublet creation is ruled out now..

    Mesh::getRelations112(connect[0], connect[1], nghs);

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
SwapQuadEdge::get_boundary_nodes_chain()
{
    vector<Edge> bndedges;
    bndedges.reserve(6);

    Edge sharededge(connect[0], connect[1]);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Vertex *vf0 = adjFaces[i]->getNodeAt(j + 0);
            Vertex *vf1 = adjFaces[i]->getNodeAt(j + 1);
            Edge edge(vf0, vf1);
            if (!edge.isSameAs(sharededge))
                bndedges.push_back(edge);
        }
    }

    Mesh::make_chain(bndedges);
    Mesh::rotate_chain(bndedges, connect[0]);
    bound_nodes = Mesh::chain_nodes(bndedges);
    assert(bound_nodes.size() == 6);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void SwapQuadEdge::update_front()
{
    size_t nSize;
    if (check_fronts)
    {
        NodeSequence vneighs;
        while (1)
        {
            int values_updated = 0;
            for (int i = 0; i < 6; i++)
            {
                vneighs = bound_nodes[i]->getRelations0();
                int minid = bound_nodes[i]->getLayerID();
                nSize   = vneighs.size();
                for (size_t j = 0; j < nSize; j++)
                    minid = min(minid, vneighs[j]->getLayerID() + 1);

                if (bound_nodes[i]->getLayerID() != minid)
                {
                    bound_nodes[i]->setLayerID(minid);
                    values_updated = 1;
                }
            }
            if (!values_updated) break;
        }
    }

}

////////////////////////////////////////////////////////////////////////////////

void
SwapQuadEdge::backup()
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
SwapQuadEdge::rollback()
{
    // We need to manual delete the relation
    connect[0]->removeRelation0(connect[1]);
    connect[1]->removeRelation0(connect[0]);

    // remove all relations by deactivating the elements ...
    mesh->deactivate(adjFaces[0]);
    mesh->deactivate(adjFaces[1]);

    // Change the connectivity of the two quads and reactivate the element...
    adjFaces[0]->setNodes(bkp_data.face1Connect);
    adjFaces[1]->setNodes(bkp_data.face2Connect);

    mesh->reactivate(adjFaces[0]);
    mesh->reactivate(adjFaces[1]);

    connect[0] = bkp_data.diagonalConnect[0];
    connect[1] = bkp_data.diagonalConnect[1];

    for (int i = 0; i < 6; i++)
    {
        Vertex *v = bound_nodes[i];
        v->setXYZCoords(bkp_data.pCoords[v]);
    }

    assert(!adjFaces[0]->isRemoved());
    assert(!adjFaces[1]->isRemoved());

    update_front();
}
////////////////////////////////////////////////////////////////////////////////

int
SwapQuadEdge::make_new_diagonal_at(int pos, bool bound_check)
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

    mesh->deactivate(adjFaces[0]);
    mesh->deactivate(adjFaces[1]);

    // Change the connectivity of the two quads.
    NodeSequence qConnect(4);
    qConnect[0] = bound_nodes[(pos + 0) % 6];
    qConnect[1] = bound_nodes[(pos + 1) % 6];
    qConnect[2] = bound_nodes[(pos + 2) % 6];
    qConnect[3] = bound_nodes[(pos + 3) % 6];
    connect[0] = qConnect[0];
    adjFaces[0]->setNodes(qConnect);
    mesh->reactivate(adjFaces[0]);

    qConnect[0] = bound_nodes[(pos + 3) % 6];
    qConnect[1] = bound_nodes[(pos + 4) % 6];
    qConnect[2] = bound_nodes[(pos + 5) % 6];
    qConnect[3] = bound_nodes[(pos + 6) % 6];
    connect[1] = qConnect[0];
    adjFaces[1]->setNodes(qConnect);
    mesh->reactivate(adjFaces[1]);

    assert(!adjFaces[0]->isRemoved());
    assert(!adjFaces[1]->isRemoved());

    int con1 = adjFaces[0]->concaveAt();
    int con2 = adjFaces[1]->concaveAt();
//  int con1 = adjFaces[0]->isSimple();
//  int con2 = adjFaces[1]->isSimple();
//  if( con1 && con2 ) return  0;

    if( con1 >= 0 || con2 >=0) {
        rollback();
        return 1;
    }
    return 0;

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
        FaceSet neighSet;
        for (int i = 0; i < 6; i++)
        {
            Vertex *v = bound_nodes[i];
            const FaceSequence &vfaces = v->getRelations2();
            size_t nSize = vfaces.size();
            for (size_t j = 0; j < nSize; j++)
                neighSet.insert(vfaces[j]);
        }
        assert(!neighSet.empty());

        FaceSet::const_iterator it;
        for (it = neighSet.begin(); it != neighSet.end(); ++it)
        {
            Face *face = *it;
            int pos = face->concaveAt();
            if (pos >= 0)
            {
                Vertex *v = face->getNodeAt(pos);
                if (!v->isBoundary())
                {
                    rollback();
                    return 1;
                }
            }
        }
    update_front();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
SwapQuadEdge::apply_reduce_degree_rule()
{
    if (build_boundary() != 0) return 1;

    NodeSequence neighs;

    neighs = connect[0]->getRelations0();
    int d1 = neighs.size();
    if (d1 < 4) return 1;

    neighs = connect[1]->getRelations0();
    int d2 = neighs.size();
    if (d2 < 4) return 1;

    int start_pos = getPosOf(connect[0]);
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

        if (SwapQuadEdge::is_topologically_valid_swap(d1, d2, d3, d4))
        {
            int err = make_new_diagonal_at(pos);
            if (!err) return 0;
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
SwapQuadEdge::apply_concave_rule()
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
SwapQuadEdge::apply_bound_rule()
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
SwapQuadEdge::apply_advance_front_rule()
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
SwapQuadEdge::apply_singlet_rule(Vertex *singlet)
{
    // The first node must be boundary and the other node must be internal..
    assert(connect[0]->isBoundary());

    if (connect[1]->isBoundary()) return 1;

    if (build_boundary() != 0) return 1;
    ////////////////////////////////////////////////////////////////////////////
    // Objective :: Swap quads common diagonal such that the resulting diagonal
    //              contains the Singlet node.
    ////////////////////////////////////////////////////////////////////////////
    assert(QuadCleanUp::isSinglet(singlet));

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
SwapQuadEdge::apply_deficient_rule(Vertex *vertex)
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
    if( d2 < 4 ) return 1;

    if (d1 == 4)
    {
        FaceSequence faces = vertex->getRelations2();
        if (faces.size() != 3) return 1;
        int layerid = connect[0]->getLayerID();
        faces = vertex->getRelations2();
        Face *f0 = firstFace;
        Mesh::getRelations112(connect[0], connect[1], faces);
        assert(faces.size() == 2);
        Face *f1 = NULL;
        if (faces[0] == f0) f1 = faces[1];
        if (faces[1] == f0) f1 = faces[0];
        assert(f1);


        int pos = f1->getPosOf(connect[0]);
        Vertex *vopp = f1->getNodeAt(pos + 2);
        if (connect[1]->getLayerID() > layerid && vopp->getLayerID() > layerid)
        {
            set_no_tags(mesh);
            connect[1]->setTag(1);
            vopp->setTag(1);
            mesh->saveAs("dbg.dat");
            Break();
            SwapQuadEdge edge(mesh, connect[1], vopp, f1);
            int err = edge.apply_deficient_rule(connect[0]);
            if (err) return 1;
        }
    }
    exit(0);

    d1 = connect[0]->getRelations2().size();
    d2 = connect[1]->getRelations2().size();

    // Don't create doublet at two nodes...
    if( d1 == 3 || d2 == 3 ) return 1;

    // Having these conditions set, now build the structure.
    if (build_boundary() != 0) return 3;

    update_front();

    //Find the position of triplet in the contour ...
    int pos = this->getPosOf(vertex);
    assert(pos == 1 || pos == 5);

    if (check_fronts)
    {
        int layerid = vertex->getLayerID();
        if (bound_nodes[ (pos + 3) % 6]->getLayerID() <= layerid) return 1;

        if( connect[0]->getLayerID() <= layerid && d1 < 5 ) return 1;
        if( connect[1]->getLayerID() <= layerid && d2 < 5 ) return 1;
    }

    // Create new quads whose common diagonal must contain the singlet.
    int err = make_new_diagonal_at(pos);

    /*
        if( err  ) {
            for (int i = 0; i < 6; i++) {
                 NodeSequence vneighs = bound_nodes[i]->getRelations0();
                 int minid = bound_nodes[i]->getLayerID();
                 for( int j = 0; j < vneighs.size(); j++)
                      minid = min( minid, vneighs[j]->getLayerID() );
                assert( bound_nodes[i]->getLayerID() == minid+1);
             }
        }
    */

    return err;
}

/////////////////////////////////////////////////////////////////////////////

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

    size_t ncount = 0;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        int pos = face->concaveAt();
        if (pos >= 0)
        {
            for (int i = 0; i < 2; i++)
            {
                Vertex *v0 = face->getNodeAt(pos + 1 + i);
                Vertex *v1 = face->getNodeAt(pos + 2 + i);
                SwapQuadEdge edge(mesh, v0, v1);
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

    if (!relexist0) mesh->clear_relations(0, 0);
    if (!relexist2) mesh->clear_relations(0, 2);

    if (mesh->count_concave_faces() == 0)
        mopt.shape_optimize(mesh);

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

