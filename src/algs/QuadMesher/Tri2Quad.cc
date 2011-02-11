#include <meshkit/Tri2Quad.h>

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int Tri2Quads::verify(Mesh *mesh, const vector<FacePair> &matching)
{
    int relexist2 = mesh->build_relations(0, 2);

    // This is the Graph matching verification. It is useful to verify the
    // matching, if done by other algorithms.

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face * f = mesh->getFaceAt(i);
        f->setVisitMark(0);
    }

    for (size_t i = 0; i < matching.size(); i++)
    {
        size_t f1 = matching[i].first;
        size_t f2 = matching[i].second;
        assert(f1 != f2);
        if (f1 < f2)
        {
            Face *f0 = mesh->getFaceAt(f1);
            Face *f1 = mesh->getFaceAt(f2);
            FaceSequence faceneighs = f0->getRelations212();
            if (find(faceneighs.begin(), faceneighs.end(), f1) == faceneighs.end())
            {
                cout << "Error: Face matching error: faces not neighbors  " << endl;
                exit(0);
                return 1;
            }
            assert(!f0->isVisited());
            assert(!f1->isVisited());
            f0->setVisitMark(1);
            f1->setVisitMark(1);
        }
    }

    for (size_t i = 0; i < numfaces; i++) {
        Face * f = mesh->getFaceAt(i);
       assert(f->isVisited());
   }

   if( !relexist2 ) mesh->clear_relations(0,2 );
   
   return 0;

}
///////////////////////////////////////////////////////////////////////////////

Mesh* Tri2Quads::collapse_matched_triangles(Mesh *trimesh, const vector<FacePair> &matching,
                                            int replace)
{
    if(Tri2Quads::verify(trimesh, matching) != 0) return NULL;

    Mesh *quadmesh = new Mesh;

    assert(matching.size());

/*
    vector<size_t> seq;
    seq.resize(2 * matching.size());
    for (size_t i = 0; i < matching.size(); i++)
    {
        size_t f1 = matching[i].first;
        size_t f2 = matching[i].second;
        seq[2 * i] = f1;
        seq[2 * i + 1] = f2;
    }
    sort(seq.begin(), seq.end());

    for (size_t i = 0; i < 2 * matching.size() - 1; i++)
    {
        if (seq[i] == seq[i + 1])
        {
            cout << "Error: Invalid Matching" << endl;
            exit(0);
        }
    }
    //end of checking ...
*/

    Face *tri1, *tri2, *quad;
    FaceSequence faces;

    quadmesh->addNodes(trimesh->getNodes());

    faces.resize(matching.size());
    for (size_t i = 0; i < matching.size(); i++)
    {
        size_t f1 = matching[i].first;
        size_t f2 = matching[i].second;
        tri1 = trimesh->getFaceAt(f1);
        tri2 = trimesh->getFaceAt(f2);
        quad = Face::create_quad(tri1, tri2, replace);
        quadmesh->addFace(quad);
        if (replace) delete tri2;
    }

    if (replace) trimesh->emptyAll();

    return quadmesh;
}
///////////////////////////////////////////////////////////////////////////////

int Tri2Quads::refine_boundary_triangle(Face *tri0)
{
    if (tri0->getSize(0) != 3)
        return 1;

    Vertex *bv0 = NULL;
    Vertex *bv1 = NULL;
    Vertex *bv2 = NULL;

    for (int i = 0; i < 3; i++)
    {
        Vertex *ev1 = tri0->getNodeAt((i + 1) % 3);
        Vertex *ev2 = tri0->getNodeAt((i + 2) % 3);
        if (ev1->isBoundary() && ev2->isBoundary())
        {
            bv0 = ev1;
            bv1 = ev2;
            bv2 = tri0->getNodeAt(i);
            break;
        }
    }

    if (bv0 == NULL || bv1 == NULL)
        return 2;

    Point3D p3d = Vertex::mid_point(bv0, bv1);

    Vertex *bound = Vertex::newObject();
    bound->setXYZCoords(p3d);

    trimesh->addNode(bound);

    NodeSequence tconnect(3);
    tconnect[0] = bv0;
    tconnect[1] = bound;
    tconnect[2] = bv2;
    tri0->setNodes(tconnect);

    tconnect[0] = bound;
    tconnect[1] = bv1;
    tconnect[2] = bv2;
    maxfaceID++;
    Face *tri1 = Face::newObject();
    tri1->setID(maxfaceID);
    tri1->setNodes(tconnect);
    trimesh->addFace(tri1);

    steinerNodes.push_back(bound);

    FacePair facepair;
    facepair.first = tri0->getID();
    facepair.second = tri1->getID();
    facematching.push_back(facepair);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::matchnodes(Vertex *child, Vertex *parent)
{
    child->setDualMate(parent);
    parent->setDualMate(child);

    child->setStatus(MeshEntity::REMOVE);
    parent->setStatus(MeshEntity::REMOVE);
}
///////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::matchnodes(BinaryNode *child, BinaryNode *parent)
{
    if (parent->isMatched() && !child->isMatched())
    {
        cout << "Warning: parent is already matched " << endl;
    }

    if (!child->isMatched() && !parent->isMatched())
        matchnodes(child->getDualNode(), parent->getDualNode());

    btree->unlinkNode(child);
    btree->unlinkNode(parent);
}

///////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::splitParent(BinaryNode *parent, BinaryNode *child1,
                            BinaryNode *child2)
{
    Face *parentface = parent->getDualNode()->getPrimalFace();
    Face *face1 = child1->getDualNode()->getPrimalFace();
    Face *face2 = child2->getDualNode()->getPrimalFace();

    NodeSequence connect(3);

    // Remove all existing vertex-face relations;
    Vertex *vertex;
    for (int i = 0; i < 3; i++)
    {
        vertex = parentface->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face1->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face2->getNodeAt(i);
        vertex->clearRelations(2);
    }

    // Rebuild vertex-face relations...
    for (int i = 0; i < 3; i++)
    {
        vertex = parentface->getNodeAt(i);
        vertex->addRelation2(parentface);

        vertex = face1->getNodeAt(i);
        vertex->addRelation2(face1);

        vertex = face2->getNodeAt(i);
        vertex->addRelation2(face2);
    }

    Vertex *steiner = parentface->getDualNode()->getClone();
    steiner->setID(parentface->getID());
    trimesh->addNode(steiner);
    steinerNodes.push_back(steiner);

    int edge1, edge2, edge3;

    edge1 = edge2 = edge3 = -1;
    FaceSequence neighs;
    for (int i = 0; i < 3; i++)
    {
        Vertex *v0 = parentface->getNodeAt((i + 1) % 3);
        Vertex *v1 = parentface->getNodeAt((i + 2) % 3);
        neighs = Mesh::getRelations112(v0, v1);

        if (neighs.size() == 1)
            edge3 = i;

        if (neighs.size() == 2)
        {
            if (find(neighs.begin(), neighs.end(), face1) != neighs.end())
                edge1 = i;
            if (find(neighs.begin(), neighs.end(), face2) != neighs.end())
                edge2 = i;
        }
    }

    Face *qface;
    Vertex *ev0, *ev1;

    // Match Child1 and One of the Split Triangle ...
    maxfaceID++;
    ev0 = parentface->getNodeAt((edge1 + 1) % 3);
    ev1 = parentface->getNodeAt((edge1 + 2) % 3);
    connect[0] = steiner;
    connect[1] = ev0;
    connect[2] = ev1;
    qface = Face::newObject();
    qface->setNodes(connect);
    qface->setID(maxfaceID);
    Vertex *dc1 = qface->getNewDualNode();
    dc1->setID(maxfaceID);
    trimesh->addFace(qface);
    steinerFaces.push_back(qface);
    matchnodes(face1->getDualNode(), dc1);
    BinaryNode *bnode1 = new BinaryNode(dc1);
    bnode1->setMatchMark(1);
    bnode1->setParent(parent);
    bnode1->addChild(child1);
    btree->addNode(bnode1);

    // Match Child2 and One of the Split Triangle ...
    maxfaceID++;
    ev0 = parentface->getNodeAt((edge2 + 1) % 3);
    ev1 = parentface->getNodeAt((edge2 + 2) % 3);
    connect[0] = steiner;
    connect[1] = ev0;
    connect[2] = ev1;
    qface = Face::newObject();
    qface->setID(maxfaceID);
    qface->setNodes(connect);
    Vertex *dc2 = qface->getNewDualNode();
    dc2->setID(maxfaceID);
    trimesh->addFace(qface);
    steinerFaces.push_back(qface);
    matchnodes(face2->getDualNode(), dc2);
    BinaryNode *bnode2 = new BinaryNode(dc2);
    bnode2->setMatchMark(1);
    bnode2->setParent(parent);
    bnode2->addChild(child2);
    btree->addNode(bnode2);

    // Now Parent have different connectivity ...
    ev0 = parentface->getNodeAt((edge3 + 1) % 3);
    ev1 = parentface->getNodeAt((edge3 + 2) % 3);
    connect[0] = steiner;
    connect[1] = ev0;
    connect[2] = ev1;
    parentface->setNodes(connect);
    Point3D p3d = parentface->getCentroid();
    Vertex *dc3 = parentface->getDualNode();
    dc3->setXYZCoords(p3d);
    parent->addChild(bnode1);
    parent->addChild(bnode2);
    modifiedFaces.push_back(parentface);

    for (int i = 0; i < 3; i++)
    {
        vertex = parentface->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face1->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face2->getNodeAt(i);
        vertex->clearRelations(2);
    }

    child1->setMatchMark(1);
    child2->setMatchMark(1);

    btree->removeNode(child1);
    btree->removeNode(child2);
}

////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::matchnode(BinaryNode* v)
{
    BinaryNode *parv = v->getParent();

    if (parv == NULL)
        return;

    int degree = parv->getDegree();

    if (parv->isRoot() && degree == 2)
    {
        BinaryNode *vsib = v->getSibling();
        splitParent(parv, v, vsib);
        return;
    }

    if (degree == 1 || degree == 2)
    {
        matchnodes(v, parv);
        return;
    }

    if ((degree == 3))
    {

        if (required_topology == ALL_QUADS)
        {
            BinaryNode *vsib = v->getSibling();
            splitParent(parv, v, vsib);
            return;
        }

        BinaryNode *vsib = v->getSibling();
        Vertex *d0 = v->getDualNode();
        Vertex *d1 = vsib->getDualNode();
        if (d0->getNumOfRelations(0) < d1->getNumOfRelations(0))
        {
            matchnodes(v, parv);
            btree->unlinkNode(vsib);
        }
        else
        {
            matchnodes(vsib, parv);
            btree->unlinkNode(v);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

bool Jaal::already_matched(const BinaryNode *node)
{
    return node->isMatched();
}
///////////////////////////////////////////////////////////////////////////////

BinaryNode* Tri2Quads::getChildofDegreeNParent(BNodeList &levelnodes,
                                               int nd)
{
    BinaryNode *currnode, *parent, *child;

    int ncount;
    BNodeList::const_iterator it;

    for (it = levelnodes.begin(); it != levelnodes.end(); ++it)
    {
        currnode = *it;
        parent = currnode->getParent();
        if (parent)
        {
            if (!parent->isMatched())
            {
                ncount = 0;
                if (parent->getParent())
                    ncount = 1;
                for (int i = 0; i < parent->getNumChildren(); i++)
                {
                    child = parent->getChild(i);
                    if (!child->isMatched())
                        ncount++;
                }
                if (ncount == nd)
                    return currnode;
            }
        }
    }

    return NULL;
}

///////////////////////////////////////////////////////////////////////////////

BinaryNode *Tri2Quads::getNextNode(BNodeList &levelnodes)
{
    BinaryNode *currnode = NULL;

    if (levelnodes.empty())
        return currnode;

    BNodeList::iterator it;
    for (it = levelnodes.begin(); it != levelnodes.end(); ++it)
    {
        currnode = *it;
        if (currnode->isMatched())
            btree->unlinkNode(currnode);
    }

    it = remove_if(levelnodes.begin(), levelnodes.end(), already_matched);
    levelnodes.erase(it, levelnodes.end());

    BinaryNode *child = NULL;

    // High Priority: parent having degree = 1;
    child = getChildofDegreeNParent(levelnodes, 1);

    if (!child)
        child = getChildofDegreeNParent(levelnodes, 2);

    // Low Priority: parent having degree = 3;
    if (!child)
        child = getChildofDegreeNParent(levelnodes, 3);

    return child;
}

////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::prunelevel(BNodeList &levelnodes)
{
    while (1)
    {
        BinaryNode *currnode = getNextNode(levelnodes);
        if (currnode == NULL) break;
        matchnode(currnode);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::percolateup()
{
    steinerNodes.clear();
    steinerFaces.clear();

    int height = btree->getHeight();
    BNodeList levelnodes;
    BNodeList::const_iterator it;

    //Reset all the Matching marks to 0;
    for (int i = 0; i < height; i++)
    {
        levelnodes = btree->getLevelNodes(height - i - 1);
        BinaryNode *currnode;
        for (it = levelnodes.begin(); it != levelnodes.end(); ++it)
        {
            currnode = *it;
            currnode->setMatchMark(0);
        }
    }

    // Start Prunning the level. At most the root will be unmatched.
    for (int i = 0; i < height; i++)
    {
        levelnodes = btree->getLevelNodes(height - i - 1);
        prunelevel(levelnodes);
    }

    size_t numfaces = trimesh->getSize(2);
    facematching.reserve(numfaces);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = trimesh->getFaceAt(i);
        Vertex *u = face->getDualNode();
        assert(u);
        const Vertex *v = u->getDualMate();
        if (v)
        {
            if (v->getID() > u->getID())
            {
                FacePair facepair;
                facepair.first = u->getID();
                facepair.second = v->getID();
                facematching.push_back(facepair);
            }
        }
    }

    // If the root is unmatched, bring it down to a leaf and then split the
    // leaf. Do this step after the triangles have been matched.
    BinaryNode *root = btree->getRoot();
    if (!root->isMatched())
    {
        cout << "Warning: Boundary Triangle modified " << endl;
        Face *rootface = root->getDualNode()->getPrimalFace();
        refine_boundary_triangle(rootface);
    }

}

///////////////////////////////////////////////////////////////////////////////

void Tri2Quads::maximum_tree_matching()
{

    // In order to insert any steiner point on the boundary triangle (at the root)
    // We should know which triangles and nodes are on the boundary. Therefore,
    // call this function to set the boundary flags. Building the relationship
    // at this stage is good as even the DualGraph construction require it.

    trimesh->build_relations(0, 2);
    trimesh->search_boundary();

    if (verbose)
        cout << " Creating Dual Graph ... " << endl;

    DualGraph *dgraph = new DualGraph;
    dgraph->build(trimesh);

    if (verbose)
        cout << " Building Binary Tree of Dual Graph ... " << endl;

    btree = new BinaryTree(dgraph);
    btree->build();

    //
    // Store the original tree. As this tree may get changed with insertion
    // of steiner points.
    //
    //btree->saveAs( "tree0");

    if (verbose)
        cout << " Tree Matching ... " << endl;

    percolateup();

    cout << " Percolate up done " << endl;

    // Since percolateup algorithm breaks the link from parent to children
    // relink all the nodes to rebuild the tree.
    //
    // btree->relinkAll();
    // btree->saveAs( "tree1");

    //
    // Percolateup will remove relations, so it is necessary to clear all the
    // relations.
    trimesh->clear_relations(0, 2);

    btree->clear();
    delete btree;
    delete dgraph;
}

///////////////////////////////////////////////////////////////////////////////

Mesh* Tri2Quads::getQuadMesh(Mesh *inmesh, int replace, int topo)
{
    // Preprocessing ....
    if (inmesh->isHomogeneous() != 3)
    {
        cout << "Warning: Input mesh is not triangular " << endl;
        return NULL;
    }

    if (!inmesh->isSimple())
    {
        cout << "Warning: Input mesh is not simple, use edmonds' algorithm " << endl;
        return NULL;
    }

    if (inmesh->getNumOfComponents() > 1)
    {
        cout << "Warning: There are multiple components in the mesh" << endl;
        cout << "         Algorithm works for single component " << endl;
        return NULL;
    }

    trimesh = inmesh;

    required_topology = topo;

    trimesh->enumerate(2);
    maxfaceID = trimesh->getSize(2) - 1;

    int euler0 = trimesh->getEulerCharacteristic();
    cout << " Input Euler # : " << euler0 << endl;

    ///////////////////////////////////////////////////////////////////////////
    // Generate Maximum Matching on a binary tree using Suneeta's Algorithm.
    // If the required topology is set to ALL_QUADS, steiner points( and new
    // faces) will be inserted in the input triangle mesh.  Please note that
    // this implementation doesn't produces "The" optimal soluation as
    // described in the original papers, and doesn't even guarantee that
    // the resulting quadrilaterals will be convex. This along with other
    // topological and geometric optimization are anyhow essential and
    // are carried out during the "Post Processing" step. Therefore, we
    // have sacrifised performance over quality in this implementation.
    // Roughly we can expect that about 4-5% steiner points are inserted in
    // most of the general cases.
    ///////////////////////////////////////////////////////////////////////////
    ticks start_tick = getticks();

    maximum_tree_matching();


    Mesh *quadmesh = collapse_matched_triangles(trimesh, facematching, replace);

    ticks end_tick = getticks();

    double elapsed_ticks = elapsed(end_tick, start_tick);

    cout << "Info: Tri->Quad Elapsed Performance Counter " << elapsed_ticks << endl;

    if (!quadmesh->isSimple())
        cout << "Warning: Quadrilateral Mesh is not simple " << endl;

    quadmesh->enumerate(2);
    if (!quadmesh->is_consistently_oriented())
    {
        cout << "Warning:Trying to make Quadrilateal Mesh consistently oriented " << endl;
        quadmesh->make_consistently_oriented();
        if (quadmesh->is_consistently_oriented())
            cout << "Info: Quadrilateral Mesh is now consistently oriented: Very good " << endl;
        else
            cout << "Alas ! Quadrilateral Mesh is still inconsistently oriented: Check manually " << endl;
    }

    quadmesh->enumerate(0);
    quadmesh->enumerate(2);

    //////////////////////////////////////////////////////////////////////////
    // Since Steiner points may be inserted in the mesh ( and new triangles).
    // Renumber all the nodes and faces for future processing.
    //////////////////////////////////////////////////////////////////////////

    trimesh->enumerate(0);
    trimesh->enumerate(2);

    return quadmesh;
}

///////////////////////////////////////////////////////////////////////////////

void Tri2Quads::match_tree_walk(BinaryTree *btree, BinaryNode *parent)
{
    //
    // Brings all the internal unmatched nodes at the leaf.
    //
    if (parent == NULL)
        return;
    if (parent->isLeaf())
        return;

    int numChildren = parent->getNumChildren();

    for (int i = 0; i < numChildren; i++)
    {
        BinaryNode *child1 = parent->getChild(i);
        if (!btree->isMatched(parent, child1))
        {
            int numGrandChildren = child1->getNumChildren();
            for (int j = 0; j < numGrandChildren; j++)
            {
                BinaryNode *child2 = child1->getChild(j);
                if (btree->isMatched(child1, child2))
                {
                    Vertex *np = parent->getDualNode();
                    assert(np);
                    Vertex *c1 = child1->getDualNode();
                    assert(c1);
                    Vertex *c2 = child2->getDualNode();
                    assert(c2);
                    matchnodes(np, c1);
                    c2->setDualMate(NULL);
                    c2->setStatus(MeshEntity::ACTIVE);
                    match_tree_walk(btree, child2);
                    return;
                }
            }
        }
    }

}

///////////////////////////////////////////////////////////////////////////////
