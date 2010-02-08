#include "mesh.h"

////////////////////////////////////////////////////////////////////

vector<int> getVertexFaceDegrees(Mesh *mesh) {
    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    int numnodes = mesh->getSize(0);

    vector<int> degree(numnodes);
    vector<Face*> neighs;
    for (int i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNode(i);
        neighs = v->getRelations2();
        degree[i] = neighs.size();
    }

    int mindegree = *min_element(degree.begin(), degree.end());
    int maxdegree = *max_element(degree.begin(), degree.end());

    cout << " Min Vertex-Face Degree : " << mindegree << endl;
    cout << " Max Vertex-Face Degree : " << maxdegree << endl;

    for (int i = mindegree; i <= maxdegree; i++) {
        int ncount = 0;
        for (size_t j = 0; j < degree.size(); j++)
            if (degree[j] == i) ncount++;
        cout << "Degree : " << i << " Count " << ncount << endl;
    }

    if (!relexist) mesh->clear_relations(0, 2);

    return degree;
}

////////////////////////////////////////////////////////////////////

Vertex* get_VertexOf_FaceDegree(Mesh *mesh, int n) {
    size_t numnodes = mesh->getSize(0);
    vector<Face*> neighs;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNode(i);
        neighs = v->getRelations2();
        if (neighs.size() == size_t(n)) return v;
    }
    return NULL;
}

////////////////////////////////////////////////////////////////////

vector<Face*> search_diamonds(Mesh *mesh) {
    int numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    vector<Face*> diamonds;
    for (int iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFace(iface);
        assert(face);
        if (face->getSize(0) == 4) {
            Vertex *v0 = face->getConnection(0);
            Vertex *v1 = face->getConnection(1);
            Vertex *v2 = face->getConnection(2);
            Vertex *v3 = face->getConnection(3);

            int d0 = v0->getRelations2().size();
            int d1 = v1->getRelations2().size();
            int d2 = v2->getRelations2().size();
            int d3 = v3->getRelations2().size();
            if ((d0 == 3 && d2 == 3) || (d1 == 3 && d3 == 3)) {
                diamonds.push_back(face);
            }
        }
    }

    if (!relexist) mesh->clear_relations(0, 2);

    cout << "Number of Diamonds " << diamonds.size() << endl;
    return diamonds;
}

///////////////////////////////////////////////////////////////////////////////

vector<Vertex*> search_doublets(Mesh *mesh) {
    int numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    vector<Vertex*> doublets;
    for (int i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNode(i);
        if (v->getRelations2().size() == 2)
            doublets.push_back(v);
    }

    if (!relexist) mesh->clear_relations(0, 2);

    cout << "Number of doublets : " << doublets.size() << endl;
    return doublets;
}


////////////////////////////////////////////////////////////////////

int face_close(Mesh *mesh, Face *face, Vertex *v0, Vertex *v2) {
    assert(face->getSize(0) == 4);
    assert(mesh->getAdjTable(0, 2));

    if (v0->isBoundary() || v2->isBoundary()) return 1;


    int vpos = -1;
    for (int i = 0; i < 4; i++) {
        if (face->getConnection(i) == v0) {
            vpos = i;
            break;
        }
    }

    assert(vpos >= 0);
    if (face->getConnection((vpos + 2) % 4) != v2) {
        cout << "Warning: Face-open requires opposite vertices " << endl;
        return 1;
    }

    Vertex *v1 = face->getConnection((vpos + 1) % 4);
    Vertex *v3 = face->getConnection((vpos + 3) % 4);

    vector<Face*> vr0 = v0->getRelations2();
    vector<Face*> vr1 = v1->getRelations2();
    vector<Face*> vr2 = v2->getRelations2();
    vector<Face*> vr3 = v3->getRelations2();
    cout << " Face close attempt " << endl;

    if (vr0.size() != 3) return 2;
    if (vr2.size() != 3) return 2;

    cout << " Face Close Up " << endl;

    Point3D p3d = Vertex::mid_point(v0, v2);

    Vertex *newvtx = Vertex::newObject();
    //   newvtx->setID(gid);
    newvtx->setXYZCoords(p3d);
    mesh->addNode(newvtx);

    v0->removeRelation2(face);
    v1->removeRelation2(face);
    v2->removeRelation2(face);
    v3->removeRelation2(face);

    for (size_t i = 0; i < vr0.size(); i++) {
        vr0[i]->replaceNode(v0, newvtx);
        if (vr0[i] != face)
            newvtx->addRelation2(vr0[i]);
    }

    for (size_t i = 0; i < vr2.size(); i++) {
        vr2[i]->replaceNode(v2, newvtx);
        if (vr0[i] != face)
            newvtx->addRelation2(vr2[i]);
    }

    v0->setRemoveMark(1);
    v2->setRemoveMark(1);
    face->setRemoveMark(1);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

Vertex* insert_doublet(Mesh *mesh, Face *face, Vertex *v0, Vertex *v2) {
    Point3D p3d = Vertex::mid_point(v0, v2);

    Vertex *doublet = Vertex::newObject();
    doublet->setXYZCoords(p3d);

    Vertex *o1 = NULL, *o2 = NULL;

    vector<Vertex*> connect = face->getConnection();
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
    newquad1->setConnection(connect);
    mesh->addFace(newquad1);

    face->setRemoveMark(1);

    connect[0] = doublet;
    connect[1] = v2;
    connect[2] = o2;
    connect[3] = v0;

    Face *newquad2 = new Face;
    newquad2->setConnection(connect);
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

Vertex* insert_doublet(Mesh *mesh, Face *face) {
    if (face->getSize(0) != 4) return NULL;
    Vertex *v0 = face->getConnection(0);
    Vertex *v2 = face->getConnection(2);
    return insert_doublet(mesh, face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

Vertex* insert_boundary_doublet(Mesh *mesh, Face *face) {
    if (!face->isBoundary()) return NULL;

    if (face->getSize(0) != 4) return NULL;

    int ncount = 0;
    for (int i = 0; i < 4; i++) {
        Vertex *v = face->getConnection(i);
        ncount += v->isBoundary();
    }

    if (ncount != 3) return NULL;

    Vertex *v0 = NULL, *v2 = NULL;
    for (int i = 0; i < 4; i++) {
        Vertex *v = face->getConnection(i);
        if (!v->isBoundary()) {
            v0 = face->getConnection((i + 0) % 4);
            v2 = face->getConnection((i + 2) % 4);
        }
    }

    return insert_doublet(mesh, face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

int diamond_collapse(Mesh *mesh, Face *face) {
    if (face->getSize(0) != 4) return 1;

    assert(mesh->getAdjTable(0, 2));

    Vertex *v0 = face->getConnection(0);
    Vertex *v1 = face->getConnection(1);
    Vertex *v2 = face->getConnection(2);
    Vertex *v3 = face->getConnection(3);

    vector<Face*> vr0 = v0->getRelations2();
    for (size_t i = 0; i < vr0.size(); i++)
        if (vr0[i]->isRemoved()) return 1;

    vector<Face*> vr1 = v1->getRelations2();
    for (size_t i = 0; i < vr1.size(); i++)
        if (vr1[i]->isRemoved()) return 1;

    vector<Face*> vr2 = v2->getRelations2();
    for (size_t i = 0; i < vr2.size(); i++)
        if (vr2[i]->isRemoved()) return 1;

    vector<Face*> vr3 = v3->getRelations2();
    for (size_t i = 0; i < vr3.size(); i++)
        if (vr3[i]->isRemoved()) return 1;

    int d0 = vr0.size();
    int d1 = vr1.size();
    int d2 = vr2.size();
    int d3 = vr3.size();

    if ((d0 == 3 && d2 == 3)) {
        return face_close(mesh, face, v0, v2);
    }

    if ((d1 == 3 && d3 == 3)) {
        return face_close(mesh, face, v1, v3);
    }

    return 1;
}

////////////////////////////////////////////////////////////////////

int remove_unconstrained_doublet(Mesh *mesh, Vertex *vertex) {
    if (vertex->isBoundary()) return 1;

    vector<Face*> neighs = vertex->getRelations2();

    if (neighs.size() != 2) return 1;

    Vertex *d1 = NULL, *d2 = NULL, *o1 = NULL, *o2 = NULL;

    vector<Vertex*> connect = neighs[0]->getConnection();
    if (connect.size() != 4) return 1;

    for (size_t i = 0; i < connect.size(); i++) {
        if (connect[i] == vertex) {
            d1 = connect[(i + 1) % 4];
            o1 = connect[(i + 2) % 4];
            d2 = connect[(i + 3) % 4];
            break;
        }
    }

    connect = neighs[1]->getConnection();
    if (connect.size() != 4) return 1;

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
    newquad->setConnection(connect);
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
////////////////////////////////////////////////////////////////////

void cleanup_diamonds(Mesh *mesh) {
    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    vector<Face*> faces = search_diamonds(mesh);

    int ncount = 0;
    for (size_t i = 0; i < faces.size(); i++) {
        int err = diamond_collapse(mesh, faces[i]);
        if (!err) ncount++;
    }

    laplacian_smoothing(mesh, 5);

    cout << "#Diamonds removed from the mesh : " << ncount << endl;
    if (ncount) cleanup_diamonds(mesh);

    if (!relexist) mesh->clear_relations(0, 2);

    mesh->prune();
    mesh->enumerate(0);
    mesh->enumerate(2);
}
////////////////////////////////////////////////////////////////////

void cleanup_doublets(Mesh *mesh) {
    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    vector<Vertex*> doublets = search_doublets(mesh);

    int ncount = 0;
    for (size_t i = 0; i < doublets.size(); i++) {
        int err = remove_unconstrained_doublet(mesh, doublets[i]);
        if (!err) ncount++;
    }

    laplacian_smoothing(mesh, 5);

    cout << "#Doublets removed from the mesh : " << ncount << endl;
    //if( ncount ) cleanup_doublets(mesh);

    if (!relexist) mesh->clear_relations(0, 2);

    mesh->prune();
    mesh->enumerate(0);
    mesh->enumerate(2);
}
////////////////////////////////////////////////////////////////////

void cleanup_boundary(Mesh *mesh) {
    mesh->search_boundary();

    int relexist = mesh->build_relations(0, 0);

    vector<Vertex*> degree2nodes;

    int numnodes = mesh->getSize(0);

    Vertex* node, *onode;
    vector<Vertex*> neighs;
    for (int i = 0; i < numnodes; i++) {
        node = mesh->getNode(i);
        if (node->isBoundary()) {
            neighs = node->getRelations0();
            if (neighs.size() == 2) {
                if (neighs[0]->isBoundary() && neighs[1]->isBoundary())
                    degree2nodes.push_back(node);
            }
        }
    }

    if (!relexist) mesh->clear_relations(0, 0);

    relexist = mesh->build_relations(0, 2);

    Face *boundface;
    vector<Face*> faceneighs;
    for (size_t i = 0; i < degree2nodes.size(); i++) {
        node = degree2nodes[i];
        faceneighs = node->getRelations2();
        assert(faceneighs.size() == 1);
        boundface = faceneighs[0];
        if (boundface->getSize(0) == 4) {
            int j = boundface->queryNodeAt(node);
            onode = boundface->getConnection((j + 2) % 4);
            if (!onode->isBoundary())
                insert_doublet(mesh, boundface, node, onode);
        }
    }
    mesh->prune();
    mesh->enumerate(0);
    mesh->enumerate(2);

    mesh->search_boundary();

    numnodes = mesh->getSize(0);
    for (int i = 0; i < numnodes; i++) {
        Vertex *boundnode = mesh->getNode(i);
        if (boundnode->isBoundary()) {
            faceneighs = boundnode->getRelations2();
            if (faceneighs.size() == 3) {

                bool b0 = faceneighs[0]->isBoundary();
                bool b1 = faceneighs[1]->isBoundary();
                bool b2 = faceneighs[2]->isBoundary();

                Face *internalface = 0;
                if (b0 + b1 + b2 == 2) {
                    if (!b0) internalface = faceneighs[0];
                    if (!b1) internalface = faceneighs[1];
                    if (!b2) internalface = faceneighs[2];
                }

                if (internalface) {
                    if (internalface->getSize(0) == 4) {
                        int j = internalface->queryNodeAt(boundnode);
                        Vertex *v0 = internalface->getConnection((j + 1) % 4);
                        Vertex *v2 = internalface->getConnection((j + 3) % 4);
                        face_close(mesh, internalface, v0, v2);
                    }
                }
            }
        }
    }

    if (!relexist) mesh->clear_relations(0, 2);

    laplacian_smoothing(mesh, 5);

    mesh->prune();
    mesh->enumerate(0);
    mesh->enumerate(2);

    cout << mesh->getSize(0) << endl;
    cout << mesh->getSize(2) << endl;
}
