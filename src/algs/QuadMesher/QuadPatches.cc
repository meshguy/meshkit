#include <meshkit/Mesh.h>

using namespace Jaal;

//
// QTrack is similar to "Motorcycle graph" proposed by Eppestein group.
// But somehow, I don't like this term.
//

#ifdef CSV
struct QTrack {
    const static int END_AT_TERMINALS = 0;
    const static int END_AT_CROSSINGS = 1;

    vector<Vertex*> sequence;

    bool operator ==(const QTrack & rhs) const {
        size_t nSize = sequence.size();
        if (nSize != rhs.sequence.size()) return 0;

        Vertex *v0src = sequence.front();
        Vertex *v0dst = sequence.back();
        Vertex *v1src = rhs.sequence.front();
        Vertex *v1dst = rhs.sequence.back();

        if (v0src == v1src && v0dst == v1dst) return 1;
        if (v0src == v1dst && v0dst == v1src) return 1;

        return 0;
    }

    bool operator<(const QTrack & rhs) const {
        return sequence.size() < rhs.sequence.size();
    }

    /////////////////////////////////////////////////////////////////////////////
    // There are two ways to advance along the track.
    // (1) Incremental:  All track propagate simulateneously and
    //                   at the intersection, only one is allowed
    //                   to proceed towards other irregular node.
    // (11) Greedy   :   Start from one vertex and complete its
    //                   track.
    // It is not clear which method is better, but "greedy" may likely
    // give very high aspect ratio quad patches. Incremental on the
    // other hand may produce many small patches..
    //
    /////////////////////////////////////////////////////////////////////////////

    int advance_single_step(int endat) {
        //////////////////////////////////////////////////////////////////////////
        //                   **********************
        //                   *         *          *
        //                   *         * Next     *
        //                   *         *          *
        //           Avoid   **********************  Avoid
        //                   *         *          *
        //                   *         * Current  *
        //                   *         *          *
        //                   *         *          *
        //                   **********************
        //                            Source
        // A Source vertex and Current edge is chosen.
        // We want to avoid two edges and want to select "Next" edge.
        //////////////////////////////////////////////////////////////////////////
        Vertex *v0, *v1, *v2, *v3, *v4;

        vector<Face*> adjFaces;
        vector<Vertex*> vneighs;
        set<Vertex*> vset;

        size_t index = sequence.size();
        v0 = sequence[index - 2];
        v1 = sequence[index - 1];
        v0->setVisitMark(1);

        if (endat == END_AT_CROSSINGS && v1->isVisited()) return 0;
        if (v1->isBoundary()) return 0;

        v1->setVisitMark(1);
        vneighs = v1->getRelations0();
        if (vneighs.size() != 4) return 0;

        adjFaces = Mesh::getRelations112(v0, v1);
        assert(adjFaces.size() == 2);
        v2 = Face::opposite_node(adjFaces[0], v0);
        v3 = Face::opposite_node(adjFaces[1], v0);

        vset.clear();
        vset.insert(vneighs[0]);
        vset.insert(vneighs[1]);
        vset.insert(vneighs[2]);
        vset.insert(vneighs[3]);
        vset.erase(v0);
        vset.erase(v2);
        vset.erase(v3);
        assert(vset.size() == 1);
        v4 = *vset.begin();
        sequence.push_back(v4);
        return 1;
    }

    void advance(int endat) {
        assert(sequence.size() == 2);

        // Starting node is always irregular ...
        vector<Face*> vfaces = sequence[0]->getRelations2();
        assert(vfaces.size() != 4);

        while (1) {
            int progress = advance_single_step(endat);
            if (!progress) break;
        }

/*
        // The path is reversible and therefore, we will give a direction
        // from the lower source node to higher destination node.
        if (sequence.front() > sequence.back())
            reverse(sequence.begin(), sequence.end());
        assert(sequence.front() < sequence.back());
*/
        // Checking the correctness..
        if (endat == END_AT_TERMINALS) {
            vfaces = sequence.front()->getRelations2();
            assert(vfaces.size() != 4);
            vfaces = sequence.back()->getRelations2();
            assert(vfaces.size() != 4);
            for (int i = 1; i < sequence.size() - 1; i++) {
                vfaces = sequence[i]->getRelations2();
                assert(vfaces.size() == 4);
            }

        }
    }

};

///////////////////////////////////////////////////////////////////////////////

struct StructuredMesh2D {

    StructuredMesh2D() {
        nx = ny = 0;
    }
    int nx, ny;

    vector<Face*> faces;
    vector<Vertex*> nodes;
    set<Vertex*> cornerNodes;

    void clearAll() {
        nodes.clear();
        faces.clear();
        neighs.clear();
        cornerNodes.clear();
        nx = 0;
        ny = 0;
    }

    bool operator<(const StructuredMesh2D & rhs) const {
        return this->getSize(2) < rhs.getSize(2);
    }

    size_t getSize(int e) const {
        if (e == 0) return nodes.size();
        if (e == 2) return faces.size();
        return 0;
    }

    int myID;
    vector<int> neighs;
};
#endif

///////////////////////////////////////////////////////////////////////////////

void build_submesh_topology(StructuredMesh2D &smesh) {
    smesh.neighs.clear();

    FaceSequence adjfaces;
    set<int> nset;
    for (size_t i = 0; i < smesh.faces.size(); i++) {
        Face *face = smesh.faces[i];
        for (int j = 0; j < 4; j++) {
            Vertex *v0 = face->getNodeAt((j + 0) % 4);
            Vertex *v1 = face->getNodeAt((j + 1) % 4);
            adjfaces = Mesh::getRelations112(v0, v1);
            for (int k = 0; k < adjfaces.size(); k++) {
                int nid = adjfaces[k]->getTag();
                nset.insert(nid);
            }
        }
    }

    nset.erase(smesh.myID);
    set<int>::const_iterator it1;
    for (it1 = nset.begin(); it1 != nset.end(); ++it1)
        smesh.neighs.push_back(*it1);
}

///////////////////////////////////////////////////////////////////////////////

Face *getSeedFace(Jaal::Mesh *mesh) {
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        Face *f = mesh->getFaceAt(i);
        if (!f->isVisited()) return f;
    }
    return NULL;
}

///////////////////////////////////////////////////////////////////////////////

size_t independent_components(Jaal::Mesh *mesh) {

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        Face *f = mesh->getFaceAt(i);
        f->setVisitMark(0);
        f->setTag(0);
    }

    size_t compid = 0;
    deque<Face*> faceQ;
    FaceSequence adjfaces;
    while (1) {
        Face *currface = getSeedFace(mesh);
        if (currface == NULL) break;
        faceQ.push_back(currface);
        while (!faceQ.empty()) {
            currface = faceQ.front();
            faceQ.pop_front();
            if (!currface->isVisited()) {
                currface->setTag(compid);
                currface->setVisitMark(1);
                for (int i = 0; i < 4; i++) {
                    Vertex *v0 = currface->getNodeAt((i + 0) % 4);
                    Vertex *v1 = currface->getNodeAt((i + 1) % 4);
                    if (v0->isVisited() && v1->isVisited()) continue;
                    adjfaces = Mesh::getRelations112(v0, v1);
                    for (int j = 0; j < adjfaces.size(); j++)
                        faceQ.push_back(adjfaces[j]);
                }
            }
        }
        compid++;
    }

    for (size_t i = 0; i < numfaces; i++) {
        Face *f = mesh->getFaceAt(i);
        assert(f->isVisited());
    }
    return compid;
}

/////////////////////////////////////////////////////////////////////////////////

int merge_submesh(StructuredMesh2D &amesh, StructuredMesh2D &bmesh) {
    if (amesh.myID == bmesh.myID) return 1;

    if (amesh.cornerNodes.size() != 4) return 2;
    if (bmesh.cornerNodes.size() != 4) return 2;

    vector<Vertex*> common;
    set_intersection(amesh.cornerNodes.begin(), amesh.cornerNodes.end(),
            bmesh.cornerNodes.begin(), bmesh.cornerNodes.end(),
            back_inserter(common));

    if (common.size() != 2) return 3;

    FaceSequence vfaces;
    for (int i = 0; i < amesh.faces.size(); i++) {
        for (int j = 0; j < 4; j++) {
            Vertex *v = amesh.faces[i]->getNodeAt(j);
            if (!v->isBoundary()) {
                vfaces = v->getRelations2();
                if (vfaces.size() != 4) return 4;
            }
        }
    }

    for (int i = 0; i < bmesh.faces.size(); i++) {
        for (int j = 0; j < 4; j++) {
            Vertex *v = bmesh.faces[i]->getNodeAt(j);
            if (!v->isBoundary()) {
                vfaces = v->getRelations2();
                if (vfaces.size() != 4) return 4;
            }
        }
    }

    for (int i = 0; i < bmesh.faces.size(); i++) {
        amesh.faces.push_back(bmesh.faces[i]);
        bmesh.faces[i]->setTag(amesh.myID);
    }

    set<Vertex*>::const_iterator it;
    for (it = bmesh.cornerNodes.begin(); it != bmesh.cornerNodes.end(); ++it)
        amesh.cornerNodes.insert(*it);

    assert(amesh.cornerNodes.size() == 6);

    amesh.cornerNodes.erase(common[0]);
    amesh.cornerNodes.erase(common[1]);

    bmesh.clearAll();

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////

StructuredMesh2D getQuadPatch(Jaal::Mesh *mesh, int compid) {
    StructuredMesh2D smesh;
    smesh.myID = compid;
    smesh.nx = 1;
    smesh.ny = 1;

    size_t numfaces = mesh->getSize(2);

    set<Face*> faceSet;
    set<Face*>::const_iterator fit;
    set<Vertex*> nodeSet;
    set<Vertex*>::const_iterator nit;

    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        if (face->getTag() == compid) {
            faceSet.insert(face);
            nodeSet.insert(face->getNodeAt(0));
            nodeSet.insert(face->getNodeAt(1));
            nodeSet.insert(face->getNodeAt(2));
            nodeSet.insert(face->getNodeAt(3));
        }
    }
    if (faceSet.empty()) return smesh;

    NodeSet boundNodes, cornerNodes;
    FaceSet boundFaces, cornerFaces;

    FaceSequence vfaces;
    for (nit = nodeSet.begin(); nit != nodeSet.end(); ++nit) {
        Vertex *vertex = *nit;
        vfaces = vertex->getRelations2();
        int ncount = 0;
        for (int j = 0; j < vfaces.size(); j++)
            if (faceSet.find(vfaces[j]) != faceSet.end()) ncount++;
        if (ncount == 1) {
            cornerNodes.insert(vertex);
            for (int j = 0; j < vfaces.size(); j++)
                if (faceSet.find(vfaces[j]) != faceSet.end()) cornerFaces.insert(vfaces[j]);
        }
        if (ncount != 4) boundNodes.insert(vertex);
    }

    for (fit = faceSet.begin(); fit != faceSet.end(); ++fit)
        smesh.faces.push_back(*fit);

    for (nit = nodeSet.begin(); nit != nodeSet.end(); ++nit)
        smesh.nodes.push_back(*nit);

    smesh.cornerNodes = cornerNodes;

    build_submesh_topology(smesh);

    return smesh;

}
/////////////////////////////////////////////////////////////////////////////////

int Mesh::search_quad_patches() 
{
    int nTopo = isHomogeneous();
    if (nTopo != 4) {
        cout << "Error: The mesh must be all Quads " << endl;
        return 1;
    }

    int relexist2 = build_relations(0, 2);
    int relexist0 = build_relations(0, 0);

    search_boundary();

    size_t numnodes = getSize(0);
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = getNodeAt(i);
        vertex->setVisitMark(0);
    }

    size_t numfaces = getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = getFaceAt(i);
        face->setVisitMark(0);
    }

    vector<QTrack> qpath;
    QTrack qp;
    qp.sequence.resize(2); // As we know the starting edge

    size_t nCount = 0;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = getNodeAt(i);
        if (!vertex->isBoundary()) {
            NodeSequence vnodes = vertex->getRelations0();
            if (vnodes.size() != 4)  {
                qp.sequence[0] = vertex;
                for (int j = 0; j < vnodes.size(); j++) {
                     qp.sequence[1] = vnodes[j];
                     nCount++;
                     qpath.push_back(qp);
                }
            }
        }
    }

    if (qpath.empty()) {
        cout << "Info: There are no irregular nodes in the mesh" << endl;
        return 0;
    } else
        cout << "# of branches spawned " << nCount << endl;

    for (size_t j = 0; j < qpath.size(); j++) 
        qpath[j].advance(0);

    sort(qpath.begin(), qpath.end());
    saveAs("b.dat");

    for (int i = 0; i < qpath.size(); i++) {
        if (qpath[i].sequence.front()->isBoundary()) continue;
        if (qpath[i].sequence.back()->isBoundary() ) continue;
        int found = 0;
        cout << " PATH SIZE " << qpath[i].sequence.size() << endl;
        for (int k = 0; k < qpath[i].sequence.size(); k++)
               cout << qpath[i].sequence[k]->getID() << " ";
        cout << endl;
        for (int j = 0; j < qpath.size(); j++) {
            if ((i != j) && (qpath[i] == qpath[j])) {
                found = 1;
                cout << "Same path " << i << " " << j << endl;
                for (int k = 0; k < qpath[j].sequence.size(); k++)
                    cout << qpath[j].sequence[k]->getID() << " ";
                cout << endl;
            }
        }
        assert(found);
    }

    exit(0);

    int numPatches = independent_components(this);

    deque<StructuredMesh2D> submesh(numPatches);

    // New we can merge some sub-meshes
    for (int i = 0; i < numPatches; i++) {
        StructuredMesh2D sm = getQuadPatch(this, i);
        submesh[i] = sm;
    }
    sort(submesh.begin(), submesh.end());

    while (1) {
        size_t nSize = submesh.size();
        cout << "#Submeshes " << submesh.size() << endl;
        size_t count_merged = 0;
        for (int i = 0; i < nSize; i++) {
            for (int j = i + 1; j < nSize; j++) {
                int err = merge_submesh(submesh[i], submesh[j]);
                if (!err) count_merged++;
            }
        }
        cout << " #of Submesh merged " << count_merged << endl;
        if (count_merged == 0) break;

        //
        // Retain only those submeshes which are not empty. ( Note
        // that when two submeshes merge, one of them become empty.
        //
        for (int i = 0; i < nSize; i++) {
            StructuredMesh2D sm = submesh.front();
            submesh.pop_front();
            if (sm.getSize(2)) submesh.push_back(sm);
        }
    }

    sort(submesh.begin(), submesh.end());
    for (int i = 0; i < submesh.size(); i++)
        cout << submesh[i].getSize(2) << " ";

    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = getNodeAt(i);
        vertex->setTag(vertex->isVisited());
    }

    if (relexist0)
        clear_relations(0, 0);

    if (relexist2)
        clear_relations(0, 2);
    return submesh.size();
}
