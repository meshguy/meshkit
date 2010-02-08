#include "mesh.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <set>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int Vertex::global_id = 0;

AttribKey  MeshEntity::maxAttribID = 0;
std::map<string, AttribKey> MeshEntity::attribKeyMap;

///////////////////////////////////////////////////////////////////////////////

AttribKey MeshEntity::getAttribKey(const string &s) {
      map<string,AttribKey>::const_iterator it;

      it =  attribKeyMap.find(s);

      if( it == attribKeyMap.end() ) return 0;

      return it->second;
}
///////////////////////////////////////////////////////////////////////////////

bool MeshEntity::hasAttribute(const string &s) {
     return getAttribKey(s);
}

///////////////////////////////////////////////////////////////////////////////

AttribKey MeshEntity::addAttribute(const string &s) {
     AttribKey key = getAttribKey(s);
     if( key == 0) attribKeyMap[s] = ++maxAttribID;
     return attribKeyMap[s];
}

void MeshEntity::removeAttribute(const string &s) {
     attribKeyMap.erase(s);
}
///////////////////////////////////////////////////////////////////////////////

Vertex* Vertex::newObject() {
    Vertex *v = new Vertex;
    v->setID(global_id);
    global_id++;
    return v;
}

///////////////////////////////////////////////////////////////////////////////

Point3D Vertex::mid_point(const Vertex *v0, const Vertex *v1) {

    Point3D p0 = v0->getXYZCoords();
    Point3D p1 = v1->getXYZCoords();

    Point3D pmid;
    pmid[0] = 0.5 * (p0[0] + p1[0]);
    pmid[1] = 0.5 * (p0[1] + p1[1]);
    pmid[2] = 0.5 * (p0[2] + p1[2]);

    return pmid;
}
///////////////////////////////////////////////////////////////////////////////

Point3D Face::getCentroid() const {
    Point3D pc;

    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;

    Point3D p3d;
    for (size_t inode = 0; inode < connect.size(); inode++) {
        Vertex *v = connect[inode];

        p3d = v->getXYZCoords();
        pc[0] += p3d[0];
        pc[1] += p3d[1];
        pc[2] += p3d[2];
    }

    pc[0] /= 3.0;
    pc[1] /= 3.0;
    pc[2] /= 3.0;

    return pc;

}

///////////////////////////////////////////////////////////////////////////////

vector<FaceType> Mesh::getRelation112(NodeType vtx0, NodeType vtx1) {
    vector<FaceType> v0faces = vtx0->getRelations2();
    vector<FaceType> v1faces = vtx1->getRelations2();

    if (v0faces.empty() || v1faces.empty()) {
        cout << "Warning: Vertex-Faces relations are empty " << endl;
    }

    vector<FaceType> faceneighs;
    set_intersection(v0faces.begin(), v0faces.end(),
            v1faces.begin(), v1faces.end(),
            back_inserter(faceneighs));

    return faceneighs;
}
///////////////////////////////////////////////////////////////////////////////

size_t Mesh::count_edges() {
    int relexist = build_relations(0, 0);

    size_t numnodes = getSize(0);

    vector<Vertex*> neighs;
    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++) {
        neighs = nodes[i]->getRelations0();
        for (size_t j = 0; j < neighs.size(); j++)
            if (nodes[i] > neighs[j]) ncount++;
    }
    if (!relexist) clear_relations(0, 0);

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::prune() {
    vector<NodeType> livenodes;

    livenodes.reserve(nodes.size());
    for (size_t i = 0; i < nodes.size(); i++) {
        if (!nodes[i]->isRemoved())
            livenodes.push_back(nodes[i]);
    }
    nodes = livenodes;

    vector<FaceType> livefaces;
    livefaces.reserve(faces.size());
    for (size_t i = 0; i < faces.size(); i++) {
        if (!faces[i]->isRemoved())
            livefaces.push_back(faces[i]);
    }
    faces = livefaces;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::enumerate(int etype) {
    int index = 0;

    Vertex *vertex;
    if (etype == 0) {
        index = 0;
        BOOST_FOREACH(vertex, nodes)
        vertex->setID(index++);
    }

    Face *face;
    if (etype == 2) {
        index = 0;
        BOOST_FOREACH(face, faces)
        face->setID(index++);
    }
}

///////////////////////////////////////////////////////////////////////////////

size_t Mesh::getBoundarySize(int d) const {
    int ncount = 0;

    if (d == 0) {
        for (size_t i = 0; i < nodes.size(); i++)
            if (!nodes[i]->isRemoved() && nodes[i]->isBoundary()) ncount++;
    }

    if (d == 2) {
        for (size_t i = 0; i < faces.size(); i++)
            if (!faces[i]->isRemoved() && faces[i]->isBoundary()) ncount++;
    }
    return ncount;
}



///////////////////////////////////////////////////////////////////////////////

NodeType Face::opposite_node(const FaceType tri, NodeType n1, NodeType n2) {
    NodeType tn0 = tri->getConnection(0);
    NodeType tn1 = tri->getConnection(1);
    NodeType tn2 = tri->getConnection(2);

    if (tn0 == n1 && tn1 == n2) return tn2;
    if (tn0 == n2 && tn1 == n1) return tn2;

    if (tn1 == n1 && tn2 == n2) return tn0;
    if (tn1 == n2 && tn2 == n1) return tn0;

    if (tn2 == n1 && tn0 == n2) return tn1;
    if (tn2 == n2 && tn0 == n1) return tn1;
    cout << " Warning: You should not come here " << endl;
    cout << " Face " << tn0 << " " << tn1 << " " << tn2 << endl;
    cout << " search for " << n1 << "  " << n2 << endl;
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

void Face::opposite_nodes(const FaceType quad, NodeType n1, NodeType n2,
        NodeType &n3, NodeType &n4) {
    NodeType qn0 = quad->getConnection(0);
    NodeType qn1 = quad->getConnection(1);
    NodeType qn2 = quad->getConnection(2);
    NodeType qn3 = quad->getConnection(3);

    if ((qn0 == n1 && qn1 == n2) || (qn0 == n2 && qn1 == n1)) {
        n3 = qn2;
        n4 = qn3;
        return;
    }

    if ((qn1 == n1 && qn2 == n2) || (qn1 == n2 && qn2 == n1)) {
        n3 = qn0;
        n4 = qn3;
        return;
    }

    if ((qn2 == n1 && qn3 == n2) || (qn2 == n2 && qn3 == n1)) {
        n3 = qn0;
        n4 = qn1;
        return;
    }

    if ((qn3 == n1 && qn0 == n2) || (qn3 == n2 && qn0 == n1)) {
        n3 = qn1;
        n4 = qn2;
        return;
    }

    cout << " Warning: You should not come here " << endl;
    cout << " search for " << n1 << "  " << n2 << endl;
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

FaceType create_quad(const FaceType t1, const FaceType t2) {
    vector<NodeType> connect;
    NodeType commonnodes[3];

    connect = t1->getConnection();

    int index = 0;
    for (int i = 0; i < 3; i++) {
        if (t2->hasNode(connect[i]))
            commonnodes[index++] = connect[i];
    }

    assert(index == 2);

    NodeType ot1 = Face::opposite_node(t1, commonnodes[0], commonnodes[1]);
    NodeType ot2 = Face::opposite_node(t2, commonnodes[0], commonnodes[1]);

    connect.resize(4);
    connect[0] = ot1;
    connect[1] = commonnodes[0];
    connect[2] = ot2;
    connect[3] = commonnodes[1];

    Face *qface = new Face;
    qface->setConnection(connect);
    return qface;
}

///////////////////////////////////////////////////////////////////////////////

void tri2quads(Mesh *mesh, vector<FacePair> &matching) {
    assert(matching.size());
    //
    // Checking match ....All the nodes must appear only once in the match.
    //
    vector<int> seq;
    seq.resize(2 * matching.size());
    for (size_t i = 0; i < matching.size(); i++) {
        int f1 = matching[i].first;
        int f2 = matching[i].second;
        seq[2 * i] = f1;
        seq[2 * i + 1] = f2;
    }
    sort(seq.begin(), seq.end());

    for (size_t i = 0; i < 2 * matching.size() - 1; i++) {
        if (seq[i] == seq[i + 1]) {
            cout << "Error: Invalid Matching" << endl;
            exit(0);
        }
    }

    //end of checking ...


    // Merge matching triangles into quad. Unmatched triangles remains in
    // the mesh.
    //
    size_t numfaces = mesh->getSize(2);

    Face *tri1, *tri2, *quad;
    for (size_t i = 0; i < numfaces; i++) {
        mesh->getFace(i)->setRemoveMark(0);
    }

    vector<Face*> faces;
    faces.resize(matching.size());

    for (size_t i = 0; i < matching.size(); i++) {
        int f1 = matching[i].first;
        int f2 = matching[i].second;
        tri1 = mesh->getFace(f1);
        tri2 = mesh->getFace(f2);
        tri1->setRemoveMark(1);
        tri2->setRemoveMark(1);
        quad = create_quad(tri1, tri2);
        mesh->addFace(quad);
    }
    mesh->prune();

    //
    // Check #of triangles and #quads in the mesh. if we are lucky then
    // number of triangles must be zero.
    //
    numfaces = mesh->getSize(2);

    int numtris = 0;
    int numquads = 0;

    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFace(i);
        if (face->getSize(0) == 3) numtris++;
        if (face->getSize(0) == 4) numquads++;
    }

    cout << " Total Number of Triangles in the mesh " << numtris << endl;
    cout << " Total Number of Quads in the mesh " << numquads << endl;

}

///////////////////////////////////////////////////////////////////////////////

void Mesh::saveAs(const string &s) {
    string filename = s + ".dat";
    ofstream ofile(filename.c_str(), ios::out);

    int numnodes = nodes.size();
    int numfaces = faces.size();

    int nn = numnodes;

    // if( hasdual ) nn += numfaces;

    ofile << nn << " " << numfaces << endl;

    for (int i = 0; i < numnodes; i++) {
        Point3D p3d = nodes[i]->getXYZCoords();
        ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

/*
    if (hasdual) {
        for (int i = 0; i < numfaces; i++) {
            Vertex *v = faces[i]->getDualNode();
            assert(v);
            Point3D p3d = v->getXYZCoords();
            ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
        }
    }
*/

    vector<NodeType> connect;
    for (int i = 0; i < numfaces; i++) {
        Face *face = faces[i];
        connect = face->getConnection();
        int nnodes = connect.size();
        ofile << nnodes << " ";
        for (int j = 0; j < nnodes; j++)
            ofile << connect[j]->getID() << " ";
        ofile << endl;
    }

    for (int i = 0; i < numfaces; i++) {
        Vertex *v = faces[i]->getDualNode();
        if (v) ofile << v->isRemoved() << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::build_relations02() {
    if (adjTable[0][2] == 1) return 1;

    clear_relations(0, 2);

    int numfaces = getSize(2);

    for (int iface = 0; iface < numfaces; iface++) {
        Face *face = getFace(iface);
        assert(face);
        for (int j = 0; j < face->getSize(0); j++) {
            Vertex *vtx = face->getConnection(j);
            vtx->addRelation2(face);
        }
    }
    adjTable[0][2] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::build_relations00() {
    if (adjTable[0][0] == 1) return 1;

    clear_relations(0, 0);

    int numfaces = getSize(2);

    for (int iface = 0; iface < numfaces; iface++) {
        Face *face = getFace(iface);
        assert(face);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++) {
            Vertex *v0 = face->getConnection(j);
            Vertex *v1 = face->getConnection((j + 1) % nnodes);
            v0->addRelation0(v1);
            v1->addRelation0(v0);
        }
    }
    adjTable[0][0] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::clear_relations(int src, int dst) {
    int numnodes = getSize(0);

    if (src == 0 && dst == 0) {
        for (int i = 0; i < numnodes; i++) {
            Vertex *vtx = getNode(i);
            vtx->clearRelations(0);
        }
        adjTable[0][0] = 0;
    }

    if (src == 0 && dst == 2) {
        for (int i = 0; i < numnodes; i++) {
            Vertex *vtx = getNode(i);
            vtx->clearRelations(2);
        }
        adjTable[0][2] = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::search_boundary() {
    int relexist = build_relations(0, 2);

    int numfaces = getSize(2);

    vector<FaceType> neighs;
    for (int iface = 0; iface < numfaces; iface++) {
        Face *face = getFace(iface);
        assert(face);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++) {
            Vertex *v0 = face->getConnection(j);
            Vertex *v1 = face->getConnection((j + 1) % nnodes);
            neighs = Mesh::getRelation112(v0, v1);
            if (neighs.size() == 1) {
                v0->setBoundaryMark(1);
                v1->setBoundaryMark(1);
                int bmark = max(1, face->getBoundaryMark());
                face->setBoundaryMark(bmark);
            }
        }
    }

    if (!relexist) clear_relations(0, 2);
}
///////////////////////////////////////////////////////////////////////////////

Mesh *struct_tri_grid(int nx, int ny)
{
    Mesh *trimesh = new Mesh;

    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);

    Point3D xyz;

    int index = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            xyz[0] = -1.0 + i*dx;
            xyz[1] = -1.0 + j*dy;
            xyz[2] = 0.0;
            Vertex *vnew = Vertex::newObject();
            vnew->setID(index++);
            vnew->setXYZCoords(xyz);
            trimesh->addNode(vnew);
        }
    }

    vector<Vertex*> connect(3);
    index = 0;
    Face *newtri;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i;
            int n1 = n0 + 1;
            int n2 = n1 + nx;
            int n3 = n0 + nx;
            connect[0] = trimesh->getNode(n0);
            connect[1] = trimesh->getNode(n1);
            connect[2] = trimesh->getNode(n2);
            newtri = new Face;
            newtri->setConnection(connect);
            trimesh->addFace(newtri);

            connect[0] = trimesh->getNode(n0);
            connect[1] = trimesh->getNode(n2);
            connect[2] = trimesh->getNode(n3);
            newtri = new Face;
            newtri->setConnection(connect);
            trimesh->addFace(newtri);
        }
    }
    return trimesh;
}

/////////////////////////////////////////////////////////////////////////////

Mesh *struct_quad_grid(int nx, int ny) {
    Mesh *quadmesh = new Mesh;

    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);

    Point3D xyz;

    int index = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            xyz[0] = -1.0 + i*dx;
            xyz[1] = -1.0 + j*dy;
            xyz[2] = 0.0;
            Vertex *vnew = Vertex::newObject();
            vnew->setID(index++);
            vnew->setXYZCoords(xyz);
            quadmesh->addNode(vnew);
        }
    }

    vector<Vertex*> connect(4);
    index = 0;
    Face *newquad;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i;
            int n1 = n0 + 1;
            int n2 = n1 + nx;
            int n3 = n0 + nx;
            connect[0] = quadmesh->getNode(n0);
            connect[1] = quadmesh->getNode(n1);
            connect[2] = quadmesh->getNode(n2);
            connect[3] = quadmesh->getNode(n3);
            newquad = new Face;
            newquad->setConnection(connect);
            quadmesh->addFace(newquad);
        }
    }
    return quadmesh;
}

////////////////////////////////////////////////////////////////////

void expand_strip(Face *prevface, Vertex *v0, Vertex *v1, list<Face*> &strip) {
    vector<Face*> neighs = Mesh::getRelation112(v0, v1);

    Vertex *vn0, *vn1; // Next-Edge Nodes;

    Face *nextface;
    if (neighs.size() == 2) {
        if (!neighs[0]->isVisited() && neighs[1]->isVisited()) {
            nextface = neighs[0];
            if (nextface->getSize(0) == 4) {
                Face::opposite_nodes(nextface, v0, v1, vn0, vn1);
                nextface->setVisitMark(1);
                strip.push_back(nextface);
                expand_strip(nextface, vn0, vn1, strip);
            }
        }
        if (neighs[0]->isVisited() && !neighs[1]->isVisited()) {
            nextface = neighs[1];
            if (nextface->getSize(0) == 4) {
                Face::opposite_nodes(nextface, v0, v1, vn0, vn1);
                nextface->setVisitMark(1);
                strip.push_back(nextface);
                expand_strip(nextface, vn0, vn1, strip);
            }
        }
    }
}
////////////////////////////////////////////////////////////////////

void get_quad_strip(Mesh *mesh, Face *rootface, vector<Face*> &strip1, vector<Face*> &strip2) {
    Vertex *v0, *v1;

    list<Face*> strip01, strip12, strip23, strip03;
    int numfaces = mesh->getSize(0);

    // Strip Starting from edge 0-1
    for (int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFace(i);
        face->setVisitMark(0);
    }

    v0 = rootface->getConnection(0);
    v1 = rootface->getConnection(1);

    rootface->setVisitMark(1);
    strip01.push_back(rootface);
    expand_strip(rootface, v0, v1, strip01);

    // Strip Starting from edge 2-3
    for (int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFace(i);
        face->setVisitMark(0);
    }

    v0 = rootface->getConnection(2);
    v1 = rootface->getConnection(3);

    rootface->setVisitMark(1);
    strip23.push_back(rootface);
    expand_strip(rootface, v0, v1, strip23);

    // Strip Starting from edge 1-2
    for (int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFace(i);
        face->setVisitMark(0);
    }

    v0 = rootface->getConnection(1);
    v1 = rootface->getConnection(2);

    rootface->setVisitMark(1);
    strip12.push_back(rootface);
    expand_strip(rootface, v0, v1, strip12);

    // Strip Starting from edge 0-3
    for (int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFace(i);
        face->setVisitMark(0);
    }

    v0 = rootface->getConnection(0);
    v1 = rootface->getConnection(3);

    rootface->setVisitMark(1);
    strip12.push_back(rootface);
    expand_strip(rootface, v0, v1, strip03);
}

///////////////////////////////////////////////////////////////////////////////

void mesh_quality(Mesh *mesh) {
    cout << "Vertex-Face Degree Distribution " << endl;
    getVertexFaceDegrees(mesh);
}

///////////////////////////////////////////////////////////////////////////////

#ifdef USE_MOAB

int get_moab_mesh(const Mesh *mesh, iMesh_Instance &imesh) {
    int status, err;
    int optlen = 0;
    char *options = NULL;

    iMesh_newMesh(options, &imesh, &err, optlen);
    assert(!err);

    iBase_EntityHandle newHandle;

    int numnodes = mesh->getSize(0);
    vector<iBase_EntityHandle> nodeHandles(numnodes);
    for (int i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNode(i);
        Point3D p = v->getXYZCoords();
        assert(v->getID() == i);
        iMesh_createVtx(imesh, p[0], p[1], p[2], &newHandle, &err);
        nodeHandles[i] = newHandle;
    }

    int numfaces = mesh->getSize(2);
    vector<iBase_EntityHandle> connect;

    for (int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFace(i);
        vector<Vertex*> facevtx = face->getConnection();
        int nnodes = facevtx.size();
        connect.resize(nnodes);
        for (int j = 0; j < nnodes; j++) {
            Vertex *v = face->getConnection(j);
            connect[j] = nodeHandles[v->getID()];
        }
        switch (nnodes) {
            case 3:
                iMesh_createEnt(imesh, iMesh_TRIANGLE, &connect[0], nnodes,
                        &newHandle, &status, &err);
                assert(!err);
                break;
            case 4:
                iMesh_createEnt(imesh, iMesh_QUADRILATERAL, &connect[0], nnodes,
                        &newHandle, &status, &err);
                assert(!err);
                break;
            default:
                iMesh_createEnt(imesh, iMesh_POLYGON, &connect[0], nnodes,
                        &newHandle, &status, &err);
                assert(!err);
                break;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void mesh_optimization(Mesh *mesh) {
    int err;

    iMesh_Instance imesh;
    get_moab_mesh(mesh, imesh);

    iBase_EntitySetHandle rootSet;
    iMesh_getRootSet(imesh, &rootSet, &err);

    /*
        const char *outfile = "moab.vtk";
        int namelen = strlen(outfile);
        iMesh_save(imesh, rootSet, outfile, NULL, &err, namelen, 0);
        exit(0);
     */

    Mesquite::MsqError ierr;
    Mesquite::Mesh* mesqmesh = new Mesquite::MsqIMesh(imesh, rootSet, iBase_ALL_TYPES, ierr, "fixed");
    assert(!ierr);

    //  Mesquite::MeshWriter::write_vtk(mesqmesh, "original.vtk", ierr);

    Mesquite::PlanarDomain domain(Mesquite::PlanarDomain::XY);

    Mesquite::LaplacianIQ laplacian_smoother;
    laplacian_smoother.run_instructions(mesqmesh, &domain, ierr);
    if (ierr) return;


    Mesquite::ShapeImprovementWrapper shape_wrapper(ierr);
    if (ierr) {
        cout << "Shape wrapper error " << ierr << endl;
        exit(2);
    }
    shape_wrapper.run_instructions(mesqmesh, &domain, ierr);

    if (ierr) {
        cout << "Error smoothing mesh " << ierr << endl;
        return;
    }
}

#endif
