#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cassert>
#include <fstream>
#include <math.h>
#include <string>

#include <vector>
#include <list>
#include <map>
#include <deque>
#include <algorithm>

#include <boost/foreach.hpp>
#include <boost/array.hpp>
#include <boost/utility.hpp>
#include <boost/any.hpp>

#ifdef USE_MOAB
#include <iMesh.h>
#include <MBInterface.hpp>
#include "SimpleArray.h"
#include <Mesquite_all_headers.hpp>
#endif

///////////////////////////////////////////////////////////////////////////////

using namespace std;

typedef boost::array<double, 3 > Point3D;
typedef boost::array<double, 3 > Vec3D;

class Vertex;
class Face;
class Edge;

typedef Vertex* NodeType;
typedef Edge* EdgeType;
typedef Face* FaceType;

typedef int AttribKey;

///////////////////////////////////////////////////////////////////////////////

struct Attribute {

    bool operator ==(const Attribute & rhs) const {
        return key == rhs.key;
    }
    AttribKey key;
    boost::any value;
};

///////////////////////////////////////////////////////////////////////////////

class MeshEntity {
public:
    static AttribKey getAttribKey(const string &s);
    static AttribKey addAttribute(const string &s);
    static bool hasAttribute(const string &);
    static void removeAttribute(const string &s);

    void setVisitMark(bool r) {
        visitMark = r;
    }

    bool isVisited() const {
        return visitMark;
    }

    void setRemoveMark(bool r) {
        removeMark = r;
    }

    bool isRemoved() const {
        return removeMark;
    }

    bool isBoundary() const {
        return boundarymark;
    }

    void setBoundaryMark(int b) {
        boundarymark = b;
    }

    int getBoundaryMark() const {
        return boundarymark;
    }

    bool isConstrained() const {
        if (boundarymark > 0 || constrained) return 1;
        return 0;
    }

    void setConstrainedMark(int b) {
        constrained = b;
    }

    void setID(int id) {
        gid = id;
    }

    int getID() const {
        return gid;
    }

    template<class T>
    void setAttribute(const string &s, const T &val) {
        AttribKey key = getAttribKey(s);
        setAttribute(key, val);
    }

    template<class T>
    void setAttribute(const AttribKey &k, const T &val) {
        Attribute property;
        property.key = k;

        list<Attribute>::iterator it;

        it = std::find(attriblist.begin(), attriblist.end(), property);

        if (it == attriblist.end()) {
            property.value = val;
            attriblist.push_back(property);
        } else
            it->value = val;
    }

    template<class T>
    int getAttribute(const string &s, T &val) const {
        AttribKey key = getAttribKey(s);
        return getAttribute(key, val);
    }

    template<class T>
    int getAttribute(const AttribKey &k, T &val) const {
        Attribute property;
        property.key = k;

        list<Attribute>::const_iterator it;

        it = std::find(attriblist.begin(), attriblist.end(), property);

        if (it == attriblist.end()) return 1;

        boost::any anyvalue = it->value;

        if (anyvalue.type() == typeid (T)) {
            val = boost::any_cast<T > (anyvalue);
            return 0;
        }
        return 1;
    }

protected:
    int gid;
    bool constrained;
    int boundarymark;
    volatile bool visitMark, removeMark;
    std::list<Attribute> attriblist;
    static AttribKey maxAttribID;
    static map<string, AttribKey> attribKeyMap;
};

///////////////////////////////////////////////////////////////////////////////

class Vertex : public MeshEntity {
    static int global_id;

public:

    static Vertex* newObject();

    Vertex() {
        visitMark = 0;
        removeMark = 0;
        boundarymark = 0;
        mate = 0;
    }

    static Point3D mid_point(const Vertex *v0, const Vertex *v1);

    void setXYZCoords(const Point3D &p) {
        xyz[0] = p[0];
        xyz[1] = p[1];
        xyz[2] = p[2];
    }

    const Point3D &getXYZCoords() const {
        return xyz;
    }

    int getDegree() const {
        int ncount = 0;
        for (size_t i = 0; i < relations00.size(); i++)
            if (!relations00[i]->isRemoved()) ncount++;
        return ncount;
    }

    void addRelation0(Vertex *vertex) {
        if (!hasRelation0(vertex)) relations00.push_back(vertex);
        sort(relations00.begin(), relations00.end());
    }

    void removeRelation0(const Vertex *vertex) {
        if (!hasRelation0(vertex)) {
            vector<NodeType>::iterator it;
            it = remove(relations00.begin(), relations00.end(), vertex);
            relations00.erase(it, relations00.end());
            sort(relations00.begin(), relations00.end());
        }
    }

    void addRelation2(Face* face) {
        if (!hasRelation2(face)) relations02.push_back(face);
        sort(relations02.begin(), relations02.end());
    }

    void removeRelation2(const Face *face) {

        if (!hasRelation2(face)) {
            vector<FaceType>::iterator it;
            it = remove(relations02.begin(), relations02.end(), face);
            relations02.erase(it, relations02.end());
            sort(relations02.begin(), relations02.end());
        }
    }

    void clearRelations(int t) {
        if (t == 0) relations00.clear();
        if (t == 2) relations02.clear();
    }

    const vector<NodeType> &getRelations0() const {
        return relations00;
    }

    const vector<FaceType> &getRelations2() const {
        return relations02;
    }

    bool hasRelation0(const Vertex* vertex) const {
        vector<NodeType>::const_iterator it;
        it = find(relations00.begin(), relations00.end(), vertex);
        if (it == relations00.end()) return 0;

        return 1;
    }

    bool hasRelation2(const Face* face) const {
        vector<FaceType>::const_iterator it;
        it = find(relations02.begin(), relations02.end(), face);
        if (it == relations02.end()) return 0;

        return 1;
    }

    void setDualMate(Vertex *v) {
        mate = v;
    }

    Vertex* getDualMate() const {
        return mate;
    }

private:
    // void * operator new( size_t size, void *);

#ifdef USE_MOAB
    iBase_EntityHandle vertexHandle;
#endif

    Vertex *mate;
    vector<NodeType> relations00; // vertex-vertex
    vector<FaceType> relations02; // vertex-face

    Point3D xyz;
};

///////////////////////////////////////////////////////////////////////////////

class Edge : public MeshEntity {
public:

    Edge() {
    }

    Edge(NodeType n1, NodeType n2) {
        assert(n1 != n2);
        connect[0] = n1;
        connect[1] = n2;
    }

    void setConnection(const vector<NodeType> &v) {
        assert(v.size() == 2);
        connect[0] = v[0];
        connect[1] = v[1];
    }

    NodeType getConnection(int id) const {
        return connect[id];
    }

private:
    NodeType connect[2];
};

///////////////////////////////////////////////////////////////////////////////

class Face : public MeshEntity {
public:

    Face() {
        removeMark = 0;
        boundarymark = 0, dualnode = 0;
    }

    static NodeType opposite_node(const FaceType tri, NodeType n1, NodeType n2);
    static void opposite_nodes(const FaceType quad, NodeType n1, NodeType n2,
            NodeType &n3, NodeType &n4);

    static int check_on_boundary(const FaceType tri);

    bool hasNode(const NodeType &vertex) const {
        if (find(connect.begin(), connect.end(), vertex) != connect.end()) return 1;
        return 0;
    }

    int queryNodeAt(const NodeType &vertex) const {
        for (size_t i = 0; i < connect.size(); i++)
            if (connect[i] == vertex) return i;

        return -1;
    }

    int replaceNode(const Vertex *oldvertex, Vertex *newvertex) {

        for (size_t i = 0; i < connect.size(); i++) {
            if (connect[i] == oldvertex) {
                connect[i] = newvertex;
                return 0;
            }
        }
        return 1;
    }

    int getSize(int etype) const {
        if (etype == 0) return connect.size();
        return 0;
    }

    void setConnection(const vector<NodeType> &v) {
        connect = v;
    }

    const vector<NodeType> &getConnection() const {
        return connect;
    }

    NodeType getConnection(int id) const {
        return connect[id];
    }

    Point3D getCentroid() const;

    void setDualNode(const NodeType n) {
        dualnode = n;
    }

    NodeType getDualNode() const {
        return dualnode;
    }

private:

#ifdef USE_MOAB
    iBase_EntityHandle faceHandle;
#endif
    NodeType dualnode;
    vector<NodeType> connect;
};

///////////////////////////////////////////////////////////////////////////////

class Mesh {
public:

    Mesh() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) adjTable[i][j] = 0;
        }
        hasdual = 0;
    }
    static vector<FaceType> getRelation112(NodeType v0, NodeType v1);

    bool hasDual() const {
        return hasdual;
    }

    void setDual(bool d) {
        hasdual = d;
    }

    int getNumOfComponents();

    void readData(const string &f);

    bool getAdjTable(int i, int j) const {
        return adjTable[i][j];
    }

    int getEulerCharacteristic() const;

    size_t getSize(int d) {
        if (d == 0) return nodes.size();
        if (d == 1) return count_edges();
        if (d == 2) return faces.size();
        return 0;
    }

    size_t getBoundarySize(int d) const;

    void addNode(NodeType v) {
        nodes.push_back(v);
    }

    NodeType getNode(size_t id) const {
        assert(id >= 0 && id < nodes.size());
        return nodes[id];
    }

    void addFace(FaceType v) {
        faces.push_back(v);
    }

    FaceType getFace(size_t id) const {
        return faces[id];
    }

    // Get rid of entities which are marked "removed" from the mesh. 
    // a la Lazy garbage collection.
    //
    void prune();

    // Renumber mesh entities starting from index = 0
    void enumerate(int etype);

    // Search the boundary of the mesh (nodes, edges, and faces).
    void search_boundary();

    int build_relations(int src, int dst) {
        if (src == 0 && dst == 0) return build_relations00();
        if (src == 0 && dst == 2) return build_relations02();
        return 1;
    }

    void clear_relations(int src, int dst);

    // Return all the edges of the primal mesh ...
    vector<Edge*> getEdges();

    // Return all the edges of dual mesh...
    vector<Edge*> getDualEdges() const;

    // Returm all the matching edges of the dual graph...
    vector<Edge*> getMatchingDuals() const;

    // Save the mesh and all its attributes ( in the simple format ).
    void saveAs(const string &s);

private:
    volatile char adjTable[4][4];
    size_t count_edges();

    bool hasdual;

    vector<NodeType> nodes;
    vector<EdgeType> edges;
    vector<FaceType> faces;

    string filename;
    std::map<int, int> global2local;
    void readNodes(const string &s);
    void readEdges(const string &s);
    void readFaces(const string &s);

    int build_relations00();
    int build_relations02();
};

///////////////////////////////////////////////////////////////////////////////
#include "dualgraph.h"
#include "binarytree.h"

///////////////////////////////////////////////////////////////////////////////
// Graph Matching operations ....
///////////////////////////////////////////////////////////////////////////////

int quadrangulate(Mesh *mesh);

void maximum_matching(BinaryTree *tree, vector<FacePair> &matching);
void maximum_matching(Mesh *m, vector<FacePair> &matching);
void maximum_matching(DualGraph *g, vector<FacePair> &matching);

void tri2quads(Mesh *mesh, vector<FacePair> &matching);

Mesh* hamiltonian_triangulation(Mesh *orgmesh);
Mesh* hamiltonian_quadrangulation(Mesh *orgmesh);

////////////////////////////////////////////////////////////////////////////////
// Mesh Cleanup operations ....
////////////////////////////////////////////////////////////////////////////////

void get_strips(Mesh *mesh, Face *face, vector<Face*> &strip1, vector<Face*> strip2);

vector<int> getVertexFaceDegrees(Mesh *mesh);
vector<Face*> search_diamonds(Mesh *mesh);
vector<Vertex*> search_doublets(Mesh *mesh);

void cleanup_diamonds(Mesh *mesh);
void cleanup_doublets(Mesh *mesh);
void cleanup_boundary(Mesh *mesh);

Vertex* insert_doublet(Mesh *mesh, Face *face);
Vertex* insert_boundary_doublet(Mesh *mesh, Face *face);
Vertex* insert_doublet(Mesh *mesh, Face *face, Vertex *v0, Vertex *v2);

////////////////////////////////////////////////////////////////////////////////
//Helper functions ....
////////////////////////////////////////////////////////////////////////////////
Mesh* readOffData(const string &s);

Mesh *struct_tri_grid(int nx, int ny);
Mesh *struct_quad_grid(int nx, int ny);

void mesh_quality(Mesh *mesh);

////////////////////////////////////////////////////////////////////////////////
// Mesh Optimization ...
////////////////////////////////////////////////////////////////////////////////

void laplacian_smoothing(Mesh *mesh, int numIters);

#ifdef USE_MOAB
int get_moab_mesh(const Mesh *mesh, iMesh_Instance &imesh);
void mesh_optimization(Mesh *mesh);
#endif


#endif
