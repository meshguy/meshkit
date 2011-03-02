#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cassert>
#include <fstream>
#include <math.h>
#include <string>
#include <values.h>

#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <algorithm>

#ifdef USE_MESQUITE
#include <Mesquite_all_headers.hpp>
#endif

#ifdef USE_VERDICT
#include <verdict.h>
#endif

#define FOR_EACH(container, iter ) for(iter = container.begin(); iter != container.end(); ++iter)

#include "basic_math.hpp"
#include "tfiblend.hpp"

#define JAAL_SUCCESS 0
#define JAAL_GEOMETRIC_FAILURE    1
#define JAAL_TOP0LOGICAL_FAILURE  100

#define Break() \
{       \
   cout << "Break in file " << __FILE__ << " at " << __LINE__ << endl; \
   getchar(); \
} \


template <class T>
T next(T x) {
    return ++x;
}

template <class T, class Distance>
T next(T x, Distance n) {
    std::advance(x, n);
    return x;
}

template <class T>
T prior(T x) {
    return --x;
}

template <class T, class Distance>
T prior(T x, Distance n) {
    std::advance(x, -n);
    return x;
}

///////////////////////////////////////////////////////////////////////////////
using namespace std;

typedef std::pair<size_t, size_t> FacePair;

namespace Jaal {

class Vertex;
class Face;
class Edge;
class Mesh;

typedef Vertex* PNode;
typedef Edge* PEdge;
typedef Face* PFace;

#define SEQUENCE_IS_VECTOR  1

typedef std::vector<PNode> NodeSequence;
typedef std::vector<PEdge> EdgeSequence;
typedef std::vector<PFace> FaceSequence;

typedef std::list<PNode> NodeList;
typedef std::list<PEdge> EdgeList;
typedef std::list<PFace> FaceList;

typedef std::set<PNode> NodeSet;
typedef std::set<PFace> FaceSet;

#define X_AXIS  0
#define Y_AXIS  1
#define Z_AXIS  2

///////////////////////////////////////////////////////////////////////////////

template <class T>
int split_stl_vector(const std::vector<T> &a, size_t pos, std::vector<T> &b, std::vector<T> &c) {
    if (pos >= a.size()) return 1;

    b.resize(pos);

    size_t index = 0;
    for (size_t i = 0; i < pos; i++)
        b[index++] = a[i];

    c.resize(a.size() - pos + 1);
    index = 0;
    for (size_t i = pos - 1; i < a.size(); i++)
        c[index++] = a[i];

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

template <class T>
int split_stl_vector(const std::deque<T> &a, size_t pos, std::deque<T> &b, std::deque<T> &c) {
    if (pos >= a.size()) return 1;

    b.resize(pos);

    size_t index = 0;
    for (size_t i = 0; i < pos; i++)
        b[index++] = a[i];

    c.resize(a.size() - pos + 1);
    index = 0;
    for (size_t i = pos - 1; i < a.size(); i++)
        c[index++] = a[i];

    return 0;
}

typedef int AttribKey;

struct Attribute {

    bool operator ==(const Attribute & rhs) const {
        return key == rhs.key;
    }
    AttribKey key;

#ifdef USE_BOOST_LIBS
    boost::any value;
#endif
};

template<class T>
struct POD : public Attribute {
    T value;
};

template<class T>
struct StlVec : public Attribute {
    std::vector<T> values;
};

template<class T>
struct StlList : public Attribute {
    std::list<T> values;
};

template<class T>
struct StlSet : public Attribute {
    std::set<T> values;
};

class AttributeManager {
public:

    AttributeManager() {
        currAttribID = 0;
    }

    AttribKey registerAttribute(const string &s);
    AttribKey getAttribKey(const string &s);

    void unRegisterAttribute(const string &s);
    bool hasAttribute(const string &s) const;

    std::vector<string> getAttributeNames();
private:
    int currAttribID;
    std::map<string, int> attribKeyMap;
};

///////////////////////////////////////////////////////////////////////////////
typedef size_t Handle_t;

/// Base class for all handle types

class Handle {
public:

    explicit Handle(Handle_t hid = 0) : hid(hid) {
    }

    /// Get the underlying index of this handle

    Handle_t getHandle() const {
        return hid;
    }

    /// The handle is valid iff the index is not equal to 0.

    bool is_valid() const {
        return hid != 0;
    }

    /// reset handle to be invalid

    void reset() {
        hid = 0;
    }

    /// reset handle to be invalid

    void invalidate() {
        hid = 0;
    }

    bool operator==(const Handle& rhs) const {
        return hid == rhs.hid;
    }

    bool operator!=(const Handle& rhs) const {
        return hid != rhs.hid;
    }

    bool operator<(const Handle& rhs) const {
        return hid < rhs.hid;
    }

    // this is to be used only by the iterators

    void increment() {
        ++hid;
    }

    void decrement() {
        --hid;
    }

private:
    Handle_t hid;
};

class MeshEntity {
public:

    typedef size_t idtype;
    static const int ACTIVE = 0;
    static const int REMOVE = 1;
    static const int INACTIVE = 2;

    MeshEntity() {
        iTag = 0;
        layerID = -1;
        visitMark = 0; // Default: Not yet visited.
        statusMark = ACTIVE; // Default: Active< not removable.
        constrained = 0; // Default: No Contrainted
        boundarymark = 0; // Default: Internal entity
    }

    void setVisitMark(bool r) {
        visitMark = r;
    }

    bool isVisited() const {
        return visitMark;
    }

    void setStatus(char a) {
        statusMark = a;
    }

    const bool isRemoved() const {
        if (statusMark == REMOVE) return 1;
        return 0;
    }

    bool isActive() const {
        if (statusMark == ACTIVE) return 1;
        return 0;
    }

    const bool isBoundary() const {
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

    void setID(idtype id) {
        gid = id;
    }

    const idtype &getID() const {
        return gid;
    }

    void setTag(int v) {
        iTag = v;
    }

    const int &getTag() const {
        return iTag;
    }

#ifdef USE_BOOST_LIBS

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
#endif

    void setAttribute(const AttribKey &k, Attribute *newattrib) {
        Attribute *oldattrib = getAttribute(k);
        if (oldattrib == NULL) {
            attributes.push_back(newattrib);
            return;
        }

        newattrib->key = k;
        for (int i = 0; i < attributes.size(); i++) {
            if (attributes[i]->key == k) {
                attributes[i] = newattrib;
            }
        }

        delete oldattrib;
    }

    Attribute* getAttribute(const AttribKey &k) {
        for (int i = 0; i < attributes.size(); i++)
            if (attributes[k]->key == k) return attributes[k];
        return NULL;
    }

    void deleteAttribute(const AttribKey &k) {
        Attribute *a = getAttribute(k);
        if (a == NULL) return;

        vector<Attribute*>::iterator it;
        it = remove(attributes.begin(), attributes.end(), a);
        attributes.erase(it, attributes.end());
        delete a;
    }

    void setLayerID(int l) {
        layerID = l;
    }

    int getLayerID() const {
        return layerID;
    }

    void clearRelations(int t);

    void addRelation0(Vertex *vertex);
    void addRelation2(Face* face);

    bool hasRelation0(const Vertex* vertex) const;
    bool hasRelation2(const Face* face) const;

    void removeRelation0(const Vertex *vertex);
    void removeRelation2(const Face *face);

    const NodeSequence &getRelations0() const;
    const FaceSequence &getRelations2() const;

    size_t getNumOfRelations(int e) {
        if (e == 0)
            return relations0.size();
        if (e == 2)
            return relations2.size();
        return 0;
    }

protected:
    idtype gid;
    int boundarymark;
    int iTag;
    int layerID;
    bool constrained;
    volatile bool visitMark;
    volatile char statusMark;

    NodeSequence relations0; // vertex-vertex
    FaceSequence relations2; // vertex-face

    std::vector<Attribute*> attributes;
};

struct EntityRemovedPred {

    bool operator() (const MeshEntity * entity) const {
        if (entity) return entity->isRemoved();
        return 0;
    }
};

///////////////////////////////////////////////////////////////////////////////

class Vertex : public MeshEntity {
    static size_t global_id;
public:
    typedef Handle_t VertexHandle;

    static PNode newObject();
    //  static AttributeManager  attribManager;

    Vertex() {
        visitMark = 0;
        statusMark = MeshEntity::ACTIVE;
        boundarymark = 0;
        mate = 0;
        primalface = 0;
    }

    static double length(const Vertex *v0, const Vertex *v1);
    static double length2(const Vertex *v0, const Vertex *v1);
    static Point3D mid_point(const Vertex *v0, const Vertex *v1, double alpha = 0.5);

    void setXYZCoords(const Point3D &p) {
        xyz[0] = p[0];
        xyz[1] = p[1];
        xyz[2] = p[2];
    }

    const Point3D &getXYZCoords() const {
        return xyz;
    }

    void setDualMate(Vertex *v) {
        mate = v;
    }

    const Vertex* getDualMate() const {
        return mate;
    }

    void setPrimalFace(Face *f) {
        primalface = f;
    }

    Face* getPrimalFace() const {
        return primalface;
    }

    Vertex* getClone() const;

    double getFeatureLength() const;
    double getFeatureAngle() const;
    int get_ideal_vertex_degree() const;

private:
    // void * operator new( size_t size, void *);

    Vertex *mate;
    Face *primalface;

    Point3D xyz;
};

////////////////////////////////////////////////////////////////////////////////

struct LowVertexDegreeCompare {

    bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
        NodeSequence neighs;
        neighs = vertex1->getRelations0();
        size_t d1 = neighs.size();
        neighs = vertex2->getRelations0();
        size_t d2 = neighs.size();
        return d1 < d2;
    }
};
////////////////////////////////////////////////////////////////////////////////

struct HighVertexDegreeCompare {

    bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
        NodeSequence neighs;
        neighs = vertex1->getRelations0();
        size_t d1 = neighs.size();
        neighs = vertex2->getRelations0();
        size_t d2 = neighs.size();
        return d1 > d2;
    }
};
////////////////////////////////////////////////////////////////////////////////

struct LowVertexTagCompare {

    bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
        int val1 = vertex1->getTag();
        int val2 = vertex2->getTag();
        return val1 < val2;
    }
};

////////////////////////////////////////////////////////////////////////////////

struct HighVertexTagCompare {

    bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
        int val1 = vertex1->getTag();
        int val2 = vertex2->getTag();
        return val1 > val2;
    }
};

////////////////////////////////////////////////////////////////////////////////

struct HighLayerCompare {

    bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
        int val1 = vertex1->getLayerID();
        int val2 = vertex2->getLayerID();
        return val1 > val2;
    }
};

///////////////////////////////////////////////////////////////////////////////

class Edge : public MeshEntity {
public:
    typedef Handle_t EdgeHandle;

    static PEdge newObject();
    //  static AttributeManager attribManager;

    Edge() {
    }

    Edge(PNode n1, PNode n2) {
        setNodes(n1, n2);
    }

    Vertex* getHashNode() const {
        return min(connect[0], connect[1]);
    }

    void setNodes(const NodeSequence &v) {
        assert(v.size() == 2);
        setNodes(v[0], v[1]);
    }

    void setNodes(Vertex *v1, Vertex *v2) {
        assert(v1 != v2);
        connect[0] = v1;
        connect[1] = v2;
    }

    const PNode &getNodeAt(int id) const {
        return connect[id];
    }

    bool isSameAs(const Edge &rhs) {
        if ((connect[0] == rhs.connect[0]) && (connect[1] == rhs.connect[1])) return 1;
        if ((connect[0] == rhs.connect[1]) && (connect[1] == rhs.connect[0])) return 1;
        return 0;
    }

    bool hasCrease() const {
        return 0;
    }

    Edge* getClone() const {
        Edge *newedge = new Edge(connect[0], connect[1]);
        return newedge;
    }

private:
    //  void * operator new( size_t size, void *);
    static std::map<string, AttribKey> attribKeyMap;

    PNode connect[2];
};

///////////////////////////////////////////////////////////////////////////////

class Face : public MeshEntity {
public:
    typedef Handle_t FaceHandle;

    static const int POLYGON = 0;
    static const int TRIANGLE = 3;
    static const int QUADRILATERAL = 4;

    static PFace newObject();
    //  static AttributeManager attribManager;

    static PFace create_quad(const PFace t1, const PFace t2, int replace = 0);

    // This function is for concave quads. Arrange vertices so that OpenGL can render them
    // correctly. When you save the mesh, probably the connectivity may change.
    static int quad_tessalate(const NodeSequence &orgNodes, NodeSequence &rotatedNodes);
    static int hexagon_2_quads(const NodeSequence &hexnodes, FaceSequence &quads, int start_from);

    static int is_3_sided_convex_loop_quad_meshable(const int *s, int *sdiv);
    static int is_5_sided_convex_loop_quad_meshable(const int *s, int *sdiv);
    static int is_cyclic_quad(const Point3D &p0, const Point3D &p1, const Point3D &p2, const Point3D &p3);

    // Centroid of Triangle element ...
    static Vertex *centroid(const Vertex *v0, const Vertex *v1, const Vertex *v2);

    // Centroid of Quad  element ...
    static Vertex *centroid(const Vertex *v0, const Vertex *v1, const Vertex *v2,
            const Vertex *v3, const Vertex *v4);

    // Distortion of a triangle element ...
    static double distortion(const Vertex *v0, const Vertex *v1, const Vertex *v2);

    // Distortion of a quad element ...
    static double distortion(const Vertex *v0, const Vertex *v1, const Vertex *v2,
            const Vertex *v3, const Vertex *v4);

    /////////////////////////////////////////////////////////////////////////////
    // Use modified Heron formula for find out the area of a triangle
    /////////////////////////////////////////////////////////////////////////////

    static double tri_area(const Point3D &p0, const Point3D &p1, const Point3D &p2);

    static Vec3D normal(const Vertex *v0, const Vertex *v1, const Vertex *v2);
    static Vec3D normal(const Point3D &p0, const Point3D &p1, const Point3D &v2);

    static void bilinear_weights(double xi, double eta, vector<double> &weight);
    static double linear_interpolation(const vector<double> &x, const vector<double> &w);

    /////////////////////////////////////////////////////////////////////////////
    // Calculate Area of Quadrilateral:
    // Example : Convex Case
    //           Coorindates  (0.0,0.0)  ( 1.0, 1.0), (2.0, 0.0), (1,0, -1.0)
    //           Result :  2*( 0.5*2.0*1.0 ) = 2.0;
    //           Concave case:
    //           Coodinates: ( 0.0, 0.0), 10, 2), (2.0, 0.0), (10, -2)
    //           Result : 2 *( 0.5*2.0* 2.0) = 4.0
    /////////////////////////////////////////////////////////////////////////////
    static double quad_area(const Point3D &p0, const Point3D &p1,
            const Point3D &p2, const Point3D &p3);

    static bool is_convex_quad(const Point3D &p0, const Point3D &p1,
            const Point3D &p2, const Point3D &p3);

    Face() {
        statusMark = 0;
        boundarymark = 0;
        visitMark = 0;
        dualnode = 0;
        partID = 0;
    }

    static PNode opposite_node(const PFace quad, const PNode n1);
    static PNode opposite_node(const PFace tri, PNode n1, PNode n2);
    static void opposite_nodes(const PFace quad, PNode n1, PNode n2,
            PNode &n3, PNode &n4);

    static int check_on_boundary(const PFace tri);

    Vertex* getHashNode() const {
        return *min_element(connect.begin(), connect.end());
    }

    int getType() const {
        if (connect.size() == 3) return Face::TRIANGLE;
        if (connect.size() == 4) return Face::QUADRILATERAL;
        return Face::POLYGON;
    }

    int invertedAt() const;

    bool isSameAs(const Face *rhs) const {
        if (rhs->getSize(0) != getSize(0)) return 0;

        for (int i = 0; i < getSize(0); i++)
            if (!rhs->hasNode(connect[i])) return 0;

        return 1;
    }

    bool hasNode(const PNode &vertex) const {
        if (find(connect.begin(), connect.end(), vertex) != connect.end()) return 1;
        return 0;
    }

    int getPosOf(const PNode vertex) const {
        for (size_t i = 0; i < connect.size(); i++)
            if (connect[i] == vertex) return i;

        return -1;
    }

    void reverse() {
        std::reverse(connect.begin(), connect.end());
    }

    int getOrientation(const Vertex *ev0, const Vertex *ev1) const {
        size_t nsize = connect.size();

        for (int i = 0; i < nsize; i++) {
            Vertex *v0 = connect[(i + 0) % nsize];
            Vertex *v1 = connect[(i + 1) % nsize];
            if (v0 == ev0 && v1 == ev1) return +1;
            if (v0 == ev1 && v1 == ev0) return -1;
        }
        return 0;
    }

    int replaceNode(const Vertex *oldvertex, Vertex *newvertex) {

        for (size_t i = 0; i < connect.size(); i++) {
            if (connect[i] == oldvertex) {
                connect[i] = newvertex;
                return 0;
            }
        }

        cout << " Warning: Requested node not replaced " << endl;
        return 1;
    }

    int getSize(int etype) const {
        if (etype == 0) return connect.size();
        return 0;
    }

    void setNodes(const NodeSequence &v) {
        for (size_t i = 0; i < v.size(); i++) {
            for (size_t j = 0; j < v.size(); j++) {
                if (i != j) assert(v[i] != v[j]);
            }
        }
        connect = v;
    }

    const NodeSequence &getNodes() const {
        return connect;
    }

    const PNode &getNodeAt(int id) const {
        assert(id < connect.size());
        return connect[id];
    }

    Point3D getCentroid() const;

    void setDualNode(const PNode n) {
        dualnode = n;
        dualnode->setPrimalFace(this);
    }

    PNode getDualNode() const {
        return dualnode;
    }

    bool hasBoundaryNode() const {
        for (int i = 0; i < connect.size(); i++)
            if (connect[i]->isBoundary()) return 1;
        return 0;
    }

    bool has_all_bound_nodes() const;
    bool has_boundary_edge() const;

    FaceSequence getRelations202();
    FaceSequence getRelations212();
    NodeSequence getRelations0();

    bool isConvex() const {
        if (connect.size() <= 3) return 1;

        if (connect.size() == 4)
            return is_convex_quad(connect[0]->getXYZCoords(),
                connect[1]->getXYZCoords(),
                connect[2]->getXYZCoords(),
                connect[3]->getXYZCoords());

        return 0;
    }

    double getAngleAt(const Vertex *v) const;

    vector<double> get_interior_angles() const;

    double getAspectRatio();

    double getArea() {
        if (connect.size() == 3) {
            return tri_area(connect[0]->getXYZCoords(),
                    connect[1]->getXYZCoords(),
                    connect[2]->getXYZCoords());
        }

        if (connect.size() == 4) {
            return quad_area(connect[0]->getXYZCoords(),
                    connect[1]->getXYZCoords(),
                    connect[2]->getXYZCoords(),
                    connect[3]->getXYZCoords());
        }
        return 0.0;
    }

    vector<Face> triangulate();

    bool isValid() const;

    const Vec3D &getNormal() const {
        return fnormal;
    }

    void setNormal(const Vec3D &n) {
        fnormal = n;
    }

    void setPartID(int id) {
        partID = id;
    }

    int getPartID() const {
        return partID;
    }

    int refine_quad14(NodeSequence &newnodes, FaceSequence &newfaces);
    int refine_quad15(NodeSequence &newnodes, FaceSequence &newfaces);

    void start_from_concave_corner();

    Face* getClone() const {
        Face *newface = Face::newObject();
        newface->setNodes(connect);
        return newface;
    }

private:
    // void * operator new( size_t size, void *);
    NodeSequence connect;
    static std::map<string, AttribKey> attribKeyMap;

    int partID;
    PNode dualnode;
    Vec3D fnormal;

    int refine_convex_quad15(NodeSequence &newnodes, FaceSequence &newfaces);
    int refine_concave_quad15(NodeSequence &newnodes, FaceSequence &newfaces);
};

///////////////////////////////////////////////////////////////////////////////

struct BoundingBox {

    double getLength(int d) {
        assert(d >= 0 && d < 3);
        return fabs(upperRightCorner[d] - lowerLeftCorner[d]);
    }

    void setLowerLeftCorner(const Point3D & p) {
        lowerLeftCorner = p;
    }

    const Point3D & getLowerLeftCorner() const {
        return lowerLeftCorner;
    }

    void setUpperRightCorner(const Point3D & p) {
        upperRightCorner = p;
    }

    const Point3D & getUpperRightCorner() const {
        return upperRightCorner;
    }

private:
    Point3D lowerLeftCorner;
    Point3D upperRightCorner;
};

inline
PNode
Vertex::newObject() {
    Vertex *v = new Vertex;
    assert(v);
    v->setID(global_id);
    global_id++;
    return v;
}

///////////////////////////////////////////////////////////////////////////////

inline PNode Vertex::getClone() const {
    Vertex *v = new Vertex;
    assert(v);
    v->setID(global_id);
    v->setXYZCoords(xyz);
    return v;
}

///////////////////////////////////////////////////////////////////////////////

inline PEdge Edge::newObject() {
    Edge *e = new Edge;
    assert(e);
    return e;
}

///////////////////////////////////////////////////////////////////////////////

inline PFace Face::newObject() {
    Face *f = new Face;
    assert(f);
    return f;
}

///////////////////////////////////////////////////////////////////////////////

#define TOPOLOGICAL_DISTANCE  0 
#define EUCLIDEAN_DISTANCE    1 
#define EUCLIDEAN_DISTANCE2   2 
#define CITY_BLOCK_DISTANCE   3

///////////////////////////////////////////////////////////////////////////////

struct MeshFilter {

    virtual bool pass(const Vertex * v) const {
        cout << "Warning: passing the base " << endl;
        return 1;
    }

    virtual bool pass(const Face * v) const {
        return 1;
    };
};

class Mesh {
public:
    static NodeSequence generate_nodes(size_t n);
    static FaceSequence generate_faces(size_t n);

    static FaceSequence getRelations112(PNode v0, PNode v1);
    static FaceSequence getRelations102(PNode v0, PNode v1);
    static int make_chain(vector<Edge> &edges);
    static int is_closed_chain(const vector<Edge> &edges);
    static int is_closeable_chain(const vector<Edge> &edges);
    static int rotate_chain(vector<Edge> &edges, Vertex* start_vertex);
    static NodeSequence boundary_chain_nodes(Vertex *v0, Vertex *v1);
    static NodeSequence chain_nodes(const vector<Edge> &e);

    Mesh() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) adjTable[i][j] = 0;
        }
        hasdual = 0;
        meshname = "unknown";
        boundary_known = 0;
    }

    ~Mesh() {
        emptyAll();
    }

    void setName(const string &s) {
        meshname = s;
    }

    string getName() const {
        return meshname;
    }

    void reserve(size_t nSize, int entity) {
#ifdef SEQUENCE_IS_VECTOR
        if (entity == 0) nodes.reserve(nSize);
        if (entity == 1) edges.reserve(nSize);
        if (entity == 2) faces.reserve(nSize);
#endif
    }

    bool hasDual() const {
        return hasdual;
    }

    void setDual(bool d) {
        hasdual = d;
    }

    int remove_nodes_attribute(const string &s);
    int remove_edges_attribute(const string &s);
    int remove_faces_attribute(const string &s);

    int isHomogeneous() const;

    int hasConvexCells();

    size_t count_convex_faces();
    size_t count_concave_faces();
    size_t count_inverted_faces();
    size_t count_irregular_nodes(int regdef);

    Mesh* deep_copy();
    Mesh* shallow_copy();

    BoundingBox getBoundingBox() const;

    int getNumOfComponents(bool stop_at_interface = 0);
    Mesh *getComponent(int id);

    int readFromFile(const string &f);

    bool getAdjTable(int i, int j) const {
        return adjTable[i][j];
    }

    int getEulerCharacteristic() {
        size_t F = getSize(2);
        size_t E = getSize(1);
        size_t V = getSize(0);

        return F - E + V;
    }

    Vertex* nearest_neighbour(const Vertex *v, double &d);

    size_t getSize(int d) {
        if (d == 0) return nodes.size();
        if (d == 1) return count_edges();
        if (d == 2) return faces.size();
        return 0;
    }

    size_t getCapacity(int d) {
#ifdef SEQUENCE_IS_VECTOR
        if (d == 0) return nodes.capacity();
        if (d == 2) return faces.capacity();
#endif
        return 0;
    }

    // If the topology of the mesh is changed, make sure to setBoundaryStatus = 0
    // so that it is again determined.

    void setBoundaryKnown(bool v) {
        boundary_known = v;
    }

    // If the boundary nodes and faces are knowm, return the value = 1. otherwise 0

    bool isBoundaryKnown() const {
        return boundary_known;
    }

    // A mesh is simple, when each edge is shared by at the most two faces which is
    // topological simple. A geometrically simple mesh will not have any crossing,
    // or overlapping edges, but this has not been checked right now.
    bool isSimple();

    // A mesh is consistently oriented when an edge shared by two faces, traversed in
    // opposite direction. Such a mesh will have proper normals ( either all inside or
    // or all outside ).
    bool is_consistently_oriented();

    //  Make the mesh consistent.
    void make_consistently_oriented();

    //  How many strongly connected components this mesh has ?
    int getNumOfConnectedComponents();

    //  Get #of boundary nodes, edges, and faces information.
    size_t getBoundarySize(int d) const;

    // Appends a node in the mesh.

    void addNode(PNode v) {
        nodes.push_back(v);
        v->setStatus(MeshEntity::ACTIVE);
    }

    // Appends a bulk of nodes in the mesh.

    void addNodes(const NodeSequence &vnodes) {
        for (size_t i = 0; i < vnodes.size(); i++)
            addNode(vnodes[i]);
    }

    void addNodes(NodeList &nlist) {
        while (!nlist.empty()) {
            PNode anode = nlist.front();
            nlist.pop_front();
            addNode(anode);
        }
    }

    int remove(PNode v);
    int deactivate(PNode v);
    int reactivate(PNode v);

    // Get the node at the specified position in the container.

    const PNode &getNodeAt(size_t id) const {
        assert(id < nodes.size());
        return nodes[id];
    }

    // Get All the nodes in the mesh after lazy prunning ( Garbage collection ).

    const NodeSequence &getNodes() {
        return nodes;
    }

    NodeList &get_nodes_list(); // Return nodes as STL List and  clear from the mesh.
    FaceList &get_faces_list(); // Returns faces as STL List and clear from the mesh

    //
    // Return irregular nodes in the domain. 
    // regular_count = 6   :  For triangle mesh;
    // regular_count = 4   :  For Quad mesh;
    // where         = 0   :  From the inner nodes only
    // where         = 1   :  From the boundary nodes only
    //
    NodeSequence get_irregular_nodes(int regular_count, int where = 0);

    const PEdge &getEdgeAt(size_t id) const {
        assert(id < edges.size());
        return edges[id];
    }

    EdgeSequence get_sharp_edges(double angle); //Specify the angle in degree.

    void addFeatureEdge(PEdge e) {
        feature_edges.push_back(e);
        //
        // Feature edge vertices are always constrained, by design. some codes
        // will break if this is modified in future.
        //
        Vertex *v;

        v = e->getNodeAt(0);
        v->setConstrainedMark(1);

        v = e->getNodeAt(1);
        v->setConstrainedMark(1);
    }

    bool hasFeatureEdge(const Edge &query_edge) const {
        for (size_t i = 0; i < feature_edges.size(); i++)
            if (feature_edges[i]->isSameAs(query_edge)) return 1;
        return 0;
    }

    // Add a face in the mesh. No duplication is checked..
    int addFace(PFace f);
    int remove(PFace f);

    int deactivate(PFace f);
    int reactivate(PFace f);

    // Add bulk of faces in the mesh. No duplication is checked..

    void addFaces(FaceSequence &vfaces) {
        for (size_t i = 0; i < vfaces.size(); i++)
            addFace(vfaces[i]);
    }

    void addFaces(FaceList &flist) {
        while (!flist.empty()) {
            PFace f = flist.front();
            flist.pop_front();
            addFace(f);
        }
    }

    // Get the face at the specified position in the container.

    const PFace &getFaceAt(size_t id) const {
        assert(id < faces.size());
        return faces[id];
    }

    bool contains(const Vertex *v) const {
        if (find(nodes.begin(), nodes.end(), v) == nodes.end()) return 0;
        return 1;
    }

    bool contains(const Face *f) const {
        if (find(faces.begin(), faces.end(), f) == faces.end()) return 0;
        return 1;
    }

    // Get rid of entities which are marked "removed" from the mesh.
    // a la Lazy garbage collection.
    //
    void prune();

    // Check if the lazy garbage collection is performed..
    bool isPruned() const;

    // Renumber mesh entities starting from index = 0
    void enumerate(int etype);

    // Search the boundary of the mesh (nodes, edges, and faces).
    int search_boundary();

    // Build entity-entity relations.

    int build_relations(int src, int dst, bool rebuild = 0) {
        if (src == 0 && dst == 0) return build_relations00(rebuild);
        if (src == 0 && dst == 2) return build_relations02(rebuild);
        return 1;
    }

    void normalize();

    // clean specified entity-entity relations.
    void clear_relations(int src, int dst);

    // Return all the edges of the primal mesh ...

    EdgeSequence &getEdges() {
        if (edges.empty()) build_edges();
        return edges;
    }

    // Return all the edges of dual mesh...
    EdgeSequence getDualEdges() const;

    // Returm all the matching edges of the dual graph...
    EdgeSequence getMatchingDuals() const;

    // Filter out face.
    FaceSequence filter(int facetype) const;

    // Save the mesh and all its attributes ( in the simple format ).
    int saveAs(const string &s);

    // Empty every thing in the Mesh, but don't delete the objects.
    void emptyAll();

    // Empty every thing in the Mesh, and also deallocate all the objects.
    void deleteAll();

    // Reverse the connection of all the faces in the mesh. This will be flip
    // the normal of each face.

    void reverse() {
        for (size_t i = 0; i < faces.size(); i++)
            faces[i]->reverse();
        ;
    }

    // Collect topological information. Ideally an internal vertex has 4 faces.
    vector<int> get_topological_statistics(int entity = 0, bool sorted = 1);

    // Collect nodes and faces in Depth First Sequence ...
    NodeSequence get_depth_first_ordered_nodes(Vertex *v = NULL, MeshFilter *mf = NULL);
    FaceSequence get_depth_first_ordered_faces(Face *f);

    // Collect nodes and faces in Breadth First Sequence ...
    NodeSequence get_breadth_first_ordered_nodes(Vertex *f = NULL, MeshFilter *mf = NULL);
    FaceSequence get_breadth_first_ordered_faces(Face *f);
    //
    // Creates waves of nodes/faces starting from the boundary. Each vertex/face
    // is assigned layerID, denoting the position in the wave fronts.
    // Used in smoothing the nodes in Advancing front style.
    //
    int setWavefront(int ofwhat);
    size_t setNodeWavefront(int layerid);
    size_t setFaceWavefront(int layerid);
    int verify_front_ordering(int ofwhat);

    NodeSet get_nearest_neighbors(const Point3D &p, int n = 1);
    NodeSet get_topological_neighbors(const NodeSequence &n, int k = 1);

    int remove_unreferenced_nodes();

    //  int check_unused_objects();

    //
    // Connect a strip of faces. Termination occurs when the face is reached
    // (1) At the boundary (2) Return back to the starting face.
    // There will be two strips per quadrilateral which are orthogonal to each
    // other. Using in the Ring Collapse Simplification
    //
    void get_quad_strips(Face *rootface, FaceSequence &strip1,
            FaceSequence &strip2);

    //
    void set_strip_markers();

    // Collect all the boundary nodes in the mesh..
    NodeSequence get_bound_nodes();

    //
    // Collect all the boundary faces in the mesh. The boundary faces could be
    // defined in two ways, bound_what flag is used for that purpose.
    // bound_what  = 0   At least one vertex is on the boundary
    //             = 1   At least one edge is on the boundary.
    //
    FaceSequence get_bound_faces(int bound_what = 1);

    //
    // Get the Aspect ratio of each face in the mesh. Aspect ratio is defined
    // as min_edge_length/max_edge_length..
    //
    vector<double> getAspectRatio(bool sorted = 1);

    //
    // Get unsigned surface area of the mesh. The calculation is for both
    // 2D and 3D surface elements...
    //
    double getSurfaceArea();

    // Check the Convexity..
    int check_convexity();

    // Collect Quadlity Statistics
    int get_quality_statistics(const string &s);

    // Get the nodes coordinates as a stream ( Used for interface with other
    // software. example Mesquite )
    const vector<double> getCoordsArray();

    // Get the node connectivity as a stream ( Used for interface with other
    // software. example Mesquite )
    const vector<size_t> getNodesArray();

    // Reset the coordinates of nodes, ( Generally it comes from optimization
    // modules. The array size must match in order to reflect changes.
    int setCoordsArray(const vector<double> &v);

    // Check if there are duplicate faces in the mesh;

    int hasDuplicates(int what);

    void getMinMaxVertexDegrees(int &mind, int &maxd);
    void getMinMaxFaceArea(double &mind, double &maxd);

    void setFacesNormal();

    int refine_quads14();
    int refine_quads15();

    int refine_quad14(Face *f);
    int refine_quad15(Face *f);

    void collect_garbage();
    int search_quad_patches();

private:
    volatile char adjTable[4][4];
    void build_edges();
    size_t count_edges();
    string meshname;

    double normalize_factor;
    double getLength(int dir) const;

    bool hasdual;
    bool boundary_known;

    // Contains all the mesh entities.
    NodeSequence nodes;
    EdgeSequence edges;
    FaceSequence faces;

    // Contains the nodes, edges and faces as STL list instead of vector.
    // This is dynamic structure to allow efficient refinement purposes.
    //
    NodeList nodelist;
    EdgeList edgelist;
    FaceList facelist;

    // Contains selected mesh entities.
    NodeSequence feature_nodes;
    EdgeSequence feature_edges;
    FaceSequence feature_faces;

    // Lazy garbage collection containers.
    NodeList garbageNodes;
    FaceList garbageFaces;

    // Build Topological relations...
    int build_relations00(bool rebuild = 0);
    int build_relations02(bool rebuild = 0);

    // Build wavefronts ...
    int setNodeWavefront();
    int setFaceWavefront();
};

int mesh_unit_tests();

struct MeshImporter {
    static const int SIMPLE_FORMAT = 0;
    static const int OFF_FORMAT = 1;
    static const int OBJ_FORMAT = 2;
    static const int VTK_FORMAT = 3;
    static const int TRIANGLE_FORMAT = 4;

    Mesh * load(const string &fname, Mesh *m = NULL) {
        mesh = m;
        if (mesh == NULL) mesh = new Mesh;

        int err = 1;
        if (fname.rfind(".vtk") != string::npos) {
            err = vtk_file(fname);
        }

        if (fname.rfind(".off") != string::npos) {
            err = off_file(fname);
        }

        if (fname.rfind(".dat") != string::npos) {
            err = simple_file(fname);
        }

        if (err) {
            err = triangle_file(fname);
        }

        if (!err) return mesh;

        return NULL;
    }

private:
    Mesh *mesh;
    std::map<int, int> global2local;
    void readTriNodes(const string & s);
    void readTriEdges(const string & s);
    void readTriFaces(const string & s);

    // Different format files ...
    int vtk_file(const string & s);
    int off_file(const string & s);
    int simple_file(const string & s);
    int triangle_file(const string & s);
};

struct MeshExporter {
    static const int SIMPLE_FORMAT = 0;
    static const int OFF_FORMAT = 1;
    static const int OBJ_FORMAT = 2;
    static const int VTK_FORMAT = 3;
    static const int TRIANGLE_FORMAT = 4;

    int saveAs(Mesh *m, const string & fname) {

        int err;
        if (fname.rfind(".vtk") != string::npos) {
            err = vtk_file(m, fname);
        }

        if (fname.rfind(".off") != string::npos) {
            err = off_file(m, fname);
        }

        if (fname.rfind(".dat") != string::npos) {
            err = simple_file(m, fname);
        }

    }

private:
    int off_file(Mesh *m, const string & s);
    int simple_file(Mesh *m, const string & s);
    int vtk_file(Mesh *m, const string & s);
};

struct MeshOptimization {
    int shape_optimize(Mesh * m);
    int untangle(Mesh * m);
    int equalize_edge_length(Mesh * m);
    int bandwidth_reduction(Mesh * m);
};

struct Delaunay {
    Delaunay();

    Delaunay(Mesh * m) {
        mesh = m;
    }

    //
    // Check if the surface triangulation is Delaunay. Works only for
    // triangle mesh. We use Jonathan Shewchuk's predicates and all his
    // ideas. For 2D check with Circum-Circle and for 2D Equitorial Sphere.
    //
    bool isDelaunay();

    //
    // If the surface triangulation is not Delaunay, use edge-flips to
    // tranform the mesh into Delaunay.
    //
    int makeDelaunay();

private:
    Mesh *mesh;
};

////////////////////////////////////////////////////////////////////////////////

struct SurfPatch {
public:

    void build() {
        search_boundary();
    }

    int count_boundary_segments() {
        return corners.size();
    }

    bool isClosed();
    bool isSimple();

    NodeSequence get_boundary_segment_nodes(int i) {
        int nsize = corners.size();
        int ib = cornerPos[i];
        int ie = cornerPos[i + 1];
        return get_bound_nodes(bound_nodes[ib], bound_nodes[ie]);
    }

    FaceSet faces;
private:
    vector<int> segSize; // How many nodes on each segment.
    vector<int> cornerPos; // Positions of the corners in the bound_nodes.
    NodeSet corners; // Corners of the Blob
    NodeSet inner_nodes; // Inner nodes ( not on the boundary ) of the blob
    NodeSequence bound_nodes; // Boundary nodes
    vector<Edge> boundary; // boundary of the blob.
    NodeSequence get_bound_nodes(const Vertex *src, const Vertex * dst);
    int search_boundary();
    int getPosOf(const Vertex * v);
    void set_boundary_segments();
    int reorient_4_sided_loop();
    void start_boundary_loop_from(Vertex * v);
};

////////////////////////////////////////////////////////////////////////////////

inline void MeshEntity::addRelation0(Vertex *vertex) {
    if (!hasRelation0(vertex)) relations0.push_back(vertex);
    sort(relations0.begin(), relations0.end());
}

inline void MeshEntity::removeRelation0(const Vertex *vertex) {
    if (hasRelation0(vertex)) {
        NodeSequence::iterator it;
        it = remove(relations0.begin(), relations0.end(), vertex);
        relations0.erase(it, relations0.end());
    }
}

inline void MeshEntity::addRelation2(Face* face) {
    if (!hasRelation2(face)) relations2.push_back(face);
    sort(relations2.begin(), relations2.end());
}

inline void MeshEntity::removeRelation2(const Face *face) {
    if (hasRelation2(face)) {
        FaceSequence::iterator it;
        it = remove(relations2.begin(), relations2.end(), face);
        relations2.erase(it, relations2.end());
    }
}

inline void MeshEntity::clearRelations(int t) {
    if (t == 0) relations0.clear();
    if (t == 2) relations2.clear();
}

inline const NodeSequence &MeshEntity::getRelations0() const {
    return relations0;
}

inline const FaceSequence &MeshEntity::getRelations2() const {
    return relations2;
}
///////////////////////////////////////////////////////////////////////////////

inline bool MeshEntity::hasRelation0(const Vertex* vertex) const {
    if (relations0.empty()) return 0;

    NodeSequence::const_iterator it;

    it = find(relations0.begin(), relations0.end(), vertex);
    if (it == relations0.end()) return 0;

    return 1;
}
///////////////////////////////////////////////////////////////////////////////

inline bool MeshEntity::hasRelation2(const Face* face) const {
    if (relations2.empty()) return 0;
    FaceSequence::const_iterator it;

    it = find(relations2.begin(), relations2.end(), face);
    if (it == relations2.end()) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

inline int Vertex::get_ideal_vertex_degree() const {
    if (!isBoundary())
        return 4;
    else
        return 3;

    if (getFeatureAngle() <= 90.0) return 1;
    if (getFeatureAngle() <= 180.0) return 2;
    if (getFeatureAngle() <= 270.0) return 3;

    return 4;
}

///////////////////////////////////////////////////////////////////////////////

inline void Face::start_from_concave_corner() {
    int pos = invertedAt();
    if (pos < 1) return;

    int nsize = getSize(0);
    NodeSequence rotated(nsize);
    for (int i = 0; i < nsize; i++)
        rotated[i] = connect[ (pos + i) % nsize ];
    connect = rotated;
}
///////////////////////////////////////////////////////////////////////////////

inline NodeSequence Mesh::generate_nodes(size_t n) {
    assert(n > 0);
    NodeSequence seq(n);
    for (size_t i = 0; i < n; i++)
        seq[i] = Vertex::newObject();
    return seq;
}
///////////////////////////////////////////////////////////////////////////////

inline FaceSequence Mesh::generate_faces(size_t n) {
    assert(n > 0);
    FaceSequence seq(n);
    for (size_t i = 0; i < n; i++)
        seq[i] = Face::newObject();
    return seq;
}
///////////////////////////////////////////////////////////////////////////////

inline int Mesh::remove(PNode v) {
    if (v == NULL) return 1;

    if (v->isRemoved()) return 0;

    if (getAdjTable(0, 2) == 0) {
        cout << "Warning: vertex-face relationship is absent: vertex not removed" << endl;
        return 1;
    }

    if (v->isActive()) {
        if (adjTable[0][2]) {
            FaceSequence vfaces = v->getRelations2();
            for (size_t i = 0; i < vfaces.size(); i++) {
                if (!vfaces[i]->isRemoved()) {
                    cout << "Warning: The vertex is currently used by face: vertex not removed" << endl;
                    return 1;
                }
            }
            v->clearRelations(2);
        }

        if (adjTable[0][0]) {
            NodeSequence vneighs = v->getRelations0();
            for (size_t i = 0; i < vneighs.size(); i++)
                vneighs[i]->removeRelation0(v);
            v->clearRelations(0);
        }
    }

    v->setStatus(MeshEntity::REMOVE);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::addFace(PFace f) {
    faces.push_back(f);

    adjTable[2][0] = 1;

    if (adjTable[0][2]) {
        for (int i = 0; i < f->getSize(0); i++) {
            Vertex *v = f->getNodeAt(i);
            v->addRelation2(f);
        }
    }

    if (adjTable[0][0]) {
        int nsize = f->getSize(0);
        for (int i = 0; i < nsize; i++) {
            Vertex *vi = f->getNodeAt((i + 0) % nsize);
            Vertex *vj = f->getNodeAt((i + 1) % nsize);
            vi->addRelation0(vj);
            vj->addRelation0(vi);
        }
    }
    f->setStatus(MeshEntity::ACTIVE);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::remove(PFace f) {
    if (f->isActive()) {
        if (adjTable[0][0]) {
            int nsize = f->getSize(0);
            for (int i = 0; i < nsize; i++) {
                Vertex *vi = f->getNodeAt((i + 0) % nsize);
                Vertex *vj = f->getNodeAt((i + 1) % nsize);
                FaceSequence vfaces = Mesh::getRelations112(vi, vj);
                if (vfaces.size() == 1) {
                    vi->removeRelation0(vj);
                    vj->removeRelation0(vi);
                }
            }
        }

        if (adjTable[0][2]) {
            for (int i = 0; i < f->getSize(0); i++) {
                Vertex *v = f->getNodeAt(i);
                v->removeRelation2(f);
            }
        }
    }

    f->setStatus(MeshEntity::REMOVE);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::deactivate(PNode vi) {
    assert(!vi->isRemoved());

    if (adjTable[0][0]) {
        NodeSequence vneighs = vi->getRelations0();
        for (size_t j = 0; j < vneighs.size(); j++) {
            Vertex *vj = vneighs[j];
            vi->removeRelation0(vj);
            vj->removeRelation0(vi);
        }
    }

    vi->setStatus(MeshEntity::INACTIVE);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

inline int Mesh::deactivate(PFace face) {
    // For the time being. but if you change the API, reimplement it.
    // The "remove" function can remove the elements immediately and
    // free the memory, but deactivate doesn't free it until garbage]
    // collection is called. Therefore, an element can be reused.
    // A deactivate element, looses all the information associated with
    // it.

    assert(!face->isRemoved());

    if (adjTable[0][2]) {
        for (int i = 0; i < face->getSize(0); i++) {
            Vertex *vertex = face->getNodeAt(i);
            vertex->removeRelation2(face);
        }
    }

    face->setStatus(MeshEntity::INACTIVE);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::reactivate(PNode vertex) {

    if (adjTable[0][0]) {
        FaceSequence vfaces = vertex->getRelations2();
        for (int i = 0; i < vfaces.size(); i++) {
            int nsize = vfaces[i]->getSize(0);
            for (int j = 0; j < nsize; j++) {
                Vertex *vi = vfaces[i]->getNodeAt((i + 0) % nsize);
                Vertex *vj = vfaces[i]->getNodeAt((i + 1) % nsize);
                vi->addRelation0(vj);
                vj->addRelation0(vi);
            }
        }
    }
    vertex->setStatus(MeshEntity::ACTIVE);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::reactivate(PFace face) {
    if (adjTable[0][0]) {
        int nsize = face->getSize(0);
        for (int i = 0; i < nsize; i++) {
            Vertex *vi = face->getNodeAt((i + 0) % nsize);
            Vertex *vj = face->getNodeAt((i + 1) % nsize);
            assert(!vi->isRemoved());
            assert(!vj->isRemoved());
            vi->addRelation0(vj);
            vj->addRelation0(vi);
        }
    }

    if (adjTable[0][2]) {
        for (int i = 0; i < face->getSize(0); i++) {
            Vertex *vi = face->getNodeAt(i);
            assert(!vi->isRemoved());
            vi->addRelation2(face);
        }
    }

    face->setStatus(MeshEntity::ACTIVE);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

Mesh *create_structured_mesh(double *origin, double *length, int *griddim, int spacedim);

NodeSequence linear_interpolation(Vertex *v0, Vertex *v1, int n);

int remesh_quad_loop(Mesh *m, NodeSequence &nodes,
        int nx, int ny,
        NodeSequence &newnodes, FaceSequence &newfaces,
        bool smooth = 0);

int remesh_quad_loop(Mesh *m,
        NodeSequence &anodes,
        NodeSequence &bnodes,
        NodeSequence &cnodes,
        NodeSequence &dnodes,
        NodeSequence &newnodes,
        FaceSequence &newfaces,
        bool smooth = 0);

int remesh_tri_loop(Mesh *m,
        NodeSequence &anodes,
        NodeSequence &bnodes,
        NodeSequence &cnodes,
        int *segments,
        NodeSequence &newnodes,
        FaceSequence &newfaces,
        bool smooth = 0);

int remesh_penta_loop(Mesh *mesh,
        NodeSequence &anodes,
        NodeSequence &bnodes,
        NodeSequence &cnodes,
        NodeSequence &dnodes,
        NodeSequence &enodes,
        int *segments,
        NodeSequence &newnodes,
        FaceSequence &newfaces,
        bool smooth = 0);

//
// QTrack is similar to "Motorcycle graph" proposed by Eppestein group.
// But somehow, I don't like this term.
//

struct QTrack {
    const static int END_AT_TERMINALS = 0;
    const static int END_AT_CROSSINGS = 1;

    Mesh *mesh;

    NodeSequence sequence;

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

    void advance(int endat);
private:
    int advance_single_step(int endat);
};

vector<QTrack> generate_quad_irregular_graph(Mesh *mesh);

///////////////////////////////////////////////////////////////////////////////

struct StructuredMesh2D {

    StructuredMesh2D() {
        nx = ny = 0;
    }
    int nx, ny;

    FaceSequence faces;
    NodeSequence nodes;

    NodeSet cornerNodes;

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

///////////////////////////////////////////////////////////////////////////////
// Set Tag Values for Visualzation and sometimes debugging
///////////////////////////////////////////////////////////////////////////////

void set_no_tags(Mesh *m);
void set_visit_tags(Mesh *m);
void set_layer_tag(Mesh *m);
void set_doublet_tag(Mesh *m);
void set_allboundnodes_tag(Mesh *m);
void set_singlet_tag(Mesh *m);
void set_boundary_tag(Mesh *m);
void set_constrained_tag(Mesh *m);
void set_feature_angle_tag(Mesh *m);
void set_convexity_tag(Mesh *m);
void set_bound1node_tag(Mesh *m);
void set_regular_node_tag(Mesh *m);
void set_inverted_tag(Mesh *m);
void set_large_area_tag(Mesh *m);
void set_tiny_area_tag(Mesh *m, double val = 1.0E-06);
void set_irregular_path_tag(Mesh *m, vector<QTrack> &qp);
void set_partition_tag(Mesh *m);

///////////////////////////////////////////////////////////////////////////////
// Graph Matching operations ....
///////////////////////////////////////////////////////////////////////////////

int quadrangulate(Mesh *mesh);

////////////////////////////////////////////////////////////////////////////////
//Helper functions ....
////////////////////////////////////////////////////////////////////////////////
Mesh *struct_tri_grid(int nx, int ny);
Mesh *struct_quad_grid(int nx, int ny);

////////////////////////////////////////////////////////////////////////////////
// Mesh Optimization ...
////////////////////////////////////////////////////////////////////////////////

struct LaplaceWeight {
    virtual double get(const Vertex *apex, const Vertex * neigh) = 0;
};

struct NoWeight : public LaplaceWeight {

    double get(const Vertex *apex, const Vertex * neigh) {
        return 1.0;
    }
};

struct LaplaceLengthWeight : public LaplaceWeight {

    double get(const Vertex *apex, const Vertex * neigh) {
        double len = Vertex::length(apex, neigh);
        return len;
    }
};

struct LaplaceAreaWeight : public LaplaceWeight {

    double get(const Vertex *apex, const Vertex * neighs) {
        FaceSequence vfaces = neighs->getRelations2();
        double area = 0.0;
        for (size_t i = 0; i < vfaces.size(); i++)
            area += vfaces[i]->getArea();
        return area;
    }
};

class LaplaceSmoothing {
public:

    LaplaceSmoothing(Mesh *m, int n = 100) {
        mesh = m;
        numIters = n;
        verbose = 0;
        lambda = 1.0;
        lapweight = NULL;
        method = 0;
        maxerror = 0;
    }

    void setMethod(int m) {
        method = m;
    }

    void setNumIterations(int n) {
        numIters = n;
    }

    int execute();

    void setWeight(LaplaceWeight *w) {
        lapweight = w;
    }

    int localized_at(const NodeSequence &q);

    int convexify();

    double getError() {
        return maxerror;
    }

private:
    Mesh *mesh;
    int global_smoothing();

    struct LVertex {
        Vertex *apex;
        map<Vertex*, vector<Vertex*> > neighs;
    };

    vector<LVertex> lnodes;
    LaplaceWeight *lapweight;
    vector<double> backupCoordsArray;
    vector<double> facearea;
    int method;
    int verbose;
    int numIters;
    double lambda;
    double maxerror;
    Point3D displace(const Point3D &src, const Point3D &dst, double r);

    double update_vertex_position(Vertex *vertex);
    double update_vertex_position(Vertex *vertex, const Point3D &p);
};

int quad_concave_tests();

} // namespace Jaal

#endif

