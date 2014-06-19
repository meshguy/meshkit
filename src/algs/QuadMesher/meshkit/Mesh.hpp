#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <cassert>
#include <fstream>
#include <math.h>
#include <string>
#include <limits>

#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <algorithm>

#include <array>

#ifdef HAVE_MESQUITE
#include <Mesquite_all_headers.hpp>
#endif

#ifdef HAVE_VERDICT
#include <verdict.h>
#endif

#include "myany.hpp"

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

///////////////////////////////////////////////////////////////////////////////

using namespace std;

typedef std::pair<size_t, size_t> FacePair;

#ifdef USE_HASHMAP
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#endif

namespace Jaal {

class Vertex;
class Face;
class Edge;
class Mesh;

typedef Vertex* PNode;
typedef Edge* PEdge;
typedef Face* PFace;

typedef std::vector<PNode> NodeSequence;
typedef std::vector<PEdge> EdgeSequence;
typedef std::vector<PFace> FaceSequence;

typedef std::list<PNode> NodeList;
typedef std::list<PEdge> EdgeList;
typedef std::list<PFace> FaceList;

typedef std::set<PNode> NodeSet;
typedef std::set<PFace> FaceSet;

////////////////////////////////////////////////////////////////////////////////

template <class T>
int split_stl_vector(const std::vector<T> &a, size_t pos, std::vector<T> &b, std::vector<T> &c)
{
     size_t nSize = a.size();
     if (pos >= nSize ) return 1;

     b.resize(pos);

     size_t index = 0;
     for (size_t i = 0; i < pos; i++)
          b[index++] = a[i];

     c.resize(nSize - pos + 1);
     index = 0;

     for (size_t i = pos - 1; i < nSize; i++)
          c[index++] = a[i];

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

template <class T>
int split_stl_vector(const std::deque<T> &a, size_t pos, std::deque<T> &b, std::deque<T> &c)
{
     size_t nSize = a.size();

     if (pos >= nSize) return 1;

     b.resize(pos);

     size_t index = 0;
     for (size_t i = 0; i < pos; i++)
          b[index++] = a[i];

     c.resize(nSize - pos + 1);
     index = 0;
     for (size_t i = pos - 1; i < nSize; i++)
          c[index++] = a[i];

     return 0;
}

////////////////////////////////////////////////////////////////////////////////

struct Attribute {
    public:
        Attribute() {}
        Attribute(const std::string &s, const myany &p) {
            name = s;
            value = p;
        }

        template <typename T> Attribute(T n) { value = n; }

        std::string name;
        myany value;
};

template<class T>
struct VecAttribute {
     std::vector<T> values;
};

template<class T, int n>
struct ArrayAttribute {     
      std::array <T,n> values;
};

template<class T>
struct ListAttribute {
     std::list<T> values;
};

template<class T>
struct SetAttribute {
     std::set<T> values;
};

///////////////////////////////////////////////////////////////////////////////
struct RelationRep {
     void clearRelations(int t);

     void addRelation(Vertex *vertex);
     void addRelation(Edge *edge);
     void addRelation(Face* face);

     int  getNumRelations( int e ) const {
          if( e == 0) return relations0.size();
          if( e == 1) return relations1.size();
          if( e == 2) return relations2.size();
          return 0;
     }

     bool hasRelation(const Vertex* vertex) const;
     bool hasRelation(const Edge* edge) const;
     bool hasRelation(const Face* face) const;

     void removeRelation(const Vertex *vertex);
     void removeRelation(const Edge *edge);
     void removeRelation(const Face *face);

     int getRelations( NodeSequence &seq, bool cyclic_ordered = 0) const {
          seq.clear();
          size_t nSize = relations0.size();
          if( nSize == 0) return 1;
          seq.resize( nSize );
          for( size_t i = 0; i< nSize; i++)
               seq[i] = relations0[i];
          return 0;
     }

     int getRelations( EdgeSequence &seq, bool cyclic_ordered = 0) const {
          seq.clear();
          size_t nSize = relations1.size();
          if( nSize == 0) return 1;
          seq.resize( nSize );
          for( size_t i = 0; i< nSize; i++)
               seq[i] = relations1[i];
          return 0;
     }


     int getRelations( FaceSequence &seq, bool cyclic_ordered = 0) const {
          seq.clear();
          size_t nSize = relations2.size();
          if( nSize == 0) return 1;
          seq.resize( nSize );
          for( size_t i = 0; i< nSize; i++)
               seq[i] = relations2[i];
          return 0;
     }

     NodeSequence relations0; // vertex-vertex
     EdgeSequence relations1; // vertex-vertex
     FaceSequence relations2; // vertex-face
};

struct AttribRep {
     void removeAll() {
          for( size_t i = 0; i < attributes.size(); i++)
               delete attributes[i];
     }
     template<class T>
     int setAttribute(const string &s, const T &val) {
          int nAttribs = attributes.size();
          for( int i = 0; i < nAttribs; i++) {
               if( attributes[i]->name == s ) {
                    attributes[i]->value = val;
                    return 0;
               }
          }
          Attribute *a = new Attribute(s,val);
          attributes.push_back( a );
          return 0;
     }

     template<class T>
     int getAttribute(const string &s, T &val) const {
          int nAttribs = attributes.size();
          for( int i = 0; i < nAttribs; i++) {
               if( attributes[i]->name == s ) {
                    val = any_cast <T>( attributes[i]->value);
                    return 0;
               }
          }
          return 1;
     }

     int hasAttribute(const string &s) const {
          int nAttribs = attributes.size();
          for( int i = 0; i < nAttribs; i++)
               if( attributes[i]->name == s ) return 1;
          return 0;
     }

     void removeAttribute( const string &s) {
     }

     std::vector<Attribute*> attributes;
};

class MeshEntity {
public:

     typedef size_t idtype;
     static const int ACTIVE   = 0;
     static const int REMOVE   = 1;
     static const int INACTIVE = 2;

     MeshEntity() {
          visitMark = 0; // Default: Not yet visited.
          statusMark = ACTIVE; // Default: Active< not removable.
          constrained = 0; // Default: No Contrainted
          boundarymark = 0; // Default: Internal entity
          relationRep  = NULL;
          attribRep    = NULL;
     }

     ~MeshEntity() {
          if( attribRep ) attribRep->removeAll();
     }

     void setVisitMark(bool r) {
          visitMark = r;
     }

     bool isVisited() const {
          return visitMark;
     }

     int getStatus() const {
          return statusMark;
     }

     bool isRemoved() const {
          if (statusMark == REMOVE) return 1;
          return 0;
     }

     bool isActive() const {
          if (statusMark == ACTIVE) return 1;
          return 0;
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

     void setID(idtype id) {
          gid = id;
     }

     const idtype &getID() const {
          return gid;
     }

     template<class T>
     int setAttribute(const string &s, const T &val) {
          if( attribRep == NULL ) attribRep = new AttribRep;
          return attribRep->setAttribute(s,val);
     }

     template<class T>
     int getAttribute(const string &s, T &val) const {
          if( attribRep )
               return attribRep->getAttribute(s,val);
          return 1;
     }

     int hasAttribute(const string &s) const {
          if( attribRep )
               return attribRep->hasAttribute(s);
          return 0;
     }

     void  removeAttribute(const string &s) const {
          if( attribRep ) attribRep->removeAttribute(s);
     }

     void clearRelations(int e) {
          if( relationRep ) relationRep->clearRelations(e);
     }

     void addRelation(Vertex *vertex) {
          if( relationRep  == NULL ) relationRep = new RelationRep;
          relationRep->addRelation(vertex);
     }

     void addRelation(Edge  *edge) {
          if( relationRep  == NULL ) relationRep = new RelationRep;
          relationRep->addRelation(edge);
     }

     void addRelation(Face* face) {
          if( relationRep  == NULL ) relationRep = new RelationRep;
          relationRep->addRelation(face);
     }

     int  getNumRelations( int e ) const {
          if( relationRep )
               return relationRep->getNumRelations(e);
          return 0;
     }

     bool hasRelation(const Vertex* vertex) const {
          if( relationRep )
               return relationRep->hasRelation( vertex );
          return 0;
     }

     bool hasRelation(const Edge* edge) const {
          if( relationRep )
               return relationRep->hasRelation( edge );
          return 0;
     }

     bool hasRelation(const Face* face) const {
          if( relationRep )
               return relationRep->hasRelation( face );
          return 0;
     }

     void removeRelation(const Vertex *vertex) {
          if( relationRep ) relationRep->removeRelation( vertex );
     }
     void removeRelation(const Face *face) {
          if( relationRep ) relationRep->removeRelation( face );
     }

     int getRelations( NodeSequence &seq, bool cyclic_ordered = 0) const {
          seq.clear();
          if( relationRep )
               return relationRep->getRelations(seq,cyclic_ordered);
          return 1;
     }

     int getRelations( EdgeSequence &seq, bool cyclic_ordered = 0) const {
          seq.clear();
          if( relationRep )
               return relationRep->getRelations(seq,cyclic_ordered);
          return 1;
     }

     int getRelations( FaceSequence &seq, bool cyclic_ordered = 0) const {
          seq.clear();
          if( relationRep )
               return relationRep->getRelations(seq,cyclic_ordered);
          return 1;
     }

     void setStatus(int a) {
          statusMark = a;
     }

protected:
     idtype gid;
     int boundarymark;
     bool constrained;
     volatile bool visitMark;
     volatile short int statusMark;
     RelationRep *relationRep;
     AttribRep   *attribRep;
};

struct EntityRemovedPred {

     bool operator() (const MeshEntity * entity) const {
          if (entity) return entity->isRemoved();
          return 0;
     }
};

///////////////////////////////////////////////////////////////////////////////

struct BaseVertex : public MeshEntity {
};

class Vertex : public BaseVertex {
     static size_t global_id;
public:

     static PNode newObject();

     Vertex() {
          visitMark = 0;
          statusMark = MeshEntity::ACTIVE;
          boundarymark = 0;
     }

     static double length(const Vertex *v0, const Vertex *v1);
     static double length2(const Vertex *v0, const Vertex *v1);
     static void   mid_point(const Vertex *v0, const Vertex *v1, Point3D &p, double alpha = 0.5);
     static Vertex* mid_node(const Vertex *v0, const Vertex *v1, double alpha = 0.5);
     static double point_orient( const Point3D &p0, const Point3D &p1, const Point3D &qpoint);

     void setXYZCoords(const Point3D &p) {
          xyz[0] = p[0];
          xyz[1] = p[1];
          xyz[2] = p[2];
     }

     const Point3D &getXYZCoords() const {
          return xyz;
     }

     Vertex* getClone() const;

     double getSpanAngle() const;
     double getFeatureLength() const;
     int get_ideal_face_degree( int n ) const;

private:
     Point3D xyz;
};

////////////////////////////////////////////////////////////////////////////////

struct LowVertexDegreeCompare {
     bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
          size_t d1 = vertex1->getNumRelations(0);
          size_t d2 = vertex2->getNumRelations(0);
          return d1 < d2;
     }
};
////////////////////////////////////////////////////////////////////////////////

struct HighVertexDegreeCompare {
     bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
          size_t d1 = vertex1->getNumRelations(0);
          size_t d2 = vertex2->getNumRelations(0);
          return d1 > d2;
     }
};
////////////////////////////////////////////////////////////////////////////////

struct LowLayerCompare {

     bool operator() (const Vertex *vertex1, const Vertex * vertex2) const {
          int val1, val2;
          vertex1->getAttribute("Layer", val1);
          vertex2->getAttribute("Layer", val2);
          return val1 < val2;
     }
};

///////////////////////////////////////////////////////////////////////////////

struct BaseEdge : public MeshEntity {
};

class Edge : public BaseEdge {
public:
     static PEdge newObject();

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

     PNode connect[2];
};

///////////////////////////////////////////////////////////////////////////////

struct BaseFace : public MeshEntity {
};

class Face : public BaseFace {
public:
     static const int POLYGON = 0;
     static const int TRIANGLE = 3;
     static const int QUADRILATERAL = 4;

     static PFace newObject();

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

     static PNode opposite_node(const PFace quad, const PNode n1);
     static PNode opposite_node(const PFace tri, PNode n1, PNode n2);
     static void opposite_nodes(const PFace quad, PNode n1, PNode n2,
                                PNode &n3, PNode &n4);

     static int check_on_boundary(const PFace tri);

     Face() {
          statusMark = 0;
          boundarymark = 0;
          visitMark = 0;
     }

     Face( const Vertex *v0, const Vertex *v1, const Vertex *v2) {
          NodeSequence seq(3);
          seq[0] = const_cast<Vertex*>(v0);
          seq[1] = const_cast<Vertex*>(v1);
          seq[2] = const_cast<Vertex*>(v2);
          setNodes(seq);
     }

     Face( const Vertex *v0, const Vertex *v1, const Vertex *v2, const Vertex *v3) {
          NodeSequence seq(4);
          seq[0] = const_cast<Vertex*>(v0);
          seq[1] = const_cast<Vertex*>(v1);
          seq[2] = const_cast<Vertex*>(v2);
          seq[3] = const_cast<Vertex*>(v3);
          setNodes(seq);
     }

     Vertex* getHashNode() const {
          return *min_element(connect.begin(), connect.end());
     }

     int getType() const {
          if (connect.size() == 3) return Face::TRIANGLE;
          if (connect.size() == 4) return Face::QUADRILATERAL;
          return Face::POLYGON;
     }

     int concaveAt() const;
     bool isSimple() const;

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

     int getPosOf(const Vertex *vertex) const {
          int nSize = connect.size();
          for (int i = 0; i < nSize; i++)
               if (connect[i] == vertex) return i;

          return -1;
     }

     void reverse() {
          std::reverse(connect.begin(), connect.end());
     }

     int getOrientation(const Vertex *ev0, const Vertex *ev1) const {
          int  nSize = connect.size();

          for (int i = 0; i < nSize; i++) {
               Vertex *v0 = connect[(i + 0) % nSize];
               Vertex *v1 = connect[(i + 1) % nSize];
               if (v0 == ev0 && v1 == ev1) return +1;
               if (v0 == ev1 && v1 == ev0) return -1;
          }
          return 0;
     }

     int replaceNode(Vertex *oldvertex, Vertex *newvertex) {

          int  nSize = connect.size();
          for (int i = 0; i < nSize; i++) {
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

     int setNodes(const NodeSequence &v) {
          int nSize = v.size();
          bool err = 0;
          for (int i = 0; i < nSize; i++) {
               for (int j = i+1; j < nSize; j++)
                    if(v[i] == v[j]) err = 1;
          }
          if( !err ) {
               connect = v;
          }
          return err;
     }

     const NodeSequence &getNodes() const {
          return connect;
     }

     const PNode &getNodeAt(size_t id) const {
          return connect[id%connect.size()];
     }

     void getAvgPos( Point3D &p) const;

     bool hasBoundaryNode() const {
          int nSize = connect.size();
          for (int i = 0; i < nSize; i++)
               if (connect[i]->isBoundary()) return 1;
          return 0;
     }

     bool has_all_bound_nodes() const;
     bool has_boundary_edge() const;

     int getRelations02( FaceSequence &seq);
     int getRelations12( FaceSequence &seq);
     int getRelations( NodeSequence &seq);

     bool isConvex() const {
          if (connect.size() <= 3) return 1;

          if (connect.size() == 4)
               return is_convex_quad(connect[0]->getXYZCoords(),
                                     connect[1]->getXYZCoords(),
                                     connect[2]->getXYZCoords(),
                                     connect[3]->getXYZCoords());

          cout << "Warning: General Polygon not supported " << endl;
          return 0;
     }

     double getAngleAt(const Vertex *v) const;

     int get_interior_angles(vector<double> &angles) const;

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

     int triangulate( vector<Face> &f) const;

     bool isValid() const;

     int refine_quad14(NodeSequence &newnodes, FaceSequence &newfaces);
     int refine_quad15(NodeSequence &newnodes, FaceSequence &newfaces);

     void start_from_concave_corner();

     Face* getClone() const {
          Face *newface = Face::newObject();
          newface->setNodes(connect);
          return newface;
     }

private:
     NodeSequence connect;

     int refine_convex_quad15(NodeSequence &newnodes,  FaceSequence &newfaces);
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
Vertex::newObject()
{
     Vertex *v = new Vertex;
     assert(v);
     v->setID(global_id);
     global_id++;
     return v;
}

///////////////////////////////////////////////////////////////////////////////

inline PNode Vertex::getClone() const
{
     Vertex *v = new Vertex;
     assert(v);
     v->setID(global_id);
     v->setXYZCoords(xyz);
     return v;
}

///////////////////////////////////////////////////////////////////////////////

inline PEdge Edge::newObject()
{
     Edge *e = new Edge;
     assert(e);
     return e;
}

///////////////////////////////////////////////////////////////////////////////

inline PFace Face::newObject()
{
     Face *f = new Face;
     assert(f);
     return f;
}

///////////////////////////////////////////////////////////////////////////////
inline double Face:: getAngleAt( const Vertex *v) const
{
     int pos = getPosOf(v);
     assert( pos >= 0);

     int nnodes = getSize(0);

     if( nnodes == 3 ) {
          const Point3D &p0 =  getNodeAt((pos+0)%nnodes)->getXYZCoords();
          const Point3D &p1 =  getNodeAt((pos+1)%nnodes)->getXYZCoords();
          const Point3D &p2 =  getNodeAt((pos+2)%nnodes)->getXYZCoords();
          return Math::getTriAngle(p0, p1, p2);
     }

     if( nnodes == 4 ) {
          const Point3D &p0 =  getNodeAt((pos+0)%nnodes)->getXYZCoords();
          const Point3D &p1 =  getNodeAt((pos+1)%nnodes)->getXYZCoords();
          const Point3D &p2 =  getNodeAt((pos+2)%nnodes)->getXYZCoords();
          const Point3D &p3 =  getNodeAt((pos+3)%nnodes)->getXYZCoords();
          double a1 = Math::getTriAngle(p0, p1, p2);
          double a2 = Math::getTriAngle(p0, p2, p3);
          return a1 + a2;
     }

     cout << "Error: Not implemented for general polygons " << endl;
     exit(0);
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

inline double Vertex :: getSpanAngle() const
{
     FaceSequence vfaces;
     this->getRelations( vfaces );
     int nSize = vfaces.size();
     double sum = 0.0;
     for( int i = 0; i <  nSize; i++) {
          Face *face = vfaces[i];
          sum += face->getAngleAt( this );
     }
     return sum;
}
////////////////////////////////////////////////////////////////////////////////

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
     }
    virtual ~MeshFilter() {
     }
        ;
};

class Mesh {
public:
     static NodeSequence generate_nodes(size_t n);
     static FaceSequence generate_faces(size_t n);

     static int getRelations112(const PNode v0, const PNode v1, FaceSequence &seq);
     static int getRelations102(const PNode v0, const PNode v1, FaceSequence &seq);
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
          collect_garbage();
     }

     void  objects_from_pool( size_t n, NodeSequence &objects);
     void  objects_from_pool( size_t n, FaceSequence &objects);

     void setName(const string &s) {
          meshname = s;
     }

     string getName() const {
          return meshname;
     }

     void reserve(size_t nSize, int entity) {
          if (entity == 0) nodes.reserve(nSize);
          if (entity == 1) edges.reserve(nSize);
          if (entity == 2) faces.reserve(nSize);
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

     int getNumComponents(bool stop_at_interface = 0);
     Mesh *getComponent(int id);

     int readFromFile(const string &f);
     bool getAdjTable(int i, int j) const {
          return adjTable[i][j];
     }
     Mesh * load(const string &fname, Mesh *m);
     Mesh *getPartMesh( int p);
     int  getNumPartitions(int e );

     int getEulerCharacteristic();

     Vertex* nearest_neighbour(const Vertex *v, double &d);

     size_t getSize(int d) {
          if (d == 0) return nodes.size();
          if (d == 1) return count_edges();
          if (d == 2) return faces.size();
          return 0;
     }

     size_t getCapacity(int d) {
          if (d == 0) return nodes.capacity();
          if (d == 2) return faces.capacity();
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
     int getNumConnectedComponents();

     //  Get #of boundary nodes, edges, and faces information.
     size_t getBoundarySize(int d) ;

     // Appends a node in the mesh.

     void addNode(PNode v) {
          nodes.push_back(v);
          v->setStatus(MeshEntity::ACTIVE);
     }

     // Appends a bulk of nodes in the mesh.

     void addNodes(const NodeSequence &vnodes) {
          size_t nSize = vnodes.size();
          for (size_t i = 0; i < nSize; i++)
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

     int getNodes(NodeSequence &seq) const {
          size_t nSize = nodes.size();

          seq.clear();
          seq.reserve( nSize );
          for( size_t i = 0; i <  nSize; i++) {
               Vertex *v = getNodeAt(i);
               if( v->isActive() ) seq.push_back(v);
          }
          return 0;
     }

     int getNodes(NodeList &l) const;
     int getFaces(FaceList &l) const;

     //
     // Return irregular nodes in the domain.
     // regular_count = 6   :  For triangle mesh;
     // regular_count = 4   :  For Quad mesh;
     // where         = 0   :  From the inner nodes only
     // where         = 1   :  From the boundary nodes only
     //
     int get_irregular_nodes(NodeSequence &seq, int regular_count, int where = 0);

     const PEdge &getEdgeAt(size_t id) const {
          assert(id < edges.size());
          if( edges[id]->getNodeAt(0)->isRemoved() || edges[id]->getNodeAt(1)->isRemoved()  )
               edges[id]->setStatus( MeshEntity::REMOVE );
          return edges[id];
     }

     int get_sharp_edges(EdgeSequence &e, double angle); //Specify the angle in degree.

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
          size_t nSize = feature_edges.size();
          for (size_t i = 0; i < nSize; i++)
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
          size_t nSize = vfaces.size();
          for (size_t i = 0; i < nSize; i++)
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
          if (src == 0 && dst == 0) return build_relations00( rebuild );
          if (src == 0 && dst == 1) return build_relations01( rebuild );
          if (src == 0 && dst == 2) return build_relations02( rebuild );
          if (src == 1 && dst == 2) return build_relations12( rebuild );
          return 1;
     }

     void normalize();

     // clean specified entity-entity relations.
     void clear_relations(int src, int dst);

     // Return all the edges of the primal mesh ...
     int getEdges( EdgeSequence &seq, int directed = 0) {
          if (edges.empty()) build_edges();
          seq = edges;
          return 0;
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
     void deleteNodes();
     void deleteEdges();
     void deleteFaces();
     void deleteAll();

     // Reverse the connection of all the faces in the mesh. This will be flip
     // the normal of each face.

     void reverse() {

          size_t  nSize = faces.size();
          for (size_t i = 0; i < nSize; i++)
               faces[i]->reverse();
     }

     // Collect topological information. Ideally an internal vertex has 4 faces.
     vector<int> get_topological_statistics(int entity = 0, bool sorted = 1);

     // Collect nodes and faces in Depth First Sequence ...
     int get_depth_first_ordered_nodes(NodeSequence &seq, Vertex *v = NULL, MeshFilter *mf = NULL);
     int get_depth_first_ordered_faces(FaceSequence &seq, Face *f = NULL);

     // Collect nodes and faces in Breadth First Sequence ...
     int get_breadth_first_ordered_nodes(NodeSequence &seq, Vertex *f = NULL, MeshFilter *mf = NULL);
     int get_breadth_first_ordered_faces( FaceSequence &seq, Face *f = NULL);
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
     int get_bound_nodes( NodeSequence &seq);

     //
     // Collect all the boundary faces in the mesh. The boundary faces could be
     // defined in two ways, bound_what flag is used for that purpose.
     // bound_what  = 0   At least one vertex is on the boundary
     //             = 1   At least one edge is on the boundary.
     //
     int get_bound_faces( FaceSequence &seq, int bound_what = 1);

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
     int getCoordsArray( vector<double> &a, vector<size_t> &l2g);

     // Get the node connectivity as a stream ( Used for interface with other
     // software. example Mesquite )
     int getNodesArray( vector<size_t> &a, vector<int> &topo);

     // Reset the coordinates of nodes, ( Generally it comes from optimization
     // modules. The array size must match in order to reflect changes.
     int setCoordsArray(const vector<double> &v, const vector<size_t> &l2g);

     // Check if there are duplicate faces in the mesh;

     int hasDuplicates(int what);

     void getMinMaxVertexDegrees(int &mind, int &maxd);
     void getMinMaxFaceArea(double &mind, double &maxd);

     void setNormals();

     //
     // Refine the mesh into "n" segments and return newnodes and newfaces
     // inserted into  the mesh structure. Old faces are reused, i.e.
     // they are deactivated and then reactivated with new connectivity...
     //
     // Requirement: At least two segments must be given...
     int refine_tri_edge( Vertex *v0, Vertex *v1, int n,
                          NodeSequence &newnodes, FaceSequence &newfaces);

     int collapse_tri_edge( Vertex *v0, Vertex *v1);

     int refine_quads14();
     int refine_quads15();

     int refine_quad14(Face *f);
     int refine_quad15(Face *f);

     void collect_garbage();
     int search_quad_patches();

     void build_edges();
private:
     volatile char adjTable[4][4];
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
     // Contains selected mesh entities.
     NodeSequence feature_nodes;
     EdgeSequence feature_edges;
     FaceSequence feature_faces;

     // Lazy garbage collection containers.
     NodeList garbageNodes;
     EdgeList garbageEdges;
     FaceList garbageFaces;

     // Build Topological relations...
     int build_relations00(bool rebuild = 0);
     int build_relations01(bool rebuild = 0);
     int build_relations02(bool rebuild = 0);
     int build_relations12(bool rebuild = 0);

     // Build wavefronts ...
     int setNodeWavefront();
     int setFaceWavefront();
};

//int mesh_unit_tests();

struct MeshImporter {
     static const int SIMPLE_FORMAT = 0;
     static const int OFF_FORMAT = 1;
     static const int OBJ_FORMAT = 2;
     static const int VTK_FORMAT = 3;
     static const int TRIANGLE_FORMAT = 4;
     static const int CUBIT_FORMAT = 5;

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

          if (fname.rfind(".cub") != string::npos) {
               err = cubit_file(fname);
          }

          if (err) {
               err = triangle_file(fname);
          }

          if( err ) return NULL;

          if( !mesh->is_consistently_oriented() )
               mesh->make_consistently_oriented();

          global2local.clear();

          return mesh;
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
     int xml_file(const string & s);
     int cubit_file(const string & s);
};

struct MeshExporter {

     static const int SIMPLE_FORMAT = 0;
     static const int OFF_FORMAT = 1;
     static const int OBJ_FORMAT = 2;
     static const int VTK_FORMAT = 3;
     static const int TRIANGLE_FORMAT = 4;
     static const int XML_FORMAT = 5;

     int saveAs(Mesh *m, const string & fname) {

          int err = 1;
          if (fname.rfind(".vtk") != string::npos) {
               err = vtk_file(m, fname);
          }

          if (fname.rfind(".off") != string::npos) {
               err = off_file(m, fname);
          }

          if (fname.rfind(".dat") != string::npos) {
               err = simple_file(m, fname);
          }

          if (fname.rfind(".xml") != string::npos) {
               err = xml_file(m, fname);
          }

          return err;
     }

private:
     int off_file(Mesh *m, const string & s);
     int simple_file(Mesh *m, const string & s);
     int vtk_file(Mesh *m, const string & s);
     int xml_file(Mesh *m, const string & s);
};

class MeshOptimization {
public:
     static const int STEEPEST_DESCENT   = 0;
     static const int QUASI_NEWTON       = 1;
     static const int TRUST_REGION       = 2;
     static const int FEASIBLE_NEWTON    = 3;
     static const int CONJUGATE_GRADIENT = 4;
     static const int LAPLACIAN          = 5;

     int shape_optimize(Mesh * m, int algo = QUASI_NEWTON ,  int numiter = 10);

     int bandwidth_reduction(Mesh * m);
private:
     Mesh *inmesh;
     int execute(Mesh * m);

     int algorithm;
     int numIter;

     vector<int>     vfixed;
     vector<size_t>  vNodes;
     vector<int>     etopo;
     vector<size_t>  l2g;
     vector<double>  vCoords;

#ifdef HAVE_MESQUITE
     Mesquite::QualityImprover *improver;
     Mesquite::ArrayMesh* jaal_to_mesquite(Mesh *m);
#endif
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

inline void RelationRep::addRelation(Vertex *vertex)
{
     if (!hasRelation(vertex)) relations0.push_back(vertex);
     sort(relations0.begin(), relations0.end());
}

inline void RelationRep::addRelation(Edge *edge)
{
     if (!hasRelation(edge)) relations1.push_back(edge);
     sort(relations1.begin(), relations1.end());
}

inline void RelationRep ::addRelation(Face* face)
{
     if (!hasRelation(face)) relations2.push_back(face);
     sort(relations2.begin(), relations2.end());
}

inline void RelationRep::removeRelation(const Vertex *vertex)
{
     if (hasRelation(vertex)) {
          NodeSequence::iterator it;
          it = remove(relations0.begin(), relations0.end(), vertex);
          relations0.erase(it, relations0.end());
     }
}

inline void RelationRep::removeRelation(const Edge *edge)
{
     if (hasRelation(edge)) {
          EdgeSequence::iterator it;
          it = remove(relations1.begin(), relations1.end(), edge);
          relations1.erase(it, relations1.end());
     }
}

inline void RelationRep ::removeRelation(const Face *face)
{
     if (hasRelation(face)) {
          FaceSequence::iterator it;
          it = remove(relations2.begin(), relations2.end(), face);
          relations2.erase(it, relations2.end());
     }
}

inline void RelationRep ::clearRelations(int t)
{
     if (t == 0) relations0.clear();
     if (t == 1) relations1.clear();
     if (t == 2) relations2.clear();
}

inline bool RelationRep ::hasRelation(const Vertex* vertex) const
{
     if (relations0.empty()) return 0;

     NodeSequence::const_iterator it;

     it = find(relations0.begin(), relations0.end(), vertex);
     if (it == relations0.end()) return 0;

     return 1;
}

inline bool RelationRep ::hasRelation(const Edge* edge) const
{
     if (relations1.empty()) return 0;

     EdgeSequence::const_iterator it;

     it = find(relations1.begin(), relations1.end(), edge);
     if (it == relations1.end()) return 0;

     return 1;
}
///////////////////////////////////////////////////////////////////////////////

inline bool RelationRep ::hasRelation(const Face* face) const
{
     if (relations2.empty()) return 0;
     FaceSequence::const_iterator it;

     it = find(relations2.begin(), relations2.end(), face);
     if (it == relations2.end()) return 0;

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

inline int Vertex::get_ideal_face_degree(int n) const
{
     if( n == 3 ) {
          if (!isBoundary()) return 6;

          if (getSpanAngle() <= 90.0 ) return 1;
          if (getSpanAngle() <= 220.0) return 3;
          if (getSpanAngle() <= 300.0) return 4;
          return 5;
     }

     if( n == 4 ) {
          if (!isBoundary()) return 4;

          if (getSpanAngle() <= 100 )  return 1;
          if (getSpanAngle() <= 220.0) return 2;
          if (getSpanAngle() <= 300.0) return 3;
          return 4;
     }

     cout << "Error: Ideal vertex degree only for Quad right now " << endl;
     exit(0);
}

///////////////////////////////////////////////////////////////////////////////

inline void Face::start_from_concave_corner()
{
     int pos = concaveAt();
     if (pos < 1) return;

     int nsize = getSize(0);
     NodeSequence rotated(nsize);
     for (int i = 0; i < nsize; i++)
          rotated[i] = connect[ (pos + i) % nsize ];
     connect = rotated;
}
///////////////////////////////////////////////////////////////////////////////

inline NodeSequence Mesh::generate_nodes(size_t n)
{
     assert(n > 0);
     NodeSequence seq(n);
     for (size_t i = 0; i < n; i++)
          seq[i] = Vertex::newObject();
     return seq;
}
///////////////////////////////////////////////////////////////////////////////

inline FaceSequence Mesh::generate_faces(size_t n)
{
     assert(n > 0);
     FaceSequence seq(n);
     for (size_t i = 0; i < n; i++)
          seq[i] = Face::newObject();
     return seq;
}
///////////////////////////////////////////////////////////////////////////////

inline int Mesh::remove(PNode v)
{
     if (v == NULL) return 1;

     if (v->isRemoved()) return 0;

     size_t nSize;

     FaceSequence vfaces;
     NodeSequence vneighs;

     if( v->isActive() ) {

          if (adjTable[0][2]) {
               v->getRelations( vfaces );
               nSize = vfaces.size();
               for (size_t i = 0; i < nSize; i++) {
                    if (!vfaces[i]->isRemoved()) {
                         cout << "Warning: The vertex is currently used by face: vertex not removed" << endl;
                         cout << "Debug : Vertex ID : " << v->getID() << endl;
                         cout << "Face:  ";
                         for( int j = 0; j < vfaces[i]->getSize(0); j++)
                              cout << vfaces[i]->getNodeAt(j)->getID() << " ";
                         cout << endl;
                         return 1;
                    }
               }
               v->clearRelations(2);
          }

          if (adjTable[0][0]) {
               v->getRelations( vneighs );
               nSize = vneighs.size();
               for (size_t i = 0; i < nSize; i++)
                    vneighs[i]->removeRelation(v);
               v->clearRelations(0);
          }
     }

     v->setStatus(MeshEntity::REMOVE);
     garbageNodes.push_back(v);

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::addFace(PFace f)
{
     faces.push_back(f);
     int nSize = f->getSize(0);

     adjTable[2][0] = 1;

     if (adjTable[0][2]) {
          for (int i = 0; i < nSize; i++) {
               Vertex *v = f->getNodeAt(i);
               v->addRelation(f);
          }
     }

     if (adjTable[1][0]) {
          for (int i = 0; i < nSize; i++) {
               Vertex *v0 = f->getNodeAt(i);
               Vertex *v1 = f->getNodeAt(i+1);
               Vertex *vi = min( v0, v1);
               Vertex *vj = max( v0, v1);
               if( !vi->hasRelation(vj) ) {
                    Edge *edge = new Edge(vi, vj);
                    edges.push_back(edge);
                    vi->addRelation(vj);
                    vj->addRelation(vi);
               }
          }
     }

     if (adjTable[0][0]) {
          for (int i = 0; i < nSize; i++) {
               Vertex *vi = f->getNodeAt((i + 0));
               Vertex *vj = f->getNodeAt((i + 1));
               vi->addRelation(vj);
               vj->addRelation(vi);
          }
     }

     f->setStatus(MeshEntity::ACTIVE);
     return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::remove(PFace f)
{
     if( f == NULL ) return 0;
     if( f->isRemoved() ) return 0;

     int nSize = f->getSize(0);

     FaceSequence vfaces;

     if( f->isActive() ) {
          if (adjTable[0][0]) {
               for (int i = 0; i < nSize; i++) {
                    Vertex *vi = f->getNodeAt(i + 0);
                    Vertex *vj = f->getNodeAt(i + 1);
                    Mesh::getRelations112(vi, vj, vfaces);
                    if (vfaces.size() == 1) {
                         vi->removeRelation(vj);
                         vj->removeRelation(vi);
                    }
               }
          }

          if (adjTable[0][2]) {
               for (int i = 0; i <  nSize; i++) {
                    Vertex *v = f->getNodeAt(i);
                    v->removeRelation(f);
               }
          }
     }
     f->setStatus(MeshEntity::REMOVE);
     garbageFaces.push_back(f);

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::deactivate(PNode vi)
{
     assert(!vi->isRemoved());

     NodeSequence vneighs;
     if (adjTable[0][0]) {
          vi->getRelations( vneighs );
          int nSize = vneighs.size();
          for (int j = 0; j <  nSize; j++) {
               Vertex *vj = vneighs[j];
               vj->removeRelation(vi);
          }
          vi->clearRelations(0);
     }

     if (adjTable[0][2]) vi->clearRelations(2);

     vi->setStatus(MeshEntity::INACTIVE);

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

inline int Mesh::deactivate(PFace face)
{
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
               vertex->removeRelation(face);
          }
     }

     face->setStatus(MeshEntity::INACTIVE);

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::reactivate(PNode vertex)
{
     FaceSequence vfaces;
     if (adjTable[0][0]) {
          vertex->getRelations( vfaces );
          int numFaces = vfaces.size();
          for (int i = 0; i < numFaces; i++) {
               int nsize = vfaces[i]->getSize(0);
               for (int j = 0; j < nsize; j++) {
                    Vertex *vi = vfaces[i]->getNodeAt(i + 0);
                    Vertex *vj = vfaces[i]->getNodeAt(i + 1);
                    vi->addRelation(vj);
                    vj->addRelation(vi);
               }
          }
     }
     vertex->setStatus(MeshEntity::ACTIVE);

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline int Mesh::reactivate(PFace face)
{
     int nsize = face->getSize(0);

     if (adjTable[0][0]) {
          for (int i = 0; i < nsize; i++) {
               Vertex *vi = face->getNodeAt(i + 0);
               Vertex *vj = face->getNodeAt(i + 1);
               assert(!vi->isRemoved());
               assert(!vj->isRemoved());
               vi->addRelation(vj);
               vj->addRelation(vi);
          }
     }

     if (adjTable[0][2]) {
          for (int i = 0; i < nsize; i++) {
               Vertex *vi = face->getNodeAt(i);
               assert(!vi->isRemoved());
               vi->addRelation(face);
          }
     }

     face->setStatus(MeshEntity::ACTIVE);
     return 0;
}

///////////////////////////////////////////////////////////////////////////////

Mesh *create_structured_mesh(double *origin, double *length, int *griddim, int spacedim);
Mesh* quad_to_tri4( Mesh *quadmesh, vector<Vertex*> &steiner);
Mesh* quad_to_tri2( Mesh *quadmesh );

void set_tfi_coords(int i, int j, int nx, int ny, vector<Vertex*> &qnodes);

void linear_interpolation(Mesh *m, Vertex *v0, Vertex *v1, int n, NodeSequence &seq);

void advancing_front_triangle_cleanup ( Mesh *mesh );
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

vector<QTrack> generate_quad_partitioning(Mesh *mesh);

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
void set_layer_tag(Mesh *m, const string &s = "Layer");
void set_inverted_tag(Mesh *m, const string &s = "Inverted" );
void set_boundary_tag(Mesh *m, const string &s = "Boundary");
void set_partition_tag(Mesh *m, const string &s = "Partition");
void set_convexity_tag(Mesh *m , const string &s = "Convexity" );
void set_regular_node_tag(Mesh *m, const string &s = "Regularity" );

void set_visit_tags(Mesh *m);
void set_allboundnodes_tag(Mesh *m);
void set_constrained_tag(Mesh *m);
void set_feature_angle_tag(Mesh *m);
void set_ideal_node_tag(Mesh *m, int elemtype);
void set_large_area_tag(Mesh *m);
void set_tiny_area_tag(Mesh *m, double val = 1.0E-06);
void set_irregular_path_tag(Mesh *m, vector<QTrack> &qp);

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

struct LaplaceNoWeight : public LaplaceWeight {

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

     FaceSequence  vfaces;
     double get(const Vertex *apex, const Vertex * neighs) {
          neighs->getRelations( vfaces );
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
          map<Vertex*, NodeSequence> neighs;
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

void plot_all_quad_quality_measures( Jaal::Mesh *mesh );
void plot_all_tri_quality_measures( Jaal::Mesh *mesh );

#endif

