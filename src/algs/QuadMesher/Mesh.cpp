#include <iomanip>

#include "Mesh.hpp"
#include "basic_math.hpp"

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

size_t Vertex::global_id = 0;

#ifdef USE_BOOST_LIBS
AttribKey MeshEntity::maxAttribID = 0;
std::map<string, AttribKey> MeshEntity::attribKeyMap;

///////////////////////////////////////////////////////////////////////////////

AttribKey
MeshEntity::getAttribKey(const string &s)
{
     map<string, AttribKey>::const_iterator it;

     it = attribKeyMap.find(s);

     if (it == attribKeyMap.end()) return 0;

     return it->second;
}
///////////////////////////////////////////////////////////////////////////////

bool
MeshEntity::hasAttribute(const string &s)
{
     return getAttribKey(s);
}

///////////////////////////////////////////////////////////////////////////////

AttribKey
MeshEntity::addAttribute(const string &s)
{
     AttribKey key = getAttribKey(s);
     if (key == 0) attribKeyMap[s] = ++maxAttribID;
     return attribKeyMap[s];
}

///////////////////////////////////////////////////////////////////////////////

void
MeshEntity::removeAttribute(const string &s)
{
     attribKeyMap.erase(s);
}
#endif

///////////////////////////////////////////////////////////////////////////////

Point3D
Vertex::mid_point(const Vertex *v0, const Vertex *v1, double alpha)
{
     const Point3D &p0 = v0->getXYZCoords();
     const Point3D &p1 = v1->getXYZCoords();

     Point3D pmid;
     pmid[0] = (1 - alpha) * p0[0] + alpha * p1[0];
     pmid[1] = (1 - alpha) * p0[1] + alpha * p1[1];
     pmid[2] = (1 - alpha) * p0[2] + alpha * p1[2];

     return pmid;
}

///////////////////////////////////////////////////////////////////////////////

double
Vertex::length(const Vertex *v0, const Vertex *v1)
{
     const Point3D &p0 = v0->getXYZCoords();
     const Point3D &p1 = v1->getXYZCoords();

     double dx = p0[0] - p1[0];
     double dy = p0[1] - p1[1];
     double dz = p0[2] - p1[2];

     return sqrt(dx * dx + dy * dy + dz * dz);
}

///////////////////////////////////////////////////////////////////////////////

double
Vertex::length2(const Vertex *v0, const Vertex *v1)
{
     const Point3D &p0 = v0->getXYZCoords();
     const Point3D &p1 = v1->getXYZCoords();

     double dx = p0[0] - p1[0];
     double dy = p0[1] - p1[1];
     double dz = p0[2] - p1[2];

     return dx * dx + dy * dy + dz * dz;
}

///////////////////////////////////////////////////////////////////////////////

int
Face::quad_tessalate(const NodeSequence &orgNodes, NodeSequence &rotatedNodes)
{
     double A013 = tri_area(orgNodes[0]->getXYZCoords(),
                            orgNodes[1]->getXYZCoords(),
                            orgNodes[3]->getXYZCoords());

     double A123 = tri_area(orgNodes[1]->getXYZCoords(),
                            orgNodes[2]->getXYZCoords(),
                            orgNodes[3]->getXYZCoords());

     double A012 = tri_area(orgNodes[0]->getXYZCoords(),
                            orgNodes[1]->getXYZCoords(),
                            orgNodes[2]->getXYZCoords());

     double A023 = tri_area(orgNodes[0]->getXYZCoords(),
                            orgNodes[2]->getXYZCoords(),
                            orgNodes[3]->getXYZCoords());

     rotatedNodes.resize(4);
     if (fabs(A013) + fabs(A123) < fabs(A012) + fabs(A023)) {
          rotatedNodes[0] = orgNodes[1];
          rotatedNodes[1] = orgNodes[2];
          rotatedNodes[2] = orgNodes[3];
          rotatedNodes[3] = orgNodes[0];
          return 1;
     }
     rotatedNodes[0] = orgNodes[0];
     rotatedNodes[1] = orgNodes[1];
     rotatedNodes[2] = orgNodes[2];
     rotatedNodes[3] = orgNodes[3];
     return 0;
}
///////////////////////////////////////////////////////////////////////////////

vector<Face>
Face::triangulate() const
{
     NodeSequence rotatedNodes(4), tconnect(3);
     vector<Face> trifaces;
     if (this->getSize(0) == 4) {
          quad_tessalate(connect, rotatedNodes);
          trifaces.resize(2);
          tconnect[0] = rotatedNodes[0];
          tconnect[1] = rotatedNodes[1];
          tconnect[2] = rotatedNodes[2];
          trifaces[0].setNodes(tconnect);

          tconnect[0] = rotatedNodes[0];
          tconnect[1] = rotatedNodes[2];
          tconnect[2] = rotatedNodes[3];
          trifaces[1].setNodes(tconnect);
     }
     return trifaces;
}
///////////////////////////////////////////////////////////////////////////////

vector<double>
Face::get_interior_angles() const
{
     NodeSequence rotatedNodes;
     map<Vertex*, double> mapangles;

     int nsize = connect.size();
     for (int i = 0; i < nsize; i++)
          mapangles[ connect[i] ] = 0.0;
     quad_tessalate(connect, rotatedNodes);

     Array<double, 3 > tangles;

     const Point3D &p0 = rotatedNodes[0]->getXYZCoords();
     const Point3D &p1 = rotatedNodes[1]->getXYZCoords();
     const Point3D &p2 = rotatedNodes[2]->getXYZCoords();
     const Point3D &p3 = rotatedNodes[3]->getXYZCoords();

     tangles = Math::getAngles(p0, p1, p2);
     mapangles[ rotatedNodes[0]] += tangles[0];
     mapangles[ rotatedNodes[1]] += tangles[1];
     mapangles[ rotatedNodes[2]] += tangles[2];

     tangles = Math::getAngles(p0, p2, p3);
     mapangles[ rotatedNodes[0]] += tangles[0];
     mapangles[ rotatedNodes[2]] += tangles[1];
     mapangles[ rotatedNodes[3]] += tangles[2];

     vector<double> vals(nsize);
     for (int i = 0; i < nsize; i++)
          vals[i] = mapangles[ connect[i] ];
     return vals;
}
///////////////////////////////////////////////////////////////////////////////

PNode
Face::opposite_node(const PFace tri, PNode n1, PNode n2)
{
     PNode tn0 = tri->getNodeAt(0);
     PNode tn1 = tri->getNodeAt(1);
     PNode tn2 = tri->getNodeAt(2);

     if (tn0 == n1 && tn1 == n2)
          return tn2;
     if (tn0 == n2 && tn1 == n1)
          return tn2;

     if (tn1 == n1 && tn2 == n2)
          return tn0;
     if (tn1 == n2 && tn2 == n1)
          return tn0;

     if (tn2 == n1 && tn0 == n2)
          return tn1;
     if (tn2 == n2 && tn0 == n1)
          return tn1;

     cout << " Warning: You should not come here " << endl;
     cout << " Face " << tn0 << " " << tn1 << " " << tn2 << endl;
     cout << " search for " << n1 << "  " << n2 << endl;
     exit(0);
}

///////////////////////////////////////////////////////////////////////////////

PNode
Face::opposite_node(const PFace face, PNode n1)
{
     int pos = face->getPosOf(n1);
     if (pos >= 0) {
          int size = face->getSize(0);
          if (size == 4)
               return face->getNodeAt((pos + 2) % 4);
     }
     return NULL;
}

///////////////////////////////////////////////////////////////////////////////

int Face::is_cyclic_quad(const Point3D &p0, const Point3D &p1, const Point3D &p2,
                         const Point3D &p3)
{
     //
     // http://en.wikipedia.org/wiki/Ptolemy%27s_theorem
     //

     double d02 = Math::length(p0, p2);
     double d13 = Math::length(p1, p3);
     double d01 = Math::length(p0, p1);
     double d12 = Math::length(p1, p2);
     double d23 = Math::length(p2, p3);
     double d30 = Math::length(p3, p0);

     double lval = d02*d13;
     double rval = d01 * d23 + d12*d30;

     if (fabs(lval - rval) < 1.0E-15) return 1;

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
Face::opposite_nodes(const PFace quad, PNode n1, PNode n2,
                     PNode &n3, PNode &n4)
{
     PNode qn0 = quad->getNodeAt(0);
     PNode qn1 = quad->getNodeAt(1);
     PNode qn2 = quad->getNodeAt(2);
     PNode qn3 = quad->getNodeAt(3);

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

PFace
Face::create_quad(const PFace t1, const PFace t2, int replace)
{
     NodeSequence connect;
     PNode commonnodes[3];

     connect = t1->getNodes();

     int index = 0;
     for (int i = 0; i < 3; i++) {
          if (t2->hasNode(connect[i]))
               commonnodes[index++] = connect[i];
     }

     if(index != 2) {
          cout << "Fatal Error:  Faces don't match for comoon edge " << endl;
          exit(0);
     }

     PNode ot1 = Face::opposite_node(t1, commonnodes[0], commonnodes[1]);
     PNode ot2 = Face::opposite_node(t2, commonnodes[0], commonnodes[1]);

     connect.resize(4);
     connect[0] = ot1;
     connect[1] = commonnodes[0];
     connect[2] = ot2;
     connect[3] = commonnodes[1];

     if (!replace) {
          Face *qface = new Face;
          qface->setNodes(connect);
          return qface;
     }

     t1->setNodes(connect);
     t2->setStatus(MeshEntity::REMOVE);
     return t1;
}

///////////////////////////////////////////////////////////////////////////////

int
Face::hexagon_2_quads(const NodeSequence &hexnodes, FaceSequence &quads, int offset)
{
#ifdef DEBUG
     if( hexnodes.size() != 6) {
          cout << "Fatal Error: Given polygon is not hex " << endl;
          exit(0);
     }
#endif

     PFace face1 = Face::newObject();
     PFace face2 = Face::newObject();

     NodeSequence connect(4);

     for (int i = 0; i < 3; i++) {
          connect[0] = hexnodes[ (i + offset + 0) % 6];
          connect[1] = hexnodes[ (i + offset + 1) % 6];
          connect[2] = hexnodes[ (i + offset + 2) % 6];
          connect[3] = hexnodes[ (i + offset + 3) % 6];
          face1->setNodes(connect);

          connect[0] = hexnodes[ (i + offset + 3) % 6];
          connect[1] = hexnodes[ (i + offset + 4) % 6];
          connect[2] = hexnodes[ (i + offset + 5) % 6];
          connect[3] = hexnodes[ (i + offset + 6) % 6];
          face2->setNodes(connect);
          if (face1->isConvex() && face2->isConvex()) {
               quads.resize(2);
               quads[0] = face1;
               quads[1] = face2;
               return 0;
          }
     }
     quads.clear();
     delete face1;
     delete face2;
     return 1;
}

///////////////////////////////////////////////////////////////////////////////

const Point3D
Face::getCentroid() const
{
     Point3D pc;

     pc[0] = 0.0;
     pc[1] = 0.0;
     pc[2] = 0.0;

     int nsize = connect.size();
     for (int inode = 0; inode < nsize; inode++) {
          Vertex *v = connect[inode];
          const Point3D &p3d = v->getXYZCoords();
          pc[0] += p3d[0];
          pc[1] += p3d[1];
          pc[2] += p3d[2];
     }

     double multby = 1.0/(double)nsize;
     pc[0] *=  multby;
     pc[1] *=  multby;
     pc[2] *=  multby;

     return pc;
}

////////////////////////////////////////////////////////////////////////////////

Vec3D
Face::normal(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
     Point3D p2p0 = Math::create_vector(p2, p0);
     Point3D p1p0 = Math::create_vector(p1, p0);
     Point3D normal = Math::cross_product(p2p0, p1p0);

     double mag = Math::magnitude(normal);
     normal[0] /= mag;
     normal[1] /= mag;
     normal[2] /= mag;

     return normal;
}

////////////////////////////////////////////////////////////////////////////////

Vec3D
Face::normal(const Vertex *v0, const Vertex *v1, const Vertex *v2)
{
     return Math::normal(v0->getXYZCoords(),
                         v1->getXYZCoords(),
                         v2->getXYZCoords());
}

////////////////////////////////////////////////////////////////////////////////

double
Face::tri_area(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
     ////////////////////////////////////////////////////////////////////////////////
     // Ref: http://mathworld.wolfram.com/HeronsFormula.html
     //
     // Heron's formular  A = sqrt(s(s-a)(s-b)(s-c)) is expensive because if require
     // three square roots for each a,b, and c.
     // Instead we will use alternate formula of "Heron" which avoids three
     // expensive square roots.
     //
     // Sorting is done to reduce the truncation errors. Very similar to Kahan's
     // Original idea.
     ////////////////////////////////////////////////////////////////////////////////

     double d[3];
     d[0] = Math::length2(p1, p2);
     d[1] = Math::length2(p2, p0);
     d[2] = Math::length2(p0, p1);

     std::sort(d, d + 3); // May be we should have optimized version than STL one

     double a2 = d[0];
     double b2 = d[1];
     double c2 = d[2];

     double area = 0.25 * sqrt(4 * a2 * b2 - (a2 + b2 - c2)*(a2 + b2 - c2));
     return area;
}

///////////////////////////////////////////////////////////////////////////////////

double
Face::quad_area(const Point3D &p0, const Point3D &p1,
                const Point3D &p2, const Point3D &p3)
{
     /////////////////////////////////////////////////////////////////////////////
     // For explanation of some amazing proofs and theorems, please refer to the
     // following site. This implementation is based on this article.
     //
     // http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#Quadrilaterals
     //
     // For Bretschneider's Formula: refer to
     // http://mathworld.wolfram.com/BretschneidersFormula.html
     // Given a general quadrilateral with sides of lengths a, b, c, and d, the area
     // is given by
     //              K = 1/4sqrt(4p^2q^2-(b^2+d^2-a^2-c^2)^2)
     //
     // It seems that the first method is better because it doesn't require
     // expensive 6 lenghts. ( 4 sides + 2 diagonal ). But a proper analysis needs
     // be done. But probably, the most amazing thing about the formula is that
     // it is valid for 3D quadrilateral and for both convex and concave.
     // But why it handles concavity, I am not quite sure.
     //
     // Chaman Singh Verma
     // 16th Feb 2010.
     //////////////////////////////////////////////////////////////////////////////

     Point3D v2v0 = Math::create_vector(p2, p0);
     Point3D v3v1 = Math::create_vector(p3, p1);
     Point3D d0d1 = Math::cross_product(v2v0, v3v1);

     double area = 0.5 * Math::magnitude(d0d1);
     return area;
}
////////////////////////////////////////////////////////////////////////////////

bool
Face::is_convex_quad(const Point3D &p0, const Point3D &p1,
                     const Point3D &p2, const Point3D &p3)
{
     double qarea = quad_area(p0, p1, p2, p3);

     double tarea1, tarea2;
     tarea1 = tri_area(p0, p1, p2);
     tarea2 = tri_area(p0, p2, p3);
     if (fabs(tarea1 + tarea2 - qarea) > 1.0E-10) return 0;

     tarea1 = tri_area(p0, p1, p3);
     tarea2 = tri_area(p1, p2, p3);
     if (fabs(tarea1 + tarea2 - qarea) > 1.0E-10) return 0;

     return 1;
}

////////////////////////////////////////////////////////////////////////////////

NodeSequence
Face::getRelations0()
{
     NodeSequence vneighs;
     set<Vertex*> vset;

     int nnodes = connect.size();

     for (int i = 0; i < nnodes; i++) {
          Vertex *vertex = connect[i];
          vneighs = vertex->getRelations0();
          int numneighs = vneighs.size();
          for (int j = 0; j < numneighs; j++)
               vset.insert(vneighs[j]);
     }

     for (int i = 0; i < nnodes; i++)
          vset.erase(connect[i]);

     NodeSequence vresult;
     if (!vset.empty()) {
          set<Vertex*>::const_iterator it;
          size_t index = 0;
          for (it = vset.begin(); it != vset.end(); ++it)
               vresult[index++] = *it;
     }
     return vresult;
}

////////////////////////////////////////////////////////////////////////////////

double
Face::getAspectRatio()
{
     int nSize = connect.size();

     double minlen = MAXDOUBLE;
     double maxlen = 0.0;

     for (int i = 0; i < nSize; i++) {
          Vertex *v0 = connect[(i + 0) % nSize];
          Vertex *v1 = connect[(i + 1) % nSize];
          double len2 = Vertex::length2(v0, v1);
          if (len2 > maxlen) maxlen = len2;
          if (len2 < minlen) minlen = len2;
     }

     return sqrt(minlen / maxlen);
}

////////////////////////////////////////////////////////////////////////////////

FaceSequence
Face::getRelations202()
{
     FaceSequence faceneighs;

     int nSize = connect.size();
     for (int i = 0; i < nSize; i++) {
          Vertex *v0 = connect[(i + 0) % nSize];
          const FaceSequence &vneighs = v0->getRelations2();
          int numneighs = vneighs.size();
          for (int  j = 0; j < numneighs; j++) {
               if (vneighs[j] != this) {
                    if (find(faceneighs.begin(), faceneighs.end(), vneighs[j])
                              == faceneighs.end())
                         faceneighs.push_back(vneighs[j]);
               }
          }
     }
     return faceneighs;
}
///////////////////////////////////////////////////////////////////////////////

FaceSequence
Face::getRelations212()
{
     FaceSequence faceneighs, edgeneighs;

     int nSize = connect.size();
     for (int i = 0; i < nSize; i++) {
          Vertex *v0 = connect[i + 0];
          Vertex *v1 = connect[i + 1];
          Mesh::getRelations112(v0, v1, edgeneighs);
          int numneighs = edgeneighs.size();
          for (int j = 0; j <  numneighs; j++) {
               if (edgeneighs[j] != this) {
                    if (find(faceneighs.begin(), faceneighs.end(), edgeneighs[j])
                              == faceneighs.end())
                         faceneighs.push_back(edgeneighs[j]);
               }
          }
     }
     return faceneighs;
}

///////////////////////////////////////////////////////////////////////////////

bool
Face::isValid() const
{
     if (isRemoved()) return 0;

     int nSize = connect.size();
     for (int i = 0; i < nSize; i++)
          if (connect[i]->isRemoved()) return 0;

     for (int i = 0; i < nSize; i++) {
          int ncount = 0;
          for (int j = 0; j < nSize; j++)
               if (connect[i] == connect[j]) ncount++;
          if (ncount != 1) return 0;
     }
     return 1;
}

///////////////////////////////////////////////////////////////////////////////

bool
Face::has_boundary_edge() const
{
     int nSize = connect.size();
     FaceSequence neighs;
     for (int i = 0; i < nSize; i++) {
          Vertex *v0 = connect[(i + 0) % nSize];
          Vertex *v1 = connect[(i + 1) % nSize];
          if (v0->isBoundary() && v1->isBoundary()) {
               Mesh::getRelations112(v0, v1, neighs);
               if (neighs.size() == 1) return 1;
          }
     }
     return 0;
}
/////////////////////////////////////////////////////////////////////////////////////

void
Face::bilinear_weights(double xi, double eta, vector<double> &weight)
{
     weight.resize(4);

#ifdef DEBUG
     assert(xi >= -1.0 && xi <= 1.0);
     assert(eta >= -1.0 && eta <= 1.0);
#endif

     double coeff = 1.0 / 4.0;

     weight[0] = coeff * (1.0 - xi)*(1.0 - eta);
     weight[1] = coeff * (1.0 + xi)*(1.0 - eta);
     weight[2] = coeff * (1.0 + xi)*(1.0 + eta);
     weight[3] = coeff * (1.0 - xi)*(1.0 + eta);
}

/////////////////////////////////////////////////////////////////////////////////////

double
Face::linear_interpolation(const vector<double> &x, const vector<double> &w)
{
     size_t numNodes = x.size();
     assert(w.size() == numNodes);

     double sum = 0.0;
     for (size_t i = 0; i < numNodes; i++) sum += x[i] * w[i];

     return sum;
}

/////////////////////////////////////////////////////////////////////////////////////

int
Face::concaveAt() const
{
     int nsize = getSize(0);
     Vertex *vertex;
     double x[5], y[5], z[5], triarea;

     Point3D xyz;

     for (int i = 0; i < nsize; i++) {
          vertex = getNodeAt(i);
          xyz = vertex->getXYZCoords();
          x[0] = xyz[0];
          y[0] = xyz[1];
          z[0] = xyz[2];

          vertex = getNodeAt((i + 1) % nsize);
          xyz = vertex->getXYZCoords();
          x[1] = xyz[0];
          y[1] = xyz[1];
          z[1] = xyz[2];

          vertex = getNodeAt((i + nsize - 1) % nsize);
          xyz = vertex->getXYZCoords();
          x[2] = xyz[0];
          y[2] = xyz[1];
          z[2] = xyz[2];

          triarea = PolygonArea3D(3, x, y, z);
          if (triarea < 0.0) return i;
     }
     return -1;
}

/////////////////////////////////////////////////////////////////////////////////////
double Vertex:: point_orient( const Point3D &p0, const Point3D &p1, const Point3D &qpoint)
{
     double x[5], y[5], z[5];

     x[0] = p0[0];
     x[1] = p1[0];
     x[2] = qpoint[0];

     y[0] = p0[1];
     y[1] = p1[1];
     y[2] = qpoint[1];

     z[0] = p0[2];
     z[1] = p1[2];
     z[2] = qpoint[2];

     return PolygonArea3D(3, x, y, z);
}
/////////////////////////////////////////////////////////////////////////////////////

bool
Face::isSimple() const
{
     int nnodes =  getSize(0);
     if( nnodes == 3 ) return 1;

     double d0, d1, d2, d3;
     if( nnodes == 4 ) {
          const Point3D &v0 = getNodeAt(0)->getXYZCoords();
          const Point3D &v1 = getNodeAt(1)->getXYZCoords();
          const Point3D &v2 = getNodeAt(2)->getXYZCoords();
          const Point3D &v3 = getNodeAt(3)->getXYZCoords();

          d0  = Vertex::point_orient( v0, v1, v2);
          d1  = Vertex::point_orient( v0, v1, v3);

          d2  = Vertex::point_orient( v2, v3, v0);
          d3  = Vertex::point_orient( v2, v3, v1);

          if( (d0*d1 < 0.0) && (d2*d3 < 0.0 ) ) return 0;

          d0  = Vertex::point_orient( v1, v2, v0);
          d1  = Vertex::point_orient( v1, v2, v3);

          d2  = Vertex::point_orient( v0, v3, v1);
          d3  = Vertex::point_orient( v0, v3, v2);

          if( (d0*d1 < 0.0) && (d2*d3 < 0.0 ) ) return 0;

          return 1;
     }

     cout << "Warning: General Polygons not supported yet " << endl;

     return 0;

}

////////////////////////////////////////////////////////////////////////////////

int
Face::refine_quad15(NodeSequence &newnodes, FaceSequence &newfaces)
{
     if (!isConvex())
          return refine_concave_quad15(newnodes, newfaces);

     return refine_convex_quad15(newnodes, newfaces);
}

////////////////////////////////////////////////////////////////////////////////

int
Face::refine_concave_quad15(NodeSequence &newnodes, FaceSequence &newfaces)
{
     newnodes.resize(4);
     newfaces.resize(5);

     NodeSequence rotatedNodes(4), tconnect(3);

     Point3D xyz;

     quad_tessalate(connect, rotatedNodes);
     xyz = Vertex::mid_point(rotatedNodes[0], rotatedNodes[2], 1.0 / 3.0);
     newnodes[0] = Vertex::newObject();
     newnodes[0]->setXYZCoords(xyz);

     xyz = Vertex::mid_point(rotatedNodes[0], rotatedNodes[2], 2.0 / 3.0);
     newnodes[2] = Vertex::newObject();
     newnodes[2]->setXYZCoords(xyz);

     Face triface;
     tconnect[0] = rotatedNodes[0];
     tconnect[1] = rotatedNodes[1];
     tconnect[2] = rotatedNodes[2];
     triface.setNodes(tconnect);
     xyz = triface.getCentroid();
     newnodes[1] = Vertex::newObject();
     newnodes[1]->setXYZCoords(xyz);

     tconnect[0] = rotatedNodes[0];
     tconnect[1] = rotatedNodes[2];
     tconnect[2] = rotatedNodes[3];
     triface.setNodes(tconnect);
     xyz = triface.getCentroid();
     newnodes[3] = Vertex::newObject();
     newnodes[3]->setXYZCoords(xyz);

     NodeSequence connect(4);

     connect[0] = rotatedNodes[0];
     connect[1] = rotatedNodes[1];
     connect[2] = newnodes[1];
     connect[3] = newnodes[0];
     newfaces[0] = Face::newObject();
     newfaces[0]->setNodes(connect);

     connect[0] = rotatedNodes[1];
     connect[1] = rotatedNodes[2];
     connect[2] = newnodes[2];
     connect[3] = newnodes[1];
     newfaces[1] = Face::newObject();
     newfaces[1]->setNodes(connect);

     connect[0] = rotatedNodes[2];
     connect[1] = rotatedNodes[3];
     connect[2] = newnodes[3];
     connect[3] = newnodes[2];
     newfaces[2] = Face::newObject();
     newfaces[2]->setNodes(connect);

     connect[0] = rotatedNodes[3];
     connect[1] = rotatedNodes[0];
     connect[2] = newnodes[0];
     connect[3] = newnodes[3];
     newfaces[3] = Face::newObject();
     newfaces[3]->setNodes(connect);

     connect[0] = newnodes[0];
     connect[1] = newnodes[1];
     connect[2] = newnodes[2];
     connect[3] = newnodes[3];
     newfaces[4] = Face::newObject();
     newfaces[4]->setNodes(connect);

     return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
Face::refine_convex_quad15(NodeSequence &newnodes, FaceSequence &newfaces)
{
     newnodes.resize(4);
     newfaces.resize(5);

     Point3D xyz;
     vector<double> weight;
     vector<double> xc(4), yc(4), zc(4);
     for (int i = 0; i < 4; i++) {
          Vertex *v = this->getNodeAt(i);
          xyz = v->getXYZCoords();
          xc[i] = xyz[0];
          yc[i] = xyz[1];
          zc[i] = xyz[2];
     }

     double dl = 2.0 / 3.0;

     bilinear_weights(-1.0 + dl, -1.0 + dl, weight);
     xyz[0] = linear_interpolation(xc, weight);
     xyz[1] = linear_interpolation(yc, weight);
     xyz[2] = linear_interpolation(zc, weight);
     newnodes[0] = Vertex::newObject();
     newnodes[0]->setXYZCoords(xyz);

     bilinear_weights(-1.0 + 2 * dl, -1.0 + dl, weight);
     xyz[0] = linear_interpolation(xc, weight);
     xyz[1] = linear_interpolation(yc, weight);
     xyz[2] = linear_interpolation(zc, weight);
     newnodes[1] = Vertex::newObject();
     newnodes[1]->setXYZCoords(xyz);

     bilinear_weights(-1.0 + 2 * dl, -1.0 + 2 * dl, weight);
     xyz[0] = linear_interpolation(xc, weight);
     xyz[1] = linear_interpolation(yc, weight);
     xyz[2] = linear_interpolation(zc, weight);
     newnodes[2] = Vertex::newObject();
     newnodes[2]->setXYZCoords(xyz);

     bilinear_weights(-1.0 + dl, -1.0 + 2.0 * dl, weight);
     xyz[0] = linear_interpolation(xc, weight);
     xyz[1] = linear_interpolation(yc, weight);
     xyz[2] = linear_interpolation(zc, weight);
     newnodes[3] = Vertex::newObject();
     newnodes[3]->setXYZCoords(xyz);

     NodeSequence connect(4), oldnodes;
     oldnodes = this->getNodes();

     connect[0] = oldnodes[0];
     connect[1] = oldnodes[1];
     connect[2] = newnodes[1];
     connect[3] = newnodes[0];
     newfaces[0] = Face::newObject();
     newfaces[0]->setNodes(connect);

     connect[0] = oldnodes[1];
     connect[1] = oldnodes[2];
     connect[2] = newnodes[2];
     connect[3] = newnodes[1];
     newfaces[1] = Face::newObject();
     newfaces[1]->setNodes(connect);

     connect[0] = oldnodes[2];
     connect[1] = oldnodes[3];
     connect[2] = newnodes[3];
     connect[3] = newnodes[2];
     newfaces[2] = Face::newObject();
     newfaces[2]->setNodes(connect);

     connect[0] = oldnodes[3];
     connect[1] = oldnodes[0];
     connect[2] = newnodes[0];
     connect[3] = newnodes[3];
     newfaces[3] = Face::newObject();
     newfaces[3]->setNodes(connect);

     connect[0] = newnodes[0];
     connect[1] = newnodes[1];
     connect[2] = newnodes[2];
     connect[3] = newnodes[3];
     newfaces[4] = Face::newObject();
     newfaces[4]->setNodes(connect);

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

Vertex *
Face::centroid(const Vertex *v0, const Vertex *v1, const Vertex *v2)
{
     Vertex *vc = Vertex::newObject();

     Point3D pc, xyz;
     pc[0] = 0.0;
     pc[1] = 0.0;
     pc[2] = 0.0;

     xyz = v0->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     xyz = v1->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     xyz = v2->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     pc[0] /= 3.0;
     pc[1] /= 3.0;
     pc[2] /= 3.0;

     vc->setXYZCoords(pc);
     return vc;
}

///////////////////////////////////////////////////////////////////////////////

Vertex *
Face::centroid(const Vertex *v0, const Vertex *v1, const Vertex *v2,
               const Vertex *v3, const Vertex *v4)
{
     Vertex *vc = Vertex::newObject();

     Point3D pc, xyz;
     pc[0] = 0.0;
     pc[1] = 0.0;
     pc[2] = 0.0;

     xyz = v0->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     xyz = v1->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     xyz = v2->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     xyz = v3->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     xyz = v4->getXYZCoords();
     pc[0] += xyz[0];
     pc[1] += xyz[1];
     pc[2] += xyz[2];

     pc[0] /= 5.0;
     pc[1] /= 5.0;
     pc[2] /= 5.0;

     vc->setXYZCoords(pc);
     return vc;
}

///////////////////////////////////////////////////////////////////////////////

int
Face::is_3_sided_convex_loop_quad_meshable(const int *segments, int *partsegments)
{
     double M[6][6], rhs[6];

     // octave:1> A = [ 0  0 -1  1  0  0;
     //                -1  0  0  0  1  0;
     //                 0 -1  0  0  0  1;
     //		       1  0  0  1  0  0;
     //                 0  1  0  0  1  0;
     //                 0  0  1  0  0  1]
     //A =
     //
     //   0   0  -1   1   0   0
     //  -1   0   0   0   1   0
     //   0  -1   0   0   0   1
     //   1   0   0   1   0   0
     //   0   1   0   0   1   0
     //   0   0   1   0   0   1
     //  octave:2> inv(A)
     //  ans =

     //  -0.5  -0.5   0.5   0.5   0.5  -0.5
     //   0.5  -0.5  -0.5  -0.5   0.5   0.5
     //  -0.5   0.5  -0.5   0.5  -0.5   0.5
     //   0.5   0.5  -0.5   0.5  -0.5   0.5
     //  -0.5   0.5   0.5   0.5   0.5  -0.5
     //   0.5  -0.5   0.5  -0.5   0.5   0.5

     //   First Row
     //  -0.5  -0.5   0.5   0.5   0.5  -0.5
     M[0][0] = -0.5;
     M[0][1] = -0.5;
     M[0][2] = 0.5;
     M[0][3] = 0.5;
     M[0][4] = 0.5;
     M[0][5] = -0.5;

     //   Second Row
     //   0.5  -0.5  -0.5  -0.5   0.5   0.5
     M[1][0] = 0.5;
     M[1][1] = -0.5;
     M[1][2] = -0.5;
     M[1][3] = -0.5;
     M[1][4] = 0.5;
     M[1][5] = 0.5;

     //  Third Row
     //  -0.5   0.5  -0.5   0.5  -0.5   0.5
     M[2][0] = -0.5;
     M[2][1] = 0.5;
     M[2][2] = -0.5;
     M[2][3] = 0.5;
     M[2][4] = -0.5;
     M[2][5] = 0.5;

     // Forth Row
     //   0.5   0.5  -0.5   0.5  -0.5   0.5
     M[3][0] = 0.5;
     M[3][1] = 0.5;
     M[3][2] = -0.5;
     M[3][3] = 0.5;
     M[3][4] = -0.5;
     M[3][5] = 0.5;

     //   Fifth Row
     //  -0.5   0.5   0.5   0.5   0.5  -0.5
     M[4][0] = -0.5;
     M[4][1] = 0.5;
     M[4][2] = 0.5;
     M[4][3] = 0.5;
     M[4][4] = 0.5;
     M[4][5] = -0.5;

     //  Sixth Row
     //   0.5  -0.5   0.5  -0.5   0.5   0.5
     M[5][0] = 0.5;
     M[5][1] = -0.5;
     M[5][2] = 0.5;
     M[5][3] = -0.5;
     M[5][4] = 0.5;
     M[5][5] = 0.5;

     rhs[0] = 0.0;
     rhs[1] = 0.0;
     rhs[2] = 0.0;

     rhs[3] = segments[0];
     rhs[4] = segments[1];
     rhs[5] = segments[2];

     vector<int> x(6);
     for (int i = 0; i < 6; i++) {
          double sum = 0.0;
          for (int j = 0; j < 6; j++)
               sum += M[i][j] * rhs[j];
          x[i] = (int) sum;
     }

     for (int i = 0; i < 6; i++)
          if (x[i] < 1) return 0;

     if (x[0] + x[3] != rhs[3]) return 0;
     partsegments[0] = x[0];
     partsegments[1] = x[3];

     if (x[1] + x[4] != rhs[4]) return 0;
     partsegments[2] = x[1];
     partsegments[3] = x[4];

     if (x[2] + x[5] != rhs[5]) return 0;
     partsegments[4] = x[2];
     partsegments[5] = x[5];

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
Face::is_5_sided_convex_loop_quad_meshable(const int *segments, int *partSegments)
{
     double M[10][10], rhs[10];

     //  Equations:
     //      b0 -a2   = 0
     //      b1 -a3   = 0
     //      b2 -a4   = 0
     //      b3 -a0   = 0
     //      b4 -a1   = 0
     //      a0 + b0  = s0
     //      a1 + b1  = s1
     //      a2 + b2  = s2
     //      a3 + b3  = s3
     //      a4 + b4  = s4

     // For more details; See Guy Bunin's paper.
     // M =
     //   a0  a1 a2  a3  a4  b0   b1 b2   b3  b4
     //   0   0  -1   0   0   1   0   0   0   0
     //   0   0   0  -1   0   0   1   0   0   0
     //   0   0   0   0  -1   0   0   1   0   0
     //  -1   0   0   0   0   0   0   0   1   0
     //   0  -1   0   0   0   0   0   0   0   1
     //   1   0   0   0   0   1   0   0   0   0
     //   0   1   0   0   0   0   1   0   0   0
     //   0   0   1   0   0   0   0   1   0   0
     //   0   0   0   1   0   0   0   0   1   0
     //   0   0   0   0   1   0   0   0   0   1

     // octave:38> inv(M)
     // ans =
     //  -0.5   0.5   0.5  -0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
     //  -0.5  -0.5   0.5   0.5  -0.5   0.5   0.5  -0.5  -0.5   0.5
     //  -0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
     //   0.5  -0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
     //   0.5   0.5  -0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
     //   0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
     //   0.5   0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
     //   0.5   0.5   0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
     //  -0.5   0.5   0.5   0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
     //  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5   0.5


     //   First Row
     //  -0.5   0.5   0.5  -0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
     M[0][0] = -0.5;
     M[0][1] = 0.5;
     M[0][2] = 0.5;
     M[0][3] = -0.5;
     M[0][4] = -0.5;
     M[0][5] = 0.5;
     M[0][6] = -0.5;
     M[0][7] = -0.5;
     M[0][8] = 0.5;
     M[0][9] = 0.5;

     //  -0.5  -0.5   0.5   0.5  -0.5   0.5   0.5  -0.5  -0.5   0.5
     M[1][0] = -0.5;
     M[1][1] = -0.5;
     M[1][2] = 0.5;
     M[1][3] = 0.5;
     M[1][4] = -0.5;
     M[1][5] = 0.5;
     M[1][6] = 0.5;
     M[1][7] = -0.5;
     M[1][8] = -0.5;
     M[1][9] = 0.5;


     //  -0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
     M[2][0] = -0.5;
     M[2][1] = -0.5;
     M[2][2] = -0.5;
     M[2][3] = 0.5;
     M[2][4] = 0.5;
     M[2][5] = 0.5;
     M[2][6] = 0.5;
     M[2][7] = 0.5;
     M[2][8] = -0.5;
     M[2][9] = -0.5;

     //   0.5  -0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
     M[3][0] = 0.5;
     M[3][1] = -0.5;
     M[3][2] = -0.5;
     M[3][3] = -0.5;
     M[3][4] = 0.5;
     M[3][5] = -0.5;
     M[3][6] = 0.5;
     M[3][7] = 0.5;
     M[3][8] = 0.5;
     M[3][9] = -0.5;

     //   0.5   0.5  -0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
     M[4][0] = 0.5;
     M[4][1] = 0.5;
     M[4][2] = -0.5;
     M[4][3] = -0.5;
     M[4][4] = -0.5;
     M[4][5] = -0.5;
     M[4][6] = -0.5;
     M[4][7] = 0.5;
     M[4][8] = 0.5;
     M[4][9] = 0.5;

     //   0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
     M[5][0] = 0.5;
     M[5][1] = -0.5;
     M[5][2] = -0.5;
     M[5][3] = 0.5;
     M[5][4] = 0.5;
     M[5][5] = 0.5;
     M[5][6] = 0.5;
     M[5][7] = 0.5;
     M[5][8] = -0.5;
     M[5][9] = -0.5;

     //   0.5   0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
     M[6][0] = 0.5;
     M[6][1] = 0.5;
     M[6][2] = -0.5;
     M[6][3] = -0.5;
     M[6][4] = 0.5;
     M[6][5] = -0.5;
     M[6][6] = 0.5;
     M[6][7] = 0.5;
     M[6][8] = 0.5;
     M[6][9] = -0.5;

     //   0.5   0.5   0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
     M[7][0] = 0.5;
     M[7][1] = 0.5;
     M[7][2] = 0.5;
     M[7][3] = -0.5;
     M[7][4] = -0.5;
     M[7][5] = -0.5;
     M[7][6] = -0.5;
     M[7][7] = 0.5;
     M[7][8] = 0.5;
     M[7][9] = 0.5;

     //  -0.5   0.5   0.5   0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
     M[8][0] = -0.5;
     M[8][1] = 0.5;
     M[8][2] = 0.5;
     M[8][3] = 0.5;
     M[8][4] = -0.5;
     M[8][5] = 0.5;
     M[8][6] = -0.5;
     M[8][7] = -0.5;
     M[8][8] = 0.5;
     M[8][9] = 0.5;

     //  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5   0.5
     M[9][0] = -0.5;
     M[9][1] = -0.5;
     M[9][2] = 0.5;
     M[9][3] = 0.5;
     M[9][4] = 0.5;
     M[9][5] = 0.5;
     M[9][6] = 0.5;
     M[9][7] = -0.5;
     M[9][8] = -0.5;
     M[9][9] = 0.5;

     rhs[0] = 0.0;
     rhs[1] = 0.0;
     rhs[2] = 0.0;
     rhs[3] = 0.0;
     rhs[4] = 0.0;

     rhs[5] = segments[0];
     rhs[6] = segments[1];
     rhs[7] = segments[2];
     rhs[8] = segments[3];
     rhs[9] = segments[4];

     vector<int> x(10);
     for (int i = 0; i < 10; i++) {
          double sum = 0.0;
          for (int j = 0; j < 10; j++)
               sum += M[i][j] * rhs[j];
          x[i] = (int) sum;
     }

     for (int i = 0; i < 10; i++)
          if (x[i] < 1) return 0;

     if (x[0] != x[8]) return 0;
     if (x[0] + x[5] != rhs[5]) return 0;
     partSegments[0] = x[0];
     partSegments[1] = x[5];

     if (x[1] != x[9]) return 0;
     if (x[1] + x[6] != rhs[6]) return 0;
     partSegments[2] = x[1];
     partSegments[3] = x[6];

     if (x[2] != x[5]) return 0;
     if (x[2] + x[7] != rhs[7]) return 0;
     partSegments[4] = x[2];
     partSegments[5] = x[7];

     if (x[3] != x[6]) return 0;
     if (x[3] + x[8] != rhs[8]) return 0;
     partSegments[6] = x[3];
     partSegments[7] = x[8];

     if (x[4] != x[7]) return 0;
     if (x[4] + x[9] != rhs[9]) return 0;
     partSegments[8] = x[4];
     partSegments[9] = x[9];

     return 1;
}

///////////////////////////////////////////////////////////////////////////////
int Mesh::refine_tri_edge( Vertex *v0, Vertex *v1, int numSegments,
                           vector<Vertex*> &newnodes, vector<Face*> &newfaces)
{
     assert( getAdjTable(0,2) );

     size_t nSize;
     newnodes.clear();
     newfaces.clear();
     assert( numSegments >= 2);

     FaceSequence efaces;
     Mesh::getRelations112( v0, v1, efaces);

     int numFaces = efaces.size();
     if( numFaces == 0) return 1;

     assert( numFaces <= 2 );

     newnodes.reserve(numSegments-1);
     newfaces.reserve(2*numSegments -numFaces ); // minuse because faces will be reused..

     double dt = 1.0/(double)numSegments;

     NodeSequence edgenodes;
     edgenodes.reserve(numSegments+1);
     edgenodes.push_back(v0);
     for( int i = 1; i < numSegments; i++) {
          Point3D xyz = Vertex::mid_point(v0, v1, i*dt);
          Vertex *vmid = Vertex::newObject();
          vmid->setXYZCoords( xyz );
          newnodes.push_back(vmid);
          edgenodes.push_back(vmid);
          addNode( vmid );
     }
     edgenodes.push_back(v1);

     int pos0, pos1, pos2;
     Vertex *ov1 =  Face::opposite_node(efaces[0], v0, v1);
     assert( ov1 );

     deactivate( efaces[0] );

     nSize = edgenodes.size();
     vector<Vertex*> conn(3);

     pos0 = efaces[0]->getPosOf(v0);
     pos1 = efaces[0]->getPosOf(v1);
     pos2 = efaces[0]->getPosOf(ov1);
     for( int i = 0; i < numSegments; i++) {
          conn[pos0] = edgenodes[i];
          conn[pos1] = edgenodes[i+1];
          conn[pos2] = ov1;
          if( i == 0) {
               efaces[0]->setNodes( conn );
               reactivate( efaces[0] );
          } else {
               Face *f = Face::newObject();
               f->setNodes( conn );
               addFace(f);
               newfaces.push_back(f);
          }
     }


     if( numFaces > 1 ) {
          Vertex *ov2 =  Face::opposite_node(efaces[1], v0, v1);
          assert( ov2 );
          deactivate( efaces[1] );

          nSize = edgenodes.size();
          pos0 = efaces[1]->getPosOf(v0);
          pos1 = efaces[1]->getPosOf(v1);
          pos2 = efaces[1]->getPosOf(ov2);
          for( int i = 0; i < numSegments; i++) {
               conn[pos0] = edgenodes[i];
               conn[pos1] = edgenodes[i+1];
               conn[pos2] = ov2;
               if( i == 0) {
                    efaces[1]->setNodes( conn );
                    reactivate( efaces[1] );
               } else {
                    Face *f = Face::newObject();
                    f->setNodes( conn );
                    addFace(f);
                    newfaces.push_back(f);
               }
          }
     }
     return 0;
}

int Mesh::collapse_tri_edge( Vertex *v0, Vertex *v1)
{
     assert( getAdjTable(0,2) );

     if( v0->isBoundary() && v1->isBoundary() ) return 1;

     FaceSequence efaces;
     Mesh::getRelations112( v0, v1, efaces);

     assert( efaces.size() == 2 );
     if( efaces[0]->isRemoved() || efaces[1]->isRemoved() ) return 1;

     Vertex *keepVertex   = NULL;
     Vertex *removeVertex = NULL;

     if( v0->isBoundary() ) {
          keepVertex   = v0;
          removeVertex = v1;
     }

     if( v1->isBoundary() ) {
          keepVertex   = v1;
          removeVertex = v0;
     }

     if( !v0->isBoundary()  &&  !v1->isBoundary() ) {
          keepVertex   = v0;
          removeVertex = v1;
          Point3D p3d = Vertex::mid_point( v0, v1);
          keepVertex->setXYZCoords( p3d );
     }

     remove( efaces[0] );
     remove( efaces[1] );

     const FaceSequence &v1faces = removeVertex->getRelations2();
     int nSize = v1faces.size();
     for( int i = 0; i < nSize; i++) {
          Face *f = v1faces[i];
          assert( !f->isRemoved() );
          int err = f->replaceNode(removeVertex, keepVertex);
          assert( !err);
          reactivate( f );
     }

     removeVertex->clearRelations(0);
     removeVertex->clearRelations(2);

     remove( removeVertex );

     return 0;

}
///////////////////////////////////////////////////////////////////////////////
int
Mesh::refine_quad15(Face *face)
{
     if( face == NULL ) return 1;
     if( face->isRemoved() ) return 2;

     NodeSequence newnodes;
     FaceSequence newfaces;
     face->refine_quad15(newnodes, newfaces);

     int numnodes =  newnodes.size();
     for (int i = 0; i < numnodes; i++)
          addNode(newnodes[i]);

     int numfaces = newfaces.size();
     for (int i = 0; i < numfaces; i++)
          addFace(newfaces[i]);

     remove(face);

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

Jaal::NodeSequence
Mesh::boundary_chain_nodes(Vertex *v0, Vertex *v1)
{
     NodeSequence bndnodes;

     FaceSequence neighs;
     Mesh::getRelations112(v0, v1, neighs);

     if (neighs.size() != 2) return bndnodes;

     vector<Edge> bndedges;

     bndedges.reserve(6);
     Edge sharededge(v0, v1);

     for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 4; j++) {
               Vertex *vf0 = neighs[i]->getNodeAt((j + 0) % 4);
               Vertex *vf1 = neighs[i]->getNodeAt((j + 1) % 4);
               Edge edge(vf0, vf1);
               if (!edge.isSameAs(sharededge))
                    bndedges.push_back(edge);
          }
     }

     Mesh::make_chain(bndedges);

#ifdef SEQUENCE_IS_VECTOR
     bndnodes.reserve(6);
#endif
     bndnodes.push_back(bndedges[0].getNodeAt(0));
     bndnodes.push_back(bndedges[0].getNodeAt(1));

     size_t nSize = bndedges.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          Vertex *v0 = bndedges[i].getNodeAt(0);
          Vertex *v1 = bndedges[i].getNodeAt(1);
          if (v0 == bndnodes[i]) {
               bndnodes.push_back(v1);
          } else if (v1 == bndnodes[i])
               bndnodes.push_back(v0);
          else {
               cout << "Error in bound edges : " << endl;
               bndnodes.clear();
               return bndnodes;
          }
     }
     return bndnodes;
}

/////////////////////////////////////////////////////////////////////////////

int
Mesh::check_convexity()
{
     size_t numfaces = getSize(2);

     size_t nconvex = 0, nconcave = 0;
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          if (face->isConvex())
               nconvex++;
          else
               nconcave++;
     }
     cout << "Info: Total Number of faces " << numfaces << " #Convex : " << nconvex << " #Concave " << nconcave << endl;
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

void
Mesh::getMinMaxFaceArea(double &minarea, double &maxarea)
{
     size_t numfaces = getSize(2);

     vector<double> area(numfaces);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          area[i] = face->getArea();
     }
     minarea = *min_element(area.begin(), area.end());
     maxarea = *max_element(area.begin(), area.end());
}


////////////////////////////////////////////////////////////////////////////////

int
Mesh::getCoordsArray( vector<double> &vcoords, vector<size_t> &l2g)
{
     size_t numnodes = getSize(0);

     vcoords.reserve(3 * numnodes);
     l2g.reserve(numnodes);

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          if( !v->isRemoved() ) {
               l2g.push_back( v->getID() );
               const Point3D &xyz = v->getXYZCoords();
               vcoords.push_back( xyz[0] );
               vcoords.push_back( xyz[1] );
               vcoords.push_back( xyz[2] );
          }
     }
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
Mesh::getNodesArray( vector<size_t> &nodearray, vector<int> &elmtopo)
{
     elmtopo.clear();
     nodearray.clear();

     size_t numfaces = getSize(2);
     size_t numnodes = getSize(0);

     int topo = isHomogeneous();
     if (topo) nodearray.reserve(topo * numfaces);

     map<size_t,size_t> nodemap;

     size_t index = 0;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          if( !v->isRemoved() )
               nodemap[v->getID()] = index++;
     }

     elmtopo.reserve( numfaces );

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          if( !face->isRemoved() ) {
               elmtopo.push_back( face->getSize(0) );
               for (int j = 0; j < face->getSize(0); j++) {
                    Vertex *v = face->getNodeAt(j);
                    int lid   = nodemap[v->getID()];
                    nodearray.push_back(lid);
               }
          }
     }
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

Mesh*
Mesh::deep_copy()
{
     std::map<Vertex*, Vertex*> vmap;

     Mesh *newmesh = new Mesh;
     size_t numnodes = getSize(0);
     newmesh->reserve(numnodes, 0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vold = getNodeAt(i);
          Vertex *vnew = Vertex::newObject();
          vnew->setXYZCoords(vold->getXYZCoords());
          vmap[vold] = vnew;
          newmesh->addNode(vnew);
     }

     size_t numfaces = getSize(2);
     newmesh->reserve(numnodes, 2);
     NodeSequence connect;
     for (size_t i = 0; i < numfaces; i++) {
          Face *fold = getFaceAt(i);
          connect  = fold->getNodes();
          int nv   = connect.size();
          for (int j = 0; j < nv; j++)
               connect[j] = vmap[connect[j]];

          Face *fnew = Face::newObject();
          fnew->setNodes(connect);
          newmesh->addFace(fnew);
     }

     return newmesh;
}

////////////////////////////////////////////////////////////////////////////////

int
Mesh::setCoordsArray(const vector<double> &vcoords, const vector<size_t> &l2g)
{
     size_t numnodes = getSize(0);

     Point3D xyz;
     int index = 0;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          if( !v->isRemoved() ) {
               assert( v->getID() == l2g[index] );
               xyz[0] = vcoords[3 * index + 0];
               xyz[1] = vcoords[3 * index + 1];
               xyz[2] = vcoords[3 * index + 2];
               v->setXYZCoords(xyz);
               index++;
          }
     }
     return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
Mesh::make_chain(vector<Edge> &boundedges)
{
     size_t nSize = boundedges.size();

     Edge edge = boundedges.front();

     list<Edge> listedges;
     for (size_t i = 1; i < nSize; i++)
          listedges.push_back(boundedges[i]);

     boundedges.clear();
     boundedges.reserve(nSize);

     Vertex *curr_vertex = edge.getNodeAt(1);
     boundedges.push_back(edge);

     list<Edge>::iterator it;

     for (size_t i = 0; i < nSize; i++) {
          for (it = listedges.begin(); it != listedges.end(); ++it) {
               edge = *it;
               Vertex *v0 = edge.getNodeAt(0);
               Vertex *v1 = edge.getNodeAt(1);
               if (v0 == curr_vertex) {
                    curr_vertex = v1;
                    boundedges.push_back(edge);
                    break;
               }
               if (v1 == curr_vertex) {
                    curr_vertex = v0;
                    Edge newedge(v1, v0);
                    boundedges.push_back(newedge);
                    break;
               }
          }
          if (it != listedges.end()) listedges.erase(it);
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::is_closeable_chain(const vector<Edge> &boundedges)
{
     std::map<Vertex*, set<Vertex*> > relations00;

     for (size_t i = 0; i < boundedges.size(); i++) {
          Vertex *v0 = boundedges[i].getNodeAt(0);
          Vertex *v1 = boundedges[i].getNodeAt(1);
          relations00[v0].insert(v1);
          relations00[v1].insert(v0);
     }

     std::map<Vertex*, set<Vertex*> > ::const_iterator it;
     for (it = relations00.begin(); it != relations00.end(); ++it) {
          Vertex *v = it->first;
          if (relations00[v].size() != 2) return 0;
     }

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::is_closed_chain(const vector<Edge> &boundedges)
{
     Vertex *first_vertex = boundedges[0].getNodeAt(0);
     Vertex *curr_vertex = boundedges[0].getNodeAt(1);

     for (size_t i = 1; i < boundedges.size(); i++) {
          if (boundedges[i].getNodeAt(0) != curr_vertex) return 0;
          curr_vertex = boundedges[i].getNodeAt(1);
     }
     if (curr_vertex != first_vertex) return 0;

     return 1;
}
///////////////////////////////////////////////////////////////////////////////

int
Mesh::rotate_chain(vector<Edge> &boundedges, Vertex *first_vertex)
{
     size_t nSize = boundedges.size();
     if (!is_closeable_chain(boundedges)) return 1;

     vector<Edge> listedges(nSize);
     int istart = 0;
     for (size_t i = 0; i < nSize; i++) {
          if (boundedges[i].getNodeAt(0) == first_vertex) istart = i;
          listedges[i] = boundedges[i];
     }

     for (size_t i = 0; i < nSize; i++)
          boundedges[i] = listedges[(i + istart) % nSize];

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

int
Mesh::getRelations102(const PNode vtx0, const PNode vtx1, FaceSequence &faceneighs)
{
     faceneighs.clear();

     assert( !vtx0->isRemoved() && !vtx1->isRemoved() );

     const FaceSequence &v0faces = vtx0->getRelations2();
     const FaceSequence &v1faces = vtx1->getRelations2();

     if (v0faces.empty() || v1faces.empty()) {
          cout << "Warning: Vertex-Faces relations are empty " << endl;
          return 1;
     }

     FaceSet vset;
     for (size_t i = 0; i < v0faces.size(); i++)
          vset.insert(v0faces[i]);

     for (size_t i = 0; i < v1faces.size(); i++)
          vset.insert(v1faces[i]);

     FaceSet::const_iterator it;

     if (vset.size()) {
          faceneighs.resize(vset.size());
          int index = 0;
          for (it = vset.begin(); it != vset.end(); ++it)
               faceneighs[index++] = *it;
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::getRelations112(const PNode vtx0, const PNode vtx1, FaceSequence &faceneighs)
{
     faceneighs.clear();

     const FaceSequence &v0faces = vtx0->getRelations2();
     const FaceSequence &v1faces = vtx1->getRelations2();

     if (v0faces.empty() || v1faces.empty()) {
          cout << "Warning: Vertex-Faces relations are empty " << endl;
          return 1;
     }

     for (size_t i = 0; i < v0faces.size(); i++)
          assert(!v0faces[i]->isRemoved());

     for (size_t i = 0; i < v0faces.size() - 1; i++)
          assert(v0faces[i] < v0faces[i + 1]);

     for (size_t i = 0; i < v1faces.size(); i++)
          assert(!v1faces[i]->isRemoved());

     for (size_t i = 0; i < v1faces.size() - 1; i++)
          assert(v1faces[i] < v1faces[i + 1]);

     set_intersection(v0faces.begin(), v0faces.end(), v1faces.begin(),
                      v1faces.end(), back_inserter(faceneighs));

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::count_edges()
{
     if (getAdjTable(1, 0)) return edges.size();

     int relexist = build_relations(0, 0);

     size_t numnodes = getSize(0);

     size_t ncount = 0;
     for (size_t i = 0; i < numnodes; i++) {
          const NodeSequence &neighs = nodes[i]->getRelations0();
          for (size_t j = 0; j < neighs.size(); j++)
               if (nodes[i] > neighs[j])
                    ncount++;
     }

     if (!relexist) clear_relations(0, 0);

     return ncount;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::build_edges()
{
     int relexist0 = build_relations(0, 0);

     size_t numnodes = getSize(0);

     edges.clear();

     for (size_t i = 0; i < numnodes; i++) {
          NodeSequence &neighs = nodes[i]->getRelations0();
          for (size_t j = 0; j < neighs.size(); j++)
               if (nodes[i] > neighs[j]) {
                    Edge *newedge = new Edge(nodes[i], neighs[j]);
                    assert(newedge);
                    edges.push_back(newedge);
               }
     }

     if (!relexist0) clear_relations(0, 0);
     adjTable[1][0] = 1;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::prune()
{
     size_t nSize;

     nSize = faces.size();
     for (size_t i = 0; i < nSize; i++) {
          Face *f = faces[i];
          if (f->isRemoved()) {
               int nnodes = f->getSize(0);
               for( int j = 0; j < nnodes; j++) {
                    Vertex *v = f->getNodeAt(j);
                    v->setStatus( MeshEntity::REMOVE );
               }
          }
     }

     nSize = faces.size();
     for (size_t i = 0; i < nSize; i++) {
          Face *f = faces[i];
          if (!f->isRemoved()) {
               int nnodes = f->getSize(0);
               for( int j = 0; j < nnodes; j++) {
                    Vertex *v = f->getNodeAt(j);
                    v->setStatus( MeshEntity::ACTIVE );
               }
          }
     }

     if (adjTable[0][0]) {
          nSize = nodes.size();
          for (size_t i = 0; i < nSize; i++) {
               Vertex *vi = nodes[i];
               if( vi->isRemoved() ) {
                    const NodeSequence &vneighs = vi->getRelations0();
                    for (size_t j = 0; j < vneighs.size(); j++)
                         vneighs[j]->removeRelation0(vi);
                    vi->clearRelations(0);
               }
          }
     }

     if (adjTable[0][2]) {
          nSize = faces.size();
          for (size_t i = 0; i < nSize; i++) {
               Face *f = faces[i];
               if( f->isRemoved() ) {
                    for (int j = 0; j < f->getSize(0); j++) {
                         Vertex *v = f->getNodeAt(j);
                         v->removeRelation2(f);
                    }
               }
          }
     }

     NodeSequence::iterator vend;
     vend = remove_if(nodes.begin(), nodes.end(), EntityRemovedPred());
     nodes.erase(vend, nodes.end());

     nSize = nodes.size();
     for (size_t i = 0; i < nSize; i++)
          assert(nodes[i]->isActive());

     EdgeSequence::iterator eend;
     eend = remove_if(edges.begin(), edges.end(), EntityRemovedPred());
     edges.erase(eend, edges.end());

     nSize = edges.size();
     for(size_t i = 0; i < nSize; i++)
          assert(edges[i]->isActive());
  

     cout << " Before Prune " << faces.size() << endl;

     FaceSequence::iterator fend;
     fend = remove_if(faces.begin(), faces.end(), EntityRemovedPred());
     faces.erase(fend, faces.end());
     nSize = faces.size();
     for (size_t i = 0; i < nSize; i++)
          assert(faces[i]->isActive());

     cout << " After Prune : " << faces.size() << endl;

     enumerate(0);
     enumerate(2);
}

///////////////////////////////////////////////////////////////////////////////

bool
Mesh::isPruned() const
{
     for (size_t i = 0; i < nodes.size(); i++)
          if (!nodes[i]->isActive()) return 0;

     for (size_t i = 0; i < faces.size(); i++)
          if (!faces[i]->isActive()) return 0;

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::collect_garbage()
{
     if( !isPruned() ) prune();

     list<Face*>::const_iterator fiter;
     for (fiter = garbageFaces.begin(); fiter != garbageFaces.end(); ++fiter) {
          Face *face = *fiter;
          assert(face);
          if (face->isRemoved()) delete face;
     }
     garbageFaces.clear();

     if( nodes.empty() ) {
          list<Vertex*>::const_iterator viter;
          for (viter = garbageNodes.begin(); viter != garbageNodes.end(); ++viter) {
               Vertex *vertex = *viter;
               assert(vertex);
               if (vertex->isRemoved()) delete vertex;
          }
          garbageNodes.clear();
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::enumerate(int etype)
{
     size_t index = 0;

     NodeSequence::const_iterator viter;
     if (etype == 0) {
          index = 0;
          for (viter = nodes.begin(); viter != nodes.end(); ++viter) {
               Vertex *vertex = *viter;
               if( !vertex->isRemoved() ) vertex->setID(index++);
          }
     }

     FaceSequence::const_iterator fiter;
     if (etype == 2) {
          index = 0;
          for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
               Face *face = *fiter;
               if( !face->isRemoved() ) face->setID(index++);
          }
     }
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::getBoundarySize(int d)
{
     if (boundary_known == 0) search_boundary();

     size_t nSize, ncount = 0;

     if (d == 0) {
          nSize = nodes.size();
          for (size_t i = 0; i < nSize; i++)
               if (!nodes[i]->isRemoved() && nodes[i]->isBoundary())
                    ncount++;
     }

     if (d == 2) {
          nSize = faces.size();
          for (size_t i = 0; i < nSize; i++)
               if (!faces[i]->isRemoved() && faces[i]->isBoundary())
                    ncount++;
     }

     return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::isHomogeneous() const
{
     int maxnodes = 0;
     size_t nSize = faces.size();
     for (size_t i = 0; i < nSize; i++) {
          if (!faces[i]->isRemoved())
               maxnodes = max(maxnodes, faces[i]->getSize(0));
     }

     return maxnodes;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::saveAs(const string &s)
{
     MeshExporter mexp;
     return mexp.saveAs(this, s);
}
///////////////////////////////////////////////////////////////////////////////

Vertex*
Mesh::nearest_neighbour(const Vertex *myself, double &mindist)
{
     Vertex* nearest;
     assert(getAdjTable(0, 0));

     mindist = MAXDOUBLE;
     nearest = NULL;

     const NodeSequence &neighs = myself->getRelations0();

     for (size_t i = 0; i < neighs.size(); i++) {
          double d = Vertex::length2(myself, neighs[i]);
          if (d < mindist) {
               mindist = d;
               nearest = neighs[i];
          }
     }

     mindist = sqrt(mindist);

     return nearest;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::build_relations02(bool rebuild)
{
     if (rebuild) clear_relations(0, 2);

     int status = adjTable[0][2];

     if (adjTable[0][2] == 1) {
          size_t numnodes = getSize(0);
          for (size_t inode = 0; inode < numnodes; inode++) {
               Vertex *vtx = getNodeAt( inode);
               if( vtx ) {
                    if( vtx->isRemoved() )
                         vtx->clearRelations(2);
                    else {
                         FaceSequence &vfaces = vtx->getRelations2();
                         FaceSequence::iterator fend;
                         fend = remove_if(vfaces.begin(), vfaces.end(), EntityRemovedPred());
                         vfaces.erase(fend, vfaces.end());
                    }
               }
          }
     }

     size_t numfaces = getSize(2);
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = getFaceAt(iface);
          if( face ) {
               if( !face->isRemoved() ) {
                    int nf = face->getSize(0);
                    for (int j = 0; j < nf; j++) {
                         Vertex *vtx = face->getNodeAt(j);
                         vtx->addRelation2(face);
                    }
               }
          }
     }

     adjTable[0][2] = 1;

     return status;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::build_relations00(bool rebuild)
{
     int status = adjTable[0][0];
     // Delete all the old stuff and rebuild new. Sometimes used for debugging purpose
     // also.
     if (rebuild) clear_relations(0, 0);

     if (adjTable[0][0] == 1) {
          size_t numnodes = getSize(0);
          for (size_t inode = 0; inode < numnodes; inode++) {
               Vertex *vtx = getNodeAt( inode);
               if( vtx ) {
                    if( vtx->isRemoved() )
                         vtx->clearRelations(0);
                    else {
                         NodeSequence &vneighs = vtx->getRelations0();
                         NodeSequence::iterator vend;
                         vend = remove_if(vneighs.begin(), vneighs.end(), EntityRemovedPred());
                         vneighs.erase(vend, vneighs.end());
                    }
               }
          }
     }

     size_t numfaces = getSize(2);
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = getFaceAt(iface);
          if( face ) {
               size_t nnodes = face->getSize(0);
               for (size_t j = 0; j < nnodes; j++) {
                    Vertex *v0 = face->getNodeAt(j);
                    Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
                    v0->addRelation0(v1);
                    v1->addRelation0(v0);
               }
          }
     }
     adjTable[0][0] = 1;

     return status;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::clear_relations(int src, int dst)
{
     size_t numnodes = getSize(0);

     if (src == 0 && dst == 0) {
          for (size_t i = 0; i < numnodes; i++) {
               Vertex *vtx = getNodeAt(i);
               vtx->clearRelations(0);
          }
          adjTable[0][0] = 0;
     }

     if (src == 0 && dst == 2) {
          for (size_t i = 0; i < numnodes; i++) {
               Vertex *vtx = getNodeAt(i);
               vtx->clearRelations(2);
          }
          adjTable[0][2] = 0;
     }
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::search_boundary()
{
     if (boundary_known == 1) return 1;

     if (!isPruned()) prune();
     int relexist = build_relations(0, 2);
     size_t numfaces = getSize(2);

     int bmark;
     FaceSequence neighs;
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = getFaceAt(iface);
          size_t nnodes = face->getSize(0);
          for (size_t j = 0; j < nnodes; j++) {
               Vertex *v0 = face->getNodeAt(j);
               Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
               Mesh::getRelations112(v0, v1, neighs);
               if (neighs.size() == 1) {
                    bmark = max(1, face->getBoundaryMark());
                    face->setBoundaryMark(bmark);

                    bmark = max(1, v0->getBoundaryMark());
                    v0->setBoundaryMark(bmark);

                    bmark = max(1, v1->getBoundaryMark());
                    v1->setBoundaryMark(bmark);
               }
          }
     }

     if (!relexist)
          clear_relations(0, 2);

     // Calculate the boundary feature angles also.
     // setFeatureAngles ();

     boundary_known = 1;

     return 0;

}
///////////////////////////////////////////////////////////////////////////////

Jaal::FaceSequence
Mesh::filter(int facetype) const
{
     FaceSequence::const_iterator it;
     size_t ncount = 0;
     for (it = faces.begin(); it != faces.end(); ++it) {
          Face *face = *it;
          if (face->getType() == facetype)
               ncount++;
     }

     FaceSequence tmpfaces;
     if (ncount) {
          tmpfaces.resize(ncount);
          size_t index = 0;
          for (it = faces.begin(); it != faces.end(); ++it) {
               Face *face = *it;
               if (face->getType() == facetype)
                    tmpfaces[index++] = face;
          }
     }

     return tmpfaces;
}

///////////////////////////////////////////////////////////////////////////////

bool
Mesh::isSimple()
{
     int simple = 1;
     int relexist = build_relations(0, 2);

     size_t numfaces = getSize(2);

     FaceSequence neighs;
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = getFaceAt(iface);
          assert(face);
          int nnodes = face->getSize(0);
          for (int j = 0; j < nnodes; j++) {
               Vertex *v0 = face->getNodeAt(j);
               Vertex *v1 = face->getNodeAt(j + 1);
               Mesh::getRelations112(v0, v1, neighs);
               if (neighs.size() > 2) {
                    simple = 0;
                    break;
               }
          }
     }

     if (!relexist)
          clear_relations(0, 2);

     return simple;

}

///////////////////////////////////////////////////////////////////////////////

size_t Mesh::count_irregular_nodes(int degree_of_regular_node)
{
     int relexist = build_relations(0, 2);
     search_boundary();

     size_t numnodes = getSize(0);
     size_t ncount = 0;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = getNodeAt(i);
          if( !vertex->isRemoved() ) {
               if (!vertex->isBoundary()) {
                    const FaceSequence &vfaces = vertex->getRelations2();
                    int nSize = vfaces.size();
                    if ( nSize != degree_of_regular_node)
                         ncount++;
               }
          }
     }

     if (!relexist)
          clear_relations(0, 2);

     return ncount;
}
///////////////////////////////////////////////////////////////////////////////

bool
Mesh::is_consistently_oriented()
{
     int consistent = 1;
     int relexist = build_relations(0, 2);

     size_t numfaces = getSize(2);

     FaceSequence neighs;
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = getFaceAt(iface);
          if( !face->isRemoved() ) {
               assert(face);
               int nnodes = face->getSize(0);
               for (int j = 0; j < nnodes; j++) {
                    Vertex *v0 = face->getNodeAt(j);
                    Vertex *v1 = face->getNodeAt(j + 1);
                    Mesh::getRelations112(v0, v1, neighs);
                    if (neighs.size() == 2) {
                         assert(!neighs[0]->isRemoved());
                         assert(!neighs[1]->isRemoved());
                         int dir1 = neighs[0]->getOrientation(v0, v1);
                         int dir2 = neighs[1]->getOrientation(v0, v1);
                         if (dir1 * dir2 == 1) {
                              cout << "Warning: Mesh is not consistently oriented " << endl;
                              cout << neighs[0]->getID() << " Face 1: ";
                              for (int k = 0; k < neighs[0]->getSize(0); k++)
                                   cout << neighs[0]->getNodeAt(k)->getID() << " ";
                              cout << endl;
                              cout << neighs[1]->getID() << " Face 2: ";
                              for (int k = 0; k < neighs[1]->getSize(0); k++)
                                   cout << neighs[1]->getNodeAt(k)->getID() << " ";
                              cout << endl;
                              consistent = 0;
                              break;
                         }
                    }
               }
               if (!consistent)
                    break;
          }
     }

     if (!relexist)
          clear_relations(0, 2);

     return consistent;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::make_consistently_oriented()
{
     int relexist2 = build_relations(0, 2);

     Face *face = NULL;
     deque<Face*> faceQ;
     FaceSequence neighs;

     size_t numfaces = getSize(2);

     for (size_t iface = 0; iface < numfaces; iface++) {
          face = getFaceAt(iface);
          face->setID(iface);
          face->setVisitMark(0);
     }

     face = getFaceAt(0);
     faceQ.push_back(face);

     while (!faceQ.empty()) {
          Face *face = faceQ.front();
          faceQ.pop_front();
          if (!face->isVisited()) {
               face->setVisitMark(1);
               int nnodes = face->getSize(0);
               for (int j = 0; j < nnodes; j++) {
                    Vertex *v0 = face->getNodeAt(j);
                    Vertex *v1 = face->getNodeAt(j + 1);
                    Mesh::getRelations112(v0, v1, neighs);
                    if (neighs.size() == 2) {
                         int dir1 = neighs[0]->getOrientation(v0, v1);
                         int dir2 = neighs[1]->getOrientation(v0, v1);
                         if (dir1 * dir2 == 1) {
                              if (!neighs[0]->isVisited() && neighs[1]->isVisited())
                                   neighs[0]->reverse();

                              if (!neighs[1]->isVisited() && neighs[0]->isVisited())
                                   neighs[1]->reverse();
                         }
                         faceQ.push_back(neighs[0]);
                         faceQ.push_back(neighs[1]);
                    }
               }
          }
     }

     for (size_t iface = 0; iface < numfaces; iface++) {
          face = getFaceAt(iface);
          if (!face->isVisited())
               cout << "Error: not visited : " << face->getID() << " "
                    << face->isVisited() << endl;
     }

     if( !relexist2) clear_relations(0, 2);
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::getNumOfComponents(bool stop_at_interface)
{
     build_relations(0, 2);

     Face *face = NULL;
     deque<Face*> faceQ;
     FaceSequence neighs;

     size_t numfaces = getSize(2);

     for (size_t iface = 0; iface < numfaces; iface++) {
          face = getFaceAt(iface);
          face->setID(iface);
          face->setVisitMark(0);
     }

     int numComponents = 0;

     while (1) {
          face = NULL;
          faceQ.clear();
          for (size_t iface = 0; iface < numfaces; iface++) {
               face = getFaceAt(iface);
               if (!face->isVisited()) {
                    face->setPartID(numComponents);
                    faceQ.push_back(face);
                    break;
               }
          }

          if (faceQ.empty())
               break;

          while (!faceQ.empty()) {
               Face *face = faceQ.front();
               faceQ.pop_front();
               if (!face->isVisited()) {
                    face->setPartID(numComponents);
                    face->setVisitMark(1);
                    int nnodes = face->getSize(0);
                    for (int j = 0; j < nnodes; j++) {
                         Vertex *v0 = face->getNodeAt(j);
                         Vertex *v1 = face->getNodeAt((j + 1) % nnodes);

                         int proceed = 1;
                         if (stop_at_interface) {
                              if (v0->isConstrained() && v1->isConstrained()) {
                                   Edge edge(v0, v1);
                                   if (hasFeatureEdge(edge)) proceed = 0;
                              }
                         }

                         if (proceed) {
                              Mesh::getRelations112(v0, v1, neighs);
                              if (neighs.size() == 2) {
                                   faceQ.push_back(neighs[0]);
                                   faceQ.push_back(neighs[1]);
                              }
                         }

                    }
               }
          } // Complete one Component
          numComponents++;
     }

     for (size_t iface = 0; iface < numfaces; iface++) {
          face = getFaceAt(iface);
          if (!face->isVisited())
               cout << "Error: not visited : " << face->getID() << " "
                    << face->isVisited() << endl;
     }

     clear_relations(0, 2);
     return numComponents;
}

///////////////////////////////////////////////////////////////////////////////

NodeList &Mesh::get_nodes_list()
{
     nodelist.clear();
     for (size_t i = 0; i < nodes.size(); i++)
          nodelist.push_back(nodes[i]);
     nodes.clear();
     return nodelist;
}

FaceList &Mesh::get_faces_list()
{
     facelist.clear();
     for (size_t i = 0; i < faces.size(); i++)
          facelist.push_back(faces[i]);
     faces.clear();
     return facelist;
}

Mesh* Mesh::getComponent(int id)
{
     Mesh *submesh = new Mesh;

     size_t numnodes = getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          v->setVisitMark(0);
     }

     size_t numfaces = getSize(2);

     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = getFaceAt(iface);
          if (face->getGroupID() == id) {
               for (int j = 0; j < face->getSize(0); j++) {
                    Vertex *v = face->getNodeAt(j);
                    if (!v->isVisited()) {
                         submesh->addNode(v);
                         v->setVisitMark(1);
                    }
               }
               submesh->addFace(face);
          }
     }
     return submesh;
}
///////////////////////////////////////////////////////////////////////////////

Mesh *
Jaal::struct_tri_grid(int nx, int ny)
{
     Mesh *trimesh = new Mesh;

     double dx = 2.0 / (nx - 1);
     double dy = 2.0 / (ny - 1);

     Point3D xyz;

     int index = 0;
     for (int j = 0; j < ny; j++) {
          for (int i = 0; i < nx; i++) {
               xyz[0] = -1.0 + i * dx;
               xyz[1] = -1.0 + j * dy;
               xyz[2] = 0.0;
               Vertex *vnew = Vertex::newObject();
               vnew->setID(index++);
               vnew->setXYZCoords(xyz);
               trimesh->addNode(vnew);
          }
     }

     NodeSequence connect(3);
     index = 0;
     Face *newtri;
     for (int j = 0; j < ny - 1; j++) {
          for (int i = 0; i < nx - 1; i++) {
               int n0 = j * nx + i;
               int n1 = n0 + 1;
               int n2 = n1 + nx;
               int n3 = n0 + nx;
               connect[0] = trimesh->getNodeAt(n0);
               connect[1] = trimesh->getNodeAt(n1);
               connect[2] = trimesh->getNodeAt(n2);
               newtri = Face::newObject();
               newtri->setNodes(connect);
               trimesh->addFace(newtri);

               connect[0] = trimesh->getNodeAt(n0);
               connect[1] = trimesh->getNodeAt(n2);
               connect[2] = trimesh->getNodeAt(n3);
               newtri = Face::newObject();
               newtri->setNodes(connect);
               trimesh->addFace(newtri);
          }
     }
     return trimesh;
}

/////////////////////////////////////////////////////////////////////////////

Mesh *
Jaal::struct_quad_grid(int nx, int ny)
{
     Mesh *quadmesh = new Mesh;

     double dx = 2.0 / (nx - 1);
     double dy = 2.0 / (ny - 1);

     Point3D xyz;

     int index = 0;
     for (int j = 0; j < ny; j++) {
          for (int i = 0; i < nx; i++) {
               xyz[0] = -1.0 + i * dx;
               xyz[1] = -1.0 + j * dy;
               xyz[2] = 0.0;
               Vertex *vnew = Vertex::newObject();
               vnew->setID(index++);
               vnew->setXYZCoords(xyz);
               quadmesh->addNode(vnew);
          }
     }

     NodeSequence connect(4);
     index = 0;
     Face *newquad;
     for (int j = 0; j < ny - 1; j++) {
          for (int i = 0; i < nx - 1; i++) {
               int n0 = j * nx + i;
               int n1 = n0 + 1;
               int n2 = n1 + nx;
               int n3 = n0 + nx;
               connect[0] = quadmesh->getNodeAt(n0);
               connect[1] = quadmesh->getNodeAt(n1);
               connect[2] = quadmesh->getNodeAt(n2);
               connect[3] = quadmesh->getNodeAt(n3);
               newquad = Face::newObject();
               newquad->setNodes(connect);
               quadmesh->addFace(newquad);
          }
     }
     return quadmesh;
}

////////////////////////////////////////////////////////////////////

void
expand_strip(Face *prevface, Vertex *v0, Vertex *v1, list<Face*> &strip)
{
     FaceSequence neighs;
     Mesh::getRelations112(v0, v1, neighs);

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

void
Mesh::get_quad_strips(Face *rootface, FaceSequence &strip1,
                      FaceSequence &strip2)
{

     Vertex *v0, *v1;

     list<Face*> strip01, strip12, strip23, strip03;
     list<Face*>::const_iterator it;
     size_t numfaces = getSize(2);

     // Strip Starting from edge 0-1
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          face->setVisitMark(0);
     }

     v0 = rootface->getNodeAt(0);
     v1 = rootface->getNodeAt(1);

     rootface->setVisitMark(1);
     strip01.push_back(rootface);
     expand_strip(rootface, v0, v1, strip01);
     for (it = strip01.begin(); it != strip01.end(); ++it)
          strip1.push_back(*it);

     // Strip Starting from edge 2-3
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          face->setVisitMark(0);
     }

     v0 = rootface->getNodeAt(2);
     v1 = rootface->getNodeAt(3);

     rootface->setVisitMark(1);
     strip23.push_back(rootface);
     expand_strip(rootface, v0, v1, strip23);
     for (it = strip23.begin(); it != strip23.end(); ++it)
          strip1.push_back(*it);

     // Strip Starting from edge 1-2
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          face->setVisitMark(0);
     }

     v0 = rootface->getNodeAt(1);
     v1 = rootface->getNodeAt(2);

     rootface->setVisitMark(1);
     strip12.push_back(rootface);
     expand_strip(rootface, v0, v1, strip12);
     for (it = strip12.begin(); it != strip12.end(); ++it)
          strip2.push_back(*it);

     // Strip Starting from edge 0-3
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          face->setVisitMark(0);
     }

     v0 = rootface->getNodeAt(0);
     v1 = rootface->getNodeAt(3);

     rootface->setVisitMark(1);
     strip12.push_back(rootface);
     expand_strip(rootface, v0, v1, strip03);

     for (it = strip03.begin(); it != strip03.end(); ++it)
          strip2.push_back(*it);
}
////////////////////////////////////////////////////////////////////////////////////////

FaceSequence
Mesh::get_bound_faces(int bound_what)
{
     int relexist2 = build_relations(0, 2);

     assert(getAdjTable(0, 2));

     search_boundary();

     size_t numfaces = getSize(2);

     set<Face*> bfaces;

     if (bound_what == 0) {
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               if (face->hasBoundaryNode())
                    bfaces.insert(face);
          }
     }

     if (bound_what == 1) {
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               if (face->has_boundary_edge())
                    bfaces.insert(face);
          }
     }

     FaceSequence result;

     size_t nSize = bfaces.size();

     if (nSize) {
          result.resize(nSize);
          set<Face*>::const_iterator it;

          size_t index = 0;
          for (it = bfaces.begin(); it != bfaces.end(); ++it)
               result[index++] = *it;
     }

     if (!relexist2)
          clear_relations(0, 2);

     return result;
}

///////////////////////////////////////////////////////////////////////////////

NodeSequence Mesh::get_bound_nodes()
{
     search_boundary();

     size_t numnodes = getSize(0);

     NodeSequence seq;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          if (v->isBoundary()) seq.push_back(v);
     }
     return seq;
}
///////////////////////////////////////////////////////////////////////////////

NodeSequence Mesh::get_irregular_nodes(int regular_count, int from_where)
{
     int relexist2 = build_relations(0, 2);
     search_boundary();

     size_t numnodes = getSize(0);
     int    nSize;

     NodeSequence seq;

     // Query from the boundary nodes ...
     if (from_where == 1) {
          for (size_t i = 0; i < numnodes; i++) {
               Vertex *v = getNodeAt(i);
               if( !v->isRemoved() ) {
                    nSize  = v->getNumRelations(2);
                    if (v->isBoundary() && nSize != regular_count) seq.push_back(v);
               }
          }
     }

     // Query from the internal nodes ...
     if (from_where == 0) {
          for (size_t i = 0; i < numnodes; i++) {
               Vertex *v = getNodeAt(i);
               if( !v->isRemoved() ) {
                    nSize  = v->getNumRelations(2);
                    if (!v->isBoundary() &&  nSize != regular_count) seq.push_back(v);
               }
          }
     }

     if (!relexist2)
          clear_relations(0, 2);

     return seq;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::set_strip_markers()
{
     FaceSequence bound_faces = get_bound_faces(1);

     size_t numfaces = getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          face->setTag(0);
     }

     FaceSequence strip1, strip2;
     int id = 0;
     for (size_t i = 0; i < bound_faces.size(); i++) {
          strip1.clear();
          strip2.clear();
          get_quad_strips(bound_faces[i], strip1, strip2);
          id++;
          for (size_t j = 0; j < strip1.size(); j++)
               strip1[j]->setTag(id);
          /*
           id++;
           for( int j = 0; j < strip2.size(); j++)
           strip2[j]->setTag( id );
           */
     }
}

///////////////////////////////////////////////////////////////////////////////

vector<int>
Mesh::get_topological_statistics(int entity, bool sorted)
{
     int relexist = build_relations(0, 2);

     assert(getAdjTable(0, 2));

     int numnodes = getSize(0);

     vector<int> degree;
     degree.reserve(numnodes);
     for (int i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          if( !v->isRemoved() )
               degree.push_back(v->getNumRelations(2) );
     }

     int mindegree = *min_element(degree.begin(), degree.end());
     int maxdegree = *max_element(degree.begin(), degree.end());

     cout << " Mesh Topological Quality : " << endl;

     cout << " ************************ " << endl;
     cout << " Degree         FaceCount " << endl;
     cout << " ************************ " << endl;

     for (int i = mindegree; i <= maxdegree; i++) {
          int ncount = 0;
          for (size_t j = 0; j < degree.size(); j++)
               if (degree[j] == i)
                    ncount++;
          cout << setw(5) << i << setw(15) << ncount << endl;
     }

     if (!relexist)
          clear_relations(0, 2);

     return degree;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::setNodeWavefront(int layerid)
{
     assert(layerid >= 0);

     int relexist = build_relations(0, 0);

     size_t numNodes = getSize(0);

     size_t ncount = 0;

     if (layerid == 0) {
          search_boundary();
          for (size_t i = 0; i < numNodes; i++) {
               Vertex *vertex = getNodeAt(i);
               if (vertex->isBoundary()) {
                    vertex->setLayerID(0);
                    ncount++;
               }
          }
     } else {
          for (size_t i = 0; i < numNodes; i++) {
               Vertex *vertex = getNodeAt(i);
               if (vertex->getLayerID() == layerid - 1) {
                    NodeSequence &vnodes = vertex->getRelations0();
                    for (size_t j = 0; j < vnodes.size(); j++) {
                         int lid = vnodes[j]->getLayerID();
                         if (lid >= 0 && lid <= layerid - 1) continue;
                         vnodes[j]->setLayerID(layerid);
                         ncount++;
                    }
               }
          }
     }

     if (!relexist)
          clear_relations(0, 0);

     return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::setNodeWavefront()
{
     int relexist = build_relations(0, 0);
     search_boundary();

     size_t numNodes = getSize(0);
     for (size_t i = 0; i < numNodes; i++) {
          Vertex *v = getNodeAt(i);
          v->setLayerID(-1);
          v->setVisitMark(0);
     }

     NodeSequence vertexQ;
     NodeSet  nextQ;
     for (size_t i = 0; i < numNodes; i++) {
          Vertex *v = getNodeAt(i);
          if (v->isBoundary() && !v->isRemoved() ) {
               v->setLayerID(0);
               v->setVisitMark(1);
               vertexQ.push_back(v);
          }
     }

     if (vertexQ.empty())
          cout << "Warning: No boundary detected " << endl;

     int layerid = 1;
     size_t nSize;
     while (!vertexQ.empty()) {
          cout << "Layer " << vertexQ.size() << endl;
          nextQ.clear();
          nSize = vertexQ.size();
          for (size_t j = 0; j < nSize; j++) {
               Vertex *currVertex = vertexQ[j];
               NodeSequence &vneighs = currVertex->getRelations0();
               for (size_t i = 0; i < vneighs.size(); i++) {
                    Vertex *vn = vneighs[i];
                    if (!vn->isVisited()) nextQ.insert(vn);
               }
          }

          nSize = nextQ.size();
          if( nSize == 0) break;

          NodeSet::const_iterator it;
          for( it = nextQ.begin(); it != nextQ.end(); ++it) {
               Vertex *currVertex = *it;
               currVertex->setLayerID(layerid);
               currVertex->setVisitMark(1);
          }

          vertexQ.clear();
          vertexQ.reserve( nSize );
          for( it = nextQ.begin(); it != nextQ.end(); ++it)
               vertexQ.push_back( *it );
          layerid++;
     }

     if (!relexist)
          clear_relations(0, 0);

     cout << "Info: Number of layers in node front : " << layerid - 1 << endl;

     return layerid - 1;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::setFaceWavefront(int layerid)
{
     assert(layerid >= 0);

     int relexist = build_relations(0, 2);

     size_t numfaces = getSize(2);

     size_t ncount = 0;

     FaceSequence vfaces;
     if (layerid == 0) {
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               int nsize = face->getSize(0);
               for (int j = 0; j < nsize; j++) {
                    Vertex *v0 = face->getNodeAt(j + 0);
                    Vertex *v1 = face->getNodeAt(j + 1);
                    Mesh::getRelations112(v0, v1, vfaces);
                    if (vfaces.size() == 1) {
                         face->setLayerID(0);
                         ncount++;
                    }
               }
          }
     } else {
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               if (face->getLayerID() == layerid - 1) {
                    vfaces = face->getRelations212();
                    for (size_t j = 0; j < vfaces.size(); j++) {
                         int lid = vfaces[j]->getLayerID();
                         if (lid >= 0 && lid < layerid) continue;
                         vfaces[j]->setLayerID(layerid);
                         ncount++;
                    }
               }
          }
     }

     if (!relexist)
          clear_relations(0, 2);

     return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::setFaceWavefront()
{
     int relexist = build_relations(0, 2);

     search_boundary();

     size_t numFaces = getSize(2);

     for (size_t i = 0; i < numFaces; i++) {
          Face *f = getFaceAt(i);
          f->setLayerID(0);
          f->setVisitMark(0);
     }

     FaceSequence faceQ, nextQ;
     for (size_t i = 0; i < numFaces; i++) {
          Face *f = getFaceAt(i);
          if( !f->isRemoved() ) {
               f->setLayerID(1);
               if (f->has_boundary_edge()) {
                    f->setLayerID(0);
                    f->setVisitMark(1);
                    faceQ.push_back(f);
               }
          }
     }

     int layerid = 1;
     size_t nSize;

     while (!faceQ.empty()) {
          nextQ.clear();
          nextQ.reserve(faceQ.size());
          nSize = faceQ.size();
          for (size_t j = 0; j < nSize; j++) {
               Face *currFace = faceQ[j];
               FaceSequence neighs = currFace->getRelations212();
               for (size_t i = 0; i < neighs.size(); i++) {
                    Face *nf = neighs[i];
                    if (!nf->isVisited()) nextQ.push_back(nf);
               }
          }

          nSize = nextQ.size();
          for (size_t i = 0; i < nSize; i++) {
               Face *f = nextQ[i];
               f->setLayerID(layerid);
               f->setVisitMark(1);
          }

          layerid++;
          faceQ = nextQ;
     }

     if (!relexist)
          clear_relations(0, 2);

     cout << "Info: Number of layers in face front : " << layerid - 1 << endl;

     return layerid - 1;
}
///////////////////////////////////////////////////////////////////////////////

int
Mesh::setWavefront(int forwhat)
{
     if (forwhat == 0)
          return setNodeWavefront();

     if (forwhat == 2)
          return setFaceWavefront();

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::verify_front_ordering(int mentity)
{
     size_t numfaces = getSize(2);
     if (mentity == 0) {
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               int nsize = face->getSize(0);
               for (int j = 0; j < nsize; j++) {
                    Vertex *v0 = face->getNodeAt((j + 0) % nsize);
                    Vertex *v1 = face->getNodeAt((j + 1) % nsize);
                    int l1 = v0->getLayerID();
                    int l2 = v1->getLayerID();
                    assert(l1 >= 0 && l2 >= 0);
                    if (abs(l2 - l1) > 1) return 1;
               }
          }
     }

     if (mentity == 2) {
          int relexist2 = build_relations(0, 2);
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               FaceSequence neighs = face->getRelations212();
               int l1  = face->getLayerID();
               int nf  = neighs.size();
               assert(l1 >= 0);
               for (int j = 0; j < nf; j++) {
                    int l2 = neighs[j]->getLayerID();
                    assert(l2 >= 0);
                    if (abs(l2 - l1) > 1) return 1;
               }
          }
          if (!relexist2) clear_relations(0, 2);
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::remove_unreferenced_nodes()
{
     size_t numnodes = getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          v->setVisitMark(0);
     }

     size_t numfaces = getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          for (int j = 0; j < face->getSize(0); j++) {
               Vertex *v = face->getNodeAt(j);
               v->setVisitMark(1);
          }
     }

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          if (!v->isVisited()) cout << " Has Unreferenced Node " << endl;
     }

     /*
        for( size_t i = 0; i < numnodes; i++) {
             Vertex *v = getNodeAt(i);
             if( !v->isVisited()) v->setRemoveMark(1);
        }
      */

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

double
Mesh::getSurfaceArea()
{
     double facearea, sumArea = 0.0;

     size_t numfaces = getSize(2);
     double minarea = MAXDOUBLE;
     double maxarea = 0.0;
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          facearea = face->getArea();
          sumArea += facearea;
          if (facearea < minarea) minarea = facearea;
          if (facearea > maxarea) maxarea = facearea;
     }

     /*
        cout << "Info:   Min face Area : " << minarea << endl;
        cout << "Info:   Max face Area : " << maxarea << endl;
      */
     return sumArea;
}
///////////////////////////////////////////////////////////////////////////////

NodeSequence
Mesh::chain_nodes(const vector<Edge> &bndedges)
{
     NodeSequence bndnodes;
     bndnodes.push_back(bndedges[0].getNodeAt(0));
     bndnodes.push_back(bndedges[0].getNodeAt(1));

     size_t nSize = bndedges.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          Vertex *v0 = bndedges[i].getNodeAt(0);
          Vertex *v1 = bndedges[i].getNodeAt(1);
          if (v0 == bndnodes[i]) {
               bndnodes.push_back(v1);
          } else if (v1 == bndnodes[i])
               bndnodes.push_back(v0);
          else {
               cout << "Error in bound edges : " << endl;
               bndnodes.clear();
               return bndnodes;
          }
     }
     return bndnodes;
}

vector<double>
Mesh::getAspectRatio(bool sorted)
{
     size_t numfaces = getSize(2);

     vector<double> quality;
     quality.resize(numfaces);

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          quality[i] = face->getAspectRatio();
     }

     if (sorted)
          std::sort(quality.begin(), quality.end());

     double minval = *min_element(quality.begin(), quality.end());
     double maxval = *max_element(quality.begin(), quality.end());

     cout << "Info: Minimum Aspect Ratio  " << minval << endl;
     cout << "Info: Maximum Aspect Ratio  " << maxval << endl;

     return quality;
}
///////////////////////////////////////////////////////////////////////////////

void
Mesh::emptyAll()
{
     nodes.clear();
     edges.clear();
     faces.clear();
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::deleteNodes()
{
     size_t nSize = nodes.size();
     for (size_t i = 0; i < nSize; i++) {
          if( !nodes[i]->isRemoved() ) delete nodes[i];
     }
     nodes.clear();
}

void Mesh::deleteEdges()
{
     size_t nSize = edges.size();
     for (size_t i = 0; i < nSize; i++) {
          if( !edges[i]->isRemoved() ) delete edges[i];
     }
     edges.clear();
}

void Mesh::deleteFaces()
{
     size_t nSize = faces.size();
     for (size_t i = 0; i < nSize; i++) {
          if( !faces[i]->isRemoved() ) delete faces[i];
     }
     faces.clear();
}

void
Mesh::deleteAll()
{
     for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
               clear_relations(i, j);

     deleteFaces();
     deleteEdges();
     deleteNodes();

     boundary_known = 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Face::refine_quad14(NodeSequence &newnodes, FaceSequence &newfaces)
{
     assert(newnodes.size() == 5);
     Point3D xyz = this->getCentroid();
     newnodes[4] = Vertex::newObject();
     newnodes[4]->setXYZCoords(xyz);

     newfaces.resize(4);
     NodeSequence qC(4);
     for (int i = 0; i < 4; i++)
          newfaces[i] = Face::newObject();

     qC[0] = connect[0];
     qC[1] = newnodes[0];
     qC[2] = newnodes[4];
     qC[3] = newnodes[3];
     newfaces[0]->setNodes(qC);

     qC[0] = newnodes[0];
     qC[1] = connect[1];
     qC[2] = newnodes[1];
     qC[3] = newnodes[4];
     newfaces[1]->setNodes(qC);

     qC[0] = newnodes[1];
     qC[1] = connect[2];
     qC[2] = newnodes[2];
     qC[3] = newnodes[4];
     newfaces[2]->setNodes(qC);

     qC[0] = newnodes[3];
     qC[1] = newnodes[4];
     qC[2] = newnodes[2];
     qC[3] = connect[3];
     newfaces[3]->setNodes(qC);
     return 0;
}
///////////////////////////////////////////////////////////////////////////////

int
Mesh::refine_quads14()
{
     build_edges();

     std::multimap<Vertex*, Edge*> edgemap;
     for (size_t i = 0; i < edges.size(); i++) {
          Vertex *v0 = edges[i]->getNodeAt(0);
          Vertex *v1 = edges[i]->getNodeAt(1);
          const Point3D &xyz = Vertex::mid_point(v0, v1);
          Vertex *vnew = Vertex::newObject();
          vnew->setXYZCoords(xyz);
          edges[i]->addRelation0(vnew);
          Vertex *vhash = edges[i]->getHashNode();
          edgemap.insert(std::make_pair(vhash, edges[i]));
          this->addNode(vnew);
     }

     NodeSequence newnodes(5);
     FaceSequence newfaces(4);

     multimap<Vertex*, Edge*>::const_iterator ilower, iupper, iter;

     size_t numfaces = this->getSize(2);
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = this->getFaceAt(iface);
          for (int j = 0; j < 4; j++) {
               Vertex *v0 = face->getNodeAt((j + 0) % 4);
               Vertex *v1 = face->getNodeAt((j + 1) % 4);
               Edge edge(v0, v1);
               Vertex *vh = edge.getHashNode();
               ilower = edgemap.lower_bound(vh);
               iupper = edgemap.upper_bound(vh);
               Vertex *vmid = NULL;
               for (iter = ilower; iter != iupper; ++iter) {
                    Edge *existing_edge = iter->second;
                    if (existing_edge->isSameAs(edge)) {
                         NodeSequence vneighs = existing_edge->getRelations0();
                         assert(vneighs.size() == 1);
                         vmid = vneighs[0];
                         break;
                    }
               }
               assert(vmid);
               newnodes[j] = vmid;
          }
          face->refine_quad14(newnodes, newfaces);
          this->addNode(newnodes[4]); // Only the center node, others already registered
          for (int j = 0; j < 4; j++)
               this->addFace(newfaces[j]);
          face->setStatus(MeshEntity::REMOVE);
     }

     prune();
     enumerate(0);
     enumerate(2);
     collect_garbage();

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

int
Mesh::refine_quads15()
{
     FaceSequence newfaces;
     NodeSequence newnodes;

     int nSize;
     size_t numfaces = getSize(2);
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = getFaceAt(iface);
          face->refine_quad15(newnodes, newfaces);

          nSize = newnodes.size();
          for (int  j = 0; j < nSize; j++)
               addNode(newnodes[j]);

          nSize = newfaces.size();
          for (int j = 0; j < nSize; j++)
               addFace(newfaces[j]);

          face->setStatus(MeshEntity::REMOVE);
     }

     prune();
     enumerate(0);
     enumerate(2);
     collect_garbage();

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

/*
int
Mesh::get_quality_statistics(const string &fname)
{
    ofstream ofile(fname.c_str(), ios::out);
    if (ofile.fail())
    {
        cout << "Warning: Cann't open file " << fname << endl;
        return 1;
    }

    // Collect Element Area informtion
    size_t numfaces = getSize(2);
    vector<double> quality(numfaces);

    double minval, maxval, stddev, avgval, medianval;

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        quality[i] = face->getArea();
    }

    sort(quality.begin(), quality.end());
    minval = quality.front();
    maxval = quality.back();

    ofile << "# ********************************************* " << endl;
    ofile << "# Measure            Element Area        " << endl;
    ofile << "# Num of Elements  " << numfaces << endl;
    ofile << "# Min              " << minval << endl;
    ofile << "# Max              " << maxval << endl;
    ofile << "# Average          " << avgval << endl;
    ofile << "# MeanVal          " << medianval << endl;
    //ofile <<  "# StdDev           " <<   stddev     << endl;
    for (size_t i = 0; i < numfaces; i++)
        ofile << quality[i] << endl;
    ofile << " ********************************************* " << endl;


    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        quality[i] = face->getAspectRatio();
    }

    sort(quality.begin(), quality.end());
    minval = quality.front();
    maxval = quality.back();

    ofile << "# Measure            Aspect Ratio      " << endl;
    ofile << "# Num of Elements  " << numfaces << endl;
    ofile << "# Min              " << minval << endl;
    ofile << "# Max              " << maxval << endl;
    ofile << "# Average          " << avgval << endl;
    ofile << "# MeanVal          " << medianval << endl;
    ofile << "# StdDev           " << stddev << endl;

    for (size_t i = 0; i < numfaces; i++)
        ofile << quality[i] << endl;
    ofile << " ********************************************* " << endl;
}
 */

///////////////////////////////////////////////////////////////////////////////
double Vertex::getFeatureLength() const
{
     if (!isBoundary()) return MAXDOUBLE;

     NodeSequence vneighs = getRelations0();

     assert(!vneighs.empty());

     double minlen = MAXDOUBLE;

     int  nSize = vneighs.size();
     for (int j = 0; j < nSize; j++) {
          if (vneighs[j]->isBoundary()) {
               const Point3D &p0 = getXYZCoords();
               const Point3D &p1 = vneighs[j]->getXYZCoords();
               minlen = min(minlen, Math::length(p0, p1));
          }
     }
     return minlen;
}

///////////////////////////////////////////////////////////////////////////////

/*
void
Mesh::setFeatureLength()
{
    int relexist0 = build_relations(0, 0);
    int relexist2 = build_relations(0, 2);

    search_boundary();

    size_t numnodes = getSize(0);
    vector<Vertex*> vneighs, bndnodes, sidenodes;

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = getNodeAt(i);
        if (vertex->isBoundary())
        {
            vneighs = vertex->getRelations0();
            double minlen = MAXDOUBLE;
            for (int j = 0; j < vneighs.size(); j++)
            {
                if (vneighs[j]->isBoundary())
                {
                    Point3D p0 = vertex->getXYZCoords();
                    Point3D p1 = vneighs[j]->getXYZCoords();
                    minlen = min(minlen, Math::length(p0, p1));
                }
            }
            vertex->setFeatureLength(minlen);
        }
    }

    if (!relexist0)
        clear_relations(0, 0);

    if (!relexist2)
        clear_relations(0, 2);
}
 */

///////////////////////////////////////////////////////////////////////////////

int
Mesh::hasDuplicates(int what)
{
     map<Vertex*, FaceSequence> mapfaces;

     size_t numfaces = getSize(2);
     assert(numfaces);

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          Vertex *vertex = face->getHashNode();
          mapfaces[vertex].push_back(face);
     }

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = getFaceAt(i);
          const FaceSequence &hashfaces = mapfaces[face->getHashNode() ];
          size_t ncount = 0;
          for (size_t j = 0; j < hashfaces.size(); j++)
               if (hashfaces[j]->isSameAs(face)) ncount++;
          if (ncount != 1) {
               cout << "Warning: Mesh has some duplicate faces " << endl;
               return 1;
          }
     }
     return 0;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::count_concave_faces()
{
     size_t numfaces = getSize(2);
     size_t ncount = 0;
     for (size_t i = 0; i < numfaces; i++) {
          Face *f = getFaceAt(i);
          if(!f->isRemoved()  &&  !f->isConvex() ) ncount++;
     }
     return ncount;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::count_inverted_faces()
{
     size_t numfaces = getSize(2);
     size_t ncount = 0;
     for (size_t i = 0; i < numfaces; i++) {
          Face *f = getFaceAt(i);
          if( !f->isRemoved()  && (f->concaveAt() >= 0)) ncount++;
     }
     return ncount;
}
///////////////////////////////////////////////////////////////////////////////

void
Mesh::normalize()
{
     size_t numnodes = nodes.size();

     double xmin, xmax, ymin, ymax, zmin, zmax;
     Point3D xyz;
     xyz = nodes[0]->getXYZCoords();

     xmin = xyz[0];
     xmax = xyz[0];
     ymin = xyz[1];
     ymax = xyz[1];
     zmin = xyz[2];
     zmax = xyz[2];

     for (size_t i = 0; i < numnodes; i++) {
          xyz = nodes[i]->getXYZCoords();
          xmin = min(xmin, xyz[0]);
          xmax = max(xmax, xyz[0]);
          ymin = min(ymin, xyz[1]);
          ymax = max(ymax, xyz[1]);
          zmin = min(zmin, xyz[2]);
          zmax = max(zmax, xyz[2]);
     }

     double xlen = fabs(xmax - xmin);
     double ylen = fabs(ymax - ymin);
     double zlen = fabs(zmax - zmin);
     double scale = max(max(xlen, ylen), zlen);

     for (size_t i = 0; i < numnodes; i++) {
          xyz = nodes[i]->getXYZCoords();
          xyz[0] /= scale;
          xyz[1] /= scale;
          xyz[2] /= scale;
          nodes[i]->setXYZCoords(xyz);
     }
}

///////////////////////////////////////////////////////////////////////////////

BoundingBox
Mesh::getBoundingBox() const
{
     BoundingBox box;

     size_t numnodes = nodes.size();

     double xmin, xmax, ymin, ymax, zmin, zmax;
     Point3D xyz;
     xyz = nodes[0]->getXYZCoords();

     xmin = xyz[0];
     xmax = xyz[0];
     ymin = xyz[1];
     ymax = xyz[1];
     zmin = xyz[2];
     zmax = xyz[2];

     for (size_t i = 0; i < numnodes; i++) {
          xyz = nodes[i]->getXYZCoords();
          xmin = min(xmin, xyz[0]);
          xmax = max(xmax, xyz[0]);
          ymin = min(ymin, xyz[1]);
          ymax = max(ymax, xyz[1]);
          zmin = min(zmin, xyz[2]);
          zmax = max(zmax, xyz[2]);
     }

     xyz[0] = xmin;
     xyz[1] = ymin;
     xyz[2] = zmin;
     box.setLowerLeftCorner(xyz);

     xyz[0] = xmax;
     xyz[1] = ymax;
     xyz[2] = zmax;
     box.setUpperRightCorner(xyz);

     return box;
}

double
Mesh::getLength(int dir) const
{
     size_t numnodes = nodes.size();

     double xmin, xmax, ymin, ymax, zmin, zmax;
     Point3D xyz;
     xyz = nodes[0]->getXYZCoords();

     xmin = xyz[0];
     xmax = xyz[0];
     ymin = xyz[1];
     ymax = xyz[1];
     zmin = xyz[2];
     zmax = xyz[2];

     for (size_t i = 0; i < numnodes; i++) {
          xyz = nodes[i]->getXYZCoords();
          xmin = min(xmin, xyz[0]);
          xmax = max(xmax, xyz[0]);
          ymin = min(ymin, xyz[1]);
          ymax = max(ymax, xyz[1]);
          zmin = min(zmin, xyz[2]);
          zmax = max(zmax, xyz[2]);
     }

     double xlen = fabs(xmax - xmin);
     double ylen = fabs(ymax - ymin);
     double zlen = fabs(zmax - zmin);

     if (dir == 0) return xlen;
     if (dir == 1) return ylen;
     if (dir == 2) return zlen;

     cout << "Warning: Invalid direction : " << dir << endl;
     return 0.0;
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::linear_interpolation(Mesh *mesh, Vertex *v0, Vertex *v1, int n, NodeSequence &newnodes)
{
     assert(n >= 2);

     const Point3D &xyz0 = v0->getXYZCoords();
     const Point3D &xyz1 = v1->getXYZCoords();

     NodeSequence poolnodes;
     mesh->objects_from_pool(n-2, poolnodes);

     newnodes.resize(n);

     newnodes[0] = v0;
     newnodes[n - 1] = v1;

     double dt = 2.0 / (double) (n - 1);

     Point3D xyzt;
     int index = 0;
     for (int i = 1; i < n - 1; i++) {
          double t = -1.0 + i*dt;
          xyzt[0] = TFI::linear_interpolation(t, xyz0[0], xyz1[0]);
          xyzt[1] = TFI::linear_interpolation(t, xyz0[1], xyz1[1]);
          xyzt[2] = TFI::linear_interpolation(t, xyz0[2], xyz1[2]);
          newnodes[i] = poolnodes[index++];
          newnodes[i]->setStatus( MeshEntity::ACTIVE);
          newnodes[i]->setXYZCoords(xyzt);
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_tfi_coords(int i, int j, int nx, int ny, vector<Vertex*> &qnodes)
{
     int offset;

     offset = 0;
     const Point3D &v00 = qnodes[offset]->getXYZCoords();

     offset = i;
     const Point3D &vr0 = qnodes[offset]->getXYZCoords();

     offset = (nx - 1);
     const Point3D &v10 = qnodes[offset]->getXYZCoords();

     offset = j*nx;
     const Point3D &v0s = qnodes[offset]->getXYZCoords();

     offset = j * nx + (nx - 1);
     const Point3D &v1s = qnodes[offset]->getXYZCoords();

     offset = (ny - 1) * nx;
     const Point3D &v01 = qnodes[offset]->getXYZCoords();

     offset = (ny - 1) * nx + i;
     const Point3D &vr1 = qnodes[offset]->getXYZCoords();

     offset = (ny - 1) * nx + (nx - 1);
     const Point3D &v11 = qnodes[offset]->getXYZCoords();

     Point3D vrs;

     double dr = 2.0 / (double) (nx - 1);
     double ds = 2.0 / (double) (ny - 1);

     double r = -1.0 + i*dr;
     double s = -1.0 + j*ds;
     for (int k = 0; k < 3; k++) {
          vrs[k] = TFI::transfinite_blend(r, s,
                                          v00[k], v10[k], v11[k], v01[k],
                                          vr0[k], v1s[k], vr1[k], v0s[k]);
     }
     offset = j * nx + i;
     qnodes[offset]->setXYZCoords(vrs);

}
///////////////////////////////////////////////////////////////////////////////

int
QTrack::advance_single_step(int endat)
{
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

     NodeSet vset;

     size_t index = sequence.size();
     v0 = sequence[index - 2];
     v1 = sequence[index - 1];
     v0->setVisitMark(1);

     if (endat == END_AT_CROSSINGS && v1->isVisited()) return 0;
     if (v1->isBoundary()) return 0;

     v1->setVisitMark(1);
     const NodeSequence &vneighs = v1->getRelations0();
     if (vneighs.size() != 4) return 0;

     FaceSequence adjFaces;
     Mesh::getRelations112(v0, v1, adjFaces);
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

///////////////////////////////////////////////////////////////////////////////

void
QTrack::advance(int endat)
{
     assert(sequence.size() == 2);

     // Starting node is always irregular ...
     assert( sequence[0]->getNumRelations(2) != 4);

     while (1) {
          int progress = advance_single_step(endat);
          if (!progress) break;
     }

#ifdef DEBUG
     Vertex *endvertex;
     // Sanity Checking ....
     if (endat == END_AT_TERMINALS) {
          endvertex = sequence.front();
          if (!endvertex->isBoundary()) {
               assert( endvertex->getNumRelations(2) != 4);
          }

          endvertex = sequence.back();
          if (!endvertex->isBoundary()) {
               assert( endvertex->getNumRelations(2) != 4);
          }

          for (size_t i = 1; i < sequence.size() - 1; i++) {
               assert(!sequence[i]->isBoundary());
               assert( endvertex->getNumRelations(2) == 4);
          }
     }
#endif

}

///////////////////////////////////////////////////////////////////////////////

vector<QTrack> Jaal::generate_quad_partitioning(Mesh *mesh)
{
     vector<QTrack> qpath;

     int nTopo = mesh->isHomogeneous();
     if (nTopo != 4) {
          cout << "Error: The mesh must be all Quads " << endl;
          return qpath;
     }

     int relexist2 = mesh->build_relations(0, 2);
     int relexist0 = mesh->build_relations(0, 0);

     mesh->search_boundary();

     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          vertex->setGroupID(0);
          vertex->setVisitMark(0);
     }

     QTrack qp;
     qp.mesh = mesh;

     int found;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if (!vertex->isBoundary() && !vertex->isRemoved() ) {
               const NodeSequence &vnodes = vertex->getRelations0();
               int numneighs = vnodes.size();
               if(numneighs != 4) {
                    for (int j = 0; j < numneighs; j++) {
                         qp.sequence.resize(2); // As we know the starting edge
                         qp.sequence[0] = vertex;
                         qp.sequence[1] = vnodes[j];
                         qp.advance(0); // Terminate at irregular nodes only..
                         found = 0;
                         for (size_t k = 0; k < qpath.size(); k++) {
                              if (qpath[k] == qp) {
                                   found = 1;
                                   break;
                              }
                         }
                         if (!found) qpath.push_back(qp);
                    }
               }
          }
     }
     sort(qpath.begin(), qpath.end());

     int partid = 1;
     for( size_t i = 0; i < qpath.size(); i++) {
          for( size_t j = 0; j < qpath[i].sequence.size(); j++)
               qpath[i].sequence[j]->setGroupID( partid );
          partid++;
     }

     size_t nSize = mesh->getSize(2);
     for( size_t i = 0; i < nSize; i++) {
          Face *f = mesh->getFaceAt(i);
          f->setGroupID(0);
     }

     deque<Face*> faceQ;
     map<int, int> faceCount;
     FaceSequence fneighs;
     partid = 1;
     while(1) {
          Face *seedface = NULL;
          for( size_t i = 0; i < nSize; i++) {
               Face *f = mesh->getFaceAt(i);
               if( !f->isRemoved()  && f->getGroupID() == 0) {
                    seedface = f;
                    break;
               }
          }
          if( seedface == NULL ) break;
          faceCount[partid] = 0;

          faceQ.clear();
          faceQ.push_back(seedface);

          while(!faceQ.empty() ) {

               Face *currface = faceQ.front();
               faceQ.pop_front();

               if( !currface->isRemoved() && currface->getGroupID() == 0 ) {
                    currface->setGroupID( partid );
                    faceCount[partid]++;
                    for( int i = 0; i < 4; i++) {
                         Vertex *v0 = currface->getNodeAt(i);
                         Vertex *v1 = currface->getNodeAt(i+1);
                         if( v0->getGroupID() == 0 || v1->getGroupID() == 0) {
                              Mesh::getRelations112( v0, v1, fneighs);
                              for( size_t j = 0; j < fneighs.size(); j++)
                                   if( fneighs[j]->getGroupID() == 0)
                                        faceQ.push_back( fneighs[j] );
                         }
                    }
               }
          }
          partid++;
     }

     cout << "#of Partitions : " << faceCount.size() << endl;
     map<int,int>::const_iterator it;
     for( it = faceCount.begin(); it != faceCount.end(); ++it)
          cout << it->first << " " << it->second << endl;


     if (!relexist2) mesh->clear_relations(0, 2);
     if (!relexist0) mesh->clear_relations(0, 0);

     return qpath;
}

////////////////////////////////////////////////////////////////////////////////

void
Jaal::set_layer_tag(Mesh *mesh)
{
     assert(1);
     //  mesh->setWavefront (0);

     size_t numnodes = mesh->getSize(0);
     int max_tag_val = -1;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          max_tag_val = max(max_tag_val, vertex->getLayerID());
     }

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          int lid = vertex->getLayerID();
          if (lid < 0)
               vertex->setTag(max_tag_val + 1);
          else
               vertex->setTag(lid);
     }

     //  mesh->setWavefront (2);
     size_t numfaces = mesh->getSize(2);

     max_tag_val = -1;
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          max_tag_val = max(max_tag_val, face->getLayerID());
     }

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          int lid = face->getLayerID();
          if (lid < 0)
               face->setTag(max_tag_val);
          else
               face->setTag(lid);
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_convexity_tag(Mesh *mesh)
{
     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          vertex->setTag(1);
     }

     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setTag(0);
          if (face->isConvex())
               face->setTag(1);
          else {
               for (int j = 0; j < face->getSize(0); j++) {
                    Vertex *v = face->getNodeAt(j);
                    v->setTag(0);
               }
          }

     }
}

void
Jaal::set_ideal_node_tag(Mesh *mesh, int elem)
{
     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          if( !v->isRemoved() ) {
               int ideal_degree = v->get_ideal_face_degree(elem);
               int curr_degree  = v->getNumRelations(2);
               v->setTag(1);
               if( curr_degree != ideal_degree ) v->setTag(0);
          }
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_boundary_tag(Mesh *mesh)
{
     int relexist = mesh->build_relations(0, 2);

     if (!mesh->isBoundaryKnown())
          mesh->search_boundary();

     size_t numnodes = mesh->getSize(0);

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if (vertex->isBoundary())
               vertex->setTag(2);
          else
               vertex->setTag(1);
     }

     size_t numfaces = mesh->getSize(2);

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setTag(1);
     }

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if (face->has_boundary_edge())
               face->setTag(0);
     }

     if (!relexist) mesh->clear_relations(0, 2);

}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_constrained_tag(Mesh *mesh)
{
     size_t numnodes = mesh->getSize(0);

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if (vertex->isConstrained())
               vertex->setTag(0);
          else
               vertex->setTag(1);
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_large_area_tag(Mesh *mesh)
{
     size_t numfaces = mesh->getSize(2);

     double sumarea = 0.0;
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          sumarea += face->getArea();
     }
     double avgarea = sumarea / (double) numfaces;
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if (face->getArea() > avgarea)
               face->setTag(0);
          else
               face->setTag(1);
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_tiny_area_tag(Mesh *mesh, double tolarea)
{
     size_t numnodes = mesh->getSize(0);
     size_t numfaces = mesh->getSize(2);

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          vertex->setTag(1);
     }

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if (face->getArea() < tolarea) {
               face->setTag(0);
               for (int j = 0; j < face->getSize(0); j++) {
                    Vertex *vertex = face->getNodeAt(j);
                    vertex->setTag(0);
               }
          } else
               face->setTag(1);
     }
}

///////////////////////////////////////////////////////////////////////////////

bool
Face::has_all_bound_nodes() const
{
     for (int i = 0; i < getSize(0); i++) {
          Vertex *v = getNodeAt(i);
          if (!v->isBoundary()) return 0;
     }
     return 1;
}
///////////////////////////////////////////////////////////////////////////////

void
Mesh::setFacesNormal()
{
     double x[100], y[100], z[100], nx, ny, nz;
     Point3D xyz;
     Vec3D normal;

     for (size_t i = 0; i < faces.size(); i++) {
          int nsize = faces[i]->getSize(0);
          for (int j = 0; j < nsize; j++) {
               Vertex *vtx = faces[i]->getNodeAt(j);
               xyz = vtx->getXYZCoords();
               x[j] = xyz[0];
               y[j] = xyz[1];
               z[j] = xyz[2];
          }
          PolygonNormal3D(nsize, x, y, z, &nx, &ny, &nz);
          normal[0] = nx;
          normal[1] = ny;
          normal[2] = nz;
          faces[i]->setNormal(normal);
     }
}

///////////////////////////////////////////////////////////////////////////////
void Mesh :: objects_from_pool( size_t n, vector<Vertex*> &objects)
{
     objects.clear();

     if( n == 0) return;

     objects.reserve( n );

     size_t ncount = 0;
     while( !garbageNodes.empty() ) {
          Vertex *v = garbageNodes.front();
          garbageNodes.pop_front();
          if( v->isRemoved() ) {
               v->setStatus( MeshEntity::INACTIVE);
               objects.push_back(v);
               ncount++;
               if( ncount == n) break;
          }
     }

     for( size_t i = ncount; i < n; i++) {
          Vertex *v = Vertex::newObject();
          v->setStatus( MeshEntity::INACTIVE);
          objects.push_back(v);
          addNode( v );
     }

     assert( objects.size() == n );
}

///////////////////////////////////////////////////////////////////////////////

void Mesh :: objects_from_pool( size_t n, vector<Face*> &objects)
{
     objects.clear();

     if( n == 0) return;

     objects.reserve( n );

     size_t ncount = 0;
     while( !garbageFaces.empty() ) {
          Face *f = garbageFaces.front();
          garbageFaces.pop_front();
          if( f->isRemoved() ) {
               f->setStatus( MeshEntity::INACTIVE);
               objects.push_back(f);
               ncount++;
               if( ncount == n) break;
          }
     }

     for( size_t i = ncount; i < n; i++) {
          Face *f = Face::newObject();
          f->setStatus( MeshEntity::INACTIVE);
          objects.push_back(f);
          addFace(f);
     }

     assert( objects.size() == n );
}
///////////////////////////////////////////////////////////////////////////////

NodeSequence Mesh::get_breadth_first_ordered_nodes(Vertex *vstart, MeshFilter *filter)
{
     assert(vstart != NULL);

     NodeSequence seq;

     int relexist0 = build_relations(0, 0);

     size_t numnodes = getSize(0);
     if (numnodes == 0) return seq;

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          v->setVisitMark(0);
          v->setLayerID(0);
     }

     if (vstart == 0) vstart = getNodeAt(0);

#ifdef SEQUENCE_IS_VECTOR
     seq.reserve(numnodes);
#endif

     list<Vertex*> vertexQ;
     vertexQ.push_back(vstart);
     NodeSequence vneighs;

     int proceed = 1;
     while (!vertexQ.empty()) {
          Vertex *curr_vertex = vertexQ.front();
          vertexQ.pop_front();
          int currlevel = curr_vertex->getLayerID();
          if (filter) {
               if (curr_vertex != vstart) proceed = filter->pass(curr_vertex);
          }
          if (!curr_vertex->isVisited()) {
               seq.push_back(curr_vertex);
               if (!proceed) break;
               curr_vertex->setVisitMark(1);
               vneighs = curr_vertex->getRelations0();
               for (size_t i = 0; i < vneighs.size(); i++) {
                    if (!vneighs[i]->isVisited()) {
                         vertexQ.push_back(vneighs[i]);
                         vneighs[i]->setLayerID(currlevel + 1);
                    }
               }
          }
     }

     if (!relexist0)
          clear_relations(0, 0);

     // Free unused memory in sequence...
     if (!seq.empty()) NodeSequence(seq).swap(seq);

     return seq;
}
///////////////////////////////////////////////////////////////////////////////

NodeSequence Mesh::get_depth_first_ordered_nodes(Vertex *vstart, MeshFilter *filter)
{
     int relexist0 = build_relations(0, 0);

     size_t numnodes = getSize(0);
     list<Vertex*> vertexQ;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          v->setVisitMark(0);
     }

     NodeSequence seq;

     if (vstart == 0) vstart = getNodeAt(0);
     vertexQ.push_back(vstart);
     NodeSequence vneighs;

     while (!vertexQ.empty()) {
          Vertex *curr_vertex = vertexQ.front();
          vertexQ.pop_front();
          if (!curr_vertex->isVisited()) {
               seq.push_back(curr_vertex);
               curr_vertex->setVisitMark(1);
               vneighs = curr_vertex->getRelations0();
               for (size_t i = 0; i < vneighs.size(); i++) {
                    if (!vneighs[i]->isVisited())
                         vertexQ.push_front(vneighs[i]);
               }
          }
     }
#ifdef DEBUG
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = getNodeAt(i);
          assert(v->isVisited());
     }
     assert(seq.size() == numnodes);
#endif

     if (!relexist0)
          clear_relations(0, 0);

     return seq;
}

///////////////////////////////////////////////////////////////////////////////

EdgeSequence Mesh::get_sharp_edges(double creaseAngle)
{
     EdgeSequence sharp_edges;
     /*
        if( edges.empty() ) build_edges();

        FaceSequence efaces;
        for( size_t i = 0; i < edges.size(): i++) {
             PEdge edge = getEdge(i);
             Vertex *v1 = edge->getNodeAt(0);
             Vertex *v1 = edge->getNodeAt(1);
             efaces = Mesh::getRelations112(v0, v1);
             if( efaces.size() == 2 ) {
                 Vec3D f1normal = efaces[0]->getNormal();
                 Vec3D f2normal = efaces[1]->getNormal();
                 double angle = Math::getAngle(fn1, fn2, ANGLE_IN_DEGREES);
                 if (angle <= 90 && angle >= creaseAngle)
                     sharp_edges.push_back(edge->getClone() );
                 else if (angle >= 90 && fabs(180 - angle) >= creaseAngle)
                     sharp_edges.push_back(edge->getClone() );
             }
         }
      */
     return sharp_edges;
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_no_tags(Mesh *mesh)
{
     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          vertex->setTag(0);
     }

     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setTag(0);
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_visit_tags(Mesh *mesh)
{
     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if (vertex->isVisited())
               vertex->setTag(0);
          else
               vertex->setTag(1);
     }

     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if (face->isVisited())
               face->setTag(2);
          else
               face->setTag(3);
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_allboundnodes_tag(Mesh *mesh)
{
     mesh->search_boundary();
     size_t numfaces = mesh->getSize(2);

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setTag(1);
          if (face->has_all_bound_nodes())
               face->setTag(0);
     }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_partition_tag(Mesh *mesh)
{
     size_t numfaces = mesh->getSize(2);

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setTag(face->getPartID());
     }
}


///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_bound1node_tag(Mesh *mesh)
{
     int rel0exist = mesh->build_relations(0, 0);

     mesh->search_boundary();

     size_t numnodes = mesh->getSize(0);

     NodeSequence neighs;

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          vertex->setTag(1);
          if (!vertex->isBoundary()) {
               neighs = vertex->getRelations0();
               int ncount = 0;
               for (size_t j = 0; j < neighs.size(); j++)
                    if (neighs[j]->isBoundary()) ncount++;
               if (ncount > 1)
                    vertex->setTag(0);
          }
     }

     if (!rel0exist)
          mesh->clear_relations(0, 0);
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_inverted_tag(Mesh *mesh)
{
     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          v->setTag(1);
     }

     size_t ncount = 0;
     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setTag(1);
          int pos = face->concaveAt();
          if (pos >= 0) {
               face->setTag(0);
               Vertex *v = face->getNodeAt(pos);
               v->setTag(0);
               ncount++;
          }
     }

}

////////////////////////////////////////////////////////////////////////////////

void
Jaal::set_irregular_path_tag(Mesh *mesh, vector<QTrack> &qpath)
{
     size_t numnodes = mesh->getSize(0);
     Vertex *v;
     for (size_t i = 0; i < numnodes; i++) {
          v = mesh->getNodeAt(i);
          v->setTag(0);
     }

     for (size_t i = 0; i < qpath.size(); i++) {
          v = qpath[i].sequence.front();
          v->setTag(1);
          v = qpath[i].sequence.back();
          v->setTag(1);

          for (size_t j = 1; j < qpath[i].sequence.size() - 1; j++) {
               v = qpath[i].sequence[j];
               cout << "CSV " << v->getID() << endl;
               v->setTag(2);
          }
     }
}

////////////////////////////////////////////////////////////////////////////////

int Jaal::SurfPatch::getPosOf(const Vertex *v)
{
     for (size_t i = 0; i < bound_nodes.size(); i++)
          if (bound_nodes[i] == v) return i;

     cout << "Error: Vertex not found on the boundary " << endl;
     exit(0);

     return -1;
}

////////////////////////////////////////////////////////////////////////////////

NodeSequence SurfPatch::get_bound_nodes(const Vertex *src, const Vertex *dst)
{
     int start_pos = getPosOf(src);
     int end_pos = getPosOf(dst);
     int nsize = bound_nodes.size();

     if (end_pos == 0) end_pos = nsize;
     assert(end_pos > start_pos);

     NodeSequence seq(end_pos - start_pos + 1);
     int index = 0;
     for (int i = start_pos; i <= end_pos; i++)
          seq[index++] = bound_nodes[i % nsize];

     return seq;
}

////////////////////////////////////////////////////////////////////////////////

int Jaal::SurfPatch::search_boundary()
{
     corners.clear();
     boundary.clear();
     Vertex *vertex;

     assert(!faces.empty());

     // We need to rebuild relations locally to identfy corners and boundary.
     FaceSet::const_iterator fiter;
     std::map<Vertex*, FaceSet> relations02;

     for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
          Face *face = *fiter;
          for (int j = 0; j < face->getSize(0); j++) {
               vertex = face->getNodeAt(j);
               relations02[vertex].insert(face);
          }
     }

     // A boundary edge must have exactly one face neighbor...
     FaceSequence faceneighs;
     for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
          Face *face = *fiter;
          int nnodes = face->getSize(0);
          for (int j = 0; j < nnodes; j++) {
               Vertex *v0 = face->getNodeAt((j + 0) % nnodes);
               Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
               faceneighs.clear();
               assert(relations02[v0].size() > 0);
               assert(relations02[v1].size() > 0);
               set_intersection(relations02[v0].begin(), relations02[v0].end(),
                                relations02[v1].begin(), relations02[v1].end(),
                                back_inserter(faceneighs));
               if (faceneighs.size() == 1) {
                    Edge newedge(v0, v1);
                    boundary.push_back(newedge);
               }
          }
     }

     // Sequence the chain and start from one of the corner...
     int err = Mesh::make_chain(boundary);
     if (err) return 2;

     //
     // Identify corners in the mesh.
     // Should we check only one vertex per edge ?
     //
     FaceSet neighs;
     int boundSize = boundary.size();
     for (int k = 0; k < boundSize; k++) {
          vertex = boundary[k].getNodeAt(0);
          neighs = relations02[vertex];
          if (neighs.size() == 1) corners.insert(vertex);

          vertex = boundary[k].getNodeAt(1);
          neighs = relations02[vertex];
          if (neighs.size() == 1) corners.insert(vertex);
     }

     // Start the chain from one of the corners.
     err = Mesh::rotate_chain(boundary, *corners.begin());
     if (err) return 3;

     // Collect the sequence of boundary nodes...,
     bound_nodes.resize( boundSize );
     for (int k = 0; k < boundSize; k++) {
          bound_nodes[k] = boundary[k].getNodeAt(0); // Only the first node.
     }

     //
     // Collect the inner nodes of the patch. These nodes will be deleted, if
     // the remesh operation is successful...
     //
     inner_nodes.clear();
     for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
          Face *face = *fiter;
          int nnodes = face->getSize(0);
          for (int j = 0; j < nnodes; j++) {
               Vertex *v = face->getNodeAt(j);
               if (find(bound_nodes.begin(), bound_nodes.end(), v) == bound_nodes.end())
                    inner_nodes.insert(v);
          }
     }

     // Split the boundary loop into segments.
     // (i.e. End of the segments are the corners identified earlier )
     set_boundary_segments();

     return 0;
}

////////////////////////////////////////////////////////////////////

void Jaal::SurfPatch::set_boundary_segments()
{
     // Although this stage will not come in this algorithm...
     if (corners.size() == 0) return;

     cornerPos.resize(corners.size() + 1);

     NodeSet::const_iterator it;
     int index = 0;
     for (it = corners.begin(); it != corners.end(); ++it) {
          cornerPos[index++] = getPosOf(*it);
     }

     cornerPos[corners.size()] = bound_nodes.size();

     sort(cornerPos.begin(), cornerPos.end());

     segSize.resize(corners.size());

     for (size_t i = 0; i < corners.size(); i++)
          segSize[i] = cornerPos[(i + 1)] - cornerPos[i] + 1;
}

////////////////////////////////////////////////////////////////////

int Jaal::SurfPatch::reorient_4_sided_loop()
{
     //Always remeshable, nothing has to be done...
     if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3])) return 0;

     //////////////////////////////////////////////////////////////////////////
     // Defination:  A four sided convex loop has four segments.
     // Objectives:  A four sided convex loop must be orietned such that
     //   1.  First segment must be smaller than the third one, because
     //       we need to create a triangle patch based at segment#3.
     //
     //   2.  If there are two choices, then the side having irregular
     //       node must be given higher priority. ( Here irregulaty means
     //       that vertex valency < 4 ).
     //
     //   3.  A side having less number of nodes on the first segment than
     //       the third is given preference.
     //
     // Pre-Conditions  :  A loop must be oriented ( CW or CCW ).
     //
     // Date: 17th Nov. 2010.
     //////////////////////////////////////////////////////////////////////////

     Vertex *start_corner = NULL;

     if (segSize[0] == segSize[2]) {
          if (min(segSize[1], segSize[3]) == 2) return 1;
          //  Either Segment 2 or 3 must be starting node
          if (segSize[1] < segSize[3])
               start_corner = bound_nodes[ cornerPos[1] ];
          else
               start_corner = bound_nodes[ cornerPos[3] ];
          start_boundary_loop_from(start_corner);
     }

     if (min(segSize[0], segSize[2]) == 2) return 1;

     if (segSize[2] < segSize[0]) {
          start_corner = bound_nodes[ cornerPos[2] ];
          start_boundary_loop_from(start_corner);
     }

     // By this stage, the loop must be reoriented correctly.
     assert(segSize[0] < segSize[2]);

     cout << " Careful to change this code " << endl;
     // Great, we found one irregular node on the first boundary...
     //   if( has_irregular_node_on_first_segment() ) return 1;

     // If the segment 2 and 3 have same size, Alas, nothing can be done.
     if (segSize[1] == segSize[3]) return 1;

     if (min(segSize[1], segSize[3]) == 2) return 1;

     if (segSize[3] < segSize[1]) {
          start_corner = bound_nodes[ cornerPos[3] ];
          start_boundary_loop_from(start_corner);
     } else {
          start_corner = bound_nodes[ cornerPos[1] ];
          start_boundary_loop_from(start_corner);
     }

     //
     // Note that we didn't check for irregular node here. So if this segment
     // has at least one irregular node, then we are lucky. Otherwise decision
     // to remesh it done based wthere remeshing will result in the reduction
     // of irregular nodes in patch.
     //
     assert(segSize[0] < segSize[2]);
     return 0;
}
////////////////////////////////////////////////////////////////////////////////

void Jaal::SurfPatch::start_boundary_loop_from(Vertex *vmid)
{
     assert(corners.find(vmid) != corners.end());

     NodeSequence::iterator middle;
     middle = find(bound_nodes.begin(), bound_nodes.end(), vmid);
     assert(middle != bound_nodes.end());

     std::rotate(bound_nodes.begin(), middle, bound_nodes.end());
     assert(bound_nodes[0] == vmid);

     set_boundary_segments();
}
///////////////////////////////////////////////////////////////////////////////

Mesh * Jaal::create_structured_mesh(double *origin, double *length,
                                    int *grid_dim, int space_dim)
{
     if (space_dim < 2 || space_dim > 3) return NULL;
     Mesh *mesh = new Mesh;

     int nx = grid_dim[0];
     int ny = grid_dim[1];
     int nz = 1;
     if (space_dim == 3) nz = grid_dim[2];

     double xorg, yorg, zorg = 0.0, dx, dy, dz = 0.0;

     xorg = origin[0];
     yorg = origin[1];

     if (space_dim == 3) zorg = origin[2];

     dx = length[0] / (double) (nx - 1);
     dy = length[1] / (double) (ny - 1);
     if (nz > 2) dz = length[2] / (double) (nz - 1);


     vector<Vertex*> nodes(nx * ny * nz);
     int offset;
     Point3D xyz;
     for (int k = 0; k < nz; k++) {
          for (int j = 0; j < ny; j++) {
               for (int i = 0; i < nx; i++) {
                    offset = k * nx * ny + j * nx + i;
                    Vertex *v = Vertex::newObject();
                    xyz[0] = xorg + i*dx;
                    xyz[1] = yorg + j*dy;
                    xyz[2] = zorg + k*dz;
                    v->setXYZCoords(xyz);
                    v->setID(offset);
                    nodes[ offset ] = v;
                    mesh->addNode(v);
               }
          }
     }

     NodeSequence connect;

     if (space_dim == 2) {
          connect.resize(4);
          for (int j = 0; j < ny - 1; j++) {
               for (int i = 0; i < nx - 1; i++) {
                    offset = j * nx + i;
                    connect[0] = nodes[offset];
                    connect[1] = nodes[offset + 1];
                    connect[2] = nodes[offset + 1 + nx];
                    connect[3] = nodes[offset + nx];
                    Face *face = Face::newObject();
                    face->setNodes(connect);
                    mesh->addFace(face);
               }
          }
     }
     return mesh;
}

////////////////////////////////////////////////////////////////////////////////

int Jaal::quad_concave_tests()
{
     NodeSequence nodes = Mesh::generate_nodes(4);
     Point3D xyz;

     xyz[0] = -1.0;
     xyz[1] = 0.0;
     xyz[2] = 0.0;
     nodes[0]->setXYZCoords(xyz);

     xyz[0] = 0.0;
     xyz[1] = -0.00001;
     xyz[2] = 0.0;
     nodes[1]->setXYZCoords(xyz);

     xyz[0] = 1.0;
     xyz[1] = 0.0;
     xyz[2] = 0.0;
     nodes[2]->setXYZCoords(xyz);

     xyz[0] = 0.0;
     xyz[1] = 1.0;
     xyz[2] = 0.0;
     nodes[3]->setXYZCoords(xyz);

     Face *face = Face::newObject();
     face->setNodes(nodes);

     assert(face->concaveAt() == -1);

     xyz[0] = -0.0;
     xyz[1] = 0.000001;
     xyz[2] = 0.0;
     nodes[1]->setXYZCoords(xyz);
     assert(face->concaveAt() == 1);

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

Mesh* Jaal::quad_to_tri4( Mesh *quadmesh, vector<Vertex*> &steiner)
{
     if( quadmesh == NULL ) return NULL;

     Mesh *trimesh = new Mesh;

     size_t numnodes = quadmesh->getSize(0);
     size_t numfaces = quadmesh->getSize(2);

     trimesh->reserve( numnodes + numfaces, 0);
     trimesh->reserve( 4*numfaces, 2);

     for( size_t i = 0; i < numnodes; i++) {
          Vertex *v = quadmesh->getNodeAt(i);
          if( !v->isRemoved() )
               trimesh->addNode(v);
     }
     steiner.clear();
     steiner.reserve( numfaces );

     Point3D p3d;
     Face  *tface;
     vector<Vertex*> connect(3);
     for( size_t i = 0; i < numfaces; i++) {
          Face *f = quadmesh->getFaceAt(i);
          if( !f->isRemoved() ) {
               int pos = f->concaveAt();
               if( pos >= 0) {
                    p3d = Vertex::mid_point( f->getNodeAt(pos), f->getNodeAt(pos+2) );
               } else {
                    p3d = f->getCentroid();
               }
               Vertex *v0 = Vertex::newObject();
               v0->setXYZCoords( p3d );
               trimesh->addNode(v0);
               steiner.push_back(v0);
               connect[0] = v0;
               for( int j = 0; j < 4; j++) {
                    connect[1] = f->getNodeAt(j);
                    connect[2] = f->getNodeAt(j+1);
                    tface  = Face::newObject();
                    tface->setNodes(connect);
                    trimesh->addFace( tface );
                    assert( tface->concaveAt() < 0);
               }
          }
     }
     return trimesh;
}

///////////////////////////////////////////////////////////////////////////////


Mesh* Jaal::quad_to_tri2( Mesh *quadmesh )
{
     if( quadmesh == NULL ) return NULL;

     Mesh *trimesh = new Mesh;

     size_t numnodes = quadmesh->getSize(0);
     size_t numfaces = quadmesh->getSize(2);

     trimesh->reserve( numnodes,  0);
     trimesh->reserve( 2*numfaces, 2);

     for( size_t i = 0; i < numnodes; i++) {
          Vertex *v = quadmesh->getNodeAt(i);
          if( !v->isRemoved() )
               trimesh->addNode(v);
     }

     Face  *tface;
     vector<Vertex*> connect(3);
     for( size_t i = 0; i < numfaces; i++) {
          Face *f = quadmesh->getFaceAt(i);
          if( !f->isRemoved() ) {
               int pos = f->concaveAt();
               if( pos >= 0) {
                    connect[0] = f->getNodeAt(pos);
                    connect[1] = f->getNodeAt(pos+1);
                    connect[2] = f->getNodeAt(pos+2);
                    tface  = Face::newObject();
                    tface->setNodes(connect);
                    trimesh->addFace( tface );

                    connect[0] = f->getNodeAt(pos);
                    connect[1] = f->getNodeAt(pos+2);
                    connect[2] = f->getNodeAt(pos+3);
                    tface  = Face::newObject();
                    tface->setNodes(connect);
                    trimesh->addFace( tface );
               } else {
                    connect[0] = f->getNodeAt(0);
                    connect[1] = f->getNodeAt(1);
                    connect[2] = f->getNodeAt(2);
                    tface  = Face::newObject();
                    tface->setNodes(connect);
                    trimesh->addFace( tface );

                    connect[0] = f->getNodeAt(0);
                    connect[1] = f->getNodeAt(2);
                    connect[2] = f->getNodeAt(3);
                    tface  = Face::newObject();
                    tface->setNodes(connect);
                    trimesh->addFace( tface );
               }
          }
     }
     return trimesh;
}

///////////////////////////////////////////////////////////////////////////////

void Jaal::advancing_front_triangle_cleanup( Mesh *mesh)
{
     int relexist2 = mesh->build_relations(0, 2);
     int relexist0 = mesh->build_relations(0, 0);

     mesh->search_boundary();

     size_t numNodes = mesh->getSize(0);
     NodeSequence currlayer;
     NodeSet   nextlayer;

     for(size_t i = 0; i < numNodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          if( v->isBoundary() ) {
               v->setLayerID(0);
               currlayer.push_back(v);
          } else
               v->setLayerID(INT_MAX);
     }

     Jaal::MeshOptimization mopt;
     size_t nSize;
     mesh->make_consistently_oriented();
     mopt.shape_optimize(mesh);

     LaplaceNoWeight lw;
     LaplaceSmoothing lapsmooth(mesh);
     lapsmooth.setWeight(&lw);
     lapsmooth.setNumIterations(100);

     int curr_layer_id = 0;
     int nirregular0, nirregular1, numNeighs;

     Vertex *v0, *v1, *ov1, *ov2, *ov3;
     NodeSequence newnodes;
     FaceSequence newfaces, edgeneighs;

     while(1) {
          cout << " PROCESS LAYER : " << curr_layer_id << endl;
          nSize = currlayer.size();
          nirregular0 = 0;
          for( size_t i = 0; i < nSize; i++) {
               Vertex *v = currlayer[i];
               if( !v->isRemoved()  ) {
                    int curr_degree  = v->getNumRelations(2);
                    int ideal_degree = v->get_ideal_face_degree(3);
                    if( curr_degree != ideal_degree ) nirregular0++;
               }
          }

          nSize = currlayer.size();
          for( size_t i = 0; i < nSize; i++) {
               Vertex *v = currlayer[i];
               if( !v->isRemoved() ) {
                    int curr_degree  = v->getNumRelations(2);
                    int ideal_degree = v->get_ideal_face_degree( 3 );
                    if( curr_degree < ideal_degree ) {
                         FaceSequence &v0faces = v->getRelations2();
                         numNeighs = v0faces.size();
                         for( int j = 0; j < numNeighs; j++) {
                              int pos = v0faces[j]->getPosOf(v);
                              v0  = v0faces[j]->getNodeAt( pos + 1 );
                              v1  = v0faces[j]->getNodeAt( pos + 2 );
                              Mesh::getRelations112( v0, v1, edgeneighs);
                              if( edgeneighs.size() == 2 ) {
                                   ov1 = Face::opposite_node(edgeneighs[0], v0, v1);
                                   ov2 = Face::opposite_node(edgeneighs[1], v0, v1);
                                   ov3 = NULL;
                                   assert( ov1 != ov2 );
                                   if( ov1 == v ) ov3 = ov2;
                                   if( ov2 == v ) ov3 = ov1;
                                   assert( ov3 );
                                   if( ov3->getLayerID() == INT_MAX) {
                                        int numSegments = ideal_degree - curr_degree + 1;
                                        mesh->refine_tri_edge( v0, v1, numSegments, newnodes, newfaces);
                                        for( size_t k = 0; k < newnodes.size(); k++)
                                             newnodes[k]->setLayerID( INT_MAX );
                                        break;
                                   }
                              }
                         }
                         curr_degree  = v->getNumRelations(2);
                         if( curr_degree != ideal_degree ) {
                              cout << "Warning: low ideal vertex degree not achieved : " << v->getID() << endl;
                              set_layer_tag(mesh);
                              mesh->saveAs("tmp.off");
                              exit(0);
                         }
                    }

                    if( curr_degree > ideal_degree ) {
                         int excess_degree = curr_degree-ideal_degree;
                         for( int k = 0; k <  excess_degree; k++) {
                              FaceSequence &v1faces = v->getRelations2();
                              numNeighs = v1faces.size();
                              for( int j = 0; j < numNeighs; j++) {
                                   int pos = v1faces[j]->getPosOf(v);
                                   v0  = v1faces[j]->getNodeAt( pos + 1 );
                                   v1  = v1faces[j]->getNodeAt( pos + 2 );
                                   if( v0->getLayerID() == INT_MAX && v1->getLayerID() == INT_MAX) {
                                        Mesh::getRelations112( v0, v1, edgeneighs);
                                        assert( edgeneighs.size() == 2 );
                                        ov1 = Face::opposite_node(edgeneighs[0], v0, v1);
                                        ov2 = Face::opposite_node(edgeneighs[1], v0, v1);
                                        ov3 = NULL;
                                        assert( ov1 != ov2 );
                                        if( ov1 == v ) ov3 = ov2;
                                        if( ov2 == v ) ov3 = ov1;
                                        assert( ov3 );
                                        if( ov3->getLayerID() == INT_MAX) {
                                             mesh->collapse_tri_edge(v0, v1);
                                             break;
                                        }
                                   }
                              }
                         }
                         curr_degree  = v->getNumRelations(2);
                         if( curr_degree != ideal_degree ) {
                              cout << "Warning: high ideal vertex degree not achieved : " << v->getID() << endl;
                              set_layer_tag(mesh);
                              mesh->saveAs( "tmp.off");
                              exit(0);
                         }
                    }
               }
          }

          nSize = currlayer.size();
          nirregular1 = 0;
          for( size_t i = 0; i < nSize; i++) {
               Vertex *v = currlayer[i];
               if( !v->isRemoved()  ) {
                    int curr_degree  = v->getNumRelations(2);
                    int ideal_degree = v->get_ideal_face_degree( 3 );
                    if( curr_degree != ideal_degree ) nirregular1++;
               }

          }
          assert( nirregular1 <= nirregular0);
//        lapsmooth.execute();
          mopt.shape_optimize(mesh);

          cout << "Layer : " << curr_layer_id << endl;
          cout << "# of Irregular nodes before swapping : " << nirregular0 << endl;
          cout << "# of Irregular nodes after swapping  : " << nirregular1 << endl;

          nextlayer.clear();
          nSize = currlayer.size();
          for( size_t i = 0; i < nSize; i++) {
               Vertex *v = currlayer[i];
               NodeSequence &vneighs = v->getRelations0();
               for( size_t k = 0; k < vneighs.size(); k++) {
                    int lid = vneighs[k]->getLayerID();
                    if( lid > curr_layer_id ) {
                         vneighs[k]->setLayerID( curr_layer_id+1 );
                         nextlayer.insert( vneighs[k] );
                    }
               }
          }
          if( nextlayer.empty() ) break;

          NodeSet::const_iterator it;
          currlayer.resize(nextlayer.size() );
          int index = 0;
          for( it = nextlayer.begin(); it != nextlayer.end(); ++it)
               currlayer[index++] = *it;
          curr_layer_id++;
     }

     vector<int>  less_than_ideal, more_than_ideal, total_ideal;

     int numLayers = curr_layer_id;
     less_than_ideal.resize( numLayers );
     more_than_ideal.resize( numLayers );
     total_ideal.resize( numLayers );

     for( int i = 0; i < numLayers; i++) {
          less_than_ideal[i] = 0;
          more_than_ideal[i] = 0;
          total_ideal[i] = 0;
     }

     numNodes = mesh->getSize(0);
     int final_irregular = 0;
     for( size_t i = 0; i < numNodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          if( !v->isRemoved()) {
               int lid   = v->getLayerID();
               int curr_degree  = v->getNumRelations(2);
               int ideal_degree = v->get_ideal_face_degree(3);
               if( curr_degree != ideal_degree ) {
                    final_irregular++;
                    if( curr_degree < ideal_degree) less_than_ideal[lid]++;
                    if( curr_degree > ideal_degree) more_than_ideal[lid]++;
               } else
                    total_ideal[lid]++;
          }
     }

     cout << " Layer   Less   More  Ideal " << endl;
     for( int i = 0; i < numLayers; i++)
          cout << i << setw(10) <<  less_than_ideal[i]
               << setw(10) <<  more_than_ideal[i]
               << setw(10) <<  total_ideal[i] << endl;
     cout << " Final # of irregular nodes : " << final_irregular << endl;

     mopt.shape_optimize(mesh);
     Jaal::set_ideal_node_tag(mesh, 3);
     mesh->saveAs("final.off");
     if (!relexist2) mesh->clear_relations(0, 2);
     if (!relexist0) mesh->clear_relations(0, 0);
}


