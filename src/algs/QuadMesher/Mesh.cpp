#include <iomanip>

#include <meshkit/Mesh.hpp>
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
    Point3D p0 = v0->getXYZCoords();
    Point3D p1 = v1->getXYZCoords();

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
    Point3D p0 = v0->getXYZCoords();
    Point3D p1 = v1->getXYZCoords();

    double dx = p0[0] - p1[0];
    double dy = p0[1] - p1[1];
    double dz = p0[2] - p1[2];

    return sqrt(dx * dx + dy * dy + dz * dz);
}

///////////////////////////////////////////////////////////////////////////////

double
Vertex::length2(const Vertex *v0, const Vertex *v1)
{

    Point3D p0 = v0->getXYZCoords();
    Point3D p1 = v1->getXYZCoords();

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
    if (fabs(A013) + fabs(A123) < fabs(A012) + fabs(A023))
    {
        cout << " Reordered the nodes " << endl;
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
Face::triangulate()
{
    NodeSequence rotatedNodes(4), tconnect(3);
    vector<Face> trifaces;
    if (this->getSize(0) == 4)
    {
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

    Point3D p0 = rotatedNodes[0]->getXYZCoords();
    Point3D p1 = rotatedNodes[1]->getXYZCoords();
    Point3D p2 = rotatedNodes[2]->getXYZCoords();
    Point3D p3 = rotatedNodes[3]->getXYZCoords();

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
    if (pos >= 0)
    {
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

void
Face::opposite_nodes(const PFace quad, PNode n1, PNode n2,
                     PNode &n3, PNode &n4)
{
    PNode qn0 = quad->getNodeAt(0);
    PNode qn1 = quad->getNodeAt(1);
    PNode qn2 = quad->getNodeAt(2);
    PNode qn3 = quad->getNodeAt(3);

    if ((qn0 == n1 && qn1 == n2) || (qn0 == n2 && qn1 == n1))
    {
        n3 = qn2;
        n4 = qn3;
        return;
    }

    if ((qn1 == n1 && qn2 == n2) || (qn1 == n2 && qn2 == n1))
    {
        n3 = qn0;
        n4 = qn3;
        return;
    }

    if ((qn2 == n1 && qn3 == n2) || (qn2 == n2 && qn3 == n1))
    {
        n3 = qn0;
        n4 = qn1;
        return;
    }

    if ((qn3 == n1 && qn0 == n2) || (qn3 == n2 && qn0 == n1))
    {
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
    for (int i = 0; i < 3; i++)
    {
        if (t2->hasNode(connect[i]))
            commonnodes[index++] = connect[i];
    }

    assert(index == 2);

    PNode ot1 = Face::opposite_node(t1, commonnodes[0], commonnodes[1]);
    PNode ot2 = Face::opposite_node(t2, commonnodes[0], commonnodes[1]);

    connect.resize(4);
    connect[0] = ot1;
    connect[1] = commonnodes[0];
    connect[2] = ot2;
    connect[3] = commonnodes[1];

    if (!replace)
    {
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
    assert(hexnodes.size() == 6);

    PFace face1 = Face::newObject();
    PFace face2 = Face::newObject();

    NodeSequence connect(4);

    for (int i = 0; i < 3; i++)
    {
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
        if (face1->isConvex() && face2->isConvex())
        {
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

Point3D
Face::getCentroid() const
{
    Point3D pc;

    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;

    Point3D p3d;
    for (size_t inode = 0; inode < connect.size(); inode++)
    {
        Vertex *v = connect[inode];
        p3d = v->getXYZCoords();
        pc[0] += p3d[0];
        pc[1] += p3d[1];
        pc[2] += p3d[2];
    }

    pc[0] /= (double) connect.size();
    pc[1] /= (double) connect.size();
    pc[2] /= (double) connect.size();

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

    for (int i = 0; i < nnodes; i++)
    {
        Vertex *vertex = connect[i];
        vneighs = vertex->getRelations0();
        for (size_t j = 0; j < vneighs.size(); j++)
            vset.insert(vneighs[j]);
    }

    for (int i = 0; i < nnodes; i++)
        vset.erase(connect[i]);

    NodeSequence vresult;
    if (!vset.empty())
    {
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

    for (int i = 0; i < nSize; i++)
    {
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
    FaceSequence faceneighs, vneighs;

    int nSize = connect.size();
    for (int i = 0; i < nSize; i++)
    {
        Vertex *v0 = connect[(i + 0) % nSize];
        vneighs = v0->getRelations2();
        for (size_t j = 0; j < vneighs.size(); j++)
        {
            if (vneighs[j] != this)
            {
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
    for (int i = 0; i < nSize; i++)
    {
        Vertex *v0 = connect[(i + 0) % nSize];
        Vertex *v1 = connect[(i + 1) % nSize];
        edgeneighs = Mesh::getRelations112(v0, v1);
        for (size_t j = 0; j < edgeneighs.size(); j++)
        {
            if (edgeneighs[j] != this)
            {
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

    for (int i = 0; i < connect.size(); i++)
        if (connect[i]->isRemoved()) return 0;

    for (int i = 0; i < connect.size(); i++)
    {
        int ncount = 0;
        for (int j = 0; j < connect.size(); j++)
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
    for (int i = 0; i < nSize; i++)
    {
        Vertex *v0 = connect[(i + 0) % nSize];
        Vertex *v1 = connect[(i + 1) % nSize];
        if (v0->isBoundary() && v1->isBoundary())
        {
            FaceSequence neighs = Mesh::getRelations112(v0, v1);
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

    assert(xi >= -1.0 && xi <= 1.0);
    assert(eta >= -1.0 && eta <= 1.0);

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
    int numNodes = x.size();
    assert(w.size() == numNodes);

    double sum = 0.0;
    for (int i = 0; i < numNodes; i++) sum += x[i] * w[i];

    return sum;
}

/////////////////////////////////////////////////////////////////////////////////////

int
Face::invertedAt() const
{
    int nsize = getSize(0);
    double x[5], y[5], z[5], triarea;
    Vertex *vertex;
    Point3D xyz;

    for (int i = 0; i < nsize; i++)
    {
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
    for (int i = 0; i < 4; i++)
    {
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
    for (int i = 0; i < 6; i++)
    {
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
    for (int i = 0; i < 10; i++)
    {
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

int
Mesh::refine_quad15(Face *face)
{
    NodeSequence newnodes;
    FaceSequence newfaces;
    face->refine_quad15(newnodes, newfaces);

    for (int i = 0; i < newnodes.size(); i++)
        addNode(newnodes[i]);

    for (int i = 0; i < newfaces.size(); i++)
        addFace(newfaces[i]);

    remove(face);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

Jaal::NodeSequence
Mesh::boundary_chain_nodes(Vertex *v0, Vertex *v1)
{
    NodeSequence bndnodes;

    FaceSequence neighs = Mesh::getRelations112(v0, v1);

    if (neighs.size() != 2) return bndnodes;

    vector<Edge> bndedges;

    bndedges.reserve(6);
    Edge sharededge(v0, v1);

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
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

    for (int i = 1; i < bndedges.size() - 1; i++)
    {
        Vertex *v0 = bndedges[i].getNodeAt(0);
        Vertex *v1 = bndedges[i].getNodeAt(1);
        if (v0 == bndnodes[i])
        {
            bndnodes.push_back(v1);
        }
        else if (v1 == bndnodes[i])
            bndnodes.push_back(v0);
        else
        {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        area[i] = face->getArea();
    }
    minarea = *min_element(area.begin(), area.end());
    maxarea = *max_element(area.begin(), area.end());
}


////////////////////////////////////////////////////////////////////////////////

const vector<double>
Mesh::getCoordsArray()
{
    size_t numnodes = getSize(0);

    vector<double> vcoords(3 * numnodes);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        Point3D xyz = v->getXYZCoords();
        vcoords[3 * i + 0] = xyz[0];
        vcoords[3 * i + 1] = xyz[1];
        vcoords[3 * i + 2] = xyz[2];
    }
    return vcoords;
}

////////////////////////////////////////////////////////////////////////////////

const vector<size_t>
Mesh::getNodesArray()
{
    size_t numfaces = getSize(2);

    vector<size_t> nodearray;

    int topo = isHomogeneous();
    if (topo) nodearray.reserve(topo * numfaces);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        for (int j = 0; j < face->getSize(0); j++)
        {
            Vertex *v = face->getNodeAt(j);
            nodearray.push_back(v->getID());
        }
    }

    return nodearray;
}

////////////////////////////////////////////////////////////////////////////////

Mesh*
Mesh::deep_copy()
{
    std::map<Vertex*, Vertex*> vmap;

    Mesh *newmesh = new Mesh;
    size_t numnodes = getSize(0);
    newmesh->reserve(numnodes, 0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vold = getNodeAt(i);
        Vertex *vnew = Vertex::newObject();
        vnew->setXYZCoords(vold->getXYZCoords());
        vmap[vold] = vnew;
        newmesh->addNode(vnew);
    }

    size_t numfaces = getSize(2);
    newmesh->reserve(numnodes, 2);
    NodeSequence connect;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *fold = getFaceAt(i);
        connect = fold->getNodes();
        for (int j = 0; j < connect.size(); j++)
            connect[j] = vmap[connect[j]];

        Face *fnew = Face::newObject();
        fnew->setNodes(connect);
        newmesh->addFace(fnew);
    }

    return newmesh;
}

////////////////////////////////////////////////////////////////////////////////

int
Mesh::setCoordsArray(const vector<double> &vcoords)
{
    size_t numnodes = getSize(0);

    if (vcoords.size() != 3 * numnodes) return 1;

    Point3D xyz;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        xyz[0] = vcoords[3 * i + 0];
        xyz[1] = vcoords[3 * i + 1];
        xyz[2] = vcoords[3 * i + 2];
        v->setXYZCoords(xyz);
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
    for (size_t i = 1; i < boundedges.size(); i++)
        listedges.push_back(boundedges[i]);
    boundedges.clear();

    boundedges.reserve(nSize);

    Vertex *first_vertex = edge.getNodeAt(0);
    Vertex *curr_vertex = edge.getNodeAt(1);

    boundedges.push_back(edge);

    list<Edge>::iterator it;

    for (size_t i = 0; i < nSize; i++)
    {
        for (it = listedges.begin(); it != listedges.end(); ++it)
        {
            edge = *it;
            Vertex *v0 = edge.getNodeAt(0);
            Vertex *v1 = edge.getNodeAt(1);
            if (v0 == curr_vertex)
            {
                curr_vertex = v1;
                boundedges.push_back(edge);
                break;
            }
            if (v1 == curr_vertex)
            {
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

    for (size_t i = 0; i < boundedges.size(); i++)
    {
        Vertex *v0 = boundedges[i].getNodeAt(0);
        Vertex *v1 = boundedges[i].getNodeAt(1);
        relations00[v0].insert(v1);
        relations00[v1].insert(v0);
    }

    std::map<Vertex*, set<Vertex*> > ::const_iterator it;
    for (it = relations00.begin(); it != relations00.end(); ++it)
    {
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

    for (size_t i = 1; i < boundedges.size(); i++)
    {
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
    for (size_t i = 0; i < nSize; i++)
    {
        if (boundedges[i].getNodeAt(0) == first_vertex) istart = i;
        listedges[i] = boundedges[i];
    }

    for (size_t i = 0; i < nSize; i++)
        boundedges[i] = listedges[(i + istart) % nSize];

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

FaceSequence
Mesh::getRelations102(PNode vtx0, PNode vtx1)
{
    FaceSequence faceneighs;

    FaceSequence v0faces = vtx0->getRelations2();
    FaceSequence v1faces = vtx1->getRelations2();

    if (v0faces.empty() || v1faces.empty())
    {
        cout << "Warning: Vertex-Faces relations are empty " << endl;
        return faceneighs;
    }

    FaceSet vset;
    for (size_t i = 0; i < v0faces.size(); i++)
        vset.insert(v0faces[i]);

    for (size_t i = 0; i < v1faces.size(); i++)
        vset.insert(v1faces[i]);

    std::set<PFace>::iterator it;

    if (vset.size())
    {
        faceneighs.resize(vset.size());
        int index = 0;
        for (it = vset.begin(); it != vset.end(); ++it)
            faceneighs[index++] = *it;
    }

    return faceneighs;
}

///////////////////////////////////////////////////////////////////////////////

FaceSequence
Mesh::getRelations112(const PNode vtx0, const PNode vtx1)
{
    FaceSequence faceneighs;
    FaceSequence v0faces = vtx0->getRelations2();
    FaceSequence v1faces = vtx1->getRelations2();

    if (v0faces.empty() || v1faces.empty())
    {
        cout << "Warning: Vertex-Faces relations are empty " << endl;
        return faceneighs;
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

    return faceneighs;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::count_edges()
{
    if (getAdjTable(1, 1)) return edges.size();

    int relexist = build_relations(0, 0);

    size_t numnodes = getSize(0);

    NodeSequence neighs;
    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        neighs = nodes[i]->getRelations0();
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

    NodeSequence neighs;

    edges.clear();

    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        neighs = nodes[i]->getRelations0();
        for (size_t j = 0; j < neighs.size(); j++)
            if (nodes[i] > neighs[j])
            {
                Edge *newedge = new Edge(nodes[i], neighs[j]);
                assert(newedge);
                edges.push_back(newedge);
            }
    }

    if (!relexist0) clear_relations(0, 0);
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::prune()
{
    if (adjTable[0][0])
    {
        NodeSequence relations0;
        for (size_t i = 0; i < nodes.size(); i++)
        {
            relations0 = nodes[i]->getRelations0();
            for (size_t j = 0; j < relations0.size(); j++)
            {
                if (relations0[j]->isRemoved())
                    nodes[i]->removeRelation0(relations0[j]);
            }
        }
    }

    if (adjTable[0][2])
    {
        FaceSequence relations2;
        for (size_t i = 0; i < nodes.size(); i++)
        {
            relations2 = nodes[i]->getRelations2();
            for (size_t j = 0; j < relations2.size(); j++)
            {
                if (relations2[j]->isRemoved())
                    nodes[i]->removeRelation2(relations2[j]);
            }
        }
    }

    for (size_t i = 0; i < nodes.size(); i++)
    {
        Vertex *v = nodes[i];
        if (v->isRemoved())
        {
            if (find(garbageNodes.begin(), garbageNodes.end(), v) == garbageNodes.end())
                garbageNodes.push_back(v);
        }
    }

    for (size_t i = 0; i < faces.size(); i++)
    {
        Face *f = faces[i];
        if (f->isRemoved())
        {
            if (find(garbageFaces.begin(), garbageFaces.end(), f) == garbageFaces.end())
                garbageFaces.push_back(f);
        }
    }

    NodeSequence::iterator vend;
    vend = remove_if(nodes.begin(), nodes.end(), EntityRemovedPred());
    nodes.erase(vend, nodes.end());

    FaceSequence::iterator fend;
    fend = remove_if(faces.begin(), faces.end(), EntityRemovedPred());
    faces.erase(fend, faces.end());

    enumerate(0);
    enumerate(2);

    for (size_t i = 0; i < nodes.size(); i++)
        assert(nodes[i]->isActive());

    for (size_t i = 0; i < faces.size(); i++)
        assert(faces[i]->isActive());
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
    list<Face*>::const_iterator fiter;
    for (fiter = garbageFaces.begin(); fiter != garbageFaces.end(); ++fiter)
    {
        Face *face = *fiter;
        assert(face);
        if (face->isRemoved()) delete face;
    }
    garbageFaces.clear();

    list<Vertex*>::const_iterator viter;
    for (viter = garbageNodes.begin(); viter != garbageNodes.end(); ++viter)
    {
        Vertex *vertex = *viter;
        assert(vertex);
        if (vertex->isRemoved()) delete vertex;
    }
    garbageNodes.clear();
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::enumerate(int etype)
{
    size_t index = 0;

    NodeSequence::const_iterator viter;
    if (etype == 0)
    {
        index = 0;
        for (viter = nodes.begin(); viter != nodes.end(); ++viter)
        {
            Vertex *vertex = *viter;
            vertex->setID(index++);
        }
    }

    FaceSequence::const_iterator fiter;
    if (etype == 2)
    {
        index = 0;
        for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
        {
            Face *face = *fiter;
            face->setID(index++);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::getBoundarySize(int d) const
{
    size_t ncount = 0;

    if (d == 0)
    {
        for (size_t i = 0; i < nodes.size(); i++)
            if (!nodes[i]->isRemoved() && nodes[i]->isBoundary())
                ncount++;
    }

    if (d == 2)
    {
        for (size_t i = 0; i < faces.size(); i++)
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
    for (size_t i = 0; i < faces.size(); i++)
    {
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

    NodeSequence neighs = myself->getRelations0();

    for (size_t i = 0; i < neighs.size(); i++)
    {
        double d = Vertex::length2(myself, neighs[i]);
        if (d < mindist)
        {
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
    if (!isPruned()) prune();

    size_t numfaces = getSize(2);

    if (getAdjTable(0, 2))
    {
        for (size_t iface = 0; iface < numfaces; iface++)
        {
            Face *face = getFaceAt(iface);
            assert(!face->isRemoved());
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *vtx = face->getNodeAt(j);
                assert(!vtx->isRemoved());
                assert(vtx->hasRelation2(face));
            }
        }
    }

    // Delete all the old stuff and rebuild new. Sometimes used for debugging purpose
    // also.
    if (rebuild) clear_relations(0, 2);

    if (adjTable[0][2] == 1) return 1;

    clear_relations(0, 2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        for (int j = 0; j < face->getSize(0); j++)
        {
            Vertex *vtx = face->getNodeAt(j);
            vtx->addRelation2(face);
        }
    }
    adjTable[0][2] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::build_relations00(bool rebuild)
{
    if (!isPruned()) prune();

    // Delete all the old stuff and rebuild new. Sometimes used for debugging purpose
    // also.
    if (rebuild) clear_relations(0, 0);

    if (adjTable[0][0] == 1) return 1;

    clear_relations(0, 0);

    size_t numfaces = getSize(2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        size_t nnodes = face->getSize(0);
        for (size_t j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            v0->addRelation0(v1);
            v1->addRelation0(v0);
        }
    }
    adjTable[0][0] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::clear_relations(int src, int dst)
{
    size_t numnodes = getSize(0);

    if (src == 0 && dst == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vtx = getNodeAt(i);
            vtx->clearRelations(0);
        }
        adjTable[0][0] = 0;
    }

    if (src == 0 && dst == 2)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
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

    size_t numnodes = getSize(0);
    size_t numfaces = getSize(2);

    int bmark;
    FaceSequence neighs;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        size_t nnodes = face->getSize(0);
        for (size_t j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            neighs = Mesh::getRelations112(v0, v1);
            if (neighs.size() == 1)
            {
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
    for (it = faces.begin(); it != faces.end(); ++it)
    {
        Face *face = *it;
        if (face->getType() == facetype)
            ncount++;
    }

    FaceSequence tmpfaces;
    if (ncount)
    {
        tmpfaces.resize(ncount);
        size_t index = 0;
        for (it = faces.begin(); it != faces.end(); ++it)
        {
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
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            neighs = Mesh::getRelations112(v0, v1);
            if (neighs.size() > 2)
            {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = getNodeAt(i);
        if (!vertex->isBoundary())
        {
            FaceSequence vfaces = vertex->getRelations2();
            if (vfaces.size() != degree_of_regular_node)
                ncount++;
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
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            neighs = Mesh::getRelations112(v0, v1);
            if (neighs.size() == 2)
            {
                assert(!neighs[0]->isRemoved());
                assert(!neighs[1]->isRemoved());
                int dir1 = neighs[0]->getOrientation(v0, v1);
                int dir2 = neighs[1]->getOrientation(v0, v1);
                if (dir1 * dir2 == 1)
                {
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

    if (!relexist)
        clear_relations(0, 2);

    return consistent;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::make_consistently_oriented()
{
    build_relations(0, 2);

    Face *face = NULL;
    deque<Face*> faceQ;
    FaceSequence neighs;

    size_t numfaces = getSize(2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        face = getFaceAt(iface);
        face->setID(iface);
        face->setVisitMark(0);
    }

    face = getFaceAt(0);
    faceQ.push_back(face);

    while (!faceQ.empty())
    {
        Face *face = faceQ.front();
        faceQ.pop_front();
        if (!face->isVisited())
        {
            face->setVisitMark(1);
            int nnodes = face->getSize(0);
            for (int j = 0; j < nnodes; j++)
            {
                Vertex *v0 = face->getNodeAt(j);
                Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
                neighs = Mesh::getRelations112(v0, v1);
                if (neighs.size() == 2)
                {
                    int dir1 = neighs[0]->getOrientation(v0, v1);
                    int dir2 = neighs[1]->getOrientation(v0, v1);
                    if (dir1 * dir2 == 1)
                    {
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

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        face = getFaceAt(iface);
        if (!face->isVisited())
            cout << "Error: not visited : " << face->getID() << " "
            << face->isVisited() << endl;
    }

    clear_relations(0, 2);
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

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        face = getFaceAt(iface);
        face->setID(iface);
        face->setVisitMark(0);
    }

    int numComponents = 0;

    while (1)
    {
        face = NULL;
        faceQ.clear();
        for (size_t iface = 0; iface < numfaces; iface++)
        {
            face = getFaceAt(iface);
            if (!face->isVisited())
            {
                face->setPartID(numComponents);
                faceQ.push_back(face);
                break;
            }
        }

        if (faceQ.empty())
            break;

        while (!faceQ.empty())
        {
            Face *face = faceQ.front();
            faceQ.pop_front();
            if (!face->isVisited())
            {
                face->setPartID(numComponents);
                face->setVisitMark(1);
                int nnodes = face->getSize(0);
                for (int j = 0; j < nnodes; j++)
                {
                    Vertex *v0 = face->getNodeAt(j);
                    Vertex *v1 = face->getNodeAt((j + 1) % nnodes);

                    int proceed = 1;
                    if (stop_at_interface)
                    {
                        if (v0->isConstrained() && v1->isConstrained())
                        {
                            Edge edge(v0, v1);
                            if (hasFeatureEdge(edge)) proceed = 0;
                        }
                    }

                    if (proceed)
                    {
                        neighs = Mesh::getRelations112(v0, v1);
                        if (neighs.size() == 2)
                        {
                            faceQ.push_back(neighs[0]);
                            faceQ.push_back(neighs[1]);
                        }
                    }

                }
            }
        } // Complete one Component
        numComponents++;
    }

    for (size_t iface = 0; iface < numfaces; iface++)
    {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setVisitMark(0);
    }

    size_t numfaces = getSize(2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        if (face->getPartID() == id)
        {
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *v = face->getNodeAt(j);
                if (!v->isVisited())
                {
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
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
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
    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
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
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
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
    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
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
    FaceSequence neighs = Mesh::getRelations112(v0, v1);

    Vertex *vn0, *vn1; // Next-Edge Nodes;

    Face *nextface;
    if (neighs.size() == 2)
    {
        if (!neighs[0]->isVisited() && neighs[1]->isVisited())
        {
            nextface = neighs[0];
            if (nextface->getSize(0) == 4)
            {
                Face::opposite_nodes(nextface, v0, v1, vn0, vn1);
                nextface->setVisitMark(1);
                strip.push_back(nextface);
                expand_strip(nextface, vn0, vn1, strip);
            }
        }
        if (neighs[0]->isVisited() && !neighs[1]->isVisited())
        {
            nextface = neighs[1];
            if (nextface->getSize(0) == 4)
            {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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

    if (bound_what == 0)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            if (face->hasBoundaryNode())
                bfaces.insert(face);
        }
    }

    if (bound_what == 1)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            if (face->has_boundary_edge())
                bfaces.insert(face);
        }
    }

    FaceSequence result;

    size_t nSize = bfaces.size();

    if (nSize)
    {
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
    for (size_t i = 0; i < numnodes; i++)
    {
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

    NodeSequence seq;
    FaceSequence vfaces;

    // Query from the boundary nodes ...
    if (from_where == 1)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *v = getNodeAt(i);
            vfaces = v->getRelations2();
            if (v->isBoundary() && vfaces.size() != regular_count) seq.push_back(v);
        }
    }

    // Query from the internal nodes ...
    if (from_where == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *v = getNodeAt(i);
            vfaces = v->getRelations2();
            if (!v->isBoundary() && vfaces.size() != regular_count) seq.push_back(v);
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        face->setTag(0);
    }

    FaceSequence strip1, strip2;
    int id = 0;
    for (size_t i = 0; i < bound_faces.size(); i++)
    {
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

    vector<int> degree(numnodes);
    FaceSequence neighs;
    for (int i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        neighs = v->getRelations2();
        degree[i] = neighs.size();
    }

    int mindegree = *min_element(degree.begin(), degree.end());
    int maxdegree = *max_element(degree.begin(), degree.end());

    cout << " Mesh Topological Quality : " << endl;

    cout << " *********************** " << endl;
    cout << " Degree        Count " << endl;
    cout << " *********************** " << endl;

    for (int i = mindegree; i <= maxdegree; i++)
    {
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

    if (layerid == 0)
    {
        search_boundary();
        for (size_t i = 0; i < numNodes; i++)
        {
            Vertex *vertex = getNodeAt(i);
            if (vertex->isBoundary())
            {
                vertex->setLayerID(0);
                ncount++;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < numNodes; i++)
        {
            Vertex *vertex = getNodeAt(i);
            if (vertex->getLayerID() == layerid - 1)
            {
                NodeSequence vnodes = vertex->getRelations0();
                for (size_t j = 0; j < vnodes.size(); j++)
                {
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
    for (size_t i = 0; i < numNodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setLayerID(-1);
        v->setVisitMark(0);
    }

    NodeSequence vertexQ, nextQ, vneighs;
    for (size_t i = 0; i < numNodes; i++)
    {
        Vertex *v = getNodeAt(i);
        if (v->isBoundary())
        {
            v->setLayerID(0);
            v->setVisitMark(1);
            vertexQ.push_back(v);
        }
    }

    if (vertexQ.empty())
    {
        cout << "Warning: No boundary detected " << endl;
    }

    int layerid = 1;
    while (!vertexQ.empty())
    {
        nextQ.clear();
#ifdef SEQUENCE_IS_VECTOR
        nextQ.reserve(vertexQ.size());
#endif
        for (size_t j = 0; j < vertexQ.size(); j++)
        {
            Vertex *currVertex = vertexQ[j];
            vneighs = currVertex->getRelations0();
            for (size_t i = 0; i < vneighs.size(); i++)
            {
                Vertex *vn = vneighs[i];
                if (!vn->isVisited()) nextQ.push_back(vn);
            }
        }

        for (size_t j = 0; j < nextQ.size(); j++)
        {
            Vertex *currVertex = nextQ[j];
            currVertex->setLayerID(layerid);
            currVertex->setVisitMark(1);
        }
        vertexQ = nextQ;
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

    if (layerid == 0)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            int nsize = face->getSize(0);
            for (int j = 0; j < nsize; j++)
            {
                Vertex *v0 = face->getNodeAt((j + 0) % nsize);
                Vertex *v1 = face->getNodeAt((j + 1) % nsize);
                FaceSequence vfaces = Mesh::getRelations112(v0, v1);
                if (vfaces.size() == 1)
                {
                    face->setLayerID(0);
                    ncount++;
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            if (face->getLayerID() == layerid - 1)
            {
                FaceSequence vfaces = face->getRelations212();
                for (size_t j = 0; j < vfaces.size(); j++)
                {
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

    for (size_t i = 0; i < numFaces; i++)
    {
        Face *f = getFaceAt(i);
        f->setLayerID(0);
        f->setVisitMark(0);
    }

    FaceSequence faceQ, nextQ;
    for (size_t i = 0; i < numFaces; i++)
    {
        Face *f = getFaceAt(i);
        f->setLayerID(1);
        assert(!f->isRemoved());
        if (f->has_boundary_edge())
        {
            f->setLayerID(0);
            f->setVisitMark(1);
            faceQ.push_back(f);
        }
    }

    int layerid = 1;
    FaceSequence neighs;

    while (!faceQ.empty())
    {
        nextQ.clear();
#ifdef SEQUENCE_IS_VECTOR
        nextQ.reserve(faceQ.size());
#endif
        for (size_t j = 0; j < faceQ.size(); j++)
        {
            Face *currFace = faceQ[j];
            neighs = currFace->getRelations212();
            for (size_t i = 0; i < neighs.size(); i++)
            {
                Face *nf = neighs[i];
                if (!nf->isVisited())
                    nextQ.push_back(nf);
            }
        }

        for (size_t i = 0; i < nextQ.size(); i++)
        {
            Face *f = nextQ[i];
            f->setLayerID(layerid);
            f->setVisitMark(1);
        }

        layerid++;
        faceQ = nextQ;
    }

    for (size_t i = 0; i < numFaces; i++)
    {
        Face *f = getFaceAt(i);
        assert(f->isVisited());
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
    if (mentity == 0)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            int nsize = face->getSize(0);
            for (int j = 0; j < nsize; j++)
            {
                Vertex *v0 = face->getNodeAt((j + 0) % nsize);
                Vertex *v1 = face->getNodeAt((j + 1) % nsize);
                int l1 = v0->getLayerID();
                int l2 = v1->getLayerID();
                assert(l1 >= 0 && l2 >= 0);
                if (abs(l2 - l1) > 1) return 1;
            }
        }
    }

    if (mentity == 2)
    {
        int relexist2 = build_relations(0, 2);
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            FaceSequence neighs = face->getRelations212();
            int l1 = face->getLayerID();
            assert(l1 >= 0);
            for (int j = 0; j < neighs.size(); j++)
            {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setVisitMark(0);
    }

    size_t numfaces = getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        for (int j = 0; j < face->getSize(0); j++)
        {
            Vertex *v = face->getNodeAt(j);
            v->setVisitMark(1);
        }
    }

    for (size_t i = 0; i < numnodes; i++)
    {
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

#ifdef REMOVE_LATER 

int
Mesh::check_unused_objects()
{
    set<Vertex*> vset;
    for (size_t i = 0; i < nodes.size(); i++)
        vset.insert(nodes[i]);

    if (vset.size() != nodes.size())
        cout << "Warning: There are some duplicated nodes in the mesh " << endl;

    size_t numfaces = getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *f = getFaceAt(i);
        if (!f->isRemoved())
        {
            for (int j = 0; j < f->getSize(0); j++)
            {
                Vertex *v = f->getNodeAt(j);
                if (v->isRemoved())
                    cout << "Goofed up: Face vertex is deleted " << endl;
                if (vset.find(v) == vset.end())
                    cout << "Warning: A face vertex is not in the node list "
                        << v->getID() << endl;
            }
        }
    }
    return 0;
}
#endif

///////////////////////////////////////////////////////////////////////////////

/*
bool
Mesh::isDelaunay()
{
    bool retval = 1;
    int relexist0 = build_relations(0, 0);

    Point3D pa, pb, pc, pd, pCenter;

    size_t numfaces = getSize(2);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        face->setTag(1);
    }
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        pa = face->getNodeAt(0)->getXYZCoords();
        pb = face->getNodeAt(1)->getXYZCoords();
        pc = face->getNodeAt(2)->getXYZCoords();
        TriCircumCenter3D(&pa[0], &pb[0], &pc[0], &pCenter[0]);
        double radius2 = Math::length2(pa, pCenter);
        NodeSequence neighs = face->getRelations0();
        for (size_t j = 0; j < neighs.size(); j++)
        {
            pd = neighs[j]->getXYZCoords();
            if (Math::length2(pd, pCenter) < radius2)
            {
                face->setTag(0);
                retval = 0;
                break;
            }
        }
    }
    if (!relexist0)
        clear_relations(0, 0);
}
 */
///////////////////////////////////////////////////////////////////////////////

double
Mesh::getSurfaceArea()
{
    double facearea, sumArea = 0.0;

    size_t numfaces = getSize(2);
    double minarea = MAXDOUBLE;
    double maxarea = 0.0;
    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (int i = 1; i < bndedges.size() - 1; i++)
    {
        Vertex *v0 = bndedges[i].getNodeAt(0);
        Vertex *v1 = bndedges[i].getNodeAt(1);
        if (v0 == bndnodes[i])
        {
            bndnodes.push_back(v1);
        }
        else if (v1 == bndnodes[i])
            bndnodes.push_back(v0);
        else
        {
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

    for (size_t i = 0; i < numfaces; i++)
    {
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

void
Mesh::deleteAll()
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            clear_relations(i, j);

    for (size_t i = 0; i < nodes.size(); i++)
        delete nodes[i];
    nodes.clear();

    for (size_t i = 0; i < edges.size(); i++)
        delete edges[i];
    edges.clear();

    for (size_t i = 0; i < faces.size(); i++)
        delete faces[i];
    faces.clear();

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
    for (size_t i = 0; i < edges.size(); i++)
    {
        Vertex *v0 = edges[i]->getNodeAt(0);
        Vertex *v1 = edges[i]->getNodeAt(1);
        Point3D xyz = Vertex::mid_point(v0, v1);
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
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = this->getFaceAt(iface);
        for (int j = 0; j < 4; j++)
        {
            Vertex *v0 = face->getNodeAt((j + 0) % 4);
            Vertex *v1 = face->getNodeAt((j + 1) % 4);
            Edge edge(v0, v1);
            Vertex *vh = edge.getHashNode();
            ilower = edgemap.lower_bound(vh);
            iupper = edgemap.upper_bound(vh);
            Vertex *vmid = NULL;
            for (iter = ilower; iter != iupper; ++iter)
            {
                Edge *existing_edge = iter->second;
                if (existing_edge->isSameAs(edge))
                {
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
    size_t numfaces = getSize(2);
    FaceSequence newfaces;
    NodeSequence newnodes;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        face->refine_quad15(newnodes, newfaces);
        for (int j = 0; j < newnodes.size(); j++)
            addNode(newnodes[j]);
        for (int j = 0; j < newfaces.size(); j++)
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
/*
double Vertex ::getFeatureAngle() const
{
  FaceSequence vfaces = getRelations2();
  assert( !vfaces.empty() ); 

  double theta = 0.0;
  for( int i = 0; i < vfaces.size(); i++) 
       theta += vfaces[i]->getAngleAt( this );

  return theta;
}
 */

///////////////////////////////////////////////////////////////////////////////

/*
void
Mesh::setFeatureAngles()
{
    int relexist0 = build_relations(0, 0);
    int relexist2 = build_relations(0, 2);

    size_t numnodes = getSize(0);
    vector<Vertex*> vneighs, bndnodes, sidenodes;

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = getNodeAt(i);
        if (vertex->isBoundary())
        {
            vneighs = vertex->getRelations0();
            bndnodes.clear();
            for (int j = 0; j < vneighs.size(); j++)
                if (vneighs[j]->isBoundary()) bndnodes.push_back(vneighs[j]);
            sidenodes.clear();
            for (int j = 0; j < bndnodes.size(); j++)
            {
                vector<Face*> vfaces = Mesh::getRelations112(vertex, bndnodes[j]);
                if (vfaces.size() == 1)
                    sidenodes.push_back(bndnodes[j]);
            }
            assert(sidenodes.size() == 2);
            Point3D p0 = vertex->getXYZCoords();
            Point3D p1 = sidenodes[0]->getXYZCoords();
            Point3D p2 = sidenodes[1]->getXYZCoords();
            Vec3D v1 = Math::create_vector(p1, p0);
            Vec3D v2 = Math::create_vector(p2, p0);
            double angle = Math::getVectorAngle(v1, v2, ANGLE_IN_DEGREES);
            vertex->setFeatureAngle(angle);
        }
    }

    if (!relexist0)
        clear_relations(0, 0);

    if (!relexist2)
        clear_relations(0, 2);
}
 */

///////////////////////////////////////////////////////////////////////////////

double Vertex::getFeatureLength() const
{
    if (!isBoundary()) return MAXDOUBLE;

    NodeSequence vneighs = getRelations0();

    assert(!vneighs.empty());

    double minlen = MAXDOUBLE;
    for (int j = 0; j < vneighs.size(); j++)
    {
        if (vneighs[j]->isBoundary())
        {
            Point3D p0 = getXYZCoords();
            Point3D p1 = vneighs[j]->getXYZCoords();
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

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        Vertex *vertex = face->getHashNode();
        mapfaces[vertex].push_back(face);
    }

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        const FaceSequence &hashfaces = mapfaces[face->getHashNode() ];
        size_t ncount = 0;
        for (size_t j = 0; j < hashfaces.size(); j++)
            if (hashfaces[j]->isSameAs(face)) ncount++;
        if (ncount != 1)
        {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *f = getFaceAt(i);
        if (!f->isConvex()) ncount++;
    }
    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::count_inverted_faces()
{
    size_t numfaces = getSize(2);
    size_t ncount = 0;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *f = getFaceAt(i);
        if (f->invertedAt() >= 0) ncount++;
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

    for (int i = 0; i < numnodes; i++)
    {
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

    for (int i = 0; i < numnodes; i++)
    {
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

    for (int i = 0; i < numnodes; i++)
    {
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

    for (int i = 0; i < numnodes; i++)
    {
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

NodeSequence
Jaal::linear_interpolation(Vertex *v0, Vertex *v1, int n)
{
    assert(n >= 2);
    NodeSequence newnodes;

    Point3D xyz0 = v0->getXYZCoords();
    Point3D xyz1 = v1->getXYZCoords();

    newnodes.resize(n);
    newnodes[0] = v0;
    newnodes[n - 1] = v1;

    double dt = 2.0 / (double) (n - 1);

    Point3D xyzt;
    for (int i = 1; i < n - 1; i++)
    {
        double t = -1.0 + i*dt;
        xyzt[0] = TFI::linear_interpolation(t, xyz0[0], xyz1[0]);
        xyzt[1] = TFI::linear_interpolation(t, xyz0[1], xyz1[1]);
        xyzt[2] = TFI::linear_interpolation(t, xyz0[2], xyz1[2]);
        newnodes[i] = Vertex::newObject();
        newnodes[i]->setXYZCoords(xyzt);
    }

    return newnodes;
}

///////////////////////////////////////////////////////////////////////////////

void
set_tfi_coords(int i, int j, int nx, int ny, vector<Vertex*> &qnodes)
{
    int offset;

    offset = 0;
    const Point3D v00 = qnodes[offset]->getXYZCoords();

    offset = i;
    const Point3D vr0 = qnodes[offset]->getXYZCoords();

    offset = (nx - 1);
    const Point3D v10 = qnodes[offset]->getXYZCoords();

    offset = j*nx;
    const Point3D v0s = qnodes[offset]->getXYZCoords();

    offset = j * nx + (nx - 1);
    const Point3D v1s = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx;
    const Point3D v01 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + i;
    const Point3D vr1 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + (nx - 1);
    const Point3D v11 = qnodes[offset]->getXYZCoords();

    Point3D vrs;

    double dr = 2.0 / (double) (nx - 1);
    double ds = 2.0 / (double) (ny - 1);

    double r = -1.0 + i*dr;
    double s = -1.0 + j*ds;
    for (int k = 0; k < 3; k++)
    {
        vrs[k] = TFI::transfinite_blend(r, s,
                                        v00[k], v10[k], v11[k], v01[k],
                                        vr0[k], v1s[k], vr1[k], v0s[k]);
    }
    offset = j * nx + i;
    qnodes[offset]->setXYZCoords(vrs);

}
///////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_quad_loop(Mesh *mesh,
                       NodeSequence &boundnodes, int nx, int ny,
                       NodeSequence &newnodes, FaceSequence &newfaces,
                       bool smooth)
{

    newnodes.clear();
    newfaces.clear();

    assert(boundnodes.size() == 2 * nx + 2 * ny - 4);
    for (size_t i = 0; i < boundnodes.size(); i++)
        assert(!boundnodes[i]->isRemoved());

    vector<Vertex*> qnodes(nx * ny);

    //
    // Put the boundary nodes on the structured mesh: The orientation is
    // Counter clockwise ( south->east->north->west );
    //

    int offset, index = 0;

    // South Side ...
    index = 0;
    for (int i = 0; i < nx; i++)
    {
        offset = i;
        qnodes[i] = boundnodes[index++];
        Point3D xyz = qnodes[i]->getXYZCoords();
    }

    // East Side
    for (int j = 1; j < ny; j++)
    {
        offset = j * nx + (nx - 1);
        qnodes[offset] = boundnodes[index++];
    }

    // North Side
    for (int i = nx - 2; i >= 0; i--)
    {
        offset = (ny - 1) * nx + i;
        qnodes[offset] = boundnodes[index++];
    }

    // West Side
    for (int j = ny - 2; j >= 1; j--)
    {
        offset = j*nx;
        qnodes[j * nx] = boundnodes[index++];
    }

    // Now all internal nodes ....
    newnodes.resize((nx - 2)*(ny - 2));

    index = 0;
    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            Vertex *v = Vertex::newObject();
            offset = j * nx + i;
            qnodes[offset] = v;
            newnodes[index++] = v;
            set_tfi_coords(i, j, nx, ny, qnodes); // Coordinates values
        }
    }

    newfaces.resize((nx - 1)*(ny - 1));
    NodeSequence qc(4);

    // Create new faces ...
    index = 0;
    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
            int offset = j * nx + i;
            qc[0] = qnodes[offset];
            qc[1] = qnodes[offset + 1];
            qc[2] = qnodes[offset + 1 + nx];
            qc[3] = qnodes[offset + nx];
            Face *face = Face::newObject();
            face->setNodes(qc);
            newfaces[index++] = face;
        }
    }

    // Update the mesh ...
    for (size_t i = 0; i < newnodes.size(); i++)
        mesh->addNode(newnodes[i]);

    for (size_t i = 0; i < newfaces.size(); i++)
        mesh->addFace(newfaces[i]);

    if (smooth)
    {
        // Perform some laplacian smoothing inside the local mesh...
        LaplaceLengthWeight lw;
        LaplaceSmoothing lapsmooth(mesh);
        lapsmooth.setWeight(&lw);
        lapsmooth.setNumIterations(10);
        lapsmooth.localized_at(newnodes);
    }

    // Check for any inversion of the element, if there is inversion,
    // undo everthing (i.e. remove new nodes and faces).
    //
    for (size_t i = 0; i < newfaces.size(); i++)
    {
        if (newfaces[i]->invertedAt() >= 0)
        {
            for (size_t i = 0; i < newfaces.size(); i++)
                mesh->remove(newfaces[i]);
            for (size_t i = 0; i < newnodes.size(); i++)
                mesh->remove(newnodes[i]);
            newnodes.clear();
            newfaces.clear();
            return 1;
        }
    }

    return 0;

}

///////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_quad_loop(Mesh *mesh,
                       NodeSequence &anodes,
                       NodeSequence &bnodes,
                       NodeSequence &cnodes,
                       NodeSequence &dnodes,
                       NodeSequence &newnodes,
                       FaceSequence &newfaces,
                       bool smooth)
{
    newnodes.clear();
    newfaces.clear();

    assert(anodes.size() == cnodes.size());
    assert(bnodes.size() == dnodes.size());

    int nx = anodes.size();
    int ny = bnodes.size();

    for (int i = 0; i < nx; i++)
    {
        assert(mesh->contains(anodes[i]));
        assert(mesh->contains(cnodes[i]));
    }

    for (int i = 0; i < ny; i++)
    {
        assert(mesh->contains(bnodes[i]));
        assert(mesh->contains(dnodes[i]));
    }

    NodeSequence boundnodes(2 * nx + 2 * ny - 4);
    int index = 0;

    // Append anodes..
    for (int i = 0; i < anodes.size(); i++)
        boundnodes[index++] = anodes[i];

    // Append bnodes ...
    if (bnodes.front() != anodes.back())
        reverse(bnodes.begin(), bnodes.end());

    assert(anodes.back() == bnodes.front());
    for (int i = 1; i < bnodes.size(); i++)
        boundnodes[index++] = bnodes[i];

    // Append cnodes ...
    if (cnodes.front() != bnodes.back())
        reverse(cnodes.begin(), cnodes.end());

    assert(bnodes.back() == cnodes.front());
    for (int i = 1; i < cnodes.size(); i++)
        boundnodes[index++] = cnodes[i];

    // Append dnodes ...
    if (dnodes.front() != cnodes.back())
        reverse(dnodes.begin(), dnodes.end());

    assert(cnodes.back() == dnodes.front());
    for (int i = 1; i < dnodes.size(); i++)
        boundnodes[index++] = dnodes[i];

    // Ensure that loop is closed ...
    assert(anodes.front() == dnodes.back());

    return remesh_quad_loop(mesh, boundnodes, nx, ny, newnodes, newfaces, smooth);
}

////////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_tri_loop(Mesh *mesh,
                      NodeSequence &anodes,
                      NodeSequence &bnodes,
                      NodeSequence &cnodes,
                      int *partition,
                      NodeSequence &newnodes,
                      FaceSequence &newfaces,
                      bool smooth)
{
    // First thing to do is to clear the existing record.
    newnodes.clear();
    newfaces.clear();

    // We need atleast three nodes on each side ...
    if (anodes.size() < 3) return 1;
    if (bnodes.size() < 3) return 1;
    if (cnodes.size() < 3) return 1;

    int segments[3], partSegments[6];

    if (partition == NULL)
    {
        segments[0] = anodes.size() - 1;
        segments[1] = bnodes.size() - 1;
        segments[2] = cnodes.size() - 1;

        if (!Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments))
            return 1;
    }
    else
    {
        for (int i = 0; i < 6; i++)
            partSegments[i] = partition[i];
    }

    int err;
    if (anodes.back() == bnodes.back())
    {
        reverse(bnodes.begin(), bnodes.end());
        swap(partSegments[2], partSegments[3]);
    }

    if (anodes.front() == cnodes.front())
    {
        reverse(cnodes.begin(), cnodes.end());
        swap(partSegments[4], partSegments[5]);
    }

    // Ensure that input is closed loop of three sides ...
    assert(anodes.front() == cnodes.back());
    assert(bnodes.front() == anodes.back());
    assert(cnodes.front() == bnodes.back());

    // Ensure that segments are valid ...
    assert(anodes.size() == partSegments[0] + partSegments[1] + 1);
    assert(bnodes.size() == partSegments[2] + partSegments[3] + 1);
    assert(cnodes.size() == partSegments[4] + partSegments[5] + 1);

    // Split Side "A" nodes into a1nodes and b1nodes.
    NodeSequence a1nodes, b1nodes;
    int a1 = partSegments[0];
    int b1 = partSegments[1];
    err = split_stl_vector(anodes, a1 + 1, a1nodes, b1nodes);
    if (err) return 1;

    // Split Side "B" nodes into a2nodes and b2nodes.
    NodeSequence a2nodes, b2nodes;
    int a2 = partSegments[2];
    int b2 = partSegments[3];
    err = split_stl_vector(bnodes, a2 + 1, a2nodes, b2nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a3nodes, b3nodes;
    int a3 = partSegments[4];
    int b3 = partSegments[5];
    err = split_stl_vector(cnodes, a3 + 1, a3nodes, b3nodes);
    if (err) return 1;

    // Splitting nodes on each side
    Vertex *ca = a1nodes.back();
    Vertex *cb = a2nodes.back();
    Vertex *cc = a3nodes.back();
    Vertex *co = Face::centroid(ca, cb, cc);

    NodeSequence oa_nodes = linear_interpolation(co, ca, a2 + 1);
    NodeSequence ob_nodes = linear_interpolation(co, cb, a3 + 1);
    NodeSequence oc_nodes = linear_interpolation(co, cc, a1 + 1);

    mesh->addNode(co);
    newnodes.push_back(co);

    for (size_t i = 1; i < oa_nodes.size() - 1; i++)
    {
        mesh->addNode(oa_nodes[i]);
        newnodes.push_back(oa_nodes[i]);
    }

    for (size_t i = 1; i < ob_nodes.size() - 1; i++)
    {
        mesh->addNode(ob_nodes[i]);
        newnodes.push_back(ob_nodes[i]);
    }

    for (size_t i = 1; i < oc_nodes.size() - 1; i++)
    {
        mesh->addNode(oc_nodes[i]);
        newnodes.push_back(oc_nodes[i]);
    }

    NodeSequence qnodes;
    FaceSequence qfaces;

    err = remesh_quad_loop(mesh, b2nodes, a3nodes, oc_nodes, ob_nodes, qnodes, qfaces, 0);
    if (!err)
    {
        assert(!qfaces.empty());
        for (size_t i = 0; i < qnodes.size(); i++)
            newnodes.push_back(qnodes[i]);
        for (size_t i = 0; i < qfaces.size(); i++)
            newfaces.push_back(qfaces[i]);
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a1nodes, oa_nodes, oc_nodes, b3nodes, qnodes, qfaces, 0);
        if (!err)
        {
            assert(!qfaces.empty());
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a2nodes, ob_nodes, oa_nodes, b1nodes, qnodes, qfaces, 0);
        if (!err)
        {
            assert(!qfaces.empty());
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (err)
    {
        for (size_t i = 0; i < newfaces.size(); i++)
            mesh->remove(newfaces[i]);
        for (size_t i = 0; i < newnodes.size(); i++)
            mesh->remove(newnodes[i]);
        newnodes.clear();
        newfaces.clear();
        return 2;
    }

    if (smooth)
    {
        LaplaceLengthWeight lw;
        LaplaceSmoothing lapsmooth(mesh);
        lapsmooth.setWeight(&lw);
        lapsmooth.setNumIterations(10);
        lapsmooth.localized_at(newnodes);
    }

    // Check for any inversion of the element, if there is inversion,
    // undo everthing (i.e. remove new nodes and faces).
    //
    for (size_t i = 0; i < newfaces.size(); i++)
    {
        if (newfaces[i]->invertedAt() >= 0)
        {
            for (size_t i = 0; i < newfaces.size(); i++)
                mesh->remove(newfaces[i]);
            for (size_t i = 0; i < newnodes.size(); i++)
                mesh->remove(newnodes[i]);
            newnodes.clear();
            newfaces.clear();
            return 1;
        }
    }


    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_penta_loop(Mesh *mesh,
                        NodeSequence &anodes,
                        NodeSequence &bnodes,
                        NodeSequence &cnodes,
                        NodeSequence &dnodes,
                        NodeSequence &enodes,
                        int *partition, NodeSequence &newnodes,
                        FaceSequence &newfaces, bool smooth)
{
    // First clear the existing records.
    newnodes.clear();
    newfaces.clear();

    int segments[5], partSegments[10];

    if (partition == NULL)
    {
        segments[0] = anodes.size() - 1;
        segments[1] = bnodes.size() - 1;
        segments[2] = cnodes.size() - 1;
        segments[3] = dnodes.size() - 1;
        segments[4] = enodes.size() - 1;
        if (!Face::is_5_sided_convex_loop_quad_meshable(segments, partSegments))
        {
            cout << "Warning: Triangle patch not quad meshable " << endl;
            return 1;
        }
    }
    else
    {
        for (int i = 0; i < 10; i++)
            partSegments[i] = partition[i];
    }

    int err;
    int index = 0;

    // Ensure that input is closed loop of three sides ...
    assert(anodes.front() == enodes.back());
    assert(bnodes.front() == anodes.back());
    assert(cnodes.front() == bnodes.back());
    assert(dnodes.front() == cnodes.back());
    assert(enodes.front() == dnodes.back());

    // Ensure that segments are valid ...
    assert(anodes.size() == partSegments[0] + partSegments[1] + 1);
    assert(bnodes.size() == partSegments[2] + partSegments[3] + 1);
    assert(cnodes.size() == partSegments[4] + partSegments[5] + 1);
    assert(dnodes.size() == partSegments[6] + partSegments[7] + 1);
    assert(enodes.size() == partSegments[8] + partSegments[9] + 1);

    // Split Side "A" nodes into a1nodes and b1nodes.
    NodeSequence a1nodes, b1nodes;
    int a1 = partSegments[0];
    int b1 = partSegments[1];
    assert(a1 > 0 && b1 > 0);
    err = split_stl_vector(anodes, a1 + 1, a1nodes, b1nodes);
    if (err) return 1;

    // Split Side "B" nodes into a2nodes and b2nodes.
    NodeSequence a2nodes, b2nodes;
    int a2 = partSegments[2];
    int b2 = partSegments[3];
    assert(a2 > 0 && b2 > 0);
    err = split_stl_vector(bnodes, a2 + 1, a2nodes, b2nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a3nodes, b3nodes;
    int a3 = partSegments[4];
    int b3 = partSegments[5];
    assert(a3 > 0 && b3 > 0);
    err = split_stl_vector(cnodes, a3 + 1, a3nodes, b3nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a4nodes, b4nodes;
    int a4 = partSegments[6];
    int b4 = partSegments[7];
    assert(a4 > 0 && b4 > 0);
    split_stl_vector(dnodes, a4 + 1, a4nodes, b4nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a5nodes, b5nodes;
    int a5 = partSegments[8];
    int b5 = partSegments[9];
    assert(a5 > 0 && b5 > 0);
    err = split_stl_vector(enodes, a5 + 1, a5nodes, b5nodes);
    if (err) return 1;

    assert(a1 == b4);
    assert(b1 == a3);

    assert(a2 == b5);
    assert(b2 == a4);

    assert(b3 == a5);

    // Splitting nodes on each side
    Vertex *ca = a1nodes.back();
    Vertex *cb = a2nodes.back();
    Vertex *cc = a3nodes.back();
    Vertex *cd = a4nodes.back();
    Vertex *ce = a5nodes.back();
    Vertex *co = Face::centroid(ca, cb, cc, cd, ce);

    NodeSequence oa_nodes = linear_interpolation(co, ca, a2 + 1);
    NodeSequence ob_nodes = linear_interpolation(co, cb, a3 + 1);
    NodeSequence oc_nodes = linear_interpolation(co, cc, a4 + 1);
    NodeSequence od_nodes = linear_interpolation(co, cd, a5 + 1);
    NodeSequence oe_nodes = linear_interpolation(co, ce, a1 + 1);

    mesh->addNode(co);
    newnodes.push_back(co);

    for (size_t i = 1; i < oa_nodes.size() - 1; i++)
    {
        mesh->addNode(oa_nodes[i]);
        newnodes.push_back(oa_nodes[i]);
    }

    for (size_t i = 1; i < ob_nodes.size() - 1; i++)
    {
        mesh->addNode(ob_nodes[i]);
        newnodes.push_back(ob_nodes[i]);
    }

    for (size_t i = 1; i < oc_nodes.size() - 1; i++)
    {
        mesh->addNode(oc_nodes[i]);
        newnodes.push_back(oc_nodes[i]);
    }

    for (size_t i = 1; i < od_nodes.size() - 1; i++)
    {
        mesh->addNode(od_nodes[i]);
        newnodes.push_back(od_nodes[i]);
    }

    for (size_t i = 1; i < oe_nodes.size() - 1; i++)
    {
        mesh->addNode(oe_nodes[i]);
        newnodes.push_back(oe_nodes[i]);
    }

    NodeSequence qnodes;
    FaceSequence qfaces;

    err = remesh_quad_loop(mesh, a1nodes, oa_nodes, oe_nodes, b5nodes, qnodes, qfaces, 0);
    if (!err)
    {
        for (size_t i = 0; i < qnodes.size(); i++)
            newnodes.push_back(qnodes[i]);
        for (size_t i = 0; i < qfaces.size(); i++)
            newfaces.push_back(qfaces[i]);
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a2nodes, ob_nodes, oa_nodes, b1nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a3nodes, oc_nodes, ob_nodes, b2nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a4nodes, od_nodes, oc_nodes, b3nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a5nodes, oe_nodes, od_nodes, b4nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (err) return 2;

    if (smooth)
    {
        LaplaceLengthWeight lw;
        LaplaceSmoothing lapsmooth(mesh);
        lapsmooth.setWeight(&lw);
        lapsmooth.setNumIterations(10);
        lapsmooth.localized_at(newnodes);
    }

    // Check for any inversion of the element, if there is inversion,
    // undo everthing (i.e. remove new nodes and faces).
    //
    for (size_t i = 0; i < newfaces.size(); i++)
    {
        if (newfaces[i]->invertedAt() >= 0)
        {
            for (size_t i = 0; i < newnodes.size(); i++)
                mesh->remove(newnodes[i]);
            for (size_t i = 0; i < newfaces.size(); i++)
                mesh->remove(newfaces[i]);
            newnodes.clear();
            newfaces.clear();
            return 1;
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

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
    FaceSequence adjFaces;
    NodeSequence vneighs;

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

///////////////////////////////////////////////////////////////////////////////

void
QTrack::advance(int endat)
{
    assert(sequence.size() == 2);

    // Starting node is always irregular ...
    FaceSequence vfaces = sequence[0]->getRelations2();
    assert(vfaces.size() != 4);

    while (1)
    {
        int progress = advance_single_step(endat);
        if (!progress) break;
    }

    Vertex *endvertex;

    // Sanity Checking ....
    if (endat == END_AT_TERMINALS)
    {
        endvertex = sequence.front();
        if (!endvertex->isBoundary())
        {
            vfaces = endvertex->getRelations2();
            assert(vfaces.size() != 4);
        }

        endvertex = sequence.back();
        if (!endvertex->isBoundary())
        {
            vfaces = endvertex->getRelations2();
            assert(vfaces.size() != 4);
        }

        for (int i = 1; i < sequence.size() - 1; i++)
        {
            vfaces = sequence[i]->getRelations2();
            assert(!sequence[i]->isBoundary());
            assert(vfaces.size() == 4);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

vector<QTrack> Jaal::generate_quad_irregular_graph(Mesh *mesh)
{
    vector<QTrack> qpath;
    int nTopo = mesh->isHomogeneous();
    if (nTopo != 4)
    {
        cout << "Error: The mesh must be all Quads " << endl;
        return qpath;
    }

    int relexist2 = mesh->build_relations(0, 2);
    int relexist0 = mesh->build_relations(0, 0);

    mesh->search_boundary();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
    }

    QTrack qp;
    qp.mesh = mesh;

    int found;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (!vertex->isBoundary())
        {
            NodeSequence vnodes = vertex->getRelations0();
            if (vnodes.size() != 4)
            {
                for (size_t j = 0; j < vnodes.size(); j++)
                {
                    qp.sequence.resize(2); // As we know the starting edge
                    qp.sequence[0] = vertex;
                    qp.sequence[1] = vnodes[j];
                    qp.advance(0); // Terminate at irregular nodes only..
                    found = 0;
                    for (int k = 0; k < qpath.size(); k++)
                    {
                        if (qpath[k] == qp)
                        {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        max_tag_val = max(max_tag_val, vertex->getLayerID());
    }

    for (size_t i = 0; i < numnodes; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        max_tag_val = max(max_tag_val, face->getLayerID());
    }

    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setTag(1);
    }

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(0);
        if (face->isConvex())
            face->setTag(1);
        else
        {
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *v = face->getNodeAt(j);
                v->setTag(0);
            }
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

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
            vertex->setTag(2);
        else
            vertex->setTag(1);
    }

    size_t numfaces = mesh->getSize(2);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(1);
    }

    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (size_t i = 0; i < numnodes; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        sumarea += face->getArea();
    }
    double avgarea = sumarea / (double) numfaces;
    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setTag(1);
    }

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if (face->getArea() < tolarea)
        {
            face->setTag(0);
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *vertex = face->getNodeAt(j);
                vertex->setTag(0);
            }
        }
        else
            face->setTag(1);
    }
}

///////////////////////////////////////////////////////////////////////////////

bool
Face::has_all_bound_nodes() const
{
    for (int i = 0; i < getSize(0); i++)
    {
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

    for (size_t i = 0; i < faces.size(); i++)
    {
        int nsize = faces[i]->getSize(0);
        for (int j = 0; j < nsize; j++)
        {
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

NodeSequence Mesh::get_breadth_first_ordered_nodes(Vertex *vstart, MeshFilter *filter)
{
    assert(vstart != NULL);

    NodeSequence seq;

    int relexist0 = build_relations(0, 0);

    size_t numnodes = getSize(0);
    if (numnodes == 0) return seq;

    for (size_t i = 0; i < numnodes; i++)
    {
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
    while (!vertexQ.empty())
    {
        Vertex *curr_vertex = vertexQ.front();
        vertexQ.pop_front();
        int currlevel = curr_vertex->getLayerID();
        if (filter)
        {
            if (curr_vertex != vstart) proceed = filter->pass(curr_vertex);
        }
        if (!curr_vertex->isVisited())
        {
            seq.push_back(curr_vertex);
            if (!proceed) break;
            curr_vertex->setVisitMark(1);
            vneighs = curr_vertex->getRelations0();
            for (size_t i = 0; i < vneighs.size(); i++)
            {
                if (!vneighs[i]->isVisited())
                {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setVisitMark(0);
    }

    NodeSequence seq;

    if (vstart == 0) vstart = getNodeAt(0);
    vertexQ.push_back(vstart);
    NodeSequence vneighs;

    while (!vertexQ.empty())
    {
        Vertex *curr_vertex = vertexQ.front();
        vertexQ.pop_front();
        if (!curr_vertex->isVisited())
        {
            seq.push_back(curr_vertex);
            curr_vertex->setVisitMark(1);
            vneighs = curr_vertex->getRelations0();
            for (size_t i = 0; i < vneighs.size(); i++)
            {
                if (!vneighs[i]->isVisited())
                    vertexQ.push_front(vneighs[i]);
            }
        }
    }
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        assert(v->isVisited());
    }
    assert(seq.size() == numnodes);

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

/*
bool
Face::isValid() const
{
    if (isRemoved()) return 0;

    for (int i = 0; i < connect.size(); i++)
        if (connect[i]->isRemoved()) return 0;

    for (int i = 0; i < connect.size(); i++)
    {
        int ncount = 0;
        for (int j = 0; j < connect.size(); j++)
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
    for (int i = 0; i < nSize; i++)
    {
        Vertex *v0 = connect[(i + 0) % nSize];
        Vertex *v1 = connect[(i + 1) % nSize];
        if (v0->isBoundary() && v1->isBoundary())
        {
            FaceSequence neighs = Mesh::getRelations112(v0, v1);
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

    assert(xi >= -1.0 && xi <= 1.0);
    assert(eta >= -1.0 && eta <= 1.0);

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
    int numNodes = x.size();
    assert(w.size() == numNodes);

    double sum = 0.0;
    for (int i = 0; i < numNodes; i++) sum += x[i] * w[i];

    return sum;
}

/////////////////////////////////////////////////////////////////////////////////////


int
Face::invertedAt() const
{
    int nsize = getSize(0);
    double x[5], y[5], z[5], triarea;
    Vertex *vertex;
    Point3D xyz;

    for (int i = 0; i < nsize; i++)
    {
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
 */

#ifdef CSV

/////////////////////////////////////////////////////////////////////////////////////

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
    for (int i = 0; i < 4; i++)
    {
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
    //		   1  0  0  1  0  0;
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
    for (int i = 0; i < 6; i++)
    {
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
    for (int i = 0; i < 10; i++)
    {
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

int
Mesh::refine_quad15(Face *face)
{
    NodeSequence newnodes;
    FaceSequence newfaces;
    face->refine_quad15(newnodes, newfaces);

    for (int i = 0; i < newnodes.size(); i++)
        addNode(newnodes[i]);

    for (int i = 0; i < newfaces.size(); i++)
        addFace(newfaces[i]);

    remove(face);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

Jaal::NodeSequence
Mesh::boundary_chain_nodes(Vertex *v0, Vertex *v1)
{
    NodeSequence bndnodes;

    FaceSequence neighs = Mesh::getRelations112(v0, v1);

    if (neighs.size() != 2) return bndnodes;

    vector<Edge> bndedges;

    bndedges.reserve(6);
    Edge sharededge(v0, v1);

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
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

    for (int i = 1; i < bndedges.size() - 1; i++)
    {
        Vertex *v0 = bndedges[i].getNodeAt(0);
        Vertex *v1 = bndedges[i].getNodeAt(1);
        if (v0 == bndnodes[i])
        {
            bndnodes.push_back(v1);
        }
        else if (v1 == bndnodes[i])
            bndnodes.push_back(v0);
        else
        {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        area[i] = face->getArea();
    }
    minarea = *min_element(area.begin(), area.end());
    maxarea = *max_element(area.begin(), area.end());
}


////////////////////////////////////////////////////////////////////////////////

const vector<double>
Mesh::getCoordsArray()
{
    size_t numnodes = getSize(0);

    vector<double> vcoords(3 * numnodes);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        Point3D xyz = v->getXYZCoords();
        vcoords[3 * i + 0] = xyz[0];
        vcoords[3 * i + 1] = xyz[1];
        vcoords[3 * i + 2] = xyz[2];
    }
    return vcoords;
}

////////////////////////////////////////////////////////////////////////////////

const vector<size_t>
Mesh::getNodesArray()
{
    size_t numfaces = getSize(2);

    vector<size_t> nodearray;

    int topo = isHomogeneous();
    if (topo) nodearray.reserve(topo * numfaces);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        for (int j = 0; j < face->getSize(0); j++)
        {
            Vertex *v = face->getNodeAt(j);
            nodearray.push_back(v->getID());
        }
    }

    return nodearray;
}

////////////////////////////////////////////////////////////////////////////////

Mesh*
Mesh::deep_copy()
{
    std::map<Vertex*, Vertex*> vmap;

    Mesh *newmesh = new Mesh;
    size_t numnodes = getSize(0);
    newmesh->reserve(numnodes, 0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vold = getNodeAt(i);
        Vertex *vnew = Vertex::newObject();
        vnew->setXYZCoords(vold->getXYZCoords());
        vmap[vold] = vnew;
        newmesh->addNode(vnew);
    }

    size_t numfaces = getSize(2);
    newmesh->reserve(numnodes, 2);
    NodeSequence connect;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *fold = getFaceAt(i);
        connect = fold->getNodes();
        for (int j = 0; j < connect.size(); j++)
            connect[j] = vmap[connect[j]];

        Face *fnew = Face::newObject();
        fnew->setNodes(connect);
        newmesh->addFace(fnew);
    }

    return newmesh;
}

////////////////////////////////////////////////////////////////////////////////

int
Mesh::setCoordsArray(const vector<double> &vcoords)
{
    size_t numnodes = getSize(0);

    if (vcoords.size() != 3 * numnodes) return 1;

    Point3D xyz;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        xyz[0] = vcoords[3 * i + 0];
        xyz[1] = vcoords[3 * i + 1];
        xyz[2] = vcoords[3 * i + 2];
        v->setXYZCoords(xyz);
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
    for (size_t i = 1; i < boundedges.size(); i++)
        listedges.push_back(boundedges[i]);
    boundedges.clear();

    boundedges.reserve(nSize);

    Vertex *first_vertex = edge.getNodeAt(0);
    Vertex *curr_vertex = edge.getNodeAt(1);

    boundedges.push_back(edge);

    list<Edge>::iterator it;

    for (size_t i = 0; i < nSize; i++)
    {
        for (it = listedges.begin(); it != listedges.end(); ++it)
        {
            edge = *it;
            Vertex *v0 = edge.getNodeAt(0);
            Vertex *v1 = edge.getNodeAt(1);
            if (v0 == curr_vertex)
            {
                curr_vertex = v1;
                boundedges.push_back(edge);
                break;
            }
            if (v1 == curr_vertex)
            {
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

    for (size_t i = 0; i < boundedges.size(); i++)
    {
        Vertex *v0 = boundedges[i].getNodeAt(0);
        Vertex *v1 = boundedges[i].getNodeAt(1);
        relations00[v0].insert(v1);
        relations00[v1].insert(v0);
    }

    std::map<Vertex*, set<Vertex*> > ::const_iterator it;
    for (it = relations00.begin(); it != relations00.end(); ++it)
    {
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

    for (size_t i = 1; i < boundedges.size(); i++)
    {
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
    for (size_t i = 0; i < nSize; i++)
    {
        if (boundedges[i].getNodeAt(0) == first_vertex) istart = i;
        listedges[i] = boundedges[i];
    }

    for (size_t i = 0; i < nSize; i++)
        boundedges[i] = listedges[(i + istart) % nSize];

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

FaceSequence
Mesh::getRelations102(PNode vtx0, PNode vtx1)
{
    FaceSequence faceneighs;

    FaceSequence v0faces = vtx0->getRelations2();
    FaceSequence v1faces = vtx1->getRelations2();

    if (v0faces.empty() || v1faces.empty())
    {
        cout << "Warning: Vertex-Faces relations are empty " << endl;
        return faceneighs;
    }

    FaceSet vset;
    for (size_t i = 0; i < v0faces.size(); i++)
        vset.insert(v0faces[i]);

    for (size_t i = 0; i < v1faces.size(); i++)
        vset.insert(v1faces[i]);

    std::set<PFace>::iterator it;

    if (vset.size())
    {
        faceneighs.resize(vset.size());
        int index = 0;
        for (it = vset.begin(); it != vset.end(); ++it)
            faceneighs[index++] = *it;
    }

    return faceneighs;
}

///////////////////////////////////////////////////////////////////////////////

FaceSequence
Mesh::getRelations112(const PNode vtx0, const PNode vtx1)
{
    FaceSequence faceneighs;
    FaceSequence v0faces = vtx0->getRelations2();
    FaceSequence v1faces = vtx1->getRelations2();

    if (v0faces.empty() || v1faces.empty())
    {
        cout << "Warning: Vertex-Faces relations are empty " << endl;
        return faceneighs;
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

    return faceneighs;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::count_edges()
{
    if (getAdjTable(1, 1)) return edges.size();

    int relexist = build_relations(0, 0);

    size_t numnodes = getSize(0);

    NodeSequence neighs;
    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        neighs = nodes[i]->getRelations0();
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

    NodeSequence neighs;

    edges.clear();

    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        neighs = nodes[i]->getRelations0();
        for (size_t j = 0; j < neighs.size(); j++)
            if (nodes[i] > neighs[j])
            {
                Edge *newedge = new Edge(nodes[i], neighs[j]);
                assert(newedge);
                edges.push_back(newedge);
            }
    }

    if (!relexist0) clear_relations(0, 0);
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::prune()
{
    if (adjTable[0][0])
    {
        NodeSequence relations0;
        for (size_t i = 0; i < nodes.size(); i++)
        {
            relations0 = nodes[i]->getRelations0();
            for (size_t j = 0; j < relations0.size(); j++)
            {
                if (relations0[j]->isRemoved())
                    nodes[i]->removeRelation0(relations0[j]);
            }
        }
    }

    if (adjTable[0][2])
    {
        FaceSequence relations2;
        for (size_t i = 0; i < nodes.size(); i++)
        {
            relations2 = nodes[i]->getRelations2();
            for (size_t j = 0; j < relations2.size(); j++)
            {
                if (relations2[j]->isRemoved())
                    nodes[i]->removeRelation2(relations2[j]);
            }
        }
    }

    for (size_t i = 0; i < nodes.size(); i++)
    {
        Vertex *v = nodes[i];
        if (v->isRemoved())
        {
            if (find(garbageNodes.begin(), garbageNodes.end(), v) == garbageNodes.end())
                garbageNodes.push_back(v);
        }
    }

    for (size_t i = 0; i < faces.size(); i++)
    {
        Face *f = faces[i];
        if (f->isRemoved())
        {
            if (find(garbageFaces.begin(), garbageFaces.end(), f) == garbageFaces.end())
                garbageFaces.push_back(f);
        }
    }

    NodeSequence::iterator vend;
    vend = remove_if(nodes.begin(), nodes.end(), EntityRemovedPred());
    nodes.erase(vend, nodes.end());

    FaceSequence::iterator fend;
    fend = remove_if(faces.begin(), faces.end(), EntityRemovedPred());
    faces.erase(fend, faces.end());

    enumerate(0);
    enumerate(2);

    for (size_t i = 0; i < nodes.size(); i++)
        assert(nodes[i]->isActive());

    for (size_t i = 0; i < faces.size(); i++)
        assert(faces[i]->isActive());
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
    list<Face*>::const_iterator fiter;
    for (fiter = garbageFaces.begin(); fiter != garbageFaces.end(); ++fiter)
    {
        Face *face = *fiter;
        assert(face);
        if (face->isRemoved()) delete face;
    }
    garbageFaces.clear();

    list<Vertex*>::const_iterator viter;
    for (viter = garbageNodes.begin(); viter != garbageNodes.end(); ++viter)
    {
        Vertex *vertex = *viter;
        assert(vertex);
        if (vertex->isRemoved()) delete vertex;
    }
    garbageNodes.clear();
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::enumerate(int etype)
{
    size_t index = 0;

    NodeSequence::const_iterator viter;
    if (etype == 0)
    {
        index = 0;
        for (viter = nodes.begin(); viter != nodes.end(); ++viter)
        {
            Vertex *vertex = *viter;
            vertex->setID(index++);
        }
    }

    FaceSequence::const_iterator fiter;
    if (etype == 2)
    {
        index = 0;
        for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
        {
            Face *face = *fiter;
            face->setID(index++);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::getBoundarySize(int d) const
{
    size_t ncount = 0;

    if (d == 0)
    {
        for (size_t i = 0; i < nodes.size(); i++)
            if (!nodes[i]->isRemoved() && nodes[i]->isBoundary())
                ncount++;
    }

    if (d == 2)
    {
        for (size_t i = 0; i < faces.size(); i++)
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
    for (size_t i = 0; i < faces.size(); i++)
    {
        if (!faces[i]->isRemoved())
            maxnodes = max(maxnodes, faces[i]->getSize(0));
    }

    return maxnodes;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::saveAs(const string &s)
{
    if (s.rfind(".off") != string::npos)
        save_off_format(s);

    if (s.rfind(".dat") != string::npos)
        save_simple_format(s);
}
///////////////////////////////////////////////////////////////////////////////

Vertex*
Mesh::nearest_neighbour(const Vertex *myself, double &mindist)
{
    Vertex* nearest;
    assert(getAdjTable(0, 0));

    mindist = MAXDOUBLE;
    nearest = NULL;

    NodeSequence neighs = myself->getRelations0();

    for (size_t i = 0; i < neighs.size(); i++)
    {
        double d = Vertex::length2(myself, neighs[i]);
        if (d < mindist)
        {
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
    if (!isPruned()) prune();

    size_t numfaces = getSize(2);

    if (getAdjTable(0, 2))
    {
        for (size_t iface = 0; iface < numfaces; iface++)
        {
            Face *face = getFaceAt(iface);
            assert(!face->isRemoved());
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *vtx = face->getNodeAt(j);
                assert(!vtx->isRemoved());
                assert(vtx->hasRelation2(face));
            }
        }
    }

    // Delete all the old stuff and rebuild new. Sometimes used for debugging purpose
    // also.
    if (rebuild) clear_relations(0, 2);

    if (adjTable[0][2] == 1) return 1;

    clear_relations(0, 2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        for (int j = 0; j < face->getSize(0); j++)
        {
            Vertex *vtx = face->getNodeAt(j);
            vtx->addRelation2(face);
        }
    }
    adjTable[0][2] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::build_relations00(bool rebuild)
{
    if (!isPruned()) prune();

    // Delete all the old stuff and rebuild new. Sometimes used for debugging purpose
    // also.
    if (rebuild) clear_relations(0, 0);

    if (adjTable[0][0] == 1) return 1;

    clear_relations(0, 0);

    size_t numfaces = getSize(2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        size_t nnodes = face->getSize(0);
        for (size_t j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            v0->addRelation0(v1);
            v1->addRelation0(v0);
        }
    }
    adjTable[0][0] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::clear_relations(int src, int dst)
{
    size_t numnodes = getSize(0);

    if (src == 0 && dst == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vtx = getNodeAt(i);
            vtx->clearRelations(0);
        }
        adjTable[0][0] = 0;
    }

    if (src == 0 && dst == 2)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
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

    size_t numnodes = getSize(0);
    size_t numfaces = getSize(2);

    int bmark;
    FaceSequence neighs;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        size_t nnodes = face->getSize(0);
        for (size_t j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            neighs = Mesh::getRelations112(v0, v1);
            if (neighs.size() == 1)
            {
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
    for (it = faces.begin(); it != faces.end(); ++it)
    {
        Face *face = *it;
        if (face->getType() == facetype)
            ncount++;
    }

    FaceSequence tmpfaces;
    if (ncount)
    {
        tmpfaces.resize(ncount);
        size_t index = 0;
        for (it = faces.begin(); it != faces.end(); ++it)
        {
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
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            neighs = Mesh::getRelations112(v0, v1);
            if (neighs.size() > 2)
            {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = getNodeAt(i);
        if (!vertex->isBoundary())
        {
            FaceSequence vfaces = vertex->getRelations2();
            if (vfaces.size() != degree_of_regular_node)
                ncount++;
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
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        assert(face);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt(j);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            neighs = Mesh::getRelations112(v0, v1);
            if (neighs.size() == 2)
            {
                assert(!neighs[0]->isRemoved());
                assert(!neighs[1]->isRemoved());
                int dir1 = neighs[0]->getOrientation(v0, v1);
                int dir2 = neighs[1]->getOrientation(v0, v1);
                if (dir1 * dir2 == 1)
                {
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

    if (!relexist)
        clear_relations(0, 2);

    return consistent;
}

///////////////////////////////////////////////////////////////////////////////

void
Mesh::make_consistently_oriented()
{
    build_relations(0, 2);

    Face *face = NULL;
    deque<Face*> faceQ;
    FaceSequence neighs;

    size_t numfaces = getSize(2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        face = getFaceAt(iface);
        face->setID(iface);
        face->setVisitMark(0);
    }

    face = getFaceAt(0);
    faceQ.push_back(face);

    while (!faceQ.empty())
    {
        Face *face = faceQ.front();
        faceQ.pop_front();
        if (!face->isVisited())
        {
            face->setVisitMark(1);
            int nnodes = face->getSize(0);
            for (int j = 0; j < nnodes; j++)
            {
                Vertex *v0 = face->getNodeAt(j);
                Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
                neighs = Mesh::getRelations112(v0, v1);
                if (neighs.size() == 2)
                {
                    int dir1 = neighs[0]->getOrientation(v0, v1);
                    int dir2 = neighs[1]->getOrientation(v0, v1);
                    if (dir1 * dir2 == 1)
                    {
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

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        face = getFaceAt(iface);
        if (!face->isVisited())
            cout << "Error: not visited : " << face->getID() << " "
            << face->isVisited() << endl;
    }

    clear_relations(0, 2);
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

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        face = getFaceAt(iface);
        face->setID(iface);
        face->setVisitMark(0);
    }

    int numComponents = 0;

    while (1)
    {
        face = NULL;
        faceQ.clear();
        for (size_t iface = 0; iface < numfaces; iface++)
        {
            face = getFaceAt(iface);
            if (!face->isVisited())
            {
                face->setPartID(numComponents);
                faceQ.push_back(face);
                break;
            }
        }

        if (faceQ.empty())
            break;

        while (!faceQ.empty())
        {
            Face *face = faceQ.front();
            faceQ.pop_front();
            if (!face->isVisited())
            {
                face->setPartID(numComponents);
                face->setVisitMark(1);
                int nnodes = face->getSize(0);
                for (int j = 0; j < nnodes; j++)
                {
                    Vertex *v0 = face->getNodeAt(j);
                    Vertex *v1 = face->getNodeAt((j + 1) % nnodes);

                    int proceed = 1;
                    if (stop_at_interface)
                    {
                        if (v0->isConstrained() && v1->isConstrained())
                        {
                            Edge edge(v0, v1);
                            if (hasFeatureEdge(edge)) proceed = 0;
                        }
                    }

                    if (proceed)
                    {
                        neighs = Mesh::getRelations112(v0, v1);
                        if (neighs.size() == 2)
                        {
                            faceQ.push_back(neighs[0]);
                            faceQ.push_back(neighs[1]);
                        }
                    }

                }
            }
        } // Complete one Component
        numComponents++;
    }

    for (size_t iface = 0; iface < numfaces; iface++)
    {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setVisitMark(0);
    }

    size_t numfaces = getSize(2);

    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        if (face->getPartID() == id)
        {
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *v = face->getNodeAt(j);
                if (!v->isVisited())
                {
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
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
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
    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
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
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
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
    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
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
    FaceSequence neighs = Mesh::getRelations112(v0, v1);

    Vertex *vn0, *vn1; // Next-Edge Nodes;

    Face *nextface;
    if (neighs.size() == 2)
    {
        if (!neighs[0]->isVisited() && neighs[1]->isVisited())
        {
            nextface = neighs[0];
            if (nextface->getSize(0) == 4)
            {
                Face::opposite_nodes(nextface, v0, v1, vn0, vn1);
                nextface->setVisitMark(1);
                strip.push_back(nextface);
                expand_strip(nextface, vn0, vn1, strip);
            }
        }
        if (neighs[0]->isVisited() && !neighs[1]->isVisited())
        {
            nextface = neighs[1];
            if (nextface->getSize(0) == 4)
            {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
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

    if (bound_what == 0)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            if (face->hasBoundaryNode())
                bfaces.insert(face);
        }
    }

    if (bound_what == 1)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            if (face->has_boundary_edge())
                bfaces.insert(face);
        }
    }

    FaceSequence result;

    size_t nSize = bfaces.size();

    if (nSize)
    {
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
    for (size_t i = 0; i < numnodes; i++)
    {
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

    NodeSequence seq;
    FaceSequence vfaces;

    // Query from the boundary nodes ...
    if (from_where == 1)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *v = getNodeAt(i);
            vfaces = v->getRelations2();
            if (v->isBoundary() && vfaces.size() != regular_count) seq.push_back(v);
        }
    }

    // Query from the internal nodes ...
    if (from_where == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *v = getNodeAt(i);
            vfaces = v->getRelations2();
            if (!v->isBoundary() && vfaces.size() != regular_count) seq.push_back(v);
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        face->setTag(0);
    }

    FaceSequence strip1, strip2;
    int id = 0;
    for (size_t i = 0; i < bound_faces.size(); i++)
    {
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

    vector<int> degree(numnodes);
    FaceSequence neighs;
    for (int i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        neighs = v->getRelations2();
        degree[i] = neighs.size();
    }

    int mindegree = *min_element(degree.begin(), degree.end());
    int maxdegree = *max_element(degree.begin(), degree.end());

    cout << " Mesh Topological Quality : " << endl;

    cout << " *********************** " << endl;
    cout << " Degree        Count " << endl;
    cout << " *********************** " << endl;

    for (int i = mindegree; i <= maxdegree; i++)
    {
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

    if (layerid == 0)
    {
        search_boundary();
        for (size_t i = 0; i < numNodes; i++)
        {
            Vertex *vertex = getNodeAt(i);
            if (vertex->isBoundary())
            {
                vertex->setLayerID(0);
                ncount++;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < numNodes; i++)
        {
            Vertex *vertex = getNodeAt(i);
            if (vertex->getLayerID() == layerid - 1)
            {
                NodeSequence vnodes = vertex->getRelations0();
                for (size_t j = 0; j < vnodes.size(); j++)
                {
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
    for (size_t i = 0; i < numNodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setLayerID(-1);
        v->setVisitMark(0);
    }

    NodeSequence vertexQ, nextQ, vneighs;
    for (size_t i = 0; i < numNodes; i++)
    {
        Vertex *v = getNodeAt(i);
        if (v->isBoundary())
        {
            v->setLayerID(0);
            v->setVisitMark(1);
            vertexQ.push_back(v);
        }
    }

    if (vertexQ.empty())
    {
        cout << "Warning: No boundary detected " << endl;
    }

    int layerid = 1;
    while (!vertexQ.empty())
    {
        nextQ.clear();
#ifdef SEQUENCE_IS_VECTOR
        nextQ.reserve(vertexQ.size());
#endif
        for (size_t j = 0; j < vertexQ.size(); j++)
        {
            Vertex *currVertex = vertexQ[j];
            vneighs = currVertex->getRelations0();
            for (size_t i = 0; i < vneighs.size(); i++)
            {
                Vertex *vn = vneighs[i];
                if (!vn->isVisited()) nextQ.push_back(vn);
            }
        }

        for (size_t j = 0; j < nextQ.size(); j++)
        {
            Vertex *currVertex = nextQ[j];
            currVertex->setLayerID(layerid);
            currVertex->setVisitMark(1);
        }
        vertexQ = nextQ;
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

    if (layerid == 0)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            int nsize = face->getSize(0);
            for (int j = 0; j < nsize; j++)
            {
                Vertex *v0 = face->getNodeAt((j + 0) % nsize);
                Vertex *v1 = face->getNodeAt((j + 1) % nsize);
                FaceSequence vfaces = Mesh::getRelations112(v0, v1);
                if (vfaces.size() == 1)
                {
                    face->setLayerID(0);
                    ncount++;
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            if (face->getLayerID() == layerid - 1)
            {
                FaceSequence vfaces = face->getRelations212();
                for (size_t j = 0; j < vfaces.size(); j++)
                {
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

    for (size_t i = 0; i < numFaces; i++)
    {
        Face *f = getFaceAt(i);
        f->setLayerID(0);
        f->setVisitMark(0);
    }

    FaceSequence faceQ, nextQ;
    for (size_t i = 0; i < numFaces; i++)
    {
        Face *f = getFaceAt(i);
        f->setLayerID(1);
        assert(!f->isRemoved());
        if (f->has_boundary_edge())
        {
            f->setLayerID(0);
            f->setVisitMark(1);
            faceQ.push_back(f);
        }
    }

    int layerid = 1;
    FaceSequence neighs;

    while (!faceQ.empty())
    {
        nextQ.clear();
#ifdef SEQUENCE_IS_VECTOR
        nextQ.reserve(faceQ.size());
#endif
        for (size_t j = 0; j < faceQ.size(); j++)
        {
            Face *currFace = faceQ[j];
            neighs = currFace->getRelations212();
            for (size_t i = 0; i < neighs.size(); i++)
            {
                Face *nf = neighs[i];
                if (!nf->isVisited())
                    nextQ.push_back(nf);
            }
        }

        for (size_t i = 0; i < nextQ.size(); i++)
        {
            Face *f = nextQ[i];
            f->setLayerID(layerid);
            f->setVisitMark(1);
        }

        layerid++;
        faceQ = nextQ;
    }

    for (size_t i = 0; i < numFaces; i++)
    {
        Face *f = getFaceAt(i);
        assert(f->isVisited());
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
    if (mentity == 0)
    {
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            int nsize = face->getSize(0);
            for (int j = 0; j < nsize; j++)
            {
                Vertex *v0 = face->getNodeAt((j + 0) % nsize);
                Vertex *v1 = face->getNodeAt((j + 1) % nsize);
                int l1 = v0->getLayerID();
                int l2 = v1->getLayerID();
                assert(l1 >= 0 && l2 >= 0);
                if (abs(l2 - l1) > 1) return 1;
            }
        }
    }

    if (mentity == 2)
    {
        int relexist2 = build_relations(0, 2);
        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = getFaceAt(i);
            FaceSequence neighs = face->getRelations212();
            int l1 = face->getLayerID();
            assert(l1 >= 0);
            for (int j = 0; j < neighs.size(); j++)
            {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setVisitMark(0);
    }

    size_t numfaces = getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        for (int j = 0; j < face->getSize(0); j++)
        {
            Vertex *v = face->getNodeAt(j);
            v->setVisitMark(1);
        }
    }

    for (size_t i = 0; i < numnodes; i++)
    {
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

#ifdef REMOVE_LATER 

int
Mesh::check_unused_objects()
{
    set<Vertex*> vset;
    for (size_t i = 0; i < nodes.size(); i++)
        vset.insert(nodes[i]);

    if (vset.size() != nodes.size())
        cout << "Warning: There are some duplicated nodes in the mesh " << endl;

    size_t numfaces = getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *f = getFaceAt(i);
        if (!f->isRemoved())
        {
            for (int j = 0; j < f->getSize(0); j++)
            {
                Vertex *v = f->getNodeAt(j);
                if (v->isRemoved())
                    cout << "Goofed up: Face vertex is deleted " << endl;
                if (vset.find(v) == vset.end())
                    cout << "Warning: A face vertex is not in the node list "
                        << v->getID() << endl;
            }
        }
    }
    return 0;
}
#endif

///////////////////////////////////////////////////////////////////////////////

/*
bool
Mesh::isDelaunay()
{
    bool retval = 1;
    int relexist0 = build_relations(0, 0);

    Point3D pa, pb, pc, pd, pCenter;

    size_t numfaces = getSize(2);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        face->setTag(1);
    }
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        pa = face->getNodeAt(0)->getXYZCoords();
        pb = face->getNodeAt(1)->getXYZCoords();
        pc = face->getNodeAt(2)->getXYZCoords();
        TriCircumCenter3D(&pa[0], &pb[0], &pc[0], &pCenter[0]);
        double radius2 = Math::length2(pa, pCenter);
        NodeSequence neighs = face->getRelations0();
        for (size_t j = 0; j < neighs.size(); j++)
        {
            pd = neighs[j]->getXYZCoords();
            if (Math::length2(pd, pCenter) < radius2)
            {
                face->setTag(0);
                retval = 0;
                break;
            }
        }
    }
    if (!relexist0)
        clear_relations(0, 0);
}
 */
///////////////////////////////////////////////////////////////////////////////

double
Mesh::getSurfaceArea()
{
    double facearea, sumArea = 0.0;

    size_t numfaces = getSize(2);
    double minarea = MAXDOUBLE;
    double maxarea = 0.0;
    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (int i = 1; i < bndedges.size() - 1; i++)
    {
        Vertex *v0 = bndedges[i].getNodeAt(0);
        Vertex *v1 = bndedges[i].getNodeAt(1);
        if (v0 == bndnodes[i])
        {
            bndnodes.push_back(v1);
        }
        else if (v1 == bndnodes[i])
            bndnodes.push_back(v0);
        else
        {
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

    for (size_t i = 0; i < numfaces; i++)
    {
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

void
Mesh::clearAll()
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            clear_relations(i, j);

    for (size_t i = 0; i < nodes.size(); i++)
        delete nodes[i];
    nodes.clear();

    for (size_t i = 0; i < edges.size(); i++)
        delete edges[i];
    edges.clear();

    for (size_t i = 0; i < faces.size(); i++)
        delete faces[i];
    faces.clear();

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
    for (size_t i = 0; i < edges.size(); i++)
    {
        Vertex *v0 = edges[i]->getNodeAt(0);
        Vertex *v1 = edges[i]->getNodeAt(1);
        Point3D xyz = Vertex::mid_point(v0, v1);
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
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = this->getFaceAt(iface);
        for (int j = 0; j < 4; j++)
        {
            Vertex *v0 = face->getNodeAt((j + 0) % 4);
            Vertex *v1 = face->getNodeAt((j + 1) % 4);
            Edge edge(v0, v1);
            Vertex *vh = edge.getHashNode();
            ilower = edgemap.lower_bound(vh);
            iupper = edgemap.upper_bound(vh);
            Vertex *vmid = NULL;
            for (iter = ilower; iter != iupper; ++iter)
            {
                Edge *existing_edge = iter->second;
                if (existing_edge->isSameAs(edge))
                {
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
    size_t numfaces = getSize(2);
    FaceSequence newfaces;
    NodeSequence newnodes;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = getFaceAt(iface);
        face->refine_quad15(newnodes, newfaces);
        for (int j = 0; j < newnodes.size(); j++)
            addNode(newnodes[j]);
        for (int j = 0; j < newfaces.size(); j++)
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
/*
double Vertex ::getFeatureAngle() const
{
  FaceSequence vfaces = getRelations2();
  assert( !vfaces.empty() ); 

  double theta = 0.0;
  for( int i = 0; i < vfaces.size(); i++) 
       theta += vfaces[i]->getAngleAt( this );

  return theta;
}
 */

///////////////////////////////////////////////////////////////////////////////

/*
void
Mesh::setFeatureAngles()
{
    int relexist0 = build_relations(0, 0);
    int relexist2 = build_relations(0, 2);

    size_t numnodes = getSize(0);
    vector<Vertex*> vneighs, bndnodes, sidenodes;

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = getNodeAt(i);
        if (vertex->isBoundary())
        {
            vneighs = vertex->getRelations0();
            bndnodes.clear();
            for (int j = 0; j < vneighs.size(); j++)
                if (vneighs[j]->isBoundary()) bndnodes.push_back(vneighs[j]);
            sidenodes.clear();
            for (int j = 0; j < bndnodes.size(); j++)
            {
                vector<Face*> vfaces = Mesh::getRelations112(vertex, bndnodes[j]);
                if (vfaces.size() == 1)
                    sidenodes.push_back(bndnodes[j]);
            }
            assert(sidenodes.size() == 2);
            Point3D p0 = vertex->getXYZCoords();
            Point3D p1 = sidenodes[0]->getXYZCoords();
            Point3D p2 = sidenodes[1]->getXYZCoords();
            Vec3D v1 = Math::create_vector(p1, p0);
            Vec3D v2 = Math::create_vector(p2, p0);
            double angle = Math::getVectorAngle(v1, v2, ANGLE_IN_DEGREES);
            vertex->setFeatureAngle(angle);
        }
    }

    if (!relexist0)
        clear_relations(0, 0);

    if (!relexist2)
        clear_relations(0, 2);
}
 */

///////////////////////////////////////////////////////////////////////////////

double Vertex::getFeatureLength() const
{
    if (!isBoundary()) return MAXDOUBLE;

    NodeSequence vneighs = getRelations0();

    assert(!vneighs.empty());

    double minlen = MAXDOUBLE;
    for (int j = 0; j < vneighs.size(); j++)
    {
        if (vneighs[j]->isBoundary())
        {
            Point3D p0 = getXYZCoords();
            Point3D p1 = vneighs[j]->getXYZCoords();
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

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        Vertex *vertex = face->getHashNode();
        mapfaces[vertex].push_back(face);
    }

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = getFaceAt(i);
        const FaceSequence &hashfaces = mapfaces[face->getHashNode() ];
        size_t ncount = 0;
        for (size_t j = 0; j < hashfaces.size(); j++)
            if (hashfaces[j]->isSameAs(face)) ncount++;
        if (ncount != 1)
        {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *f = getFaceAt(i);
        if (!f->isConvex()) ncount++;
    }
    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

size_t
Mesh::count_inverted_faces()
{
    size_t numfaces = getSize(2);
    size_t ncount = 0;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *f = getFaceAt(i);
        if (f->invertedAt() >= 0) ncount++;
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

    for (int i = 0; i < numnodes; i++)
    {
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

    for (int i = 0; i < numnodes; i++)
    {
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

    for (int i = 0; i < numnodes; i++)
    {
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

    for (int i = 0; i < numnodes; i++)
    {
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

NodeSequence
Jaal::linear_interpolation(Vertex *v0, Vertex *v1, int n)
{
    assert(n >= 2);
    NodeSequence newnodes;

    Point3D xyz0 = v0->getXYZCoords();
    Point3D xyz1 = v1->getXYZCoords();

    newnodes.resize(n);
    newnodes[0] = v0;
    newnodes[n - 1] = v1;

    double dt = 2.0 / (double) (n - 1);

    Point3D xyzt;
    for (int i = 1; i < n - 1; i++)
    {
        double t = -1.0 + i*dt;
        xyzt[0] = TFI::linear_interpolation(t, xyz0[0], xyz1[0]);
        xyzt[1] = TFI::linear_interpolation(t, xyz0[1], xyz1[1]);
        xyzt[2] = TFI::linear_interpolation(t, xyz0[2], xyz1[2]);
        newnodes[i] = Vertex::newObject();
        newnodes[i]->setXYZCoords(xyzt);
    }

    return newnodes;
}

///////////////////////////////////////////////////////////////////////////////

void
set_tfi_coords(int i, int j, int nx, int ny, vector<Vertex*> &qnodes)
{
    int offset;

    offset = 0;
    const Point3D v00 = qnodes[offset]->getXYZCoords();

    offset = i;
    const Point3D vr0 = qnodes[offset]->getXYZCoords();

    offset = (nx - 1);
    const Point3D v10 = qnodes[offset]->getXYZCoords();

    offset = j*nx;
    const Point3D v0s = qnodes[offset]->getXYZCoords();

    offset = j * nx + (nx - 1);
    const Point3D v1s = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx;
    const Point3D v01 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + i;
    const Point3D vr1 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + (nx - 1);
    const Point3D v11 = qnodes[offset]->getXYZCoords();

    Point3D vrs;

    double dr = 2.0 / (double) (nx - 1);
    double ds = 2.0 / (double) (ny - 1);

    double r = -1.0 + i*dr;
    double s = -1.0 + j*ds;
    for (int k = 0; k < 3; k++)
    {
        vrs[k] = TFI::transfinite_blend(r, s,
                                        v00[k], v10[k], v11[k], v01[k],
                                        vr0[k], v1s[k], vr1[k], v0s[k]);
    }
    offset = j * nx + i;
    qnodes[offset]->setXYZCoords(vrs);

}
///////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_quad_loop(Mesh *mesh,
                       NodeSequence &boundnodes, int nx, int ny,
                       NodeSequence &newnodes, FaceSequence &newfaces,
                       bool smooth)
{

    newnodes.clear();
    newfaces.clear();

    assert(boundnodes.size() == 2 * nx + 2 * ny - 4);
    for (size_t i = 0; i < boundnodes.size(); i++)
        assert(!boundnodes[i]->isRemoved());

    vector<Vertex*> qnodes(nx * ny);

    //
    // Put the boundary nodes on the structured mesh: The orientation is
    // Counter clockwise ( south->east->north->west );
    //

    int offset, index = 0;

    // South Side ...
    index = 0;
    for (int i = 0; i < nx; i++)
    {
        offset = i;
        qnodes[i] = boundnodes[index++];
        Point3D xyz = qnodes[i]->getXYZCoords();
    }

    // East Side
    for (int j = 1; j < ny; j++)
    {
        offset = j * nx + (nx - 1);
        qnodes[offset] = boundnodes[index++];
    }

    // North Side
    for (int i = nx - 2; i >= 0; i--)
    {
        offset = (ny - 1) * nx + i;
        qnodes[offset] = boundnodes[index++];
    }

    // West Side
    for (int j = ny - 2; j >= 1; j--)
    {
        offset = j*nx;
        qnodes[j * nx] = boundnodes[index++];
    }

    // Now all internal nodes ....
    newnodes.resize((nx - 2)*(ny - 2));

    index = 0;
    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            Vertex *v = Vertex::newObject();
            offset = j * nx + i;
            qnodes[offset] = v;
            newnodes[index++] = v;
            set_tfi_coords(i, j, nx, ny, qnodes); // Coordinates values
        }
    }

    newfaces.resize((nx - 1)*(ny - 1));
    NodeSequence qc(4);

    // Create new faces ...
    index = 0;
    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
            int offset = j * nx + i;
            qc[0] = qnodes[offset];
            qc[1] = qnodes[offset + 1];
            qc[2] = qnodes[offset + 1 + nx];
            qc[3] = qnodes[offset + nx];
            Face *face = Face::newObject();
            face->setNodes(qc);
            newfaces[index++] = face;
        }
    }

    // Update the mesh ...
    for (size_t i = 0; i < newnodes.size(); i++)
        mesh->addNode(newnodes[i]);

    for (size_t i = 0; i < newfaces.size(); i++)
        mesh->addFace(newfaces[i]);

    if (smooth)
    {
        // Perform some laplacian smoothing inside the local mesh...
        LaplaceLengthWeight lw;
        LaplaceSmoothing lapsmooth(mesh);
        lapsmooth.setWeight(&lw);
        lapsmooth.setNumIterations(10);
        lapsmooth.localized_at(newnodes);
    }

    // Check for any inversion of the element, if there is inversion,
    // undo everthing (i.e. remove new nodes and faces).
    //
    for (size_t i = 0; i < newfaces.size(); i++)
    {
        if (newfaces[i]->invertedAt() >= 0)
        {
            for (size_t i = 0; i < newfaces.size(); i++)
                mesh->remove(newfaces[i]);
            for (size_t i = 0; i < newnodes.size(); i++)
                mesh->remove(newnodes[i]);
            newnodes.clear();
            newfaces.clear();
            return 1;
        }
    }

    return 0;

}

///////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_quad_loop(Mesh *mesh,
                       NodeSequence &anodes,
                       NodeSequence &bnodes,
                       NodeSequence &cnodes,
                       NodeSequence &dnodes,
                       NodeSequence &newnodes,
                       FaceSequence &newfaces,
                       bool smooth)
{
    newnodes.clear();
    newfaces.clear();

    assert(anodes.size() == cnodes.size());
    assert(bnodes.size() == dnodes.size());

    int nx = anodes.size();
    int ny = bnodes.size();

    for (int i = 0; i < nx; i++)
    {
        assert(mesh->contains(anodes[i]));
        assert(mesh->contains(cnodes[i]));
    }

    for (int i = 0; i < ny; i++)
    {
        assert(mesh->contains(bnodes[i]));
        assert(mesh->contains(dnodes[i]));
    }

    NodeSequence boundnodes(2 * nx + 2 * ny - 4);
    int index = 0;

    // Append anodes..
    for (int i = 0; i < anodes.size(); i++)
        boundnodes[index++] = anodes[i];

    // Append bnodes ...
    if (bnodes.front() != anodes.back())
        reverse(bnodes.begin(), bnodes.end());

    assert(anodes.back() == bnodes.front());
    for (int i = 1; i < bnodes.size(); i++)
        boundnodes[index++] = bnodes[i];

    // Append cnodes ...
    if (cnodes.front() != bnodes.back())
        reverse(cnodes.begin(), cnodes.end());

    assert(bnodes.back() == cnodes.front());
    for (int i = 1; i < cnodes.size(); i++)
        boundnodes[index++] = cnodes[i];

    // Append dnodes ...
    if (dnodes.front() != cnodes.back())
        reverse(dnodes.begin(), dnodes.end());

    assert(cnodes.back() == dnodes.front());
    for (int i = 1; i < dnodes.size(); i++)
        boundnodes[index++] = dnodes[i];

    // Ensure that loop is closed ...
    assert(anodes.front() == dnodes.back());

    return remesh_quad_loop(mesh, boundnodes, nx, ny, newnodes, newfaces, smooth);
}

////////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_tri_loop(Mesh *mesh,
                      NodeSequence &anodes,
                      NodeSequence &bnodes,
                      NodeSequence &cnodes,
                      int *partition,
                      NodeSequence &newnodes,
                      FaceSequence &newfaces,
                      bool smooth)
{
    // First thing to do is to clear the existing record.
    newnodes.clear();
    newfaces.clear();

    // We need atleast three nodes on each side ...
    if (anodes.size() < 3) return 1;
    if (bnodes.size() < 3) return 1;
    if (cnodes.size() < 3) return 1;

    int segments[3], partSegments[6];

    if (partition == NULL)
    {
        segments[0] = anodes.size() - 1;
        segments[1] = bnodes.size() - 1;
        segments[2] = cnodes.size() - 1;

        if (!Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments))
            return 1;
    }
    else
    {
        for (int i = 0; i < 6; i++)
            partSegments[i] = partition[i];
    }

    int err;
    if (anodes.back() == bnodes.back())
    {
        reverse(bnodes.begin(), bnodes.end());
        swap(partSegments[2], partSegments[3]);
    }

    if (anodes.front() == cnodes.front())
    {
        reverse(cnodes.begin(), cnodes.end());
        swap(partSegments[4], partSegments[5]);
    }

    // Ensure that input is closed loop of three sides ...
    assert(anodes.front() == cnodes.back());
    assert(bnodes.front() == anodes.back());
    assert(cnodes.front() == bnodes.back());

    // Ensure that segments are valid ...
    assert(anodes.size() == partSegments[0] + partSegments[1] + 1);
    assert(bnodes.size() == partSegments[2] + partSegments[3] + 1);
    assert(cnodes.size() == partSegments[4] + partSegments[5] + 1);

    // Split Side "A" nodes into a1nodes and b1nodes.
    NodeSequence a1nodes, b1nodes;
    int a1 = partSegments[0];
    int b1 = partSegments[1];
    err = split_stl_vector(anodes, a1 + 1, a1nodes, b1nodes);
    if (err) return 1;

    // Split Side "B" nodes into a2nodes and b2nodes.
    NodeSequence a2nodes, b2nodes;
    int a2 = partSegments[2];
    int b2 = partSegments[3];
    err = split_stl_vector(bnodes, a2 + 1, a2nodes, b2nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a3nodes, b3nodes;
    int a3 = partSegments[4];
    int b3 = partSegments[5];
    err = split_stl_vector(cnodes, a3 + 1, a3nodes, b3nodes);
    if (err) return 1;

    // Splitting nodes on each side
    Vertex *ca = a1nodes.back();
    Vertex *cb = a2nodes.back();
    Vertex *cc = a3nodes.back();
    Vertex *co = Face::centroid(ca, cb, cc);

    NodeSequence oa_nodes = linear_interpolation(co, ca, a2 + 1);
    NodeSequence ob_nodes = linear_interpolation(co, cb, a3 + 1);
    NodeSequence oc_nodes = linear_interpolation(co, cc, a1 + 1);

    mesh->addNode(co);
    newnodes.push_back(co);

    for (size_t i = 1; i < oa_nodes.size() - 1; i++)
    {
        mesh->addNode(oa_nodes[i]);
        newnodes.push_back(oa_nodes[i]);
    }

    for (size_t i = 1; i < ob_nodes.size() - 1; i++)
    {
        mesh->addNode(ob_nodes[i]);
        newnodes.push_back(ob_nodes[i]);
    }

    for (size_t i = 1; i < oc_nodes.size() - 1; i++)
    {
        mesh->addNode(oc_nodes[i]);
        newnodes.push_back(oc_nodes[i]);
    }

    NodeSequence qnodes;
    FaceSequence qfaces;

    err = remesh_quad_loop(mesh, b2nodes, a3nodes, oc_nodes, ob_nodes, qnodes, qfaces, 0);
    if (!err)
    {
        assert(!qfaces.empty());
        for (size_t i = 0; i < qnodes.size(); i++)
            newnodes.push_back(qnodes[i]);
        for (size_t i = 0; i < qfaces.size(); i++)
            newfaces.push_back(qfaces[i]);
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a1nodes, oa_nodes, oc_nodes, b3nodes, qnodes, qfaces, 0);
        if (!err)
        {
            assert(!qfaces.empty());
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a2nodes, ob_nodes, oa_nodes, b1nodes, qnodes, qfaces, 0);
        if (!err)
        {
            assert(!qfaces.empty());
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (err)
    {
        for (size_t i = 0; i < newfaces.size(); i++)
            mesh->remove(newfaces[i]);
        for (size_t i = 0; i < newnodes.size(); i++)
            mesh->remove(newnodes[i]);
        newnodes.clear();
        newfaces.clear();
        return 2;
    }

    if (smooth)
    {
        LaplaceLengthWeight lw;
        LaplaceSmoothing lapsmooth(mesh);
        lapsmooth.setWeight(&lw);
        lapsmooth.setNumIterations(10);
        lapsmooth.localized_at(newnodes);
    }

    // Check for any inversion of the element, if there is inversion,
    // undo everthing (i.e. remove new nodes and faces).
    //
    for (size_t i = 0; i < newfaces.size(); i++)
    {
        if (newfaces[i]->invertedAt() >= 0)
        {
            for (size_t i = 0; i < newfaces.size(); i++)
                mesh->remove(newfaces[i]);
            for (size_t i = 0; i < newnodes.size(); i++)
                mesh->remove(newnodes[i]);
            newnodes.clear();
            newfaces.clear();
            return 1;
        }
    }


    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
Jaal::remesh_penta_loop(Mesh *mesh,
                        NodeSequence &anodes,
                        NodeSequence &bnodes,
                        NodeSequence &cnodes,
                        NodeSequence &dnodes,
                        NodeSequence &enodes,
                        int *partition, NodeSequence &newnodes,
                        FaceSequence &newfaces, bool smooth)
{
    // First clear the existing records.
    newnodes.clear();
    newfaces.clear();

    int segments[5], partSegments[10];

    if (partition == NULL)
    {
        segments[0] = anodes.size() - 1;
        segments[1] = bnodes.size() - 1;
        segments[2] = cnodes.size() - 1;
        segments[3] = dnodes.size() - 1;
        segments[4] = enodes.size() - 1;
        if (!Face::is_5_sided_convex_loop_quad_meshable(segments, partSegments))
        {
            cout << "Warning: Triangle patch not quad meshable " << endl;
            return 1;
        }
    }
    else
    {
        for (int i = 0; i < 10; i++)
            partSegments[i] = partition[i];
    }

    int err;
    int index = 0;

    // Ensure that input is closed loop of three sides ...
    assert(anodes.front() == enodes.back());
    assert(bnodes.front() == anodes.back());
    assert(cnodes.front() == bnodes.back());
    assert(dnodes.front() == cnodes.back());
    assert(enodes.front() == dnodes.back());

    // Ensure that segments are valid ...
    assert(anodes.size() == partSegments[0] + partSegments[1] + 1);
    assert(bnodes.size() == partSegments[2] + partSegments[3] + 1);
    assert(cnodes.size() == partSegments[4] + partSegments[5] + 1);
    assert(dnodes.size() == partSegments[6] + partSegments[7] + 1);
    assert(enodes.size() == partSegments[8] + partSegments[9] + 1);

    // Split Side "A" nodes into a1nodes and b1nodes.
    NodeSequence a1nodes, b1nodes;
    int a1 = partSegments[0];
    int b1 = partSegments[1];
    assert(a1 > 0 && b1 > 0);
    err = split_stl_vector(anodes, a1 + 1, a1nodes, b1nodes);
    if (err) return 1;

    // Split Side "B" nodes into a2nodes and b2nodes.
    NodeSequence a2nodes, b2nodes;
    int a2 = partSegments[2];
    int b2 = partSegments[3];
    assert(a2 > 0 && b2 > 0);
    err = split_stl_vector(bnodes, a2 + 1, a2nodes, b2nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a3nodes, b3nodes;
    int a3 = partSegments[4];
    int b3 = partSegments[5];
    assert(a3 > 0 && b3 > 0);
    err = split_stl_vector(cnodes, a3 + 1, a3nodes, b3nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a4nodes, b4nodes;
    int a4 = partSegments[6];
    int b4 = partSegments[7];
    assert(a4 > 0 && b4 > 0);
    split_stl_vector(dnodes, a4 + 1, a4nodes, b4nodes);
    if (err) return 1;

    // Split Side "C" nodes into a3nodes and b3nodes.
    NodeSequence a5nodes, b5nodes;
    int a5 = partSegments[8];
    int b5 = partSegments[9];
    assert(a5 > 0 && b5 > 0);
    err = split_stl_vector(enodes, a5 + 1, a5nodes, b5nodes);
    if (err) return 1;

    assert(a1 == b4);
    assert(b1 == a3);

    assert(a2 == b5);
    assert(b2 == a4);

    assert(b3 == a5);

    // Splitting nodes on each side
    Vertex *ca = a1nodes.back();
    Vertex *cb = a2nodes.back();
    Vertex *cc = a3nodes.back();
    Vertex *cd = a4nodes.back();
    Vertex *ce = a5nodes.back();
    Vertex *co = Face::centroid(ca, cb, cc, cd, ce);

    NodeSequence oa_nodes = linear_interpolation(co, ca, a2 + 1);
    NodeSequence ob_nodes = linear_interpolation(co, cb, a3 + 1);
    NodeSequence oc_nodes = linear_interpolation(co, cc, a4 + 1);
    NodeSequence od_nodes = linear_interpolation(co, cd, a5 + 1);
    NodeSequence oe_nodes = linear_interpolation(co, ce, a1 + 1);

    mesh->addNode(co);
    newnodes.push_back(co);

    for (size_t i = 1; i < oa_nodes.size() - 1; i++)
    {
        mesh->addNode(oa_nodes[i]);
        newnodes.push_back(oa_nodes[i]);
    }

    for (size_t i = 1; i < ob_nodes.size() - 1; i++)
    {
        mesh->addNode(ob_nodes[i]);
        newnodes.push_back(ob_nodes[i]);
    }

    for (size_t i = 1; i < oc_nodes.size() - 1; i++)
    {
        mesh->addNode(oc_nodes[i]);
        newnodes.push_back(oc_nodes[i]);
    }

    for (size_t i = 1; i < od_nodes.size() - 1; i++)
    {
        mesh->addNode(od_nodes[i]);
        newnodes.push_back(od_nodes[i]);
    }

    for (size_t i = 1; i < oe_nodes.size() - 1; i++)
    {
        mesh->addNode(oe_nodes[i]);
        newnodes.push_back(oe_nodes[i]);
    }

    NodeSequence qnodes;
    FaceSequence qfaces;

    err = remesh_quad_loop(mesh, a1nodes, oa_nodes, oe_nodes, b5nodes, qnodes, qfaces, 0);
    if (!err)
    {
        for (size_t i = 0; i < qnodes.size(); i++)
            newnodes.push_back(qnodes[i]);
        for (size_t i = 0; i < qfaces.size(); i++)
            newfaces.push_back(qfaces[i]);
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a2nodes, ob_nodes, oa_nodes, b1nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a3nodes, oc_nodes, ob_nodes, b2nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a4nodes, od_nodes, oc_nodes, b3nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a5nodes, oe_nodes, od_nodes, b4nodes, qnodes, qfaces, 0);
        if (!err)
        {
            for (size_t i = 0; i < qnodes.size(); i++)
                newnodes.push_back(qnodes[i]);
            for (size_t i = 0; i < qfaces.size(); i++)
                newfaces.push_back(qfaces[i]);
        }
    }

    if (err) return 2;

    if (smooth)
    {
        LaplaceLengthWeight lw;
        LaplaceSmoothing lapsmooth(mesh);
        lapsmooth.setWeight(&lw);
        lapsmooth.setNumIterations(10);
        lapsmooth.localized_at(newnodes);
    }

    // Check for any inversion of the element, if there is inversion,
    // undo everthing (i.e. remove new nodes and faces).
    //
    for (size_t i = 0; i < newfaces.size(); i++)
    {
        if (newfaces[i]->invertedAt() >= 0)
        {
            for (size_t i = 0; i < newnodes.size(); i++)
                mesh->remove(newnodes[i]);
            for (size_t i = 0; i < newfaces.size(); i++)
                mesh->remove(newfaces[i]);
            newnodes.clear();
            newfaces.clear();
            return 1;
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

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
    FaceSequence adjFaces;
    NodeSequence vneighs;

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

///////////////////////////////////////////////////////////////////////////////

void
QTrack::advance(int endat)
{
    assert(sequence.size() == 2);

    // Starting node is always irregular ...
    FaceSequence vfaces = sequence[0]->getRelations2();
    assert(vfaces.size() != 4);

    while (1)
    {
        int progress = advance_single_step(endat);
        if (!progress) break;
    }

    Vertex *endvertex;

    // Sanity Checking ....
    if (endat == END_AT_TERMINALS)
    {
        endvertex = sequence.front();
        if (!endvertex->isBoundary())
        {
            vfaces = endvertex->getRelations2();
            assert(vfaces.size() != 4);
        }

        endvertex = sequence.back();
        if (!endvertex->isBoundary())
        {
            vfaces = endvertex->getRelations2();
            assert(vfaces.size() != 4);
        }

        for (int i = 1; i < sequence.size() - 1; i++)
        {
            vfaces = sequence[i]->getRelations2();
            assert(!sequence[i]->isBoundary());
            assert(vfaces.size() == 4);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

vector<QTrack> Jaal::generate_quad_irregular_graph(Mesh *mesh)
{
    vector<QTrack> qpath;
    int nTopo = mesh->isHomogeneous();
    if (nTopo != 4)
    {
        cout << "Error: The mesh must be all Quads " << endl;
        return qpath;
    }

    int relexist2 = mesh->build_relations(0, 2);
    int relexist0 = mesh->build_relations(0, 0);

    mesh->search_boundary();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
    }

    QTrack qp;
    qp.mesh = mesh;

    int found;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (!vertex->isBoundary())
        {
            NodeSequence vnodes = vertex->getRelations0();
            if (vnodes.size() != 4)
            {
                for (int j = 0; j < vnodes.size(); j++)
                {
                    qp.sequence.resize(2); // As we know the starting edge
                    qp.sequence[0] = vertex;
                    qp.sequence[1] = vnodes[j];
                    qp.advance(0); // Terminate at irregular nodes only..
                    found = 0;
                    for (int k = 0; k < qpath.size(); k++)
                    {
                        if (qpath[k] == qp)
                        {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        max_tag_val = max(max_tag_val, vertex->getLayerID());
    }

    for (size_t i = 0; i < numnodes; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        max_tag_val = max(max_tag_val, face->getLayerID());
    }

    for (size_t i = 0; i < numfaces; i++)
    {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setTag(1);
    }

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(0);
        if (face->isConvex())
            face->setTag(1);
        else
        {
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *v = face->getNodeAt(j);
                v->setTag(0);
            }
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

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
            vertex->setTag(2);
        else
            vertex->setTag(1);
    }

    size_t numfaces = mesh->getSize(2);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(1);
    }

    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (size_t i = 0; i < numnodes; i++)
    {
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
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        sumarea += face->getArea();
    }
    double avgarea = sumarea / (double) numfaces;
    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setTag(1);
    }

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if (face->getArea() < tolarea)
        {
            face->setTag(0);
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *vertex = face->getNodeAt(j);
                vertex->setTag(0);
            }
        }
        else
            face->setTag(1);
    }
}

///////////////////////////////////////////////////////////////////////////////

bool
Face::hasAllBoundNodes() const
{
    for (int i = 0; i < getSize(0); i++)
    {
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

    for (size_t i = 0; i < faces.size(); i++)
    {
        int nsize = faces[i]->getSize(0);
        for (int j = 0; j < nsize; j++)
        {
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

NodeSequence Mesh::get_breadth_first_ordered_nodes(Vertex *vstart, MeshFilter *filter)
{
    assert(vstart != NULL);

    NodeSequence seq;

    int relexist0 = build_relations(0, 0);

    size_t numnodes = getSize(0);
    if (numnodes == 0) return seq;

    for (size_t i = 0; i < numnodes; i++)
    {
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
    while (!vertexQ.empty())
    {
        Vertex *curr_vertex = vertexQ.front();
        vertexQ.pop_front();
        int currlevel = curr_vertex->getLayerID();
        if (filter)
        {
            if (curr_vertex != vstart) proceed = filter->pass(curr_vertex);
        }
        if (!curr_vertex->isVisited())
        {
            seq.push_back(curr_vertex);
            if (!proceed) break;
            curr_vertex->setVisitMark(1);
            vneighs = curr_vertex->getRelations0();
            for (size_t i = 0; i < vneighs.size(); i++)
            {
                if (!vneighs[i]->isVisited())
                {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        v->setVisitMark(0);
    }

    NodeSequence seq;

    if (vstart == 0) vstart = getNodeAt(0);
    vertexQ.push_back(vstart);
    NodeSequence vneighs;

    while (!vertexQ.empty())
    {
        Vertex *curr_vertex = vertexQ.front();
        vertexQ.pop_front();
        if (!curr_vertex->isVisited())
        {
            seq.push_back(curr_vertex);
            curr_vertex->setVisitMark(1);
            vneighs = curr_vertex->getRelations0();
            for (size_t i = 0; i < vneighs.size(); i++)
            {
                if (!vneighs[i]->isVisited())
                    vertexQ.push_front(vneighs[i]);
            }
        }
    }
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        assert(v->isVisited());
    }
    assert(seq.size() == numnodes);

    if (!relexist0)
        clear_relations(0, 0);

    return seq;
}
#endif

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_no_tags(Mesh *mesh)
{
    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setTag(0);
    }

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(0);
    }
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_visit_tags(Mesh *mesh)
{
    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isVisited())
            vertex->setTag(0);
        else
            vertex->setTag(1);
    }

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if (face->isVisited())
            face->setTag(2);
        else
            face->setTag(3);
    }
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_allboundnodes_tag(Mesh *mesh)
{
    mesh->search_boundary();
    size_t numfaces = mesh->getSize(2);

    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (size_t i = 0; i < numfaces; i++)
    {
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

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setTag(1);
        if (!vertex->isBoundary())
        {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        v->setTag(1);
    }

    size_t ncount = 0;
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(1);
        int pos = face->invertedAt();
        if (pos >= 0)
        {
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
    for (size_t i = 0; i < numnodes; i++)
    {
        v = mesh->getNodeAt(i);
        v->setTag(0);
    }

    for (size_t i = 0; i < qpath.size(); i++)
    {
        v = qpath[i].sequence.front();
        v->setTag(1);
        v = qpath[i].sequence.back();
        v->setTag(1);

        for (size_t j = 1; j < qpath[i].sequence.size() - 1; j++)
        {
            v = qpath[i].sequence[j];
            v->setTag(2);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int Jaal::SurfPatch::getPosOf(const Vertex *v)
{
    int pos = -1;
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
    set<Face*>::const_iterator fiter;
    std::map<Vertex*, FaceSet> relations02;

    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        Face *face = *fiter;
        for (int j = 0; j < face->getSize(0); j++)
        {
            vertex = face->getNodeAt(j);
            relations02[vertex].insert(face);
        }
    }

    // A boundary edge must have exactly one face neighbor...
    FaceSequence faceneighs;
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        Face *face = *fiter;
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++)
        {
            Vertex *v0 = face->getNodeAt((j + 0) % nnodes);
            Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
            faceneighs.clear();
            assert(relations02[v0].size() > 0);
            assert(relations02[v1].size() > 0);
            set_intersection(relations02[v0].begin(), relations02[v0].end(),
                             relations02[v1].begin(), relations02[v1].end(),
                             back_inserter(faceneighs));
            if (faceneighs.size() == 1)
            {
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
    for (int k = 0; k < boundary.size(); k++)
    {
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
    bound_nodes.resize(boundary.size());
    for (int k = 0; k < boundary.size(); k++)
    {
        bound_nodes[k] = boundary[k].getNodeAt(0); // Only the first node.
    }

    //
    // Collect the inner nodes of the patch. These nodes will be deleted, if
    // the remesh operation is successful...
    //
    inner_nodes.clear();
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        Face *face = *fiter;
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++)
        {
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

    set<Vertex*>::const_iterator it;
    int index = 0;
    for (it = corners.begin(); it != corners.end(); ++it)
    {
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

    if (segSize[0] == segSize[2])
    {
        if (min(segSize[1], segSize[3]) == 2) return 1;
        //  Either Segment 2 or 3 must be starting node
        if (segSize[1] < segSize[3])
            start_corner = bound_nodes[ cornerPos[1] ];
        else
            start_corner = bound_nodes[ cornerPos[3] ];
        start_boundary_loop_from(start_corner);
    }

    if (min(segSize[0], segSize[2]) == 2) return 1;

    if (segSize[2] < segSize[0])
    {
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

    if (segSize[3] < segSize[1])
    {
        start_corner = bound_nodes[ cornerPos[3] ];
        start_boundary_loop_from(start_corner);
    }
    else
    {
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
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
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

    if (space_dim == 2)
    {
        connect.resize(4);
        for (int j = 0; j < ny - 1; j++)
        {
            for (int i = 0; i < nx - 1; i++)
            {
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

    assert(face->invertedAt() == -1);

    xyz[0] = -0.0;
    xyz[1] = 0.000001;
    xyz[2] = 0.0;
    nodes[1]->setXYZCoords(xyz);
    assert(face->invertedAt() == 1);

    return 0;
}

