#include "Mesh.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <set>

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

size_t Vertex::global_id = 0;

#ifdef USE_BOOST_LIBS
AttribKey MeshEntity::maxAttribID = 0;
std::map<string, AttribKey> MeshEntity::attribKeyMap;

///////////////////////////////////////////////////////////////////////////////

AttribKey MeshEntity::getAttribKey(const string &s)
{
  map<string,AttribKey>::const_iterator it;

  it = attribKeyMap.find(s);

  if( it == attribKeyMap.end() ) return 0;

  return it->second;
}
///////////////////////////////////////////////////////////////////////////////

bool MeshEntity::hasAttribute(const string &s)
{
  return getAttribKey(s);
}

///////////////////////////////////////////////////////////////////////////////

AttribKey MeshEntity::addAttribute(const string &s)
{
  AttribKey key = getAttribKey(s);
  if( key == 0) attribKeyMap[s] = ++maxAttribID;
  return attribKeyMap[s];
}

void MeshEntity::removeAttribute(const string &s)
{
  attribKeyMap.erase(s);
}
#endif

///////////////////////////////////////////////////////////////////////////////

NodeType Vertex::newObject()
{
  Vertex *v = new Vertex;
  assert(v);
  v->setID(global_id);
  global_id++;
  return v;
}

///////////////////////////////////////////////////////////////////////////////

NodeType Vertex::getClone() const
{
  Vertex *v = new Vertex;
  assert(v);
  v->setID(global_id);
  v->setXYZCoords(xyz);
  return v;
}

///////////////////////////////////////////////////////////////////////////////

EdgeType Edge::newObject()
{
  Edge *e = new Edge;
  assert(e);
  return e;
}

FaceType Face::newObject()
{
  Face *f = new Face;
  assert(f);
  return f;
}

///////////////////////////////////////////////////////////////////////////////

Point3D Vertex::mid_point(const Vertex *v0, const Vertex *v1, double alpha)
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

double Vertex::length(const Vertex *v0, const Vertex *v1)
{

  Point3D p0 = v0->getXYZCoords();
  Point3D p1 = v1->getXYZCoords();

  double dx = p0[0] - p1[0];
  double dy = p0[1] - p1[1];
  double dz = p0[2] - p1[2];

  return sqrt(dx * dx + dy * dy + dz * dz);
}

///////////////////////////////////////////////////////////////////////////////

double Vertex::length2(const Vertex *v0, const Vertex *v1)
{

  Point3D p0 = v0->getXYZCoords();
  Point3D p1 = v1->getXYZCoords();

  double dx = p0[0] - p1[0];
  double dy = p0[1] - p1[1];
  double dz = p0[2] - p1[2];

  return dx * dx + dy * dy + dz * dz;
}

///////////////////////////////////////////////////////////////////////////////

Point3D Face::getCentroid() const
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

  pc[0] /= 3.0;
  pc[1] /= 3.0;
  pc[2] /= 3.0;

  return pc;
}

////////////////////////////////////////////////////////////////////////////////

Vec3D Face :: normal( const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
  Point3D p2p0   = Math::create_vector( p2, p0);
  Point3D p1p0   = Math::create_vector( p1, p0);
  Point3D normal = Math::cross_product( p2p0, p1p0 ); 

  double mag = Math::magnitude( normal );
  normal[0] /= mag;
  normal[1] /= mag;
  normal[2] /= mag;

  return normal;
}

////////////////////////////////////////////////////////////////////////////////

Vec3D Face :: normal( const Vertex *v0, const Vertex *v1, const Vertex *v2)
{
  return Math::normal( v0->getXYZCoords(),
                       v1->getXYZCoords(),
		       v2->getXYZCoords() );
}

////////////////////////////////////////////////////////////////////////////////

double Face :: tri_area( const Point3D &p0, const Point3D &p1, const Point3D &p2)
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
   d[0] = Math::length2(p1,p2);
   d[1] = Math::length2(p2,p0);
   d[2] = Math::length2(p0,p1);

   std::sort( d, d + 3 ); // May be we should have optimized version than STL one 

   double a2 = d[0];
   double b2 = d[1];
   double c2 = d[2];

   double area = 0.25*sqrt(4*a2*b2 - (a2+b2-c2)*(a2+b2-c2) );
   return area;
}

///////////////////////////////////////////////////////////////////////////////////

double Face :: quad_area( const Point3D &p0, const Point3D &p1, 
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

  Point3D v2v0   = Math::create_vector( p2, p0);
  Point3D v3v1   = Math::create_vector( p3, p1);
  Point3D d0d1   = Math::cross_product( v2v0, v3v1 ); 

  double area = 0.5*Math::magnitude( d0d1 );
  return area;
}
////////////////////////////////////////////////////////////////////////////////
bool Face :: is_convex_quad(const Point3D &p0, const Point3D &p1, 
                            const Point3D &p2, const Point3D &p3)
{
  double qarea  = quad_area(p0, p1, p2, p3);

  double tarea1, tarea2;
  tarea1 = tri_area( p0, p1, p2);
  tarea2 = tri_area( p0, p2, p3);
  if( fabs(tarea1 + tarea2- qarea) > 1.0E-10 ) return 0;

  tarea1 = tri_area( p0, p1, p3);
  tarea2 = tri_area( p1, p2, p3);
  if( fabs(tarea1 + tarea2 - qarea) > 1.0E-10 ) return 0;

  return 1;
}

////////////////////////////////////////////////////////////////////////////////
vector<Vertex*> Face :: getRelations0()
{
   vector<Vertex*> vneighs;
   set<Vertex*>    vset;

   for( int i = 0; i < connect.size(); i++) {
        Vertex *vertex = connect[i];
	vneighs = vertex->getRelations0();
	for( size_t j = 0; j < vneighs.size(); j++) 
	     vset.insert( vneighs[j] );
   }

   for( int i = 0; i < connect.size(); i++) 
        vset.erase( connect[i] );

   vector<Vertex*> vresult;
   if( !vset.empty() ) {
        set<Vertex*>::const_iterator it;
	size_t index = 0;
	for( it = vset.begin(); it != vset.end(); ++it)
	     vresult[index++] = *it;
   }
   return vresult;
}
////////////////////////////////////////////////////////////////////////////////

int Mesh :: check_convexity()
{
   size_t numfaces = getSize(2);

   int itag;
   for( size_t i = 0; i < numfaces; i++) {
        Face *face = getFaceAt(i);
	itag = 0;
	if( face->getSize(0) == 4 ) 
	    itag = Face::is_convex_quad( face->getNodeAt(0)->getXYZCoords(),
	                           face->getNodeAt(1)->getXYZCoords(),
	                           face->getNodeAt(2)->getXYZCoords(),
	                           face->getNodeAt(3)->getXYZCoords() );
        if( itag == 0) 
	    face->setTag(1);
        else
	    face->setTag(2);

   }
   return 0;
}
////////////////////////////////////////////////////////////////////////////////

double Face :: getAspectRatio() 
{
    int nSize = connect.size();

    double minlen = MAXDOUBLE;
    double maxlen = 0.0;

    for( int i = 0; i < nSize; i++) {
         Vertex *v0  = connect[(i+0)%nSize];
         Vertex *v1  = connect[(i+1)%nSize];
	 double len2 = Vertex::length2( v0, v1 );
	 if( len2  > maxlen ) maxlen = len2;
	 if( len2  < minlen ) minlen = len2;
    }

    return sqrt( minlen/maxlen );
}

///////////////////////////////////////////////////////////////////////////////
int Mesh :: make_chain( vector<Edge> &boundedges )
{
   size_t nSize = boundedges.size();

   Edge edge = boundedges.front();

   list<Edge>  listedges;
   for( int i = 1; i < boundedges.size(); i++) 
        listedges.push_back( boundedges[i] );
   boundedges.clear();

   boundedges.reserve( nSize );
         
   Vertex *first_vertex = edge.getNodeAt(0);
   Vertex *curr_vertex  = edge.getNodeAt(1);

   boundedges.push_back(edge);

   list<Edge>::iterator it;

   for( int i = 0; i < nSize; i++) {
       for( it = listedges.begin(); it != listedges.end(); ++it) 
       {
            edge = *it;
	    Vertex *v0 = edge.getNodeAt(0);
	    Vertex *v1 = edge.getNodeAt(1);
	    if( v0 == curr_vertex ) {
	        curr_vertex =  v1;
	        boundedges.push_back( edge );
	        break;
            }
	    if( v1 == curr_vertex ) {
	        curr_vertex =  v0;
		Edge newedge(v1,v0);
	        boundedges.push_back( newedge );
	        break;
            }
       }
       if( it != listedges.end() ) listedges.erase(it);
   }

   return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh :: is_closeable_chain( const vector<Edge> &boundedges ) 
{
   std::map<Vertex*, set<Vertex*> > relations00;

   for( size_t i  = 0; i < boundedges.size(); i++) {
        Vertex *v0 = boundedges[i].getNodeAt(0);
        Vertex *v1 = boundedges[i].getNodeAt(1);
	relations00[v0].insert(v1);
	relations00[v1].insert(v0);
   }

   std::map<Vertex*, set<Vertex*> > :: const_iterator it;
   for( it = relations00.begin(); it != relations00.end(); ++it){
        Vertex *v = it->first;
	if( relations00[v].size() != 2 ) return 0;
   }

   return 1;
}

///////////////////////////////////////////////////////////////////////////////

int  Mesh :: is_closed_chain( const vector<Edge> &boundedges ) 
{
   Vertex *first_vertex = boundedges[0].getNodeAt(0);
   Vertex *curr_vertex  = boundedges[0].getNodeAt(1);

   for ( size_t i = 1; i < boundedges.size(); i++) {
         if( boundedges[i].getNodeAt(0) != curr_vertex) return 0;
	 curr_vertex = boundedges[i].getNodeAt(1);
   }
   if( curr_vertex != first_vertex ) return 0;

   return 1;
}
///////////////////////////////////////////////////////////////////////////////
int Mesh :: rotate_chain( vector<Edge> &boundedges, Vertex *first_vertex)
{
   size_t nSize = boundedges.size();
   if( !is_closeable_chain(boundedges) ) return 1;

   vector<Edge>  listedges(nSize);
   int istart = 0;
   for( size_t i = 0; i < nSize; i++) {
        if( boundedges[i].getNodeAt(0) == first_vertex ) istart = i;
	listedges[i] = boundedges[i];
   }

   for( size_t i = 0; i < nSize; i++)
	boundedges[i] = listedges[(i+ istart)%nSize];

   return 0;
}
///////////////////////////////////////////////////////////////////////////////


vector<FaceType> Mesh::getRelations102(NodeType vtx0, NodeType vtx1)
{
  vector<FaceType> v0faces = vtx0->getRelations2();
  vector<FaceType> v1faces = vtx1->getRelations2();

  if (v0faces.empty() || v1faces.empty())
  {
    cout << "Warning: Vertex-Faces relations are empty " << endl;
  }

  set<FaceType> vset;
  for( int i = 0; i < v0faces.size(); i++) 
       vset.insert( v0faces[i] );

  for( int i = 0; i < v1faces.size(); i++) 
       vset.insert( v1faces[i] );

  std::set<FaceType>::iterator it;

  vector<FaceType> faceneighs;
  if( vset.size() ) {
     faceneighs.resize( vset.size() );
     int index = 0;
     for( it = vset.begin(); it != vset.end(); ++it)
          faceneighs[index++] = *it;
  }

  return faceneighs;
}

///////////////////////////////////////////////////////////////////////////////

vector<FaceType> Mesh::getRelations112(NodeType vtx0, NodeType vtx1)
{
  vector<FaceType> v0faces = vtx0->getRelations2();
  vector<FaceType> v1faces = vtx1->getRelations2();

  if (v0faces.empty() || v1faces.empty())
  {
    cout << "Warning: Vertex-Faces relations are empty " << endl;
  }

  for( size_t i = 0; i < v0faces.size()-1; i++) {
       assert( v0faces[i] < v0faces[i+1] );
       assert( !v0faces[i]->isRemoved() );
  }

  for( size_t i = 0; i < v1faces.size()-1; i++) {
       assert( v1faces[i] < v1faces[i+1] );
       assert( !v1faces[i]->isRemoved() );
  }

  vector<FaceType> faceneighs;
  set_intersection(v0faces.begin(), v0faces.end(), v1faces.begin(),
      v1faces.end(), back_inserter(faceneighs));

  return faceneighs;
}

///////////////////////////////////////////////////////////////////////////////

vector<FaceType> Face::getRelations202()
{
  vector<FaceType> faceneighs, vneighs;

  int nSize = connect.size();
  for (int i = 0; i < nSize; i++)
  {
    Vertex *v0 = connect[(i + 0) % nSize];
    vneighs = v0->getRelations2();
    for (int j = 0; j < vneighs.size(); j++)
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

vector<FaceType> Face::getRelations212()
{
  vector<FaceType> faceneighs, edgeneighs;

  int nSize = connect.size();
  for (int i = 0; i < nSize; i++)
  {
    Vertex *v0 = connect[(i + 0) % nSize];
    Vertex *v1 = connect[(i + 1) % nSize];
    edgeneighs = Mesh::getRelations112(v0, v1);
    for (int j = 0; j < edgeneighs.size(); j++)
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

size_t Mesh::count_edges()
{
  int relexist = build_relations(0, 0);

  size_t numnodes = getSize(0);

  vector<Vertex*> neighs;
  size_t ncount = 0;
  for (size_t i = 0; i < numnodes; i++)
  {
    neighs = nodes[i]->getRelations0();
    for (size_t j = 0; j < neighs.size(); j++)
      if (nodes[i] > neighs[j])
        ncount++;
  }

  if (!relexist)
    clear_relations(0, 0);

  return ncount;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::prune()
{
  NodeContainer livenodes;

  livenodes.reserve(nodes.size());
  for (size_t i = 0; i < nodes.size(); i++)
  {
    if (!nodes[i]->isRemoved())
      livenodes.push_back(nodes[i]);
    else
      delete nodes[i];
  }
  nodes = livenodes;

  FaceContainer livefaces;
  livefaces.reserve(faces.size());
  for (size_t i = 0; i < faces.size(); i++)
  {
    if (!faces[i]->isRemoved())
      livefaces.push_back(faces[i]);
    else
      delete faces[i];
  }
  faces = livefaces;
}

///////////////////////////////////////////////////////////////////////////////

bool Mesh::isPruned() const 
{
   for( size_t i = 0; i < nodes.size(); i++) 
        if( nodes[i]->isRemoved() ) return 0;

   for( size_t i = 0; i < faces.size(); i++) 
        if( faces[i]->isRemoved() ) return 0;

   return 1;
}
///////////////////////////////////////////////////////////////////////////////

void Mesh::enumerate(int etype)
{
  size_t index = 0;

  NodeContainer::const_iterator viter;
  if (etype == 0)
  {
    index = 0;
    for (viter = nodes.begin(); viter != nodes.end(); ++viter)
    {
      Vertex *vertex = *viter;
      vertex->setID(index++);
    }
  }

  FaceContainer::const_iterator fiter;
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

size_t Mesh::getBoundarySize(int d) const
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
int Mesh::isHomogeneous() const
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

NodeType Face::opposite_node(const FaceType tri, NodeType n1, NodeType n2)
{
  NodeType tn0 = tri->getNodeAt(0);
  NodeType tn1 = tri->getNodeAt(1);
  NodeType tn2 = tri->getNodeAt(2);

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

void Face::opposite_nodes(const FaceType quad, NodeType n1, NodeType n2,
    NodeType &n3, NodeType &n4)
{
  NodeType qn0 = quad->getNodeAt(0);
  NodeType qn1 = quad->getNodeAt(1);
  NodeType qn2 = quad->getNodeAt(2);
  NodeType qn3 = quad->getNodeAt(3);

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

FaceType Face::create_quad(const FaceType t1, const FaceType t2)
{
  vector<NodeType> connect;
  NodeType commonnodes[3];

  connect = t1->getNodes();

  int index = 0;
  for (int i = 0; i < 3; i++)
  {
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
  qface->setNodes(connect);
  return qface;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::saveAs(const string &s)
{
  string filename = s + ".dat";
  ofstream ofile(filename.c_str(), ios::out);

  size_t numnodes = nodes.size();
  size_t numfaces = faces.size();

  size_t nn = numnodes;

  ofile << nn << " " << numfaces << endl;

  for (size_t i = 0; i < numnodes; i++)
  {
    Point3D p3d = nodes[i]->getXYZCoords();
    ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
  }

  vector<NodeType> connect;
  for (size_t i = 0; i < numfaces; i++)
  {
    Face *face = faces[i];
    connect = face->getNodes();
    int nnodes = connect.size();
    ofile << nnodes << " ";
    for (int j = 0; j < nnodes; j++)
    {
      int vid = connect[j]->getID();
      if (vid >= numnodes)
      {
        cout << "Vertex indexing out of range " << vid << endl;
        exit(0);
      }
      ofile << vid << " ";
    }
    ofile << endl;
  }

  for (size_t i = 0; i < numnodes; i++)
  {
    Vertex *vertex = nodes[i];
    ofile << vertex->getTag() << " ";
  }
  ofile << endl;

  for (size_t i = 0; i < numfaces; i++)
  {
    Face *face = faces[i];
    ofile << face->getTag() << " ";
  }
}

///////////////////////////////////////////////////////////////////////////////
int Mesh::nearest_neighbour(const Vertex *myself, Vertex* &nearest,
    double &mindist)
{
  assert(getAdjTable(0, 0));

  mindist = MAXDOUBLE;
  nearest = NULL;

  vector<Vertex*> neighs = myself->getRelations0();

  for (int i = 0; i < neighs.size(); i++)
  {
    double d = Vertex::length2(myself, neighs[i]);
    if (d < mindist)
    {
      mindist = d;
      nearest = neighs[i];
    }
  }

  mindist = sqrt(mindist);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::build_relations02()
{
  if (adjTable[0][2] == 1)
    return 1;

  clear_relations(0, 2);

  size_t numfaces = getSize(2);

  for (size_t iface = 0; iface < numfaces; iface++)
  {
    Face *face = getFaceAt(iface);
    assert(face);
    for (size_t j = 0; j < face->getSize(0); j++)
    {
      Vertex *vtx = face->getNodeAt(j);
      vtx->addRelation2(face);
    }
  }
  adjTable[0][2] = 1;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::build_relations00()
{
  if (adjTable[0][0] == 1)
    return 1;

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

void Mesh::clear_relations(int src, int dst)
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

int Mesh::search_boundary()
{
  if( boundary_status == 1 ) return 1;

  if( !isPruned() ) prune();

  int relexist = build_relations(0, 2);

  size_t numnodes = getSize(0);
  for (size_t inode = 0; inode < numnodes; inode++) 
  {
       Vertex *v = getNodeAt(inode);
       v->setTag(0);
  }

  size_t numfaces = getSize(2);
  vector<FaceType> neighs;
  for (size_t iface = 0; iface < numfaces; iface++)
  {
    Face *face = getFaceAt(iface);
    assert(face);
    face->setTag(0);
    size_t nnodes = face->getSize(0);
    for (size_t j = 0; j < nnodes; j++)
    {
      Vertex *v0 = face->getNodeAt(j);
      Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
      neighs = Mesh::getRelations112(v0, v1);
      if (neighs.size() == 1)
      {
        v0->setBoundaryMark(1);
        v1->setBoundaryMark(1);
        int bmark = max(1, face->getBoundaryMark());
        face->setBoundaryMark(bmark);
        face->setTag(1);
        v0->setTag(1);
        v1->setTag(1);
      }
    }
  }

  if (!relexist)
    clear_relations(0, 2);

  boundary_status = 1;

  return 0;

}
///////////////////////////////////////////////////////////////////////////////
Mesh::FaceContainer Mesh::filter(int facetype) const
{
  FaceContainer::const_iterator it;
  size_t ncount = 0;
  for (it = faces.begin(); it != faces.end(); ++it)
  {
    Face *face = *it;
    if (face->getType() == facetype)
      ncount++;
  }

  FaceContainer tmpfaces;
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
bool Mesh::isSimple()
{
  int simple = 1;
  int relexist = build_relations(0, 2);

  size_t numfaces = getSize(2);

  vector<FaceType> neighs;
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
bool Mesh::isConsistentlyOriented()
{
  int consistent = 1;
  int relexist = build_relations(0, 2);

  size_t numfaces = getSize(2);

  vector<FaceType> neighs;
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
        int dir1 = neighs[0]->getOrientation(v0, v1);
        int dir2 = neighs[1]->getOrientation(v0, v1);
        if (dir1 * dir2 == 1)
        {
          cout << "Warning: Mesh is not consistently oriented " << endl;
          cout << "Face 1: ";
          for (int k = 0; k < neighs[0]->getSize(0); k++)
            cout << neighs[0]->getNodeAt(k)->getID() << " ";
          cout << endl;
          cout << "Face 2: ";
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

void Mesh::makeConsistentlyOriented()
{
  build_relations(0, 2);

  Face *face = NULL;
  deque<Face*> faceQ;
  vector<FaceType> neighs;

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

int Mesh::getNumOfComponents()
{
  build_relations(0, 2);

  Face *face = NULL;
  deque<Face*> faceQ;
  vector<FaceType> neighs;

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
        faceQ.push_back(face);
        break;
      }
    }

    if (faceQ.empty())
      break;

    numComponents++;
    while (!faceQ.empty())
    {
      Face *face = faceQ.front();
      faceQ.pop_front();
      if (!face->isVisited())
      {
        face->setVisitMark(1);
        int nnodes = face->getSize(0);
        assert(nnodes == 4);
        for (int j = 0; j < nnodes; j++)
        {
          Vertex *v0 = face->getNodeAt(j);
          Vertex *v1 = face->getNodeAt((j + 1) % nnodes);
          neighs = Mesh::getRelations112(v0, v1);
          if (neighs.size() == 2)
          {
            faceQ.push_back(neighs[0]);
            faceQ.push_back(neighs[1]);
          }
        }
      }
    } // Complete one Component
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
Mesh *struct_tri_grid(int nx, int ny)
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

  vector<Vertex*> connect(3);
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
      newtri = new Face;
      newtri->setNodes(connect);
      trimesh->addFace(newtri);

      connect[0] = trimesh->getNodeAt(n0);
      connect[1] = trimesh->getNodeAt(n2);
      connect[2] = trimesh->getNodeAt(n3);
      newtri = new Face;
      newtri->setNodes(connect);
      trimesh->addFace(newtri);
    }
  }
  return trimesh;
}

/////////////////////////////////////////////////////////////////////////////

Mesh * Jaal::struct_quad_grid(int nx, int ny)
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

  vector<Vertex*> connect(4);
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
      newquad = new Face;
      newquad->setNodes(connect);
      quadmesh->addFace(newquad);
    }
  }
  return quadmesh;
}

////////////////////////////////////////////////////////////////////

void expand_strip(Face *prevface, Vertex *v0, Vertex *v1, list<Face*> &strip)
{
  vector<Face*> neighs = Mesh::getRelations112(v0, v1);

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

void Mesh::get_quad_strips(Face *rootface, vector<Face*> &strip1,
    vector<Face*> &strip2)
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

vector<Face*> Mesh::get_bound_faces(int bound_what)
{
  int relexist = build_relations(0, 2);

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
      if (face->hasBoundaryEdge())
        bfaces.insert(face);
    }
  }

  vector<Face*> result;

  size_t nSize = bfaces.size();

  if (nSize)
  {
    result.resize(nSize);
    set<Face*>::const_iterator it;

    size_t index = 0;
    for (it = bfaces.begin(); it != bfaces.end(); ++it)
      result[index++] = *it;
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::set_strip_markers()
{
  vector<Face*> bound_faces = get_bound_faces(1);

  size_t numfaces = getSize(2);
  for (size_t i = 0; i < numfaces; i++)
  {
    Face *face = getFaceAt(i);
    face->setTag(0);
  }

  vector<Face*> strip1, strip2;
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

iBase_EntityHandle Mesh::get_MOAB_Handle(iMesh_Instance imesh, Vertex *vertex)
{
  int err;

  iBase_EntityHandle newHandle = vertex->get_MOAB_Handle();
  if (newHandle)
    return newHandle;

  Point3D p = vertex->getXYZCoords();
  iMesh_createVtx(imesh, p[0], p[1], p[2], &newHandle, &err);
  assert(!err);
  vertex->set_MOAB_Handle(newHandle);

  return newHandle;
}

///////////////////////////////////////////////////////////////////////////////
iBase_EntityHandle Mesh::get_MOAB_Handle(iMesh_Instance imesh, Face *face)
{
  int status, err;

  // Not a good way. Shouldn't delete it if the connectivity is unchanged.

  iBase_EntityHandle newHandle = face->get_MOAB_Handle();
  if (newHandle) 
      iMesh_deleteEnt(imesh, newHandle, &err);

  vector<iBase_EntityHandle> connect;

  int nnodes = face->getSize(0);
  connect.resize(nnodes);

  for (int j = 0; j < nnodes; j++)
  {
    Vertex *v = face->getNodeAt(j);
    connect[j] = get_MOAB_Handle(imesh, v);
  }

  switch (nnodes) {
  case 3:
    iMesh_createEnt(imesh, iMesh_TRIANGLE, &connect[0], nnodes, &newHandle,
        &status, &err);
    assert(!err);
    break;
  case 4:
    iMesh_createEnt(imesh, iMesh_QUADRILATERAL, &connect[0], nnodes,
        &newHandle, &status, &err);
    assert(!err);
    break;
  default:
    iMesh_createEnt(imesh, iMesh_POLYGON, &connect[0], nnodes, &newHandle,
        &status, &err);
    assert(!err);
    break;
  }

  face->set_MOAB_Handle(newHandle);

  return newHandle;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::toMOAB(iMesh_Instance &imesh, iBase_EntitySetHandle entitySet)
{
  int  err;

  if (imesh == 0)
  {
    int optlen = 0;
    char *options = NULL;
    iMesh_newMesh(options, &imesh, &err, optlen);
    assert(!err);
  }

  search_boundary();

  const char *tagname = "fixed";
  int namelen = strlen(tagname);

  iBase_TagHandle idtag;

  iMesh_getTagHandle(imesh, tagname, &idtag, &err, namelen);

  if( err ) 
      iMesh_createTag(imesh, tagname, 1, iBase_INTEGER, &idtag, &err, namelen);

  assert( !err );

  iBase_EntityHandle newHandle;

  size_t numnodes = getSize(0);
  for (size_t i = 0; i < numnodes; i++)
  {
    Vertex *v = getNodeAt(i);
    newHandle = get_MOAB_Handle(imesh, v);
    assert( newHandle );
    if (entitySet)
      iMesh_addEntToSet(imesh, newHandle, entitySet, &err);
    int bmark = v->getBoundaryMark();
    iMesh_setIntData(imesh, newHandle, idtag, bmark, &err);
    assert( !err );
  }

  size_t numfaces = getSize(2);
  for (size_t i = 0; i < numfaces; i++)
  {
    Face *f = getFaceAt(i);
    newHandle = get_MOAB_Handle(imesh, f);
    if (entitySet)
      iMesh_addEntToSet(imesh, newHandle, entitySet, &err);
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

int Mesh::fromMOAB(iMesh_Instance imesh, iBase_EntitySetHandle entitySet)
{
  int err;
  int numNodes, numFaces;

  if (entitySet == 0)
    iMesh_getRootSet(imesh, &entitySet, &err);

  iMesh_getNumOfType(imesh, entitySet, iBase_VERTEX, &numNodes, &err);
  assert(!err);

  if (numNodes == 0)
  {
    cout << "Warning: There are no nodes in iMesh " << endl;
    return 1;
  }

  reserve(numNodes, 0);

  SimpleArray<iBase_EntityHandle> nodeHandles;
  iMesh_getEntities(imesh, entitySet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
      ARRAY_INOUT(nodeHandles), &err);

  typedef map<iBase_EntityHandle, Vertex*> Moab2JaalNodeMap;
  Moab2JaalNodeMap moab2jaalNodes;

  Point3D p3d;
  double x, y, z;
  for (int i = 0; i < numNodes; i++)
  {
    Vertex *vtx = Vertex::newObject();
    moab2jaalNodes[nodeHandles[i]] = vtx;
    iMesh_getVtxCoord(imesh, nodeHandles[i], &x, &y, &z, &err);
    p3d[0] = x;
    p3d[1] = y;
    p3d[2] = z;
    vtx->setXYZCoords(p3d);
    vtx->setID(i);
    vtx->set_MOAB_Handle(nodeHandles[i]);
    addNode(vtx);
  }

  iMesh_getNumOfType(imesh, entitySet, iBase_FACE, &numFaces, &err);

  if (numNodes == 0)
  {
    cout << "Warning: There are no faces in iMesh " << endl;
    return 1;
  }

  reserve(numFaces, 2);

  vector<Vertex*> connect(3);
  SimpleArray<iBase_EntityHandle> tfaceHandles, qfaceHandles, facenodes;

  iMesh_getEntities(imesh, entitySet, iBase_FACE, iMesh_TRIANGLE, ARRAY_INOUT(
      tfaceHandles), &err);

  size_t numTris = tfaceHandles.size();
  if (numTris)
  {
    connect.resize(3);
    for (size_t i = 0; i < numTris; i++)
    {
      iMesh_getEntAdj(imesh, tfaceHandles[i], iBase_VERTEX, ARRAY_INOUT(
          facenodes), &err);
      for (int j = 0; j < 3; j++)
        connect[j] = moab2jaalNodes[facenodes[j]];
      Face *face = new Face;
      face->setNodes(connect);
      face->set_MOAB_Handle(tfaceHandles[i]);
      addFace(face);
    }
    facenodes.clear();
  }

  iMesh_getEntities(imesh, entitySet, iBase_FACE, iMesh_QUADRILATERAL,
      ARRAY_INOUT(qfaceHandles), &err);

  size_t numQuads = qfaceHandles.size();
  if (numQuads)
  {
    connect.resize(4);
    for (int i = 0; i < numQuads; i++)
    {
      iMesh_getEntAdj(imesh, qfaceHandles[i], iBase_VERTEX, ARRAY_INOUT(
          facenodes), &err);
      for (int j = 0; j < 4; j++)
        connect[j] = moab2jaalNodes[facenodes[j]];
      Face *face = new Face;
      face->setNodes(connect);
      face->set_MOAB_Handle(qfaceHandles[i]);
      addFace(face);
    }
  }

}

/////////////////////////////////////////////////////////////////////////////////////

vector<int> Mesh::get_topological_statistics( int entity, bool sorted )
{
  int relexist = build_relations(0, 2);

  assert(getAdjTable(0, 2));

  int numnodes = getSize(0);

  vector<int> degree(numnodes);
  vector<Face*> neighs;
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

int Mesh::setNodeWavefront()
{
  int relexist = build_relations(0, 0);

  size_t numNodes = getSize(0);

  for (size_t i = 0; i < numNodes; i++)
  {
    Vertex *v = getNodeAt(i);
    v->setLayerID(-1);
    v->setVisitMark(0);
  }

  deque<Vertex*> vertexQ;
  for (size_t i = 0; i < numNodes; i++)
  {
    Vertex *v = getNodeAt(i);
    if (v->isBoundary())
    {
      v->setLayerID(0);
      vertexQ.push_back(v);
    }
  }

  vector<Vertex*> neighs;
  while (vertexQ.empty())
  {
    Vertex *currVertex = vertexQ.front();
    vertexQ.pop_front();
    if (currVertex->isVisited())
    {
      int layerid = currVertex->getLayerID();
      neighs = currVertex->getRelations0();
      for (size_t i = 0; i < neighs.size(); i++)
      {
        Vertex *vn = neighs[i];
        if (!vn->isVisited())
        {
          vertexQ.push_back(vn);
          vn->setLayerID(layerid + 1);
        }
      }
      currVertex->setVisitMark(1);
    }
  }

  for (size_t i = 0; i < numNodes; i++)
  {
    Vertex *v = getNodeAt(i);
    assert(v->isVisited());
  }

  if (!relexist)
    clear_relations(0, 0);

  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int Mesh::setFaceWavefront()
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


  deque<Face*> faceQ, nextQ;
  for (size_t i = 0; i < numFaces; i++)
  {
    Face *f = getFaceAt(i);
    if (f->hasBoundaryEdge())
    {
      f->setLayerID(1);
      faceQ.push_back(f);
      /*
      if( f->getTag() == 0 ) {
           cout << f->getNodeAt(0)->isBoundary() << " "
           << f->getNodeAt(1)->isBoundary() << " "
           << f->getNodeAt(2)->isBoundary() << " "
           << f->getNodeAt(3)->isBoundary() << endl;
          f->setTag(5);
          f->setTag(1);
      }
      */
    } else
      assert( f->getTag() == 0 );
  }


  vector<Face*> neighs;
  int layerid = 1;
  while (!faceQ.empty())
  {
    for (size_t i = 0; i < faceQ.size(); i++)
      faceQ[i]->setVisitMark(1);

    nextQ.clear();

    for (size_t j = 0; j < faceQ.size(); j++)
    {
      Face *currFace = faceQ[j];
      neighs = currFace->getRelations212();
      for (size_t i = 0; i < neighs.size(); i++)
      {
        Face *vn = neighs[i];
        if (!vn->isVisited())
          nextQ.push_back(vn);
      }
    }

    layerid++;
    for (size_t i = 0; i < nextQ.size(); i++)
      nextQ[i]->setLayerID(layerid);

    faceQ = nextQ;
  }

  for (size_t i = 0; i < numFaces; i++)
  {
    Face *f = getFaceAt(i);
    assert(f->isVisited());
    int itag = f->getLayerID();
    f->setTag(itag);
  }

  if (!relexist)
    clear_relations(0, 2);

  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int Mesh::setWavefront(int forwhat)
{
  if (forwhat == 0)
    return setNodeWavefront();
  if (forwhat == 2)
    return setFaceWavefront();

  return 1;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::check_unused_objects()
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

///////////////////////////////////////////////////////////////////////////////
bool Mesh::isDelaunay() 
{
  bool retval = 1;
  int relexist = build_relations(0,0);

  Point3D  pa, pb, pc, pd, pCenter;

  size_t numfaces = getSize(2);

  for (size_t i = 0; i < numfaces; i++)
  {
    Face *face = getFaceAt(i);
    face->setTag(1);
  }
   for( size_t i = 0; i < numfaces; i++) {
       Face *face = getFaceAt(i);
       pa = face->getNodeAt(0)->getXYZCoords();
       pb = face->getNodeAt(1)->getXYZCoords();
       pc = face->getNodeAt(2)->getXYZCoords();
       TriCircumCenter3D( &pa[0], &pb[0], &pc[0], &pCenter[0]);
       double radius2 = Math::length2( pa, pCenter );
       vector<Vertex*> neighs = face->getRelations0();
       for( size_t j = 0; j < neighs.size(); j++) {
            pd = neighs[j]->getXYZCoords();
            if( Math::length2( pd, pCenter ) < radius2) {
	        face->setTag(0);
	        retval = 0;
		break;
            }
        }
    }


   if (!relexist)
       clear_relations(0, 0);
}
///////////////////////////////////////////////////////////////////////////////

double Mesh :: getSurfaceArea()
{
   double facearea, sumArea = 0.0;

   size_t numfaces = getSize(2);
   double minarea = MAXDOUBLE;
   double maxarea = 0.0;
   for( size_t i = 0; i < numfaces; i++) {
       Face *face = getFaceAt(i);
       facearea   = face->getArea();
       sumArea    += facearea;
       if( facearea < minarea) minarea = facearea;
       if( facearea > maxarea) maxarea = facearea;
   }

   cout << "Info:   Min face Area : " << minarea << endl;
   cout << "Info:   Max face Area : " << maxarea << endl;

   return sumArea;
}
///////////////////////////////////////////////////////////////////////////////

vector<int> Jaal::getVertexFaceDegrees(iMesh_Instance &imesh)
{
  Mesh *jmesh = new Mesh;
  jmesh->fromMOAB(imesh);
  vector<int> quality = jmesh->get_topological_statistics();
  delete jmesh;
  return quality;
}

///////////////////////////////////////////////////////////////////////////////
vector<float> Mesh :: getAspectRatio( bool sorted )
{
 size_t numfaces = getSize(2);

 vector<float> quality;
 quality.resize(numfaces );

 for( size_t i = 0; i < numfaces; i++) {
      Face *face = getFaceAt(i);
      quality[i] = face->getAspectRatio();
 }

 if( sorted )
     std::sort( quality.begin(), quality.end() );

  float minval = *min_element( quality.begin(), quality.end() );
  float maxval = *max_element( quality.begin(), quality.end() );

  cout << "Info: Minimum Aspect Ratio  " << minval << endl; 
  cout << "Info: Maximum Aspect Ratio  " << maxval << endl; 

 return  quality;
}
///////////////////////////////////////////////////////////////////////////////
void Mesh :: emptyAll()
{
    nodes.clear();
    edges.clear();
    faces.clear();
}
///////////////////////////////////////////////////////////////////////////////

void Mesh :: clearAll()
{
    for( size_t i = 0; i < nodes.size(); i++)
    delete nodes[i];
    nodes.clear();

    for( size_t i = 0; i < edges.size(); i++)
    delete edges[i];
    edges.clear();

    for( size_t i = 0; i < faces.size(); i++)
    delete faces[i];
    faces.clear();
}
///////////////////////////////////////////////////////////////////////////////
int Mesh::get_quality_statistics( const string &fname)
{
  ofstream ofile(fname.c_str(), ios::out);
  if( ofile.fail() ) {
      cout << "Warning: Cann't open file " << fname << endl;
      return 1;
  }

  // Collect Element Area informtion
  size_t numfaces = getSize(2);
  vector<double>  quality(numfaces);

  double minval, maxval, stddev, avgval, medianval;

  for( size_t i = 0; i < numfaces; i++) {
       Face *face = getFaceAt(i);
       quality[i] = face->getArea();
  }

  sort( quality.begin(), quality.end() );
  minval = quality.front();
  maxval = quality.back();

  ofile <<  "# ********************************************* " << endl;
  ofile <<  "# Measure            Element Area        " << endl;
  ofile <<  "# Num of Elements  " <<   numfaces   << endl;
  ofile <<  "# Min              " <<   minval     << endl;
  ofile <<  "# Max              " <<   maxval     << endl;
  ofile <<  "# Average          " <<   avgval     << endl;
  ofile <<  "# MeanVal          " <<   medianval  << endl;
  ofile <<  "# StdDev           " <<   stddev     << endl;
  for( size_t i = 0; i < numfaces; i++)
       ofile << quality[i] << endl;
  ofile <<  " ********************************************* " << endl;


  for( size_t i = 0; i < numfaces; i++) {
       Face *face = getFaceAt(i);
       quality[i] = face->getAspectRatio();
  }

  sort( quality.begin(), quality.end() );
  minval = quality.front();
  maxval = quality.back();

  ofile <<  "# Measure            Aspect Ratio      " << endl;
  ofile <<  "# Num of Elements  " <<   numfaces   << endl;
  ofile <<  "# Min              " <<   minval     << endl;
  ofile <<  "# Max              " <<   maxval     << endl;
  ofile <<  "# Average          " <<   avgval     << endl;
  ofile <<  "# MeanVal          " <<   medianval  << endl;
  ofile <<  "# StdDev           " <<   stddev     << endl;

  for( size_t i = 0; i < numfaces; i++)
       ofile << quality[i] << endl;
  ofile <<  " ********************************************* " << endl;

}
///////////////////////////////////////////////////////////////////////////////

#ifdef USE_MESQUITE

using namespace Mesquite;
using Mesquite::MsqError;

int run_global_smoother( Mesquite::Mesh* mesh, MsqError& err )
{
  double OF_value = 0.0001;

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::SUM, err);
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;

  // creates the feas newt optimization procedures
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func, true );
  pass1->use_global_patch();
  if (err) return 1;

  QualityAssessor stop_qa( mean_ratio );

  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  if (err) return 1;
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(mesh, err);
  cout << " Error " << err << endl;
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "feasible-newton-result.vtk", err);
  if (err) return 1;
  cout << "Wrote \"feasible-newton-result.vtk\"" << endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int run_local_smoother( Mesquite::Mesh* mesh, MsqError& err )
{
  double OF_value = 0.0001;

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::SUM, err);
  if (err) return 1;

  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;

  // creates the smart laplacian optimization procedures
  SmartLaplacianSmoother* pass1 = new SmartLaplacianSmoother( obj_func );

  QualityAssessor stop_qa( mean_ratio );

  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // adds 1 pass of pass1 to mesh_set
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(mesh, err);
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "smart-laplacian-result.vtk", err);
  if (err) return 1;
  cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Jaal::mesh_shape_optimization(iMesh_Instance imesh)
{
  int err;

  iBase_EntitySetHandle rootSet;
  iMesh_getRootSet(imesh, &rootSet, &err);

  Mesquite::MsqError ierr;
  Mesquite::Mesh* mesqmesh = new Mesquite::MsqIMesh(imesh, rootSet, iBase_FACE, ierr, "fixed");
  assert(!ierr);

  Mesquite::PlanarDomain domain(Mesquite::PlanarDomain::XY);

/*
  Mesquite::LaplacianIQ laplacian_smoother;
  laplacian_smoother.run_instructions(mesqmesh, &domain, ierr);
  if (ierr) return 1;
  run_global_smoother( mesqmesh, ierr );
*/

  cout << " Improvment :" << endl;
  Mesquite::ShapeImprovementWrapper shape_wrapper(ierr);
  if (ierr)
  {
    cout << "Shape wrapper error " << ierr << endl;
    exit(2);
  }
  shape_wrapper.run_instructions(mesqmesh, &domain, ierr);

  if (ierr)
  {
    cout << "Error smoothing mesh " << ierr << endl;
    return 1;
  }

  return 0;

}
#endif

///////////////////////////////////////////////////////////////////////////////
