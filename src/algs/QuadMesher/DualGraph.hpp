///////////////////////////////////////////////////////////////////////////////
//
//  Description:  Construct Dual Graph of from a Mesh. The dual of a triangle
//  ************  is its centroid.
//  Dual Graph can be represented as Node-Node Adjacency or through explcitly 
//  creating and storing the edges. When Adjacency is set to (1,0)  Edges are created 
//  and stored. When Adjacency is set to (0,0)  Node-Node relationships are stored.
//
//  Chaman Singh Verma
//  Department of Computer Sciences,
//  University of Wisconsin Madison.
//  Date 15th Jan 2010.
//
//  License:  Free to download and modify. 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef DUALGRAPH_H
#define DUALGRAPH_H

#include <meshkit/Mesh.hpp>

BEGIN_JAAL_NAMESPACE

class DualGraph
{
public:
  static PNode getNewDualNode(Face *f);

  DualGraph()
  {
    adjTable[0][0] = 1;
  }

  ~DualGraph()
  {
    clear();
  }

  void setAdjTable(int src, int dst)
  {
    if (src == 0 && dst == 0)
    {
      if (adjTable[1][0])
        node_node_adjacency_rep();
      return;
    }

    if (src == 1 && dst == 0)
    {
      if (adjTable[0][0])
        edge_rep();
      return;
    }
  }

  size_t getSize(int d) const
  {
    if (d == 0)
      return nodes.size();
    if (d == 1)
      return edges.size();
    return 0;
  }

  int build(Mesh *m);

  const PNode &getNode(int i) const
  {
    return nodes[i];
  }

  const PEdge &getEdge(int i) const
  {
    assert(adjTable[1][0]);
    return edges[i];
  }

  void saveAs(const string &s)
  {
    writeNodes(s);
    writeEdges(s);
  }

  void clear()
  {
    for (size_t i = 0; i < nodes.size(); i++)
      delete nodes[i];
    nodes.clear();

    for (size_t i = 0; i < edges.size(); i++)
      delete edges[i];
    edges.clear();
  }

private:
  char adjTable[4][4];

  Mesh *mesh;

  void node_node_adjacency_rep();
  void edge_rep();

  NodeSequence  nodes;
  EdgeSequence  edges;

  void addEdge(const PNode n1, const PNode n2);
  void getCentroid(int faceid, Point3D &p);
  void writeNodes(const string &s);
  void writeEdges(const string &s);
};

inline
void DualGraph::addEdge(const PNode v1, const PNode v2)
{
  v1->addRelation0(v2);
  v2->addRelation0(v1);
}

inline 
PNode DualGraph::getNewDualNode(Face *face) 
{
  PNode  dualnode = Vertex::newObject();
  Point3D p3d = face->getCentroid();
  dualnode->setXYZCoords(p3d);

  face->setDualNode( dualnode );
  dualnode->setPrimalFace(face);
  return dualnode;
}


END_JAAL_NAMESPACE

#endif
