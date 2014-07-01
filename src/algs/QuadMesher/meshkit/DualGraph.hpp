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

#include "meshkit/Mesh.hpp"

namespace Jaal {

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
      nodes.clear();
      edges.clear();
  }

  void deleteAll()
  {
    size_t nSize;

    nSize = nodes.size();
    for (size_t i = 0; i < nSize; i++)
      delete nodes[i];
    nodes.clear();

    nSize = edges.size();
    for (size_t i = 0; i < nSize; i++)
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
  v1->addRelation(v2);
  v2->addRelation(v1);
}

inline 
PNode DualGraph::getNewDualNode(Face *face) 
{
  PNode dualnode = Vertex::newObject();
  Point3D p3d;
  face->getAvgPos( p3d );
  dualnode->setXYZCoords(p3d);

  face->setAttribute("DualNode", dualnode);
  dualnode->setAttribute("PrimalFace", face);

  return dualnode;
}


} // namespace Jaal

#endif
