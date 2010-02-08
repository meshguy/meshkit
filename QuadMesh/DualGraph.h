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

#include "Mesh.h"

BEGIN_JAAL_NAMESPACE

class DualGraph
{
public:

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

  NodeType getNode(int i)
  {
    return nodes[i];
  }

  EdgeType getEdge(int i) const
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

  vector<NodeType> nodes;
  vector<EdgeType> edges;

  void addEdge(const NodeType n1, const NodeType n2);
  void getCentroid(int faceid, Point3D &p);
  void writeNodes(const string &s);
  void writeEdges(const string &s);
};

END_JAAL_NAMESPACE

#endif