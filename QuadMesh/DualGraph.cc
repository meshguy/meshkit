#include "DualGraph.h"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void DualGraph::addEdge(const NodeType v1, const NodeType v2)
{
  v1->addRelation0(v2);
  v2->addRelation0(v1);
}

///////////////////////////////////////////////////////////////////////////////

int DualGraph::build(Mesh *m)
{
  mesh = m;
  nodes.clear();
  edges.clear();

  size_t numnodes = mesh->getSize(0);
  size_t numfaces = mesh->getSize(2);

  int relexist = mesh->build_relations(0, 2);

  Vertex *dualvtx;

  nodes.resize(numfaces);

  Point3D p3d;
  Face *face;
  for (size_t iface = 0; iface < numfaces; iface++)
  {
    face = mesh->getFace(iface);
    assert(face);
    face->setID(iface);
    dualvtx = face->getNewDualNode();
    dualvtx->setID(iface);
    nodes[iface] = dualvtx;
  }

  NodeType v0, v1, dv0, dv1;

  vector<FaceType> neighs;
  for (size_t iface = 0; iface < numfaces; iface++)
  {
    face = mesh->getFace(iface);
    int nc = face->getSize(0);
    for (int j = 0; j < nc; j++)
    {
      v0 = face->getConnection((j + 0) % nc);
      v1 = face->getConnection((j + 1) % nc);
      neighs = Mesh::getRelations112(v0, v1);
      if (neighs.size() == 2)
      {
        dv0 = neighs[0]->getDualNode();
        dv1 = neighs[1]->getDualNode();
        addEdge(dv0, dv1);
      } else
      {
        dv0 = neighs[0]->getDualNode();
        dv0->setBoundaryMark(1);
      }
    }
  }

  if (!relexist)
    mesh->clear_relations(0, 2);


  adjTable[0][0] = 1;
  mesh->setDual(1);
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

void DualGraph::node_node_adjacency_rep()
{
  assert( adjTable[1][0]);

  int numedges = edges.size();

  for (int i = 0; i < numedges; i++)
  {
    Edge *edge = edges[i];
    Vertex *v0 = edge->getConnection(0);
    Vertex *v1 = edge->getConnection(1);
    addEdge(v0, v1);
    delete edge;
  }
  edges.clear();

  adjTable[1][0] = 0;
  adjTable[0][0] = 1;
}
///////////////////////////////////////////////////////////////////////////////

void DualGraph::edge_rep()
{
  vector<NodeType> vneighs, eConnect(2);

  assert( adjTable[0][0]);

  edges.clear();

  size_t numfaces = mesh->getSize(2);
  edges.reserve(2 * numfaces);
  for (size_t iface = 0; iface < numfaces; iface++)
  {
    Vertex *dualvtx = nodes[iface];
    vneighs = dualvtx->getRelations0();
    for (size_t i = 0; i < vneighs.size(); i++)
    {
      Edge *newedge = new Edge(dualvtx, vneighs[i]);
      edges.push_back(newedge);
    }
    dualvtx->clearRelations(0);
  }

  adjTable[1][0] = 1;
  adjTable[0][0] = 0;
}
///////////////////////////////////////////////////////////////////////////////

void DualGraph::writeNodes(const string &fname)
{
  string filename = fname + ".node";

  ofstream ofile(filename.c_str(), ios::out);
  if (ofile.fail())
  {
    cout << "Warning: cann't open node file " << filename << endl;
    return;
  }

  int numnodes = nodes.size();
  ofile << numnodes << " 2 0 0" << endl;

  for (int i = 0; i < numnodes; i++)
  {
    Point3D p3d = nodes[i]->getXYZCoords();
    ofile << i << " " << p3d[0] << " " << p3d[1] << " " << p3d[1] << endl;
  }
}

///////////////////////////////////////////////////////////////////////////////

void DualGraph::writeEdges(const string &fname)
{

  string filename = fname + ".edge";
  ofstream ofile(filename.c_str(), ios::out);
  if (ofile.fail())
  {
    cout << "Warning: cann't open node file " << filename << endl;
    return;
  }

  edge_rep();

  int numedges = edges.size();
  ofile << numedges << " 0 " << endl;

  for (int i = 0; i < numedges; i++)
  {
    int id0 = edges[i]->getConnection(0)->getID();
    int id1 = edges[i]->getConnection(1)->getID();
    ofile << i << " " << id0 << " " << id1 << endl;
  }

}

///////////////////////////////////////////////////////////////////////////////