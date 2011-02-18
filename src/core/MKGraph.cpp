#include "meshkit/MKGraph.hpp"
#include "meshkit/MeshOp.hpp"
#include "lemon/core.h"
#include "lemon/adaptors.h"
#include "lemon/connectivity.h"

namespace MeshKit 
{
    
MKGraph::MKGraph() 
        : mkGraph(), nodeMap(mkGraph, NULL), rootNode(NULL), leafNode(NULL)
{
}
    
MKGraph::~MKGraph() 
{
}

void MKGraph::clear_graph() 
{
    // get all the non-root, non-leaf nodes
  std::vector<lemon::ListDigraph::Node> nodes;
  for (lemon::ListDigraph::NodeIt nit(mkGraph); nit != lemon::INVALID; ++nit) 
    if (nit != rootNode->get_node() && nit != leafNode->get_node())
      nodes.push_back(nit);
  
    // now delete all those nodes
  for (std::vector<lemon::ListDigraph::Node>::iterator vit = nodes.begin(); vit != nodes.end(); vit++) {
    if (nodeMap[*vit]) delete nodeMap[*vit];
    else mkGraph.erase(*vit);
  }
  
    // restore an edge between the root and leaf
  if (!mkGraph.valid(lemon::ArcLookUp<lemon::ListDigraph>(mkGraph)(rootNode->get_node(), leafNode->get_node())))
    mkGraph.addArc(rootNode->get_node(), leafNode->get_node());
}

void MKGraph::print_graph() 
{
  for (lemon::ListDigraph::NodeIt nit(mkGraph); nit != lemon::INVALID; ++nit) {
    std::cout << "Node: ";
    if (nit == rootNode->get_node()) std::cout << "root node" << std::endl;
    else if (nit == leafNode->get_node()) std::cout << "leaf node" << std::endl;
    else if (nodeMap[nit]) std::cout << nodeMap[nit]->get_name() << std::endl;
    else std::cout << "(no MeshOp)" << std::endl;
  }
}

    //! Get the MeshOp corresponding to a graph node
MeshOp *MKGraph::get_meshop(lemon::ListDigraph::Node node) const
{
  return dynamic_cast<MeshOp*>(nodeMap[node]);
}
  
GraphNode *MKGraph::other_node(lemon::ListDigraph::Arc arc, GraphNode *node) const 
{
  lemon::ListDigraph::Node src = mkGraph.source(arc);
  if (src != node->get_node()) return get_node(src);
  else return get_node(mkGraph.target(arc));
}

MeshOp *MKGraph::find_meshop(std::string op_name) const
{
  GraphNode *node = find_node(op_name);
  if (node) return dynamic_cast<MeshOp*>(node);
}
    
GraphNode *MKGraph::find_node(std::string op_name) const
{
    // run BFS on forward graph
  lemon::Bfs<lemon::ListDigraph> bfs(mkGraph);
  bfs.init();
  bfs.addSource(rootNode->get_node());
  while (!bfs.emptyQueue()) {
    lemon::ListDigraph::Node nd = bfs.processNextNode();
    assert(nd != lemon::INVALID && (nodeMap[nd] || (nd == leafNode->get_node() || nd == rootNode->get_node())));
    if (nodeMap[nd] && nodeMap[nd]->get_name() == op_name) return nodeMap[nd];
  }
  return NULL;
}

void MKGraph::insert_node(GraphNode *inserted, GraphNode *before) 
{
    // get the corresponding Lemon nodes
  lemon::ListDigraph::Node linserted = inserted->get_node(), 
      lbefore = before->get_node();
  
    // if inserted is a leaf node (i.e. is connected to leafNode), disconnect from that
  if (mkGraph.target(inserted->out_arcs()) == leafNode->get_node())
    mkGraph.erase(inserted->out_arcs());

    // if before is a root node (i.e. is connected to rootNode), also disconnect that
  if (mkGraph.source(before->in_arcs()) == rootNode->get_node())
    mkGraph.erase(before->in_arcs());
    
    // now link inserted to the root, and to before
  mkGraph.addArc(rootNode->get_node(), inserted->get_node());
  mkGraph.addArc(inserted->get_node(), before->get_node());
}

void MKGraph::add_arc(GraphNode *source, GraphNode *target) 
{
    // add an arc from one node to another, e.g. to add a dependency between them
    // get the corresponding Lemon nodes
  lemon::ListDigraph::Node lsource = source->get_node(), 
      ltarget = target->get_node();
  
    // if inserted is a leaf node (i.e. is connected to leafNode), disconnect from that
  if (mkGraph.target(source->out_arcs()) == leafNode->get_node())
    mkGraph.erase(source->out_arcs());

    // if before is a root node (i.e. is connected to rootNode), also disconnect that
  if (mkGraph.source(target->in_arcs()) == rootNode->get_node())
    mkGraph.erase(target->in_arcs());
    
    // now link them
  mkGraph.addArc(lsource, ltarget);
}

//! Run setup on the graph
void MKGraph::setup() 
{
    // run BFS on reversed graph
  lemon::ReverseDigraph<lemon::ListDigraph> rg(mkGraph);
  lemon::Bfs<lemon::ReverseDigraph<lemon::ListDigraph> > rbfs1(rg);
  rbfs1.init();
  rbfs1.addSource(leafNode->get_node());
  while (!rbfs1.emptyQueue()) {
    lemon::ListDigraph::Node nd = rbfs1.processNextNode();
    assert(nd != lemon::INVALID && (nodeMap[nd] || (nd == leafNode->get_node() || nd == rootNode->get_node())));
    if (nodeMap[nd]) nodeMap[nd]->setup_called(false);
  }

  bool called_one;
  do {
    called_one = false;
    lemon::Bfs<lemon::ReverseDigraph<lemon::ListDigraph> > rbfs2(rg);
    rbfs2.init();
    rbfs2.addSource(leafNode->get_node());
    while (!rbfs2.emptyQueue()) {
      lemon::ListDigraph::Node nd = rbfs2.processNextNode();
      assert(nd != lemon::INVALID && (nodeMap[nd] || (nd == leafNode->get_node() || nd == rootNode->get_node())));
      if (nodeMap[nd] && !nodeMap[nd]->setup_called()) {
        nodeMap[nd]->setup_this();
        nodeMap[nd]->setup_called(true);
        called_one = true;
      }
    }
  }
  while (called_one);
}

//! Run execute on the graph
void MKGraph::execute() 
{

  typedef lemon::IterableIntMap<lemon::ListDigraph, lemon::ListDigraph::Node> topomap;

  // Run execute_this on all nodes in topological order
  topomap topo_levels( mkGraph ); 
  lemon::topologicalSort( mkGraph, topo_levels );

  for( unsigned int i = 0; i<topo_levels.size(); ++i){
    for( topomap::ItemIt j(topo_levels, i); j != lemon::INVALID; ++j ){
      GraphNode* gn = nodeMap[ j ];
      assert( gn );
      gn->execute_called(false);
    }
  }
  
  for( unsigned int i = 0; i<topo_levels.size(); ++i){

    for( topomap::ItemIt j(topo_levels, i); j != lemon::INVALID; ++j ){

      GraphNode* gn = nodeMap[ j ];
      assert( gn );
      if (!gn->execute_called()) {
        gn->execute_this();
        gn->execute_called(true);
      }
    }
  }

}

} // namespace MeshKit
