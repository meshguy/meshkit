#ifndef MKGRAPH_HPP
#define MKGRAPH_HPP

/** \file MKGraph.hpp
 */
#include "lemon/list_graph.h"
#include "lemon/bfs.h"

namespace MeshKit {

class MeshOp;
class GraphNode;
    
/** \class MKGraph MKGraph.hpp "meshkit/MKGraph.hpp"
 * \brief Class encapsulating manipulations of the meshing operation graph
 *
 * Mesh generation can be phrased as a graph of operations, where each node in the graph is a
 * distinct operation (some of them generating mesh, others just operating on mesh) and each arc
 * a dependency between one operation and another.  The graph in MeshKit has single root and leaf
 * node, just as a convenience for starting up forward and reverse traversals of the graph.
 * This class also defines the typedefs used for certain graph data structures like nodes and arcs.
 */

class MKGraph
{
public:

    /** \brief Bare constructor
     */
  MKGraph();
  
    //! destructor
  ~MKGraph();

    //! Get the graph
  lemon::ListDigraph &get_graph();

    //! Get the (const) graph
  const lemon::ListDigraph &get_graph() const;

    //! Get root node
  GraphNode *root_node() const;
  
    //! Get leaf node
  GraphNode *leaf_node() const;

    //! Get access to the node map
  lemon::ListDigraph::NodeMap<GraphNode*> &node_map();

    //! Get access to the (const) node map
  const lemon::ListDigraph::NodeMap<GraphNode*> &node_map() const;

    //! Get the MeshOp corresponding to a graph node
  MeshOp *get_meshop(lemon::ListDigraph::Node node) const;
  
    //! Get the GraphNode corresponding to a graph node
  GraphNode *get_node(lemon::ListDigraph::Node node) const;
  
    //! Get the GraphNode corresponding to the source node of this arc
  GraphNode *source(lemon::ListDigraph::Arc arc) const;
  
    //! Get the GraphNode corresponding to the target node of this arc
  GraphNode *target(lemon::ListDigraph::Arc arc) const;

    /** \brief Get the GraphNode at the other end of an arc from that specified
     *
     * \param arc The graph arc being queried
     * \param node The other node
     * \return The node on the other end of the specified arc
     */
  GraphNode *other_node(lemon::ListDigraph::Arc arc, GraphNode *node) const;
  
    /** \brief Find an existing GraphNode in the graph, starting from the root
     * \param op_name GraphNode name being requested
     * \return Pointer to GraphNode found, NULL if not found
     */
  GraphNode *find_node(std::string op_name) const;

    /** \brief Find an existing MeshOp in the graph, starting from the root
     * \param op_name GraphNode name being requested
     * \return Pointer to MeshOp found, NULL if not found
     */
  MeshOp *find_meshop(std::string op_name) const;

    /** \brief Insert a node in the graph just before the specified node
     * \param inserted Node being inserted
     * \param before Node before which new node is inserted
     */
  void insert_node(GraphNode *inserted, GraphNode *before);

    /** \brief add a graph edge from one node to another
     * \param source Source node
     * \param target Target node
     */
  void add_arc(GraphNode *source, GraphNode *target);
  
    //! Run setup on the graph
  virtual void setup();

    //! Run execute on the graph
  virtual void execute();

    //! Run both setup and execute on the graph
  virtual void setup_and_execute();

protected:
    //! The GraphNode graph
  lemon::ListDigraph mkGraph;

    //! Root and leaf nodes of graph
  GraphNode *rootNode, *leafNode;

    //! Map from graph nodes to GraphNodes
  lemon::ListDigraph::NodeMap<GraphNode*> nodeMap;

private:

    /** \brief Copy constructor, not implemented
     *
     * Lemon does not allow copy construction of graphs, so disallow it in MK too
     */
  MKGraph(const MKGraph &);
  
};

    //! Get (const) graph
inline const lemon::ListDigraph &MKGraph::get_graph() const
{
  return mkGraph;
}

    //! Get graph
inline lemon::ListDigraph &MKGraph::get_graph()
{
  return mkGraph;
}

    //! Get root node
inline GraphNode *MKGraph::root_node() const
{
  return rootNode;
}

    //! Get leaf node
inline GraphNode *MKGraph::leaf_node() const
{
  return leafNode;
}
  
    //! Get access to the node map
inline const lemon::ListDigraph::NodeMap<GraphNode*> &MKGraph::node_map() const
{
  return nodeMap;
}
  
    //! Get access to the node map
inline lemon::ListDigraph::NodeMap<GraphNode*> &MKGraph::node_map() 
{
  return nodeMap;
}
  
    //! Get the GraphNode corresponding to a graph node
inline GraphNode *MKGraph::get_node(lemon::ListDigraph::Node node) const
{
  return nodeMap[node];
}
  
inline GraphNode *MKGraph::source(lemon::ListDigraph::Arc arc) const
{
  return nodeMap[mkGraph.source(arc)];
}

inline GraphNode *MKGraph::target(lemon::ListDigraph::Arc arc) const
{
  return nodeMap[mkGraph.target(arc)];
}
  
    //! Run both setup and execute on the graph
inline void MKGraph::setup_and_execute() 
{
  setup();
  
  execute();
}

} // namespace MeshKit

#endif

  
