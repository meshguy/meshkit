#ifndef GRAPHNODE_HPP
#define GRAPHNODE_HPP

#include "lemon/list_graph.h"
#include <string>

namespace MeshKit {

class MKGraph;
    
/** \class GraphNode GraphNode.hpp "meshkit/GraphNode.hpp"
 * \brief A parent class for anything that will be a graph node
 *
 * This class encapsulates operations on graph nodes, like inserting/removing
 * in graph, inserting new parents, etc.  This class HasA Lemon node object 
 * (which is really just an index into the graph node vector).  Instances of
 * this class keep a pointer to a MKGraph object, which is really usually just
 * the MKCore instance, where the graph itself is stored.
 */
class GraphNode
{
public:

    //! Friend class, to allow more direct access

  friend class MKGraph;
  
    //! Copy constructor
  GraphNode(const GraphNode &graph_node);
  
    //! Bare constructor
  GraphNode(MKGraph *graph);

    //! Destructor
  virtual ~GraphNode();
  
    //! Get the associated MKGraph object
  MKGraph *get_graph() const;

    //! Get the graph node corresponding to this GraphNode
  lemon::ListDigraph::Node get_node() const;

    //! Return an iterator over incoming graph edges
  lemon::ListDigraph::InArcIt in_arcs() const;
  
    //! Return an iterator over outgoing graph edges
  lemon::ListDigraph::OutArcIt out_arcs() const;

    /** \brief Return the GraphNode at the other end of a connected graph edge
     * \param arc Edge being queried
     * \return GraphNode corresponding to the other graph node
     */
  GraphNode *other_node(lemon::ListDigraph::Arc arc);

    //! Get node name
  virtual std::string get_name() const;

    //! Set node name
  virtual void set_name(std::string new_name);

    //! Setup function, called in reverse order before execute
  virtual void setup();

    //! Execute function, called in forward order after setup
  virtual void execute();
  
    //! Pure virtual, derived class must define
  virtual void setup_this();

    //! Pure virtual, derived class must define
  virtual void execute_this();

protected:
    //! MKGraph associated with this GraphNode
  MKGraph *mkGraph;

    //! The graph node associated with this GraphNode
  lemon::ListDigraph::Node graphNode;

    //! Local name for this node
  std::string nodeName;

private:
};

inline GraphNode::GraphNode(const GraphNode &node) 
{
    // create a node for this meshop
  graphNode = get_graph()->get_graph().addNode();

    // this is a copy; set in/out edges from copied node
  for (lemon::ListDigraph::InArcIt ait(get_graph()->get_graph(), node.get_node());
       ait != lemon::INVALID; ++ait) 
    get_graph()->get_graph().addArc(get_graph()->source(ait)->get_node(), graphNode);
  for (lemon::ListDigraph::OutArcIt ait(get_graph()->get_graph(), node.get_node());
       ait != lemon::INVALID; ++ait) 
    get_graph()->get_graph().addArc(graphNode, get_graph()->target(ait)->get_node());
}
    
inline GraphNode::GraphNode(MKGraph *graph) 
{
    // create a Lemon node 
  graphNode = graph->get_graph().addNode();

    // set the nodemap value for this meshop
  graph->node_map()[graphNode] = this;
  
    // might not be root/leaf nodes yet, if this is a root or leaf node being created
  if (graph->root_node() && graph->leaf_node()) {
    graph->get_graph().addArc(graph->root_node()->get_node(), graphNode);
    graph->get_graph().addArc(graphNode, graph->leaf_node()->get_node());
  }
}

inline GraphNode::~GraphNode() 
{
  get_graph()->get_graph().erase(graphNode);
}

inline std::string GraphNode::get_name() const
{
  return nodeName;
}

inline void GraphNode::set_name(std::string new_name)
{
  nodeName = new_name;
}

inline lemon::ListDigraph::InArcIt GraphNode::in_arcs() const
{
  return lemon::ListDigraph::InArcIt(get_graph()->get_graph(), graphNode);
}

inline lemon::ListDigraph::OutArcIt GraphNode::out_arcs() const
{
  return lemon::ListDigraph::OutArcIt(get_graph()->get_graph(), graphNode);
}

inline MKGraph *GraphNode::get_graph() const
{
  return mkGraph;
}

    //! Get the graph node corresponding to this GraphNode
inline lemon::ListDigraph::Node GraphNode::get_node() const
{
  return graphNode;
}

inline GraphNode *GraphNode::other_node(lemon::ListDigraph::Arc arc)
{
  return get_graph()->other_node(arc, this);
}

    //! Setup function, called in reverse order before execute
inline void GraphNode::setup() 
{
  setup_this();

  for (lemon::ListDigraph::InArcIt ait = in_arcs(); ait != lemon::INVALID; ++ait) 
    mkGraph->source(ait)->setup();
}
  
    //! Execute function, called in forward order after setup
inline void GraphNode::execute() 
{
  execute_this();

  for (lemon::ListDigraph::OutArcIt ait = out_arcs(); ait != lemon::INVALID; ++ait) 
    mkGraph->target(ait)->execute();
}

    //! Pure virtual, derived class must define
inline void GraphNode::setup_this() 
{}

    //! Pure virtual, derived class must define
inline void GraphNode::execute_this() 
{}

} // namespace MeshKit

#endif

  
