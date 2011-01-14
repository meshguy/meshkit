#ifndef GRAPH_NODE_HPP
#define GRAPH_NODE_HPP

#include <vector>

namespace MeshKit {

/** \class GraphNode GraphNode.hpp "meshkit/GraphNode.hpp"
 * \brief A node in the mesh operation graph
 *
 * This class encapsulates the functionality for placing nodes in the operation graph,
 * for setup/execution of the graph, and for setup/execution of individual nodes.
 * All meshing operations in MeshKit are derived from this class.  The MKCore instance
 * also derives from this class, and serves as both the single root node and the single leaf
 * node of the graph.
 *
 */
class GraphNode
{
public:

    //! Copy constructor
  GraphNode(const GraphNode &mesh_op);
  
    //! Constructor
  GraphNode() {}
  
    //! Destructor
  virtual ~GraphNode();
  
    //! Return parent meshops
  std::vector<GraphNode*> &parents();
      
    //! Return parent meshops
  virtual const std::vector<GraphNode*> &parents() const;
      
    //! Return child meshops
  virtual std::vector<GraphNode*> &children();
      
    //! Return child meshops
  virtual const std::vector<GraphNode*> &children() const;
      
    //! Add parent meshop
  virtual void add_parent(GraphNode *par);

    //! Remove parent meshop
  virtual void remove_parent(GraphNode *par);

    //! Add child meshop
  virtual void add_child(GraphNode *child);

    //! Remove child meshop
  virtual void remove_child(GraphNode *child);

    //! Traverse downward to find the root/leaf of the graph
  virtual GraphNode *get_root();

    //! Return whether I'm the root
  virtual bool is_root();
  
    //! Setup function, called in reverse order before execute
  virtual void setup();

    //! Execute function, called in forward order after setup
  virtual void execute();
  
    //! Virtual, so that if derived class doesn't define, this will be a no-op
  virtual void setup_this();

    //! Virtual, so that if derived class doesn't define, this will be a no-op
  virtual void execute_this();

protected:
    //! Parent GraphNode's
  std::vector<GraphNode*> nodeParents;

    //! Child GraphNode's
  std::vector<GraphNode*> nodeChildren;

private:

};

inline std::vector<GraphNode*> &GraphNode::parents() 
{
  return nodeParents;
}

inline const std::vector<GraphNode*> &GraphNode::parents() const
{
  return nodeParents;
}

inline std::vector<GraphNode*> &GraphNode::children() 
{
  return nodeChildren;
}
    
inline const std::vector<GraphNode*> &GraphNode::children() const
{
  return nodeChildren;
}
    
inline void GraphNode::setup_this()
{
}
    
inline void GraphNode::execute_this()
{
}
    
}

#endif

  
