#include <assert.h>

#include "meshkit/GraphNode.hpp"

namespace MeshKit 
{
    
GraphNode::GraphNode(const GraphNode &node)
{
  for (std::vector<GraphNode*>::const_iterator vit = node.parents().begin(); vit != node.parents().end(); vit++)
    (*vit)->add_child(this);
  
  for (std::vector<GraphNode*>::const_iterator vit = node.children().begin(); vit != node.children().end(); vit++)
    add_child(*vit);
}
  
GraphNode::~GraphNode()
{
  if (!nodeParents.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<GraphNode*>::reverse_iterator vit = nodeParents.rbegin(); vit != nodeParents.rend(); vit++)
      (*vit)->remove_child(this);
  }
  
  if (!nodeChildren.empty()) {
      // remove in reverse order so iterator isn't invalidated
    for (std::vector<GraphNode*>::reverse_iterator vit = nodeChildren.rbegin(); vit != nodeChildren.rend(); vit++)
      remove_child(*vit);
  }
}

void GraphNode::add_parent(GraphNode *par)
{
    // the child in a parent/child pair just takes care of its own parent list,
    // doesn't do any chekcing of its children, that's done by the parent being
    // added/removed
  assert(("Either initializing mkCore" && par->is_core() && par == this) ||
         ("Or par isn't already a parent" && !par->is_core() &&
          std::find(nodeParents.begin(), nodeParents.end(), par) == nodeParents.end()));
  
  nodeParents.push_back(par); 
}

void GraphNode::remove_parent(GraphNode *par)
{
    // the child in a parent/child pair just takes care of its own parent list,
    // doesn't do any chekcing of its children, that's done by the parent being
    // added/removed
  std::vector<GraphNode*>::iterator vit = std::find(nodeParents.begin(), nodeParents.end(), par);
  if (vit != nodeParents.end()) nodeParents.erase(vit);
  else throw Error(MK_FAILURE, "Parent operation not in parent list.");
}

void GraphNode::add_child(GraphNode *child)
{
  assert(("Either initializing mkCore" && child->is_core() && child == this) ||
         ("Or child isn't already a child" &&
          std::find(nodeChildren.begin(), nodeChildren.end(), child) == nodeChildren.end()));

    // take care of me first; if no children, will no longer be a leaf, so remove as parent of mkCore
  if (nodeChildren.size() == 1 && nodeChildren[0]->is_core()) {
    nodeChildren[0]->remove_parent(this);
    child->add_child(nodeChildren[0]);
    nodeChildren[0] = child;
  }
  else
    nodeChildren.push_back(child);

    // add myself as a parent of child
  child->add_parent(this);
}

void GraphNode::remove_child(GraphNode *child)
{
  std::vector<GraphNode*>::iterator vit = std::find(nodeChildren.begin(), nodeChildren.end(), child);
  if (vit != nodeChildren.end()) nodeChildren.erase(vit);
  else throw Error(MK_FAILURE, "Child operation not in child list.");

    // if no more children, I'm now a leaf, add as parent of mkCore
  if (nodeChildren.empty()) {
    GraphNode *root = child->get_root();
    assert(root);
    root->add_parent(this);
    nodeChildren.push_back(root);
  }

    // remove myself as parent of child
  child->remove_parent(this);
}

GraphNode *GraphNode::get_root() 
{
  assert("GraphNode should never have no children" && !nodeChildren.empty());
  GraphNode *node = nodeChildren[0];
  while (!node->is_root()) {
    assert("Child shouldn't be empty enther." && !node->nodeChildren.empty());
    node = node->nodeChildren[0];
  }
  
  return node;
}

bool GraphNode::is_root() 
{
  return false;
}

  //! Setup function, called in reverse order before execute
void GraphNode::setup() 
{
  for (std::vector<GraphNode*>::iterator vit = nodeParents.begin(); vit != nodeParents.end(); vit++)
    (*vit)->setup();
    
  setup_this();
}
  
  //! Execute function, called in forward order after setup
void GraphNode::execute() 
{
  execute_this();

  for (std::vector<GraphNode*>::iterator vit = nodeChildren.begin(); vit != nodeChildren.end(); vit++)
    (*vit)->execute();
}

}

  
