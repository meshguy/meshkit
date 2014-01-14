///////////////////////////////////////////////////////////////////////////////
//  Description: Builds a binary tree from a dual graph of triangualated mesh.
//  Chaman Singh Verma
//  University of Wisconsin Madison, USA
//  Date 15th Jan 2010
//
//  License: Free to distribute and modify.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef BTREE_H
#define BTREE_H

#include "Mesh.hpp"
#include "DualGraph.hpp"
#include "ObjectPool.hpp"

#ifdef USE_HASHMAP
#include <tr1/unordered_map>
#endif

//#include <ext/slist>
//typedef __gnu_cxx::slist<BinaryNode*> BNodeList;

namespace Jaal {

class BinaryNode;

typedef std::list<BinaryNode*> BNodeList;

class BinaryNode {
public:
     BinaryNode()
     {}

     BinaryNode(PNode n) {
          dualNode = n;
          levelID = -1;
          parent = NULL;
     }

     void setLevelID(int l) {
          levelID = l;
     }

     int getLevelID() const {
          return levelID;
     }

     bool isMatched() const {
          return dualNode->isRemoved();
     }

     void setMatchMark(char r) {
          dualNode->setStatus(r);
     }

     bool isRoot() const {
          if (parent == NULL)
               return 1;
          return 0;
     }

     bool isLeaf() const {
          if (getNumChildren() == 0)
               return 1;
          return 0;
     }

     int getDegree() const {
          int ncount = 0;
          if (parent)
               ncount = 1;
          ncount += getNumChildren();
          return ncount;
     }

     BinaryNode * getSibling() const {
          if (parent) {
               for (int i = 0; i < parent->getNumChildren(); i++) {
                    BinaryNode *child = parent->getChild(i);
                    if (child != this)
                         return child;
               }
          }
          return NULL;
     }

     void setParent(BinaryNode * p) {
          parent = p;
     }

     BinaryNode * getParent() const {
          return parent;
     }

     int getNumChildren() const {
          int ncount = 0;
          for (size_t i = 0; i < children.size(); i++)
               if (children[i].active)
                    ncount++;
          return ncount;
     }

     void addChild(BinaryNode * c) {
          ActiveNode activenode(c);
          if (find(children.begin(), children.end(), activenode) == children.end())
               children.push_back(activenode);
     }

     void removeChild(BinaryNode * child) {
          ActiveNode activenode(child);
          vector<ActiveNode>::iterator it;
          it = remove(children.begin(), children.end(), activenode);
          children.erase(it, children.end());
     }

     void unlinkChild(BinaryNode * child) {
          for (size_t i = 0; i < children.size(); i++) {
               if (children[i].node == child) {
                    children[i].active = 0;
                    return;
               }
          }
     }

     void relinkChild(BinaryNode * child) {
          for (size_t i = 0; i < children.size(); i++) {
               if (children[i].node == child) {
                    children[i].active = 1;
                    return;
               }
          }
     }

     void relinkAll() {
          for (size_t i = 0; i < children.size(); i++)
               children[i].active = 1;
     }

     BinaryNode * getChild(int cid) const {
          int index = -1;
          for (size_t i = 0; i < children.size(); i++) {
               if (children[i].active) {
                    index++;
                    if (index == cid)
                         return children[i].node;
               }
          }
          return NULL;
     }

     int getID() const {
          return dualNode->getID();
     }

     Vertex * getDualNode() const {
          return dualNode;
     }

private:
     struct ActiveNode {
          ActiveNode(BinaryNode *n) {
               active = 1;
               node = n;
          }
          bool active;
          bool operator ==(const ActiveNode &rhs) const {
               return node == rhs.node;
          }
          BinaryNode *node;
     };
     int levelID;
     PNode dualNode;
     BinaryNode* parent;

     vector<ActiveNode> children;
};

class BinaryTree {
public:

     const static int BREADTH_FIRST_TREE = 0;
     const static int DEPTH_FIRST_TREE   = 1;

     BinaryTree(DualGraph *g) {
          dgraph = g;
          treetype = BREADTH_FIRST_TREE;
     }

     ~BinaryTree() {
          clear();
     }

     void setTreeType(int t) {
          treetype = t;
     }

     void build(BinaryNode *r = NULL);

     BinaryNode* getRoot() const {
          return root;
     }

     // Remove a given node from the tree. Unlink child and parents.
     //
     void removeNode(BinaryNode *node) {
          if (node) {
               BinaryNode *parv = node->getParent();
               if (parv)
                    parv->removeChild(node);
          }
     }

     void addNode(BinaryNode *tnode) {
          tnodemap[tnode->getDualNode()] = tnode;
     }

     void unlinkNode(BinaryNode *node) {
          if (node) {
               BinaryNode *parv = node->getParent();
               if (parv)
                    parv->unlinkChild(node);
          }
     }

     bool isMatched(const BinaryNode *u, const BinaryNode *v) const {
          Vertex *du = u->getDualNode();
          Vertex *dv = v->getDualNode();
          Vertex *umate = NULL;
          Vertex *vmate = NULL;
          du->getAttribute("DualMate", umate);
          dv->getAttribute("DualMate", vmate);
          if (umate == dv && vmate == du) return 1;
/*
          if (du->getDualMate() == dv && dv->getDualMate() == du)
               return 1;
*/
          return 0;
     }

     // Give nodes at a given level. Root is at level = 0,
     const BNodeList &getLevelNodes(int l) const;

     // Total Height of the tree...
     int getHeight() const;

     // How many nodes in the tree.
     size_t getSize() const;

     // Clear all the nodes created in the tree.
     void clear();

     void deleteAll();

     void relinkAll();

     // Save the tree in GraphViz data format...
     void saveAs(const string &s);

private:
     DualGraph *dgraph;
     BinaryNode *root;
     BNodeList emptylist;
     ObjectPool<BinaryNode> pool;

     int treetype;

#ifdef USE_HASHMAP
     typedef std::tr1::unordered_map<Vertex*, BinaryNode*>  TNodeMap;
#else
     typedef std::map<Vertex*, BinaryNode*> TNodeMap;
#endif
     TNodeMap  tnodemap;

     std::map<int, BNodeList> levelnodes;

     void bfs_traverse(BinaryNode *parent, BNodeList &nextnodes);
     void bfs_traverse(BinaryNode *parent);

     void dfs_traverse(BinaryNode *parent);
     void dfs_traverse(BinaryNode *parent, BNodeList &nextnodes);
};

} // namespace Jaal

#endif
