////////////////////////////////////////////////////////////////////////////////
//            Jaal:: Triangle-to-Quad Transformation Code 
//            ********************************************
//  Description:
//  Given a triangulated mesh, convert into All-Quads or Quad-Dominated Mesh.
//  using Maximum Tree Matching Algorithm.
//
//  Maximum Cardinality using Edmond's Algorithm may give perfect matching
//  but the algorithm is too expensive for large dataset and often perfect
//  matching is not required.
//
//  If the user requires All-Quad Meshing then all the unmatched triangles
//  using Tree Matching Algorithms are split and steiner points are created.
//  In most of the cases, number of Steiner points are less than 5% of the
//  total triangles present in the original mesh.
//
//  Chaman Singh Verma
//  University of Wisconsin Madison, USA,
//  Date: 15th Jan 2010.
//  
//  License: Free to download and modify. 

//  For more information, please follow the link
//
//  http://pages.cs.wisc.edu/~csverma/CS899_09/JQuad.html
//
////////////////////////////////////////////////////////////////////////////////

#ifndef Tri2Quad_H
#define Tri2Quad_H

#include "meshkit/Mesh.hpp"
#include "meshkit/DualGraph.hpp"
#include "meshkit/BinaryTree.hpp"
#include "meshkit/cycle.hpp"         // performance counter.

namespace Jaal {

class Tri2Quads
{
public:

  const static int ALL_QUADS = 0;
  const static int QUAD_DOMINATED = 1;

  static int verify(Mesh *m, const vector<FacePair> &matching);
  static Mesh* collapse_matched_triangles(Mesh *mesh, const vector<FacePair> &matching, int replace = 0);


  Tri2Quads()
  {
    trimesh = NULL;
    btree = NULL;
    verbose = 1;
    required_topology = ALL_QUADS;
  }

  const vector<FacePair> &getMaximumDualMatching();

  Mesh* getQuadMesh(Mesh *tmesh, int replace = 0, int topo = ALL_QUADS);

  NodeSequence getSteinerNodes()  const;
  FaceSequence getInsertedFaces() const;
  FaceSequence getModifiedFaces() const;

private:
  Mesh *trimesh; // Input mesh.

  struct LVertex : public Vertex
  {
      LVertex( Vertex *v ) { vertex = v; }
      Vertex *vertex;
      Vertex  *mate;
      Face    *dual;
  };

  FaceSequence steinerFaces, modifiedFaces;
  NodeSequence steinerNodes;

  BinaryTree *btree;

  int required_topology;
  bool verbose;
  size_t maxfaceID;

  vector<FacePair> facematching;

  void splitParent(Face *parent, Face *child1, Face *child2);
  void splitParent(BinaryNode *p, BinaryNode *c1, BinaryNode *c2);

  int refine_boundary_triangle(Face *face);

  void percolateup();

  void matchnode(BinaryNode *v);
  void matchnodes(BinaryNode *child, BinaryNode *parent);
  void matchnodes(Vertex *child, Vertex *parent);

  BinaryNode* getChildofDegreeNParent(BNodeList &levelnodes, int nd);

  BinaryNode *getNextNode(BNodeList &levelnodes);
  void prunelevel(BNodeList &levelnodes);
  void maximum_tree_matching();

  void match_tree_walk(BinaryTree *tree, BinaryNode *u);
};

bool has_same_dual(const BinaryNode *nd1, const BinaryNode *nd2);

inline 
bool already_matched(const BinaryNode *node)
{
    return node->isMatched();
}

inline
void Tri2Quads::matchnodes(Vertex *child, Vertex *parent)
{
    child->setAttribute("DualMate", parent);
    parent->setAttribute("DualMate", child);
    child->setStatus(MeshEntity::REMOVE);
    parent->setStatus(MeshEntity::REMOVE);
}
///////////////////////////////////////////////////////////////////////////////////

inline 
void Tri2Quads::matchnodes(BinaryNode *child, BinaryNode *parent)
{
    if (parent->isMatched() && !child->isMatched())
    {
        cout << "Warning: parent is already matched " << endl;
    }

    if (!child->isMatched() && !parent->isMatched())
        matchnodes(child->getDualNode(), parent->getDualNode());

    btree->unlinkNode(child);
    btree->unlinkNode(parent);
}


} // namespace Jaal


#endif
