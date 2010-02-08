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

#include "Mesh.h"
#include "DualGraph.h"
#include "BinaryTree.h"

BEGIN_JAAL_NAMESPACE

class Tri2Quads
{
public:

  const static int ALL_QUADS = 0;
  const static int QUAD_DOMINATED = 1;

  static int collapse_matched_triangles(Mesh *mesh, vector<FacePair> &matching,
      int replace = 0, Mesh *qmesh = NULL);

  Tri2Quads()
  {
    trimesh = NULL;
    btree = NULL;
    verbose = 1;
    required_topology = ALL_QUADS;
  }

  vector<FacePair> getMaximumDualMatching();

#ifdef USE_MOAB
  int getQuadMesh( iMesh_Instance tmesh, iBase_EntitySetHandle eset = 0,
      bool replace = 0,
      iBase_EntitySetHandle qmesh = 0, int topo = ALL_QUADS );
#endif

  int getQuadMesh(Mesh *tmesh, bool replace = 0, Mesh *qmesh = NULL, int topo =
      ALL_QUADS);

  vector<Vertex*> getSteinerNodes() const;
  vector<Face*> getInsertedFaces() const;
  vector<Face*> getModifiedFaces() const;

private:
  Mesh *trimesh; // Input mesh.
  Mesh *quadmesh;

  vector<Face*> steinerFaces, modifiedFaces;
  vector<Vertex*> steinerNodes;

  BinaryTree *btree;

  int required_topology;
  bool verbose;
  size_t maxfaceID;

  vector<FacePair> facematching;

  void splitParent(Face *parent, Face *child1, Face *child2);
  void splitParent(BinaryNode *p, BinaryNode *c1, BinaryNode *c2);

  int quadrangulate_boundary_triangle(Face *face);

  void percolateup();

  void matchnode(BinaryNode *v);
  void matchnodes(BinaryNode *child, BinaryNode *parent);
  void matchnodes(Vertex *child, Vertex *parent);

  BinaryNode* getChildofDegreeNParent(list<BinaryNode*> &levelnodes, int nd);

  BinaryNode *getNextNode(list<BinaryNode*> &levelnodes);
  void prunelevel(list<BinaryNode*> &levelnodes);
  void maximum_tree_matching();

  void match_tree_walk(BinaryTree *tree, BinaryNode *u);
};

bool has_same_dual(const BinaryNode *nd1, const BinaryNode *nd2);
bool already_matched(const BinaryNode *node);

END_JAAL_NAMESPACE

#endif
