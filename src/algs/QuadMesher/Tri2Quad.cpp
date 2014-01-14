#include "meshkit/Tri2Quad.hpp"
#include "meshkit/StopWatch.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int Tri2Quads::verify(Mesh *mesh, const vector<FacePair> &matching)
{
     int relexist2 = mesh->build_relations(0, 2);

     // This is the Graph matching verification. It is useful to verify the
     // matching, if done by other algorithms.

     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face * f = mesh->getFaceAt(i);
          f->setVisitMark(0);
     }

     FaceSequence faceneighs;
     size_t nsize = matching.size();
     for (size_t i = 0; i < nsize; i++) {
          size_t f1 = matching[i].first;
          size_t f2 = matching[i].second;
          assert(f1 != f2);
          if (f1 < f2) {
               Face *f0 = mesh->getFaceAt(f1);
               Face *f1 = mesh->getFaceAt(f2);
               f0->getRelations12( faceneighs );
               if (find(faceneighs.begin(), faceneighs.end(), f1) == faceneighs.end()) {
                    cout << "Error: Face matching error: faces not neighbors  " << endl;
                    exit(0);
                    return 1;
               }
               assert(!f0->isVisited());
               assert(!f1->isVisited());
               f0->setVisitMark(1);
               f1->setVisitMark(1);
          }
     }

     for (size_t i = 0; i < numfaces; i++) {
          Face * f = mesh->getFaceAt(i);
          assert(f->isVisited());
     }

     if( !relexist2 ) mesh->clear_relations(0,2 );

     return 0;

}
///////////////////////////////////////////////////////////////////////////////

Mesh* Tri2Quads::collapse_matched_triangles(Mesh *trimesh, const vector<FacePair> &matching,
          int replace)
{
#ifdef DEBUG
     if(Tri2Quads::verify(trimesh, matching) != 0) return NULL;
     cout << " Verification Done " << endl;
#endif

     Mesh *quadmesh = new Mesh;

     size_t nsize = matching.size();

     assert( nsize );

     Face *tri1, *tri2, *quad;

     NodeSequence  tnodes;
     trimesh->getNodes( tnodes );
     quadmesh->addNodes(tnodes );
     quadmesh->reserve( nsize, 2 );

     for (size_t i = 0; i < nsize; i++) {
          size_t f1 = matching[i].first;
          size_t f2 = matching[i].second;
          tri1 = trimesh->getFaceAt(f1);
          tri2 = trimesh->getFaceAt(f2);
          quad = Face::create_quad(tri1, tri2, replace);
          quadmesh->addFace(quad);
          if (replace) delete tri2;
     }

     if (replace) trimesh->emptyAll();

     return quadmesh;
}
///////////////////////////////////////////////////////////////////////////////

int Tri2Quads::refine_boundary_triangle(Face *tri0)
{
     if (tri0->getSize(0) != 3)
          return 1;

     Vertex *bv0 = NULL;
     Vertex *bv1 = NULL;
     Vertex *bv2 = NULL;

     for (int i = 0; i < 3; i++) {
          Vertex *ev1 = tri0->getNodeAt(i + 1);
          Vertex *ev2 = tri0->getNodeAt(i + 2);
          if (ev1->isBoundary() && ev2->isBoundary()) {
               bv0 = ev1;
               bv1 = ev2;
               bv2 = tri0->getNodeAt(i);
               break;
          }
     }

     if (bv0 == NULL || bv1 == NULL)
          return 2;

     Point3D p3d;
     Vertex::mid_point(bv0, bv1, p3d);

     Vertex *bound = Vertex::newObject();
     bound->setXYZCoords(p3d);

     trimesh->addNode(bound);

     NodeSequence tconnect(3);
     tconnect[0] = bv0;
     tconnect[1] = bound;
     tconnect[2] = bv2;
     tri0->setNodes(tconnect);

     tconnect[0] = bound;
     tconnect[1] = bv1;
     tconnect[2] = bv2;
     maxfaceID++;
     Face *tri1 = Face::newObject();
     tri1->setID(maxfaceID);
     tri1->setNodes(tconnect);
     trimesh->addFace(tri1);

     steinerNodes.push_back(bound);

     FacePair facepair;
     facepair.first = tri0->getID();
     facepair.second = tri1->getID();
     facematching.push_back(facepair);

     return 0;
}

///////////////////////////////////////////////////////////////////////////////////


void Tri2Quads::splitParent(BinaryNode *parent, BinaryNode *child1,
                            BinaryNode *child2)
{
     Vertex* dnode = NULL;
     dnode = parent->getDualNode();

     Face *parentface = NULL;
     dnode->getAttribute("PrimalFace", parentface);

     dnode = child1->getDualNode();
     Face *face1 = NULL;
     dnode->getAttribute("PrimalFace", face1);

     dnode = child2->getDualNode();
     Face *face2 = NULL;
     dnode->getAttribute("PrimalFace", face2);

     NodeSequence connect(3);

     // Remove all existing vertex-face relations;
     Vertex *vertex;
     for (int i = 0; i < 3; i++) {
          vertex = parentface->getNodeAt(i);
          vertex->clearRelations(2);

          vertex = face1->getNodeAt(i);
          vertex->clearRelations(2);

          vertex = face2->getNodeAt(i);
          vertex->clearRelations(2);
     }

     // Rebuild vertex-face relations...
     for (int i = 0; i < 3; i++) {
          vertex = parentface->getNodeAt(i);
          vertex->addRelation(parentface);

          vertex = face1->getNodeAt(i);
          vertex->addRelation(face1);

          vertex = face2->getNodeAt(i);
          vertex->addRelation(face2);
     }

     dnode = NULL;
     parentface->getAttribute("DualNode", dnode);
     Vertex *steiner = dnode->getClone();
     steiner->setID(parentface->getID());
     trimesh->addNode(steiner);
     steinerNodes.push_back(steiner);

     int edge1, edge2, edge3;

     edge1 = edge2 = edge3 = -1;
     FaceSequence neighs;
     for (int i = 0; i < 3; i++) {
          Vertex *v0 = parentface->getNodeAt(i + 1);
          Vertex *v1 = parentface->getNodeAt(i + 2);
          Mesh::getRelations112(v0, v1, neighs);

          if (neighs.size() == 1)
               edge3 = i;

          if (neighs.size() == 2) {
               if (find(neighs.begin(), neighs.end(), face1) != neighs.end())
                    edge1 = i;
               if (find(neighs.begin(), neighs.end(), face2) != neighs.end())
                    edge2 = i;
          }
     }

     Face *qface;
     Vertex *ev0, *ev1;

     // Match Child1 and One of the Split Triangle ...
     maxfaceID++;
     ev0 = parentface->getNodeAt(edge1 + 1);
     ev1 = parentface->getNodeAt(edge1 + 2);
     connect[0] = steiner;
     connect[1] = ev0;
     connect[2] = ev1;
     qface = Face::newObject();
     qface->setNodes(connect);
     qface->setID(maxfaceID);
     Vertex *dc1 = DualGraph::getNewDualNode( qface );
     dc1->setID(maxfaceID);
     trimesh->addFace(qface);
     steinerFaces.push_back(qface);
     dnode = NULL;
     face1->getAttribute("DualNode", dnode);
     matchnodes( dnode, dc1);
     BinaryNode *bnode1 = new BinaryNode(dc1);
     bnode1->setMatchMark(1);
     bnode1->setParent(parent);
     bnode1->addChild(child1);
     btree->addNode(bnode1);

     // Match Child2 and One of the Split Triangle ...
     maxfaceID++;
     ev0 = parentface->getNodeAt(edge2 + 1);
     ev1 = parentface->getNodeAt(edge2 + 2);
     connect[0] = steiner;
     connect[1] = ev0;
     connect[2] = ev1;
     qface = Face::newObject();
     qface->setID(maxfaceID);
     qface->setNodes(connect);
     Vertex *dc2 = DualGraph::getNewDualNode( qface );
     dc2->setID(maxfaceID);
     trimesh->addFace(qface);
     steinerFaces.push_back(qface);
     dnode = NULL;
     face2->getAttribute( "DualNode", dnode);
     matchnodes( dnode, dc2);

     BinaryNode *bnode2 = new BinaryNode(dc2);
     bnode2->setMatchMark(1);
     bnode2->setParent(parent);
     bnode2->addChild(child2);
     btree->addNode(bnode2);

     // Now Parent have different connectivity ...
     ev0 = parentface->getNodeAt(edge3 + 1);
     ev1 = parentface->getNodeAt(edge3 + 2);
     connect[0] = steiner;
     connect[1] = ev0;
     connect[2] = ev1;
     parentface->setNodes(connect);
     Point3D p3d;
     parentface->getAvgPos( p3d );

     Vertex *dc3 = NULL;
     parentface->getAttribute("DualNode", dc3);
     dc3->setXYZCoords(p3d);
     parent->addChild(bnode1);
     parent->addChild(bnode2);
     modifiedFaces.push_back(parentface);

     for (int i = 0; i < 3; i++) {
          vertex = parentface->getNodeAt(i);
          vertex->clearRelations(2);

          vertex = face1->getNodeAt(i);
          vertex->clearRelations(2);

          vertex = face2->getNodeAt(i);
          vertex->clearRelations(2);
     }

     child1->setMatchMark(1);
     child2->setMatchMark(1);

     btree->removeNode(child1);
     btree->removeNode(child2);
}

////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::matchnode(BinaryNode* v)
{
     BinaryNode *parv = v->getParent();

     if (parv == NULL)
          return;

     int degree = parv->getDegree();

     if (parv->isRoot() && degree == 2) {
          BinaryNode *vsib = v->getSibling();
          splitParent(parv, v, vsib);
          return;
     }

     if (degree == 1 || degree == 2) {
          matchnodes(v, parv);
          return;
     }

     if (degree == 3) {

          if (required_topology == ALL_QUADS) {
               BinaryNode *vsib = v->getSibling();
               splitParent(parv, v, vsib);
               return;
          }

          BinaryNode *vsib = v->getSibling();
          Vertex *d0 = v->getDualNode();
          Vertex *d1 = vsib->getDualNode();
          if (d0->getNumRelations(0) < d1->getNumRelations(0)) {
               matchnodes(v, parv);
               btree->unlinkNode(vsib);
          } else {
               matchnodes(vsib, parv);
               btree->unlinkNode(v);
          }
     }
}

///////////////////////////////////////////////////////////////////////////////

BinaryNode* Tri2Quads::getChildofDegreeNParent(BNodeList &levelnodes,
          int nd)
{
     BinaryNode *currnode, *parent, *child;

     int ncount;
     BNodeList::const_iterator it;

     for (it = levelnodes.begin(); it != levelnodes.end(); ++it) {
          currnode = *it;
          parent = currnode->getParent();
          if (parent) {
               if (!parent->isMatched()) {
                    ncount = 0;
                    if (parent->getParent())
                         ncount = 1;
                    for (int i = 0; i < parent->getNumChildren(); i++) {
                         child = parent->getChild(i);
                         if (!child->isMatched())
                              ncount++;
                    }
                    if (ncount == nd)
                         return currnode;
               }
          }
     }

     return NULL;
}

///////////////////////////////////////////////////////////////////////////////

BinaryNode *Tri2Quads::getNextNode(BNodeList &levelnodes)
{
     BinaryNode *currnode = NULL;

     if (levelnodes.empty())
          return currnode;

     BNodeList::iterator it;
     for (it = levelnodes.begin(); it != levelnodes.end(); ++it) {
          currnode = *it;
          if (currnode->isMatched())
               btree->unlinkNode(currnode);
     }

     it = remove_if(levelnodes.begin(), levelnodes.end(), already_matched);
     levelnodes.erase(it, levelnodes.end());

     BinaryNode *child = NULL;

     // High Priority: parent having degree = 1;
     child = getChildofDegreeNParent(levelnodes, 1);

     if (!child)
          child = getChildofDegreeNParent(levelnodes, 2);

     // Low Priority: parent having degree = 3;
     if (!child)
          child = getChildofDegreeNParent(levelnodes, 3);

     return child;
}

////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::prunelevel(BNodeList &levelnodes)
{
     while (1) {
          BinaryNode *currnode = getNextNode(levelnodes);
          if (currnode == NULL) break;
          matchnode(currnode);
     }
}

////////////////////////////////////////////////////////////////////////////////

void Tri2Quads::percolateup()
{
     steinerNodes.clear();
     steinerFaces.clear();

     int height = btree->getHeight();
     BNodeList levelnodes;
     BNodeList::const_iterator it;

     //Reset all the Matching marks to 0;
     for (int i = 0; i < height; i++) {
          levelnodes = btree->getLevelNodes(height - i - 1);
          BinaryNode *currnode;
          for (it = levelnodes.begin(); it != levelnodes.end(); ++it) {
               currnode = *it;
               currnode->setMatchMark(0);
          }
     }

     // Start Prunning the level. At most the root will be unmatched.
     for (int i = 0; i < height; i++) {
          levelnodes = btree->getLevelNodes(height - i - 1);
          prunelevel(levelnodes);
     }

     size_t numfaces = trimesh->getSize(2);
     facematching.reserve(numfaces);

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = trimesh->getFaceAt(i);
          Vertex *u = NULL;
          face->getAttribute("DualNode", u);
           
          assert(u);
          Vertex *v = NULL;
          u->getAttribute("DualMate", v);
          if (v) {
               if (v->getID() > u->getID()) {
                    FacePair facepair;
                    facepair.first = u->getID();
                    facepair.second = v->getID();
                    facematching.push_back(facepair);
               }
          }
     }

     // If the root is unmatched, bring it down to a leaf and then split the
     // leaf. Do this step after the triangles have been matched.
     BinaryNode *root = btree->getRoot();
     if (!root->isMatched()) {
#ifdef VERBOSE
          cout << "Warning: Boundary Triangle modified " << endl;
#endif
          Vertex *dnode = root->getDualNode();
          Face *rootface = NULL;
          dnode->getAttribute("PrimalFace", rootface);
          refine_boundary_triangle(rootface);
     }

}

///////////////////////////////////////////////////////////////////////////////

void Tri2Quads::maximum_tree_matching()
{
     // In order to insert any steiner point on the boundary triangle (at the root)
     // We should know which triangles and nodes are on the boundary. Therefore,
     // call this function to set the boundary flags. Building the relationship
     // at this stage is good as even the DualGraph construction require it.

     trimesh->build_relations(0, 2);
     trimesh->search_boundary();

#ifdef VERBOSE
     cout << " Creating Dual Graph ... " << endl;
#endif

     DualGraph *dgraph = new DualGraph;
     dgraph->build(trimesh);

     if (verbose)
          cout << " Building Binary Tree of Dual Graph ... " << endl;

     btree = new BinaryTree(dgraph);
     btree->build();
     btree->saveAs("btree");

#ifdef VERBOSE
     cout << " Tree Matching ... " << endl;
#endif

     percolateup();

     // Percolateup will remove relations, so it is necessary to clear all the
     // relations.
     trimesh->clear_relations(0, 2);

     btree->deleteAll();
     dgraph->deleteAll();

     delete btree;
     delete dgraph;
}

///////////////////////////////////////////////////////////////////////////////

Mesh* Tri2Quads::getQuadMesh(Mesh *inmesh, int replace, int topo)
{

#ifdef DEBUG
     if (inmesh->isHomogeneous() != 3) {
          cout << "Warning: Input mesh is not triangular " << endl;
          return NULL;
     }

     if (!inmesh->isSimple()) {
          cout << "Warning: Input mesh is not simple, use edmonds' algorithm " << endl;
          return NULL;
     }

     if (inmesh->getNumOfComponents() > 1) {
          cout << "Warning: There are multiple components in the mesh" << endl;
          cout << "         Algorithm works for single component " << endl;
          return NULL;
     }

     int euler0 = trimesh->getEulerCharacteristic();
     cout << " Input Euler # : " << euler0 << endl;
     cout << inmesh->saveAs( "Check.dat");
#endif

     trimesh = inmesh;

     required_topology = topo;

     trimesh->enumerate(2);
     maxfaceID = trimesh->getSize(2) - 1;

     ///////////////////////////////////////////////////////////////////////////
     // Generate Maximum Matching on a binary tree using Suneeta's Algorithm.
     // If the required topology is set to ALL_QUADS, steiner points( and new
     // faces) will be inserted in the input triangle mesh.  Please note that
     // this implementation doesn't produces "The" optimal soluation as
     // described in the original papers, and doesn't even guarantee that
     // the resulting quadrilaterals will be convex. This along with other
     // topological and geometric optimization are anyhow essential and
     // are carried out during the "Post Processing" step. Therefore, we
     // have sacrifised performance over quality in this implementation.
     // Roughly we can expect that about 4-5% steiner points are inserted in
     // most of the general cases.
     ///////////////////////////////////////////////////////////////////////////
     StopWatch swatch;
     swatch.start();

     maximum_tree_matching();
     Mesh *quadmesh = collapse_matched_triangles(trimesh, facematching, replace);
     swatch.stop();
     cout << "Info: Tri->Quad Elapsed Time " << swatch.getSeconds() << endl;

     if( quadmesh ) {
          if (!quadmesh->isSimple())
               cout << "Warning: Quadrilateral Mesh is not simple " << endl;

          quadmesh->enumerate(2);
          if (!quadmesh->is_consistently_oriented()) {
               cout << "Warning:Trying to make Quadrilateal Mesh consistently oriented " << endl;
               quadmesh->make_consistently_oriented();
               if (quadmesh->is_consistently_oriented())
                    cout << "Info: Quadrilateral Mesh is now consistently oriented: Very good " << endl;
               else
                    cout << "Alas ! Quadrilateral Mesh is still inconsistently oriented: Check manually " << endl;
          }

          quadmesh->enumerate(0);
          quadmesh->enumerate(2);
     }

     //////////////////////////////////////////////////////////////////////////
     // Since Steiner points may be inserted in the mesh ( and new triangles).
     // Renumber all the nodes and faces for future processing.
     //////////////////////////////////////////////////////////////////////////

     trimesh->enumerate(0);
     trimesh->enumerate(2);

     return quadmesh;
}

///////////////////////////////////////////////////////////////////////////////

void Tri2Quads::match_tree_walk(BinaryTree *btree, BinaryNode *parent)
{
     //
     // Brings all the internal unmatched nodes at the leaf.
     //
     if (parent == NULL)
          return;
     if (parent->isLeaf())
          return;

     int numChildren = parent->getNumChildren();

     for (int i = 0; i < numChildren; i++) {
          BinaryNode *child1 = parent->getChild(i);
          if (!btree->isMatched(parent, child1)) {
               int numGrandChildren = child1->getNumChildren();
               for (int j = 0; j < numGrandChildren; j++) {
                    BinaryNode *child2 = child1->getChild(j);
                    if (btree->isMatched(child1, child2)) {
                         Vertex *np = parent->getDualNode();
                         assert(np);
                         Vertex *c1 = child1->getDualNode();
                         assert(c1);
                         Vertex *c2 = child2->getDualNode();
                         assert(c2);
                         matchnodes(np, c1);
                         c2->setAttribute("DualMate", 0);
                         c2->setStatus(MeshEntity::ACTIVE);
                         match_tree_walk(btree, child2);
                         return;
                    }
               }
          }
     }

}

///////////////////////////////////////////////////////////////////////////////
