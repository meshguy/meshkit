#include "meshkit/QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_doublet_tag(Mesh *mesh, const string &name)
{
     int relexist2 = mesh->build_relations(0, 2);

     mesh->search_boundary();

     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if( vertex->isActive() ) {
               if (QuadCleanUp::isDoublet(vertex)) {
                    vertex->setAttribute(name, 1);
               } else
                    vertex->setAttribute(name, 0);
          }
     }

     if (!relexist2)
          mesh->clear_relations(0, 2);
}

//////////////////////////////////////////////////////////////////////

bool
Doublet::isSafe() const
{
     assert(vertex);

     if (vertex->isRemoved()) return 0;

     FaceSequence apexfaces;
     vertex->getRelations( apexfaces );

     int nSize = apexfaces.size();
     for (int i = 0; i < nSize; i++) {
          Face *face = apexfaces[i];
          if (face->isRemoved()) return 0;
          if (face->isVisited()) return 0;
     }
     return 1;
}

////////////////////////////////////////////////////////////////////

void
Doublet::makeShield()
{
     FaceSequence apexfaces;
     vertex->getRelations( apexfaces );

     size_t nSize = apexfaces.size();
     for (size_t i = 0; i < nSize; i++) {
          Face *face = apexfaces[i];
          face->setVisitMark(1);
     }
}

///////////////////////////////////////////////////////////////////////////////

vector<Doublet>
QuadCleanUp::search_interior_doublets()
{
     //
     ///////////////////////////////////////////////////////////////////////////
     // An interior doublet is a vertex, which is shared by two face neighbours.
     // They are undesirables in the quadmesh as it would mean the angle is 180
     // between some adjacent edges...
     //
     ///////////////////////////////////////////////////////////////////////////

     size_t numnodes = mesh->getSize(0);

     mesh->search_boundary();

     assert(mesh->getAdjTable(0, 2));

     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setVisitMark(0);
     }

     vDoublets.clear();

     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if( !vertex->isRemoved() ) {
               if (isDoublet(vertex)) {
                    Doublet newdoublet(mesh, vertex);
                    if (newdoublet.isSafe()) {
                         newdoublet.makeShield();
                         vDoublets.push_back(newdoublet);
                    }
               }
          }
     }

     return vDoublets;
}

//////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face, Vertex *v0, Vertex *v2)
{
     //Create new vertex at the center of (v0,v2)
     Point3D p3d;
     Vertex::mid_point(v0, v2, p3d);

     Vertex *doublet = Vertex::newObject();
     doublet->setXYZCoords(p3d);

     int pos = face->getPosOf( v0 );
     if( pos < 0) return NULL;

     assert( v2 == face->getNodeAt( (pos+2) ) );

     Vertex *o1 = face->getNodeAt( pos+1 );
     Vertex *o2 = face->getNodeAt( pos+3 );

     //  Creating a doublet in the mesh changes:
     //  (1)  insert new node
     //  (2)  one old face is removed
     //  (3)  two new faces inserted.
     //

     NodeSequence connect(4);

     mesh->addNode(doublet);
     connect[0] = doublet;
     connect[1] = v0;
     connect[2] = o1;
     connect[3] = v2;

     Face *newquad1 = Face::newObject();
     newquad1->setNodes(connect);
     mesh->addFace(newquad1);

     connect[0] = doublet;
     connect[1] = v2;
     connect[2] = o2;
     connect[3] = v0;

     Face *newquad2 = Face::newObject();
     newquad2->setNodes(connect);
     mesh->addFace(newquad2);

     mesh->remove( face );

     return doublet;
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face)
{
     if (face->getSize(0) != 4) return NULL;
     Vertex *v0 = face->getNodeAt(0);
     Vertex *v2 = face->getNodeAt(2);
     return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_boundary_doublet(Face *face)
{
     ////////////////////////////////////////////////////////////////////////////
     //
     //                                  X  ( Internal Vertex)
     //                               .      .
     //	                           .           .
     //                           .               .
     //                         .                   .
     //               ********X...........X...........X************************
     //        Boundary                 Singlet                Boundary
     //
     //   There is one quad on the boundary with three nodes on the boundary and
     //   one internal nodes. In order to remove the singlet node, we artifically
     //   create one doublet between the internal node and the singlet node.
     ////////////////////////////////////////////////////////////////////////////

     if (!face->isBoundary())
          return NULL;

     if (face->getSize(0) != 4)
          return NULL;

     int ncount = 0;
     for (int i = 0; i < 4; i++) {
          Vertex *v = face->getNodeAt(i);
          ncount += v->isBoundary();
     }

     if (ncount != 3)
          return NULL;

     Vertex *v0 = NULL, *v2 = NULL;
     for (int i = 0; i < 4; i++) {
          Vertex *v = face->getNodeAt(i);
          if (!v->isBoundary()) {
               v0 = face->getNodeAt((i + 0) % 4);
               v2 = face->getNodeAt((i + 2) % 4);
          }
     }

     return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

int
Doublet::remove()
{
     if (vertex->isRemoved() || vertex->isBoundary()) return 1;

     FaceSequence neighs;
     vertex->getRelations( neighs );
     assert(neighs.size() > 0);

     if (neighs.size() != 2) return 1;

     Face *face0 = neighs[0];
     Face *face1 = neighs[1];

     assert(face0 != face1 );

     if (face0->isRemoved() || face1->isRemoved() )
          return 1;

     Vertex *d1 = NULL, *d2 = NULL, *o1 = NULL, *o2 = NULL;

     NodeSequence connect = face0->getNodes();
     if (connect.size() != 4)
          return 2;

     int nC = connect.size();
     for (int i = 0; i < nC; i++) {
          if (connect[i] == vertex) {
               d1 = connect[(i + 1) % 4];
               o1 = connect[(i + 2) % 4];
               d2 = connect[(i + 3) % 4];
               break;
          }
     }

     connect = face1->getNodes();
     if (connect.size() != 4)
          return 2;

     nC = connect.size();
     for (int i = 0; i < nC; i++) {
          if (connect[i] == vertex) {
               o2 = connect[(i + 2) % 4];
               break;
          }
     }

     assert(d1);
     assert(d2);
     assert(o1);
     assert(o2);
     assert(o1 != o2);

     // Change the connectivity of face0 and face1 is removed.
     mesh->deactivate( face0 );
     mesh->deactivate( face1 );

     connect[0] = d1;
     connect[1] = o1;
     connect[2] = d2;
     connect[3] = o2;

     // Reuse one of the face
     face0->setNodes( connect );
     mesh->reactivate( face0 );

     // Other face and doublet die.
     mesh->remove( face1  );
     mesh->remove( vertex );

     return 0;
}

///////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_doublets_once()
{
     search_interior_doublets();

     int ncount = 0;
     size_t nSize = vDoublets.size();

     for (size_t i = 0; i < nSize; i++) {
          int err = vDoublets[i].remove();
          if (!err) ncount++;
     }


     return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_interior_doublets()
{
     StopWatch swatch;
     swatch.start();

     int relexist2 = mesh->build_relations(0, 2);

     mesh->search_boundary();

     // It is possible that removal of doublet may create singlets on the boundary
     // so, it is better to call doublets first and then call to singlet removal next.
     int ncount1 = 0;
     while (1) {
          int ncount2 = remove_doublets_once();
          if (ncount2 == 0) break;
          ncount1 += ncount2;
     }

     swatch.stop();

     if (!relexist2)
          mesh->clear_relations(0, 2);

     cout << "Info: #Doublets removed : " << ncount1 << endl;
     cout << "Info:  Doublet removal execution time : " << swatch.getSeconds() << endl;

     return ncount1;
}

///////////////////////////////////////////////////////////////////////////////

