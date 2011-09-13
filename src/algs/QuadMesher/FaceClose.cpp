#include "QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_diamond_tag(Mesh *mesh, const string &aname)
{
     QuadCleanUp qClean(mesh);
     vector<Diamond> diamonds = qClean.search_diamonds();

     cout << "#Dianonds " << diamonds.size() << endl;

     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if( face->isActive() )
               face->setAttribute( aname, 0);
     }

     set<Face*>::const_iterator it;
     for (size_t i = 0; i < diamonds.size(); i++) {
          Face *face = diamonds[i].face;
          if( face->isActive() )
               face->setAttribute( aname, 1);
     }
}

////////////////////////////////////////////////////////////////////
bool
FaceClose::isSafe() const
{
     if (!face->isActive()) return 0;

     if (!vertex0->isActive()) return 0;
     if (!vertex2->isActive()) return 0;

     if (vertex0->isBoundary()) return 0;
     if (vertex2->isBoundary()) return 0;

     int vpos = face->getPosOf(vertex0);
     assert(vpos >= 0);

     if (face->getNodeAt(vpos + 2) != vertex2) {
          cout << "Warning: Face-open requires opposite vertices " << endl;
          cout << "Debug  : Face is : " << face->getNodeAt(0)->getID() << " "
               << face->getNodeAt(1)->getID() << " "
               << face->getNodeAt(2)->getID() << " "
               << face->getNodeAt(3)->getID() << endl;
          cout << "Opposite ends are " << vertex0->getID() << " " << vertex2->getID() << endl;
          return 1;
     }

     assert(face->hasNode(vertex0));
     assert(face->hasNode(vertex2));

     FaceSequence neighs;
     vertex0->getRelations( neighs );
     int nSize = neighs.size();
     for (int  j = 0; j < nSize; j++) {
          if (neighs[j] != face) {
               int val0 = neighs[j]->hasNode(vertex0);
               int val1 = neighs[j]->hasNode(vertex2);
               if (val0 + val1 == 2) return 0;
          }
     }

     vertex2->getRelations( neighs );
     nSize = neighs.size();
     for (int j = 0; j < nSize; j++) {
          if (neighs[j] != face) {
               int val0 = neighs[j]->hasNode(vertex0);
               int val1 = neighs[j]->hasNode(vertex2);
               if (val0 + val1 == 2) return 0;
          }
     }

     return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
FaceClose::build()
{
     if (!isSafe()) return 1;

     assert(face->getSize(0) == 4);

     Vertex *v0 = vertex0;
     Vertex *v2 = vertex2;

     Point3D p0 = v0->getXYZCoords();
     Point3D p2 = v2->getXYZCoords();

     Point3D p3d;
     Vertex::mid_point(v0, v2, p3d);

     // Both nodes will merge at the same point...
     vertex0->setXYZCoords(p3d);
     vertex2->setXYZCoords(p3d);

     int pass = 1;

     FaceSequence vfaces;
     if( pass ) {
          vertex0->getRelations( vfaces );
          for( size_t i = 0; i < vfaces.size(); i++) {
               if( vfaces[i] !=  face ) {
                    if (vfaces[i]->concaveAt() >= 0) {
                         pass = 0;
                         break;
                    }
               }
          }
     }

     if( pass ) {
          vertex2->getRelations( vfaces );
          for( size_t i = 0; i < vfaces.size(); i++) {
               if( vfaces[i] !=  face ) {
                    if (vfaces[i]->concaveAt() >= 0) {
                         pass = 0;
                         break;
                    }
               }
          }
     }

     if( !pass ) {
          vertex0->setXYZCoords( p0 );
          vertex2->setXYZCoords( p2 );
          return 1;
     }

     replacedNode = Vertex::newObject();
     replacedNode->setXYZCoords(p3d);

     return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
FaceClose::commit()
{
     assert(mesh->getAdjTable(0, 2));

     if (replacedNode == NULL) return 1;

     mesh->addNode(replacedNode);

     NodeSequence fnodes;
     FaceSequence vfaces;
     vertex0->getRelations( vfaces );
     size_t nSize = vfaces.size();

     for (size_t i = 0; i < nSize; i++) {
          if (vfaces[i] != face) {
               fnodes = vfaces[i]->getNodes();
               int nv = fnodes.size();
               for( int j = 0; j < nv; j++)  {
                    if( fnodes[j] == vertex0) {
                         fnodes[j] = replacedNode;
                    }
               }
               mesh->deactivate(vfaces[i]);  // remove all existing relationships.
               vfaces[i]->setNodes( fnodes );
               mesh->reactivate(vfaces[i]);  // rebuilds new relationships.
          }
     }

     vertex2->getRelations( vfaces );
     nSize = vfaces.size();
     for (size_t i = 0; i < vfaces.size(); i++) {
          if (vfaces[i] != face) {
               fnodes = vfaces[i]->getNodes();
               int nv = fnodes.size();
               for( int j = 0; j < nv; j++)   {
                    if( fnodes[j] == vertex2) {
                         fnodes[j] = replacedNode;
                    }
               }
               mesh->deactivate(vfaces[i]); // remove all existing relatioships
               vfaces[i]->setNodes( fnodes );
               mesh->reactivate(vfaces[i]); // rebuilds new relationships..
          }
     }

     // Two nodes and face go away from the mesh..
     mesh->remove(face);
     mesh->remove(vertex0);
     mesh->remove(vertex2);

     replacedNode = NULL; // So that destructor can delete if this is not used.

     return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool
QuadCleanUp::isDiamond(Face *face, int &pos, int type)
{
     pos = -1;
     const Vertex *v0 = face->getNodeAt(0);
     const Vertex *v1 = face->getNodeAt(1);
     const Vertex *v2 = face->getNodeAt(2);
     const Vertex *v3 = face->getNodeAt(3);

     int d0 = v0->getNumRelations(2);
     int d1 = v1->getNumRelations(2);
     int d2 = v2->getNumRelations(2);
     int d3 = v3->getNumRelations(2);

     // Boundary Cases ...
     if (v0->isBoundary() || v2->isBoundary()) {
          if( d0 <= v0->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;
          if( d2 <= v2->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;
          if (!v1->isBoundary() && !v3->isBoundary()) {
               if (d1 == 3 && d3 == 3) {
                    pos = 1;
                    return 1;
               }
          }
     }

     if (v1->isBoundary() || v3->isBoundary()) {
          if( d1 <= v1->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;
          if( d3 <= v3->get_ideal_face_degree( Face::QUADRILATERAL ) ) return 0;

          if (!v0->isBoundary() && !v2->isBoundary()) {
               if ((d0 == 3 && d2 == 3)) {
                    pos = 0;
                    return 1;
               }
          }
     }

     if (v0->isBoundary()) return 0;
     if (v1->isBoundary()) return 0;
     if (v2->isBoundary()) return 0;
     if (v3->isBoundary()) return 0;

     if( type == 33 ) {
          if ((d0 == 3 && d2 == 3)) {
               if( d1 < 4 || d3 < 4)  return 0;
               pos = 0;
               return 1;
          }

          if (d1 == 3 && d3 == 3) {
               if( d0 < 4 || d2 < 4)  return 0;
               pos = 1;
               return 1;
          }
     }

     if( type == 34 ) {
          if ( (d0 == 3 && d2 == 4) || ( d0 == 4 && d2 == 3)) {
               if( d1 < 4 || d3 < 4)  return 0;
               pos = 0;
               return 1;
          }

          if ( (d1 == 4 && d3 == 3) || (d1 == 3 && d3 == 4)) {
               if( d0 < 4 || d2 < 4)  return 0;
               pos = 1;
               return 1;
          }
     }

     if( type == 55 ) {
          if (d1 == 5 && d3 == 5) {
               if ((d0 == 3 && d2 == 4) || (d0 == 4 && d2 == 3)) {
                    pos = 0;
                    return 1;
               }
          }

          if (d0 == 5 && d2 == 5) {
               if ((d1 == 3 && d3 == 4) || (d1 == 4 && d3 == 3)) {
                    pos = 1;
                    return 1;
               }
          }
     }

     return 0;
}

////////////////////////////////////////////////////////////////////////////////

vector<Diamond>
QuadCleanUp::search_diamonds(int type)
{
     if (!mesh->isPruned()) mesh->prune();

     size_t numfaces = mesh->getSize(2);

     mesh->search_boundary();

     int relexist = mesh->build_relations(0, 2);

     int pos;

     assert(mesh->isBoundaryKnown());

     vDiamonds.clear();
     for (size_t iface = 0; iface < numfaces; iface++) {
          Face *face = mesh->getFaceAt(iface);
          if (isDiamond(face, pos, type)) {
               Diamond diamond(mesh, face, pos);
               vDiamonds.push_back(diamond);
          }
     }

     if (!relexist)
          mesh->clear_relations(0, 2);

     if (vDiamonds.size())
          cout << "# of Diamonds : " << vDiamonds.size() << endl;

     return vDiamonds;
}

////////////////////////////////////////////////////////////////////

int
Diamond::build()
{
     faceclose = new FaceClose(mesh, face, vertex0, vertex2);
     return faceclose->build();
}

///////////////////////////////////////////////////////////////////////////////

int
Diamond::commit()
{
     if (faceclose) {
          int err = faceclose->commit();
          delete faceclose;
          faceclose = NULL;
          return err;
     }

     return 1;
}
///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds_once()
{
     // Essential step because only the faces are updated.
     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setVisitMark(0);
     }

     int err, pos, ncount = 0;

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if(face->isActive() ) {
               if (isDiamond(face, pos)) {
                    Diamond diamond(mesh, face, pos);
                    err = diamond.build();
                    if (!err) {
                         err = diamond.commit();
                         if (!err) ncount++;
                    }
               }
          }
     }

     return ncount;
}

///////////////////////////////////////////////////////////////////////////////
int
QuadCleanUp::degree_5_dominated()
{
     cout << "Info: Removing type 34 Diamonds ... " << endl;
     int rel0exist = mesh->build_relations(0, 0); // Need for laplace smoothing
     int rel2exist = mesh->build_relations(0, 2);

     mesh->search_boundary();

     // Essential step because only the faces are updated.
     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setVisitMark(0);
     }

     int err, pos, nfound = 0, ncommit = 0;

     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          if( face->isActive() ) {
               if (isDiamond(face, pos, 34)) {
                    nfound++;
                    Diamond diamond(mesh, face, pos);
                    err = diamond.build();
                    if (!err) {
                         err = diamond.commit();
                         if (!err) ncommit++;
                    }
               }
          }
     }
     cout << "#Diamonds found " << nfound << " and removed : " << ncommit << endl;

     if (!rel0exist) mesh->clear_relations(0, 0);
     if (!rel2exist) mesh->clear_relations(0, 2);

     return ncommit;
}

///////////////////////////////////////////////////////////////////////////////


int
QuadCleanUp::remove_diamonds_in_layer(int layerid)
{

#ifdef CHANGE
     // Essential step because only the faces are updated only.
     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setVisitMark(0);
     }

     int err, pos, ncount = 0;

     size_t faceid = 0;
     while (1) {
          if (faceid >= mesh->getSize(2)) break;
          Face *face = mesh->getFaceAt(faceid);
          if (face->getLayerID() == layerid) {
               if (isDiamond(face, pos)) {
                    Diamond diamond(mesh, face, pos);
                    err = diamond.build();
                    if (!err) {
                         err = diamond.commit();
                         if (!err) ncount++;
                    }
               }
          }
          faceid++;
     }

     vDiamonds.clear();

     return ncount;
#endif
     return 0;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds()
{
     int rel0exist = mesh->build_relations(0, 0); // Need for laplace smoothing
     int rel2exist = mesh->build_relations(0, 2);

     mesh->search_boundary();

     int ncount1 = 0;
     while (1) {
          int ncount = remove_diamonds_once();
          if (ncount == 0) break;
          ncount1 += ncount;
     }

     if (!rel0exist) mesh->clear_relations(0, 0);
     if (!rel2exist) mesh->clear_relations(0, 2);

     cout << "Info: Diamonds removed ... " << ncount1 << endl;
   


     return ncount1;
}

////////////////////////////////////////////////////////////////////

