#include "meshkit/QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_singlet_tag(Mesh *mesh, const string &aname )
{
     size_t numnodes = mesh->getSize(0);

     int relexist = mesh->build_relations(0, 2);

     if (!mesh->isBoundaryKnown())
        mesh->search_boundary();

     size_t ncount = 0;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if( vertex->isActive() ) {
               if (QuadCleanUp::isSinglet(vertex)) {
                    vertex->setAttribute(aname, 1);
                    ncount++;
               } else
                    vertex->setAttribute(aname, 0);
          }
     }
     cout << "#of Singlets detected: " << ncount << endl;

     if (!relexist)
          mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////

vector<Singlet>
QuadCleanUp::search_boundary_singlets()
{
     vSinglets.clear();
     //
     // A boundary singlet is a vertex which is shared by only one face.
     // They are undesirables in the quad mesh as that would mean large
     // angle on some of the edges..
     //
     // For the flat singlet ( angle closer to 180 degree ). it is easy
     // to remove the neighbouring quad from the mesh.
     //

     int relexist = mesh->build_relations(0, 2);

     mesh->search_boundary();

     assert(mesh->getAdjTable(0, 2));

     size_t numfaces = mesh->getSize(2);
     for (size_t i = 0; i < numfaces; i++) {
          Face *face = mesh->getFaceAt(i);
          face->setVisitMark(0);
     }

     size_t numnodes = mesh->getSize(0);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *v = mesh->getNodeAt(i);
          if( v->isActive() )  {
               if (isSinglet(v)) {
                    Singlet newsinglet(mesh, v);
                    vSinglets.push_back(newsinglet);
               }
          }

     }

     if (!relexist)
          mesh->clear_relations(0, 2);

     if (vSinglets.size())
          cout << "Info: Number of Singlets detected " << vSinglets.size() << endl;

     return vSinglets;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_boundary_singlets_once()
{

     if (vSinglets.empty())
          search_boundary_singlets();

     int ncount = 0;
     size_t nSize = vSinglets.size();
     for (size_t i = 0; i <  nSize; i++) {
          int err = vSinglets[i].remove();
          if (!err) {
               ncount++;
          }
     }

     vSinglets.clear();

     return ncount;
}
//////////////////////////////////////////////////////////////////////////
int Singlet::remove()
{
     FaceSequence vfaces;
     vertex->getRelations( vfaces );
     if (vfaces.size() > 1) return 1;

     return mesh->refine_quad15(vfaces[0]);
}

//////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_boundary_singlets()
{
     int relexist = mesh->build_relations(0, 2);

     int ncount = 0;
     while (1) {
          size_t nremoved = remove_boundary_singlets_once();
          if (nremoved == 0) break;
          ncount += nremoved;
     }

     if (!relexist) mesh->clear_relations(0, 2);

     return ncount;
}

////////////////////////////////////////////////////////////////////

void
Singlet::clear()
{
     size_t nSize = newNodes.size();
     for (size_t j = 0; j < nSize; j++)
          delete newNodes[j];
     newNodes.clear();

     nSize = newFaces.size();
     for (size_t j = 0; j < nSize; j++)
          delete newFaces[j];

     newFaces.clear();

}

///////////////////////////////////////////////////////////////////////////////

int
Singlet::commit()
{

     for (size_t i = 0; i < oldFaces.size(); i++) {
          if (oldFaces[i]->isRemoved()) {
               clear();
               return 1;
          }
     }

     if (!active) return 2;

     for (size_t i = 0; i < oldNodes.size(); i++)
          mesh->remove( oldNodes[i] );
     oldNodes.clear();

     for (size_t i = 0; i < oldFaces.size(); i++)
          mesh->remove(oldFaces[i] );
     oldFaces.clear();

     for (size_t i = 0; i < newNodes.size(); i++)
          mesh->addNode(newNodes[i]);
     newNodes.clear();

     for (size_t i = 0; i < newFaces.size(); i++)
          mesh->addFace(newFaces[i]);
     newFaces.clear();

     return 0;

}

///////////////////////////////////////////////////////////////////////////////

