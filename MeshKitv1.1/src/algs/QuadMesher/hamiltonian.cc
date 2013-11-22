#include "mesh.h"
#include "QuadCleanUp.h"

///////////////////////////////////////////////////////////////////////////////

Mesh* hamiltonian_quadrangulation( Mesh *orgmesh )
{
  Mesh *trimesh = hamiltonian_triangulation( orgmesh );

  int numnodes = trimesh->getSize(0);
  int numfaces = trimesh->getSize(2);

  Face *face, *newquad;
  Vertex *ot1, *ot2;

  Mesh *quadmesh = new Mesh;

  for( int i = 0; i < numnodes; i++) 
  {
       Vertex *vtx = trimesh->getNode(i);
       quadmesh->addNode(vtx);
  }

  int relexist = trimesh->build_relations(0,2);

  for( int iface = 0; iface < numfaces; iface++) 
  {
      face = trimesh->getFace(iface);
      face->setRemoveMark(0);
  }

  vector<FaceType> neighs;
  vector<NodeType> connect(4);

  int visitmark;

  for( int iface = 0; iface < numfaces; iface++) 
  {
       face = trimesh->getFace(iface);

       NodeType v0 = face->getConnection(0);
       NodeType v1 = face->getConnection(1);
       NodeType v2 = face->getConnection(2);

       visitmark = v0->isVisited() + v1->isVisited();
       if( visitmark != 2)
       {
           neighs = Mesh::getRelation112( v0, v1 );
           if( neighs.size() == 2 ) {
	       if( !neighs[0]->isRemoved() && !neighs[1]->isRemoved() ) {
                    ot1 = Face::opposite_node( neighs[0], v0, v1);
                    ot2 = Face::opposite_node( neighs[1], v0, v1);
	            newquad =  new Face;
		    connect[0] = ot1;
		    connect[1] = v0;
		    connect[2] = ot2;
		    connect[3] = v1;
		    newquad->setConnection(connect);
		    quadmesh->addFace(newquad);
		    neighs[0]->setRemoveMark(1);
		    neighs[1]->setRemoveMark(1);
               }
	   }
       }

       visitmark = v1->isVisited() + v2->isVisited();
       if( visitmark != 2)
       {
           neighs = Mesh::getRelation112( v1, v2 );
           if( neighs.size() == 2 ) {
	       if( !neighs[0]->isRemoved() && !neighs[1]->isRemoved() ) {
                    ot1 = Face::opposite_node( neighs[0], v1, v2);
                    ot2 = Face::opposite_node( neighs[1], v1, v2);
	            newquad =  new Face;
		    connect[0] = ot1;
		    connect[1] = v1;
		    connect[2] = ot2;
		    connect[3] = v2;
		    newquad->setConnection(connect);
		    quadmesh->addFace(newquad);
		    neighs[0]->setRemoveMark(1);
		    neighs[1]->setRemoveMark(1);
               }
	   }
       }

       visitmark = v2->isVisited() + v0->isVisited();
       if( visitmark != 2)
       {
           neighs = Mesh::getRelation112( v2, v0 );
           if( neighs.size() == 2 ) {
	       if( !neighs[0]->isRemoved() && !neighs[1]->isRemoved() ) {
                    ot1 = Face::opposite_node( neighs[0], v2, v0);
                    ot2 = Face::opposite_node( neighs[1], v2, v0);
	            newquad =  new Face;
		    connect[0] = ot1;
		    connect[1] = v2;
		    connect[2] = ot2;
		    connect[3] = v0;
		    newquad->setConnection(connect);
		    quadmesh->addFace(newquad);
		    neighs[0]->setRemoveMark(1);
		    neighs[1]->setRemoveMark(1);
               }
	   }
       }

  }

  for( int iface = 0; iface < numfaces; iface++) 
  {
       face = trimesh->getFace(iface);
       if( !face->isRemoved() ) quadmesh->addFace(face);
  }

  if( !relexist) trimesh->clear_relations(0,2);

  numfaces = quadmesh->getSize(2);
  int numtris = 0, numquads = 0;
  for( int i = 0; i < numfaces; i++) {
       int nnodes = quadmesh->getFace(i)->getSize(0);
       if( nnodes == 3 ) numtris++;
       if( nnodes == 4 ) numquads++;
  }
  quadmesh->enumerate(2);

  cout << "# Triangles : " << numtris  << endl;
  cout << "# Quads     : " << numquads << endl;

  QuadCleanUp quadclean(quadmesh);

  quadclean.search_diamonds();
  quadclean.search_doublets();
       
  return quadmesh;
 
}

///////////////////////////////////////////////////////////////////////////////

Mesh* hamiltonian_triangulation( Mesh *orgmesh )
{

  int numnodes = orgmesh->getSize(0);
  int numfaces = orgmesh->getSize(2);

  Mesh *newtrimesh = new Mesh;

  for( int i = 0; i < numnodes; i++) 
  {
       Vertex *vtx = orgmesh->getNode(i);
       vtx->setVisitMark(0);
       newtrimesh->addNode(vtx);
  }


  Face *face, *newtri;
  Point3D p3d;
  Vertex *dualvtx;

  int relexist = orgmesh->build_relations(0,2);

  for( int iface = 0; iface < numfaces; iface++) 
  {
      face = orgmesh->getFace(iface); assert( face );
      face->setID(iface);
      p3d = face->getCentroid();
      dualvtx = Vertex::newObject();
      dualvtx->setID( numnodes + iface );
      dualvtx->setVisitMark(1);
      dualvtx->setXYZCoords(p3d);
      face->setDualNode( dualvtx );
      newtrimesh->addNode(dualvtx);
  }


  vector<FaceType> neighs;
  vector<NodeType> connect(3);

  for( int iface = 0; iface < numfaces; iface++) 
  {
       face = orgmesh->getFace(iface);

       NodeType v0 = face->getConnection(0);
       NodeType v1 = face->getConnection(1);
       NodeType v2 = face->getConnection(2);

       neighs = Mesh::getRelation112( v0, v1 );
       if( neighs.size() == 2 ) 
       {
           if( min( neighs[0], neighs[1] ) == face ) {
	       newtri = new Face;
	       connect[0] = neighs[0]->getDualNode();
	       connect[1] = neighs[1]->getDualNode();
	       connect[2] = v0;
	       newtri->setConnection(connect);
	       newtrimesh->addFace(newtri);

	       newtri = new Face;
	       connect[0] = neighs[0]->getDualNode();
	       connect[1] = neighs[1]->getDualNode();
	       connect[2] = v1;
	       newtri->setConnection(connect);
	       newtrimesh->addFace(newtri);
           }
       } else {
	   newtri = new Face;
	   connect[0] = neighs[0]->getDualNode();
	   connect[1] = v0;
	   connect[2] = v1;
	   newtri->setConnection(connect);
	   newtrimesh->addFace(newtri);
       }


       neighs = Mesh::getRelation112( v1, v2 );
       if( neighs.size() == 2 )
       {
           if( min( neighs[0], neighs[1] ) == face ) {
	       newtri = new Face;
	       connect[0] = neighs[0]->getDualNode();
	       connect[1] = neighs[1]->getDualNode();
	       connect[2] = v1;
	       newtri->setConnection(connect);
	       newtrimesh->addFace(newtri);

	       newtri = new Face;
	       connect[0] = neighs[0]->getDualNode();
	       connect[1] = neighs[1]->getDualNode();
	       connect[2] = v2;
	       newtri->setConnection(connect);
	       newtrimesh->addFace(newtri);
           }
      } else {
	   newtri = new Face;
	   connect[0] = neighs[0]->getDualNode();
	   connect[1] = v1;
	   connect[2] = v2;
	   newtri->setConnection(connect);
	   newtrimesh->addFace(newtri);
      }

      neighs = Mesh::getRelation112( v2, v0 );
      if( neighs.size() == 2 )
       {
           if( min( neighs[0], neighs[1] ) == face ) {
	       newtri = new Face;
	       connect[0] = neighs[0]->getDualNode();
	       connect[1] = neighs[1]->getDualNode();
	       connect[2] = v2;
	       newtri->setConnection(connect);
	       newtrimesh->addFace(newtri);

	       newtri = new Face;
	       connect[0] = neighs[0]->getDualNode();
	       connect[1] = neighs[1]->getDualNode();
	       connect[2] = v0;
	       newtri->setConnection(connect);
	       newtrimesh->addFace(newtri);
           }
      } else {
	   newtri = new Face;
	   connect[0] = neighs[0]->getDualNode();
	   connect[1] = v2;
	   connect[2] = v0;
	   newtri->setConnection(connect);
	   newtrimesh->addFace(newtri);
      }
    }

    if( !relexist) orgmesh->clear_relations(0,2);

    return newtrimesh;
}
