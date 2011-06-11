#include "Mesh.hpp"
#include "MeshRefine2D.hpp"

Mesh *SimpleTri2Quads( Mesh *mesh)
{
   assert( mesh->isHomogeneous() == 3 );

   int numNodes = mesh->getSize(0);
   int numFaces = mesh->getSize(2);

   vector<Edge*> edges = mesh->getEdges();
   int numEdges = edges.size();

   Mesh *quadmesh = new Mesh;

   quadmesh->reserve( numNodes + numEdges + numFaces, 0 );
   quadmesh->reserve( 3*numFaces, 2 );

   map<Edge*, Vertex*> edge_node;
   map<Vertex*,map<Vertex*,Edge*> >  node_edges;
   map<Vertex*,Vertex*>  tri_node_dup;

   size_t index = 0;
   for( size_t i = 0; i < numNodes; i++) {
        Vertex *vt = mesh->getNodeAt(i);
        Vertex *vc = new Vertex;
        vc->setID( index++);
        Point3D p3d=  vt->getXYZCoords();
        vc->setXYZCoords( p3d  );
        quadmesh->addNode ( vc );
        tri_node_dup[vt] = vc;
   }

   for( int i = 0; i < numEdges; i++) {
        Edge *edge = mesh->getEdgeAt(i);
        Vertex *v0 = edge->getNodeAt(0);
        Vertex *v1 = edge->getNodeAt(1);
        Point3D p3d = Vertex::mid_point(v0, v1);

        Vertex *vc = new Vertex;
        vc->setXYZCoords( p3d );
        vc->setID( index++);
        quadmesh->addNode ( vc );
        edge_node[edge] = vc;
        node_edges[v0][v1] = edge;
        node_edges[v1][v0] = edge;
   }

   Vertex *enodes[3];
   vector<Vertex*> qconn(4);
   for( int i = 0; i < numFaces; i++) {
        Face *tri = mesh->getFaceAt(i);
        enodes[0] = NULL;
        enodes[1] = NULL;
        enodes[2] = NULL;
        for( int j = 0; j < 3; j++) {
             Vertex *v0 =  tri->getNodeAt( j+1 );
             Vertex *v1 =  tri->getNodeAt( j+2 );
             Edge   *eg =  node_edges[v0][v1];
             enodes[j]  =  edge_node[eg];
        }

        assert( enodes[0] ) ;
        assert( enodes[1] ) ;
        assert( enodes[2] ) ;

        Vertex *vc = new Vertex;
        Point3D p3d = tri->getCentroid();
        vc->setXYZCoords( p3d );
        vc->setID( index++);
        quadmesh->addNode ( vc );

        Vertex *vt0 = tri_node_dup[tri->getNodeAt(0)]; assert( vt0 );
        Vertex *vt1 = tri_node_dup[tri->getNodeAt(1)]; assert( vt1 );
        Vertex *vt2 = tri_node_dup[tri->getNodeAt(2)]; assert( vt2 );

        Face *quad0 = new Face( vt0, enodes[2], vc, enodes[1] );
        Face *quad1 = new Face( vt1, enodes[0], vc, enodes[2] );
        Face *quad2 = new Face( vt2, enodes[1], vc, enodes[0] );

        quadmesh->addFace ( quad0 );
        quadmesh->addFace ( quad1 );
        quadmesh->addFace ( quad2 );
   }
   return quadmesh;
}

#ifdef MAIN_TEST
int main( int argc, char **argv )
{
   if( argc != 3 ) {
       cout << "Usage: executable <infile> <outfile> " << endl;
       return 1;
   }

   Mesh *trimesh =  new Mesh;
   trimesh->readFromFile( argv[1] );
   trimesh->get_topological_statistics();

   Mesh *quadmesh = SimpleTri2Quads( trimesh );

   quadmesh->get_topological_statistics();

   quadmesh->saveAs( argv[2] );

   delete trimesh;
   delete quadmesh;

   return 0;
}
#endif

