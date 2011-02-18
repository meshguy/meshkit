#include <meshkit/Mesh.h>
#include <meshkit/MeshRefine2D.h>

#include <meshkit/QuadCleanUp.h>
#include <meshkit/DijkstraShortestPath.h>
#include <meshkit/ObjectPool.h>

extern int QuadPatches(Jaal::Mesh *mesh);

using namespace Jaal;


void usage()
{
    cout << "Usage: Executable -i in_meshfile -o out_meshfile -c cleanOp " << endl;

    cout << " *****************************************************************" << endl;
    cout << " Option :   Mesh Cleanup Operation " << endl;
    cout << " *****************************************************************" << endl;
    cout << " 0      :   Report Mesh Quality " << endl;
    cout << " 1      :   Remove interior doublets  " << endl;
    cout << " 2      :   Remove boundary singlets  " << endl;
    cout << " 3      :   Remove diamonds " << endl;
    cout << " 4      :   Vertex Degree Reduction " << endl;
    cout << " 5      :   Laplace Smoothing (No Weight)    " << endl;
    cout << " 6      :   Laplace Smoothing (Area Weight)  " << endl;
    cout << " 7      :   Laplace Smoothing (Edge Length Weight)  " << endl;
    cout << " 8      :   Advancing front Edge Swapping " << endl;
    cout << " 9      :   Shape Optimization  " << endl;
    cout << " 10     :   Reverse Elements Connectivity  " << endl;
    cout << " 11     :   Refine QuadMesh ( Scheme 14 )  " << endl;
    cout << " 12     :   Refine QuadMesh ( Scheme 15 )  " << endl;
    cout << " 13     :   Swap Concave Faces  " << endl;
    cout << " 14     :   Refine Degree 3 Faces  " << endl;
    cout << " 15     :   Search Structured Submesh " << endl;
    cout << " 16     :   Regularization with remeshing " << endl;
    cout << " 17     :   Generate Quad-Irregular(Motorcycle) Graph " << endl;
    cout << " 18     :   Everything automatic " << endl;
}

////////////////////////////////////////////////////////////////////////////////
void angle_tests()
{
#ifdef USE_VERDICT
   double xyz[4][3];

   xyz[0][0] = -1.0;
   xyz[0][1] =  0.0;
   xyz[0][2] =  0.0;

   xyz[1][0] =  0.0;
   xyz[1][1] =  0.0;
   xyz[1][2] =  0.0;

   xyz[2][0] =  1.0;
   xyz[2][1] =  0.0;
   xyz[2][2] =  0.0;

   for( int i = 0; i < 100; i++) {
     xyz[3][0] =  0.0;
     xyz[3][1] =  i;
     xyz[3][2] =  0.0;
     cout <<  v_quad_minimum_angle(4, xyz)  << endl;
     Break();
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////

int Jaal :: mesh_unit_tests()
{
  Jaal::quad_concave_tests();

  double origin[] = {0.0, 0.0, 0.0};
  double length[] = {1.0, 1.0, 0.0};
  int  grid_dim[] = {5, 5, 1};

  Mesh *mesh = Jaal::create_structured_mesh(origin, length, grid_dim, 2);

  assert( mesh->getSize(0) == 25 );
  assert( mesh->getSize(2) == 16 );

  assert( !mesh->getAdjTable(0,0) );
  assert( !mesh->getAdjTable(0,1) );
  assert( !mesh->getAdjTable(0,2) );
  assert( !mesh->getAdjTable(0,3) );

  assert( !mesh->getAdjTable(1,0) );
  assert( !mesh->getAdjTable(1,1) );
  assert( !mesh->getAdjTable(1,2) );
  assert( !mesh->getAdjTable(1,3) );

  assert(  mesh->getAdjTable(2,0) );  // Only this is active at this time..
  assert( !mesh->getAdjTable(2,1) );
  assert( !mesh->getAdjTable(2,2) );
  assert( !mesh->getAdjTable(2,3) );

  assert( !mesh->getAdjTable(3,0) );
  assert( !mesh->getAdjTable(3,1) );
  assert( !mesh->getAdjTable(3,2) );
  assert( !mesh->getAdjTable(3,3) );

  assert( mesh->getAdjTable(1,0) == 0 ); 
  EdgeSequence edges = mesh->getEdges(); 
  assert( edges.size() == 40 );
  assert( mesh->getAdjTable(1,0) == 1 );  // Because you didn't make edges persistent

  edges = mesh->getEdges(); 
  assert( mesh->getAdjTable(1,0) == 0 );  // Now the edges are not persistent

  Face *face = mesh->getFaceAt(0);
  Vertex *vertex;
  vertex = mesh->getNodeAt(0);

  mesh->build_relations(0,2);
 
  mesh->remove( face );
//  mesh->remove( vertex );
  mesh->prune();

  cout << mesh->getSize(0) << endl;
  exit(0);

}
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
  Jaal::Mesh *mesh = new Jaal::Mesh;
  Jaal::MeshOptimization mopt;

/*
  double origin[]   = { 0.0, 0.0, 0.0};
  double length[]   = { 1.0, 1.0, 1.0};
  int    gridim[]   = { 5, 5, 5};
  mesh = Jaal::create_structured_mesh(origin, length, gridim, 2 );
*/

  string infilename, outfilename;

  int iopt, numiters = 100, topo_improve = 1;
  int cleanup_op = 0;

  while( (iopt = getopt( argc, argv, "hgtc:i:o:") ) != -1)
  {
    switch(iopt)
    {
      case 'c':
             cleanup_op = atoi( optarg );
             break;
      case 'h':
             usage();
             break;
      case 'i':
            infilename = optarg;   // Input QuadMesh File
	    break;
      case 'o':
            outfilename = optarg; // Output QuadMesh File 
	    break;
      case 'l':
            numiters  = atoi( optarg ); // Number of iterations for Laplacian Smoothing.
	    break;
      case 't':
            topo_improve  = atoi( optarg ); // Should we do topological Improvement: Default( Yes );
	    break;
      default:
            cout << "Usage: Executable t [0 1]  -l #num  -i in_meshfile -o out_meshfile -c cleanOp " << endl;
	    break;
     }
  }

  if( infilename.empty() ) {
      cout <<"Warning: No input file specified " << endl;
      usage();
      return 1 ;
  }

  if( outfilename.empty() ) {
      cout <<"Warning: No output file specified " << endl;
      usage();
      return 2;
  }

  mesh->readFromFile( infilename );
  exit(0);

  size_t ninvert  =  mesh->count_inverted_faces();
  size_t numfaces =  mesh->getSize(2);

  if( ninvert > 0.5*numfaces ) 
       mesh->reverse();

  LaplaceSmoothing lapsmooth(mesh);
  LaplaceWeight *lapweight = NULL;

  QuadCleanUp qClean(mesh);

  mesh->get_topological_statistics();

  vector<QTrack>  qpath;

  switch( cleanup_op) 
  {
    case 0:
         qClean.report();
         break;
    case 1:
         qClean.remove_interior_doublets();
         break;
    case 2:
         qClean.remove_boundary_singlets();
         break;
    case 3:
         qClean.remove_diamonds();
         break;
    case 4:
         qClean.vertex_degree_reduction();
         break;
    case 5:
         lapsmooth.setMethod(0);
         lapweight = new NoWeight();
         lapsmooth.setWeight(lapweight);
         lapsmooth.execute();
         break;
    case 6:
         lapsmooth.setMethod(0);
         lapweight = new LaplaceAreaWeight();
         lapsmooth.setWeight(lapweight);
         lapsmooth.setNumIterations(100);
         lapsmooth.execute();
         break;
    case 7:
         lapsmooth.setMethod(0);
         lapweight = new LaplaceLengthWeight();
         lapsmooth.setWeight(lapweight);
         lapsmooth.setNumIterations(100);
         lapsmooth.execute();
         break;
    case 8:
         qClean.advancing_front_edges_swap();
         break;
    case 9:
         cout << " Hello " << endl;
         mopt.shape_optimize( mesh );
         break;
    case 10:
         mesh->reverse();
         break;
    case 11:
         mesh->refine_quads14();
         break;
    case 12:
         mesh->refine_quads15();
         break;
    case 13:
         qClean.swap_concave_faces();
         break;
    case 14:
         qClean.refine_degree3_faces();
         break;
    case 15:
         mesh->search_quad_patches();
         break;
    case 16:
         qClean.remesh_defective_patches();
         break;
    case 17:
         qpath = Jaal::generate_quad_irregular_graph(mesh);
         Jaal::set_irregular_path_tag(mesh, qpath);
         break;
    case 18:
         qClean.automatic();
         break;
    }


  if( cleanup_op ) {
     mesh->get_topological_statistics();
     double minarea, maxarea;
     mesh->getMinMaxFaceArea( minarea, maxarea );

     if( minarea < 0.0 )
         cout << "Warning: There are negative area faces in the mesh " << endl;
     cout << "Min Max face areas : " << minarea << " " << maxarea << endl;
     cout << "# of faces : " << mesh->getSize(2) << endl;
     cout << " Consistency : " << mesh->is_consistently_oriented() << endl;
     cout << "# of Inverted Faces : " << mesh->count_inverted_faces() << endl;
     cout << "# of Concave Faces  : " << mesh->count_concave_faces() << endl;
    // Jaal::set_large_area_tag(mesh);
    // Jaal::set_boundary_tag(mesh);
    // Jaal::set_tiny_area_tag(mesh);
  
    // Jaal::set_layer_tag(mesh);
    // Jaal::set_constrained_tag(mesh);
    // Jaal::set_doublet_tag(mesh);
    // Jaal::set_bridge_tag(mesh);
    // Jaal::set_diamond_tag(mesh);
    // Jaal::set_singlet_tag(mesh);
    // Jaal::set_regular_node_tag(mesh);
  }

  cout << " Saving Mesh " << outfilename << endl;
  mesh->saveAs( outfilename);

  if( lapweight ) delete lapweight;

  if( mesh ) mesh->deleteAll();

  delete mesh;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//                  Garbage Stuff ...
////////////////////////////////////////////////////////////////////////////////

/*
  double origin[]   = { 0.0, 0.0, 0.0};
  double length[]   = { 1.0, 1.0, 1.0};
  int    gridim[]   = { 5, 5, 5};
  mesh = Jaal::create_structured_mesh(origin, length, gridim, 2 );
  DijkstraShortestPath djk(mesh);
  NodeSequence sq = djk.getPath(mesh->getNodeAt(4));

  NodeSequence sq = mesh->get_breadth_first_ordered_nodes(mesh->getNodeAt(12), 2);
  for( int i = 0; i < sq.size(); i++)  {
       cout << i <<  " " << sq[i]->getID() << " Level  " << sq[i]->getLayerID() << endl;
  }
  mesh->saveAs( "dbg.dat");
  exit(0);
*/
/*
  int nx = 15;
  int ny = 15;
  vector<Vertex*> anodes(nx), bnodes(ny), cnodes(nx), dnodes(ny);
  Point3D xyz;
  double dx = 1.0/(double)(nx-1);
  double dy = 1.0/(double)(ny-1);
   
  for( int i = 0; i < nx; i++) {
      xyz[0] = i*dx;
      xyz[1] = 0.0;
      xyz[2] = 0.0;
      Vertex *v = Vertex::newObject();
      v->setXYZCoords(xyz);
      anodes[i] = v; 
      mesh->addNode( v );
  }

  bnodes[0] =  anodes[nx-1];
  for( int j = 1; j < ny; j++) {
      xyz[0] = 1.0;
      xyz[1] = j*dy;
      xyz[2] = 0.0;
      Vertex *v = Vertex::newObject();
      v->setXYZCoords(xyz);
      bnodes[j] = v; 
      mesh->addNode( v );
  }

  cnodes[0] = anodes[0];
  for( int i = 1; i < nx-1; i++) {
      xyz[0] = i*dx;
      xyz[1] = i*dx;
      xyz[2] = 0.0;
      Vertex *v = Vertex::newObject();
      v->setXYZCoords(xyz);
      cnodes[i] = v; 
      mesh->addNode( v );
  }
  cnodes[nx-1] =  bnodes[ny-1];
*/

/*
  dnodes[0] =  cnodes[nx-1];
  for( int j = 1; j < ny-1; j++) {
      xyz[0] = 0.0;
      xyz[1] = 1.0-j*dy;
      xyz[2] = 0.0;
      Vertex *v = Vertex::newObject();
      v->setXYZCoords(xyz);
      dnodes[j] = v; 
      mesh->addNode( v );
  }
  dnodes[ny-1] = anodes[0];
*/

/*
  remesh_tri_patch( mesh, anodes, bnodes, cnodes);
  mesh->saveAs("dbg.dat");
  cout << " Saved " << endl;
  exit(0);
*/
   
