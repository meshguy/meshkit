#include "meshkit/Mesh.hpp"
#include "meshkit/MeshRefine2D.hpp"
#include "meshkit/QuadCleanUp.hpp"

#include "meshkit/DijkstraShortestPath.hpp"
#include "meshkit/ObjectPool.hpp"

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
     cout << " 18     :   Generate Quad-to-Tri4 " << endl;
     cout << " 19     :   Generate Quad-to-Tri2 " << endl;
     cout << " 20     :   Shift irregular nodes inside domain " << endl;
     cout << " 21     :   Everything automatic " << endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
     Jaal::Mesh *mesh = new Jaal::Mesh;
     Jaal::MeshOptimization mopt;

     /*
          double origin[]   = { 0.0, 0.0, 0.0};
          double length[]   = { 1.0, 1.0, 1.0};
          int    gridim[]   = { 6, 7, 2};
          mesh = Jaal::create_structured_mesh(origin, length, gridim, 2 );
          cout << mesh->getSize(1) << endl;
          Face *f = mesh->getFaceAt(0);
          mesh->remove(f);

          f = mesh->getFaceAt(1);
          mesh->remove(f);

          cout << mesh->getSize(1) << endl;
          mesh->saveAs( "tmp.off");
          exit(0);
     */

     string infilename, outfilename;

     int iopt, numiters = 100, topo_improve = 1;
     int cleanup_op = 0;

     while( (iopt = getopt( argc, argv, "hgtc:i:o:") ) != -1) {
          switch(iopt) {
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

     size_t ninvert  =  mesh->count_inverted_faces();
     size_t numfaces =  mesh->getSize(2);
     size_t numBound =  mesh->getBoundarySize(0);
     size_t nireg0   =  mesh->count_irregular_nodes(4);

     cout << "# of irregular nodes before cleanup : " << nireg0 << endl;

     if( ninvert > 0.5*numfaces )
          mesh->reverse();

     LaplaceSmoothing lapsmooth(mesh);
     LaplaceWeight *lapweight = NULL;

     QuadCleanUp qClean(mesh);

     vector<QTrack>  qpath;
     vector<Vertex*> steiner;
     Mesh *q2t;
     int  algo, numiter;

     StopWatch swatch;
     swatch.start();

     switch( cleanup_op) {
     case 0:
          qClean.report();
          exit(0);
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
          lapweight = new LaplaceNoWeight();
          lapsmooth.setWeight(lapweight);
          lapsmooth.setNumIterations(100);
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
//         qClean.advancing_front_edges_swap();
          break;
     case 9:
          /*
                    cout << "Choose algorithm : " << endl;
                    cout << "    Steepest Descent       : 0 " << endl;
                    cout << "    Quasi Newton(default)  : 1  " << endl;
                    cout << "    Trust Region           : 2 " << endl;
                    cout << "    Feasible Newton        : 3 " << endl;
                    cout << "    Laplacian              : 4 " << endl;
                    cin  >> algo;
                    cout << "Give number of iterations " << endl;
                    cin  >> numiter;
                    mopt.shape_optimize( mesh, algo, numiter );
          */
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
          qpath = Jaal::generate_quad_partitioning(mesh);
//      Jaal::set_irregular_path_tag(mesh, qpath);
          break;
     case 18:
          q2t = Jaal::quad_to_tri4( mesh, steiner);
          q2t->saveAs("tmesh.off");
          break;
     case 19:
          q2t = Jaal::quad_to_tri2( mesh );
          mopt.shape_optimize( q2t );
          q2t->saveAs("tmesh.off");
          break;
     case 20:
          qClean.shift_irregular_nodes();
          break;
     case 21:
          qClean.automatic();
          break;
     }
     swatch.stop();
     cout << "CleanUp time : " << swatch.getSeconds() << endl;

     if( cleanup_op ) {

          cout << "# Nodes           : " << mesh->getSize(0) << endl;
          cout << "# Faces           : " << mesh->getSize(2) << endl;
          cout << "# Inverted Faces  : " << mesh->count_inverted_faces() << endl;
          cout << "# Concave Faces   : " << mesh->count_concave_faces() << endl;
          cout << "# Irregular nodes : " << mesh->count_irregular_nodes(4) << endl;
//        cout << "Mesh Consistency  : " << mesh->is_consistently_oriented() << endl;

          // Jaal::set_large_area_tag(mesh);
          // Jaal::set_boundary_tag(mesh);
          // Jaal::set_tiny_area_tag(mesh);

          // Jaal::set_layer_tag(mesh);
          // Jaal::set_constrained_tag(mesh);
          // Jaal::set_doublet_tag(mesh);
          // Jaal::set_bridge_tag(mesh);
          // Jaal::set_singlet_tag(mesh);
          // Jaal::set_regular_node_tag(mesh);
     }

     mesh->collect_garbage();
//   Jaal::set_diamond_tag(mesh);

     cout << " Saving Mesh " << outfilename << endl;
     mesh->saveAs( outfilename);

     mesh->get_topological_statistics();
     plot_all_quad_quality_measures( mesh );

     assert( numBound == mesh->getBoundarySize(0) );

     if( lapweight ) delete lapweight;

     if( mesh ) mesh->deleteAll();

     delete mesh;

     return 0;
}

