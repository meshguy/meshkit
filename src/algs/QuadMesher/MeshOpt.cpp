#include "Mesh.hpp"
#include <ostream>

#ifdef HAVE_MESQUITE
#include "Mesquite_all_headers.hpp"
using namespace Mesquite;
#endif

///////////////////////////////////////////////////////////////////////////////

int Jaal::MeshOptimization::shape_optimize(Jaal::Mesh *mesh, int algo, int niter)
{
     algorithm = algo;
     numIter   = niter;

#ifdef HAVE_MESQUITE
     int topo = mesh->isHomogeneous();

     if( topo == 3 ) {
          execute(mesh);
          return 0;
     }

     if( topo == 4 ) {
          if( mesh->count_concave_faces() ) {
               cout << " Optimizing triangle mesh of the quadrilateral mesh: " << endl;
               vector<Vertex*> steiner;
               Mesh *trimesh = Jaal::quad_to_tri4( mesh, steiner);
               execute(trimesh);
               trimesh->deleteFaces();
               for( size_t i = 0; i < steiner.size(); i++)
                    delete steiner[i];
               delete trimesh;
          } else {
               cout << " Optimizing Quadrilateral mesh directly: " << endl;
               execute(mesh);
          }

          return 0;
     }
#endif
     return 1;
}

#ifdef HAVE_MESQUITE
///////////////////////////////////////////////////////////////////////////////
Mesquite::ArrayMesh* Jaal::MeshOptimization::jaal_to_mesquite( Jaal::Mesh *mesh )
{
     MsqError err;

     unsigned long int numnodes = mesh->getSize(0);
     unsigned long int numfaces = mesh->getSize(2);

     //
     //////////////////////////////////////////////////////////////////////////////
     // Mesquite works only when the faces are convex, and our mesh may contain some
     // truly bad elements. Therefore, we will skip those faces and mark them fix.
     // Hopefully, later by some other operations ( smoothing, decimation, refinement
     // etc.) those bad elements may get removed.

     // By ignoring some of the faces, the mesh may get disconnected. ( I need to test
     // such cases ).
     //
     // The mesh must be doublets free. ( I found out this after some pains ).
     //////////////////////////////////////////////////////////////////////////////

     mesh->search_boundary();

     size_t index = 0;
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if( vertex->isActive() ) vertex->setID( index++ );
     }
     size_t noptnodes = index;

     vfixed.resize(noptnodes);
     for (size_t i = 0; i < numnodes; i++) {
          Vertex *vertex = mesh->getNodeAt(i);
          if( vertex->isActive() ) {
               int lid = vertex->getID();
               if (vertex->isBoundary())
                    vfixed[lid] = 1;
               else
                    vfixed[lid] = 0;
          }
     }

     size_t noptfaces = 0;
     for( size_t i = 0; i < numfaces; i++) {
          Face *f = mesh->getFaceAt(i);
          if( f->isActive() ) noptfaces++;
     }

     mesh->getCoordsArray( vCoords, l2g );
     mesh->getNodesArray( vNodes, etopo );

     Mesquite::ArrayMesh *optmesh = NULL;

     int topo = mesh->isHomogeneous();

     assert( noptnodes == vCoords.size()/3 );
     assert( noptnodes == vfixed.size());

     if (topo == 4)
          optmesh = new Mesquite::ArrayMesh(3, noptnodes, &vCoords[0], &vfixed[0],
                                            noptfaces, Mesquite::QUADRILATERAL, &vNodes[0]);
     if (topo == 3) {
          optmesh = new Mesquite::ArrayMesh(3, noptnodes, &vCoords[0], &vfixed[0],
                                            noptfaces, Mesquite::TRIANGLE, &vNodes[0]);
     }

     return optmesh;
}

#endif

///////////////////////////////////////////////////////////////////////////////

int Jaal::MeshOptimization::execute(Jaal::Mesh *mesh)
{
#ifdef HAVE_MESQUITE
     MsqError err;

     Mesquite::ArrayMesh *optmesh = NULL;
     optmesh = Jaal::MeshOptimization::jaal_to_mesquite( mesh);
     if( optmesh == NULL ) return 1;

     Vector3D normal(0, 0, 1);
     Vector3D point (0, 0, 0);
     PlanarDomain mesh_plane(normal, point);

     // creates a mean ratio quality metric ...
     IdealWeightInverseMeanRatio mesh_quality;

//  IdealWeightMeanRatio mesh_quality;
//  ConditionNumberQualityMetric mesh_quality;
//  EdgeLengthQualityMetric mesh_quality;

     // sets the objective function template
     LPtoPTemplate obj_func(&mesh_quality, 2, err);

     // creates the optimization procedures

     TerminationCriterion tc_outer;

     TerminationCriterion tc_inner;
     tc_inner.add_absolute_gradient_L2_norm(1e-4);
     tc_inner.write_iterations("opt.dat", err);

     SteepestDescent  *sp = NULL;
     QuasiNewton      *qn = NULL;
     TrustRegion      *tr = NULL;
     FeasibleNewton   *fn = NULL;
     ConjugateGradient *cg = NULL;
     SmartLaplacianSmoother *lp = NULL;

     improver  = NULL;

     switch (algorithm) {
     case STEEPEST_DESCENT:
          sp = new SteepestDescent( &obj_func);   // Fastest but poor convergence in the tail.
          sp->use_global_patch();
          improver = sp;
          tc_inner.add_iteration_limit(numIter);
          improver->set_inner_termination_criterion(&tc_inner);
          break;
     case QUASI_NEWTON:
          qn = new QuasiNewton( &obj_func);   // Fastest but poor convergence in the tail.
          qn->use_global_patch();
//        qn->use_element_on_vertex_patch();
          improver = qn;
          tc_inner.add_iteration_limit(numIter);
          improver->set_inner_termination_criterion(&tc_inner);
          break;
     case TRUST_REGION:
          tr = new TrustRegion( &obj_func);   // Fastest but poor convergence in the tail.
          tr->use_global_patch();
          improver = tr;
          tc_inner.add_iteration_limit(numIter);
          improver->set_inner_termination_criterion(&tc_inner);
          break;
     case FEASIBLE_NEWTON:
          fn = new FeasibleNewton( &obj_func);   // Fastest but poor convergence in the tail.
          fn->use_global_patch();
          improver = fn;
          tc_inner.add_iteration_limit(numIter);
          improver->set_inner_termination_criterion(&tc_inner);
          break;
     case CONJUGATE_GRADIENT:
          cg = new ConjugateGradient(&obj_func);
          cg->use_global_patch();
          improver = cg;
          tc_inner.add_iteration_limit(numIter);
          improver->set_inner_termination_criterion(&tc_inner);
          break;
     case LAPLACIAN:
          lp = new SmartLaplacianSmoother( &obj_func);   // Fastest but poor convergence in the tail.
          improver = lp;
          tc_outer.add_iteration_limit(numIter);
          improver->set_outer_termination_criterion(&tc_outer);
          break;
     default:
          cout << "Warning: Invalid optimization algorithm selected: "<< endl;
     }

     if( improver ) {

          // creates a quality assessor
          QualityAssessor qa(&mesh_quality);

          // creates an instruction queue
          InstructionQueue queue;
          queue.add_quality_assessor(&qa, err);
          queue.set_master_quality_improver(improver, err);
          queue.add_quality_assessor(&qa, err);

          // do optimization of the mesh_set
          queue.run_instructions(optmesh, &mesh_plane, err);

          cout << "# of iterations " << tc_inner.get_iteration_count() << endl;

          if (err) {
               std::cout << err << std::endl;
               return 2;
          }

          mesh->setCoordsArray(vCoords, l2g);
          delete improver;
     }

     delete optmesh;

     return 0;
#endif

     return 1;
}

///////////////////////////////////////////////////////////////////////////////
