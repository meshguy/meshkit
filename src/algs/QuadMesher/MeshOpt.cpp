#include "Mesh.hpp"
#include <ostream>

#ifdef HAVE_MESQUITE
#include "Mesquite_all_headers.hpp"
using namespace Mesquite;
#endif

///////////////////////////////////////////////////////////////////////////////

int Jaal::MeshOptimization::shape_optimize(Jaal::Mesh *mesh)
{
#ifdef HAVE_MESQUITE
    int topo = mesh->isHomogeneous();

    if( topo == 3 )
    {
        shape_tri_optimize(mesh);
        return 0;
    }

    if( topo == 4 )
    {
        vector<Vertex*> steiner;
        Mesh *trimesh = Jaal::quad_to_tri4( mesh, steiner);
//      Mesh *trimesh = Jaal::quad_to_tri2( mesh );
        assert( trimesh );
        shape_tri_optimize(trimesh);
        trimesh->deleteFaces();
        for( size_t i = 0; i < steiner.size(); i++)
            delete steiner[i];
        delete trimesh;
        return 0;
    }
#endif
    return 1;
}

#ifdef HAVE_MESQUITE

///////////////////////////////////////////////////////////////////////////////
Mesquite::ArrayMesh* Jaal::MeshOptimization::jaal_to_mesquite( Jaal::Mesh *mesh)
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if( !vertex->isRemoved() )
            vertex->setID( index++ );
    }

    size_t numActiveNodes = index;
    vector<int> vfixed;
    vfixed.resize(numActiveNodes);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if( !vertex->isRemoved() )
        {
            int lid = vertex->getID();
            if (vertex->isBoundary())
                vfixed[lid] = 1;
            else
                vfixed[lid] = 0;
        }
    }

    size_t noptfaces = 0;
    for( size_t i = 0; i < numfaces; i++)
    {
        Face *f = mesh->getFaceAt(i);
        if( !f->isRemoved() ) noptfaces++;
    }

    vector<size_t>  vNodes;
    vector<double>  vCoords;
    vector<size_t>  l2g;
    vector<int>     etopo;

    mesh->getCoordsArray( vCoords, l2g );
    mesh->getNodesArray( vNodes, etopo );

    Mesquite::ArrayMesh *optmesh = NULL;

    int topo = mesh->isHomogeneous();

    numnodes = vCoords.size()/3;
    assert( vfixed.size() == numnodes );

    if (topo == 4)
        optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                          noptfaces, Mesquite::QUADRILATERAL, &vNodes[0]);
    if (topo == 3)
    {
        optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                          noptfaces, Mesquite::TRIANGLE, &vNodes[0]);
    }

    return optmesh;
}

#endif

///////////////////////////////////////////////////////////////////////////////

int Jaal::MeshOptimization::shape_tri_optimize(Jaal::Mesh *mesh)
{
#ifdef HAVE_MESQUITE
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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if( !vertex->isRemoved() )
            vertex->setID( index++ );
    }

    size_t numActiveNodes = index;
    vector<int> vfixed;
    vfixed.resize(numActiveNodes);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if( !vertex->isRemoved() )
        {
            int lid = vertex->getID();
            if (vertex->isBoundary())
                vfixed[lid] = 1;
            else
                vfixed[lid] = 0;
        }
    }

    size_t noptfaces = 0;
    for( size_t i = 0; i < numfaces; i++)
    {
        Face *f = mesh->getFaceAt(i);
        if( !f->isRemoved() ) noptfaces++;
    }

    vector<size_t>  vNodes;
    vector<double>  vCoords;
    vector<size_t>  l2g;
    vector<int>     etopo;

    mesh->getCoordsArray( vCoords, l2g );
    mesh->getNodesArray( vNodes, etopo );

    Mesquite::ArrayMesh *optmesh = NULL;

    int topo = mesh->isHomogeneous();

    numnodes = vCoords.size()/3;
    assert( vfixed.size() == numnodes );

    if (topo == 4)
        optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                          noptfaces, Mesquite::QUADRILATERAL, &vNodes[0]);
    if (topo == 3)
    {
        optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                          noptfaces, Mesquite::TRIANGLE, &vNodes[0]);
    }

    if (optmesh == NULL) return 1;

    Vector3D normal(0, 0, 1);
    Vector3D point(0, 0, 0);
    PlanarDomain mesh_plane(normal, point);

    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio mesh_quality( err );
//  IdealWeightMeanRatio mesh_quality();
//  ConditionNumberQualityMetric mesh_quality();

    // sets the objective function template
    LPtoPTemplate obj_func(&mesh_quality, 2, err);

    // creates the optimization procedures
//  ConjugateGradient* improver = new ConjugateGradient( obj_func, err );


//  SteepestDescent improver( &obj_func);   // Fastest but poor convergence in the tail.
    QuasiNewton     improver( &obj_func);   // Best Quality ( two time slower than SD )
//  TrustRegion     improver( &obj_func);
//  FeasibleNewton  improver( &obj_func);
//  SmartLaplacianSmoother improver(&obj_func); // Slightly better than useless.

    //performs optimization globally
   improver.use_global_patch();

    // creates a termination criterion and
    // add it to the optimization procedure
    // outer loop: default behavior: 1 iteration
    // inner loop: stop if gradient norm < eps
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_gradient_L2_norm(1e-4);
    tc_inner.add_iteration_limit(100);
    improver.set_inner_termination_criterion(&tc_inner);

    // creates a quality assessor
    QualityAssessor m_ratio_qa(&mesh_quality);

    // creates an instruction queue
    InstructionQueue queue;
    queue.add_quality_assessor(&m_ratio_qa, err);
    queue.set_master_quality_improver(&improver, err);
    queue.add_quality_assessor(&m_ratio_qa, err);

    // do optimization of the mesh_set
    queue.run_instructions(optmesh, &mesh_plane, err);

    cout << "# of iterations " << tc_inner.get_iteration_count() << endl;

    if (err)
    {
        std::cout << err << std::endl;
        return 2;
    }

    mesh->setCoordsArray(vCoords, l2g);
    delete optmesh;

    return 0;
#endif

    return 1;
}

///////////////////////////////////////////////////////////////////////////////


int Jaal::MeshOptimization::untangle(Mesh *mesh)
{
#ifdef HAVE_MESQUITE
    /*
        MsqError err;
        unsigned long int numnodes = mesh->getSize(0);
        unsigned long int numfaces = mesh->getSize(2);

        mesh->search_boundary();

        vector<int> vfixed(numnodes);
        for (size_t i = 0; i < numnodes; i++) {
            Vertex *vertex = mesh->getNodeAt(i);
            assert(vertex->getID() == i);
            vfixed[i] = 0;
            if (vertex->isBoundary()) {
                vfixed[i] = 1;
            }
        }

        vector<unsigned long int> vNodes;
        vNodes.reserve(4 * numfaces);

        for (size_t i = 0; i < numfaces; i++) {
            Face *face = mesh->getFaceAt(i);
            for (int j = 0; j < face->getSize(0); j++) {
                Vertex *v = face->getNodeAt(j);
                vNodes.push_back(v->getID());
            }
        }

        vector<double> vCoords;
        vector<size_t> l2g;
        mesh->getCoordsArray( vCoords, l2g );

        Mesquite::ArrayMesh *optmesh = NULL;

        int topo = mesh->isHomogeneous();

        if (topo == 4)
            optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                              numfaces, Mesquite::QUADRILATERAL, &vNodes[0]);
        if (topo == 3)
            optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                              numfaces, Mesquite::TRIANGLE, &vNodes[0]);

        if (optmesh == NULL) return 1;

        // Set Domain Constraint
        Vector3D pnt(0, 0, 0);
        Vector3D s_norm(0, 0, 1);
        PlanarDomain msq_geom(s_norm, pnt);

        // creates an intruction queue
        InstructionQueue queue1;

        // creates a mean ratio quality metric ...
        ConditionNumberQualityMetric shape_metric;
        UntangleBetaQualityMetric untangle(2);
        Randomize pass0(.05);
        // ... and builds an objective function with it
        //LInfTemplate* obj_func = new LInfTemplate(shape_metric);
        LInfTemplate obj_func(&untangle);
        LPtoPTemplate obj_func2(&shape_metric, 2, err);
        if (err) return 1;
        // creates the steepest descent optimization procedures
        ConjugateGradient improver(&obj_func, err);
        if (err) return 1;

        //SteepestDescent* pass2 = new SteepestDescent( obj_func2 );
        ConjugateGradient pass2(&obj_func2, err);
        if (err) return 1;
        pass2.use_element_on_vertex_patch();
        if (err) return 1;
        pass2.use_global_patch();
        if (err) return 1;
        QualityAssessor stop_qa = QualityAssessor(&shape_metric);
        QualityAssessor stop_qa2 = QualityAssessor(&shape_metric);

        stop_qa.add_quality_assessment(&untangle);
        // **************Set stopping criterion**************
        //untangle beta should be 0 when untangled
        TerminationCriterion sc1;
        sc1.add_relative_quality_improvement(0.000001);
        TerminationCriterion sc3;
        sc3.add_iteration_limit(10);
        TerminationCriterion sc_rand;
        sc_rand.add_iteration_limit(1);

        //StoppingCriterion sc1(&stop_qa,-1.0,.0000001);
        //StoppingCriterion sc3(&stop_qa2,.9,1.00000001);
        //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,10);
        //StoppingCriterion sc_rand(StoppingCriterion::NUMBER_OF_PASSES,1);
        //either until untangled or 10 iterations
        pass0.set_outer_termination_criterion(&sc_rand);
        improver.set_outer_termination_criterion(&sc1);
        pass2.set_inner_termination_criterion(&sc3);

        // adds 1 pass of improver to mesh_set1
        queue1.add_quality_assessor(&stop_qa, err);
        if (err) return 1;
        //queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
        //queue1.add_preconditioner(improver,err);MSQ_CHKERR(err);
        //queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
        queue1.set_master_quality_improver(&improver, err);
        if (err) return 1;
        queue1.add_quality_assessor(&stop_qa2, err);
        if (err) return 1;

        // launches optimization on mesh_set1
        queue1.run_instructions(optmesh, &msq_geom, err);
        if (err) return 1;

        mesh->setCoordsArray(vCoords);
        return 0;
    */
#endif
    return 1;
}

///////////////////////////////////////////////////////////////////////////////



#ifdef HAVE_MESQUITE
int
run_global_smoother(Mesquite::Mesh* mesh, MsqError& err)
{
    double OF_value = 0.0001;

    // creates an intruction queue
    InstructionQueue queue1;

    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
    if (err) return 1;

    mean_ratio->set_averaging_method(QualityMetric::SUM, err);
    if (err) return 1;

    // ... and builds an objective function with it
    LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
    if (err) return 1;

    // creates the feas newt optimization procedures
    FeasibleNewton* improver = new FeasibleNewton(obj_func, true);
    improver->use_global_patch();
    if (err) return 1;

    QualityAssessor stop_qa(mean_ratio);

    // **************Set stopping criterion****************
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_vertex_movement(OF_value);
    if (err) return 1;

    TerminationCriterion tc_outer;
    tc_outer.add_iteration_limit(1);
    improver->set_inner_termination_criterion(&tc_inner);
    improver->set_outer_termination_criterion(&tc_outer);

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // adds 1 pass of improver to mesh_set1
    queue1.set_master_quality_improver(improver, err);
    if (err) return 1;

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // launches optimization on mesh_set
    queue1.run_instructions(mesh, err);
    cout << " Error " << err << endl;
    if (err) return 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
run_local_smoother(Mesquite::Mesh* mesh, MsqError& err)
{
    double OF_value = 0.0001;

    // creates an intruction queue
    InstructionQueue queue1;

    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
    if (err) return 1;

    mean_ratio->set_averaging_method(QualityMetric::SUM, err);
    if (err) return 1;

    // ... and builds an objective function with it
    LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
    if (err) return 1;

    // creates the smart laplacian optimization procedures
    SmartLaplacianSmoother* improver = new SmartLaplacianSmoother(obj_func);

    QualityAssessor stop_qa(mean_ratio);

    // **************Set stopping criterion****************
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_vertex_movement(OF_value);

    TerminationCriterion tc_outer;
    tc_outer.add_iteration_limit(1);

    improver->set_inner_termination_criterion(&tc_inner);
    improver->set_outer_termination_criterion(&tc_outer);

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // adds 1 pass of improver to mesh_set
    queue1.set_master_quality_improver(improver, err);
    if (err) return 1;

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // launches optimization on mesh_set
    queue1.run_instructions(mesh, err);
    if (err) return 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
#endif
