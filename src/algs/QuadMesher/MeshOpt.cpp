#include <meshkit/Mesh.hpp>
#include <ostream>

#ifdef HAVE_MESQUITE
#include "Mesquite_all_headers.hpp"
using namespace Mesquite;

using namespace Jaal;

int Jaal::MeshOptimization::shape_optimize(Mesh *mesh)
{

    MsqError err;

    unsigned long int numnodes = mesh->getSize(0);
    unsigned long int numfaces = mesh->getSize(2);

    int relexist2 = mesh->build_relations(0, 2);

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

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    vector<int> vfixed(numnodes);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vfixed[i] = 0;
        if (vertex->isBoundary())
        {
            vfixed[i] = 1;
        }
        else
        {
            FaceSequence vfaces = vertex->getRelations2();
            if (vfaces.size() == 2)
            {
                vfaces[0]->setVisitMark(1);
                vfaces[1]->setVisitMark(1);
            }
        }
    }

    if (!relexist2) mesh->clear_relations(0, 2);

    size_t noptfaces = 0;

    vector<unsigned long int> vNodes;
    double coords[4][3];
    Point3D xyz;
    vNodes.reserve(4 * numfaces);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
//      if (face->invertedAt() >= 0) face->setVisitMark(1);
        int nsize = face->getSize(0);
        for( int j = 0; j <  nsize; j++) {
             xyz = face->getNodeAt(j)->getXYZCoords();
             coords[j][0] =  xyz[0];
             coords[j][1] =  xyz[1];
             coords[j][2] =  xyz[2];
        }
 
        if( nsize == 4 ) 
            if (v_quad_maximum_angle(nsize, coords) > 179) face->setVisitMark(1);

        if (face->isVisited())
        {
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *vertex = mesh->getNodeAt(j);
                vfixed[ vertex->getID() ] = 1;
            }
        }
        else
        {
            noptfaces++;
            for (int j = 0; j < face->getSize(0); j++)
            {
                Vertex *v = face->getNodeAt(j);
                vNodes.push_back(v->getID());
            }
        }
    }

    vector<double> vCoords = mesh->getCoordsArray();

    Mesquite::ArrayMesh *optmesh = NULL;

    int topo = mesh->isHomogeneous();

    if (topo == 4)
        optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                          noptfaces, Mesquite::QUADRILATERAL, &vNodes[0]);
    if (topo == 3) {
        cout << " Hello " << endl;
        optmesh = new Mesquite::ArrayMesh(3, numnodes, &vCoords[0], &vfixed[0],
                                          noptfaces, Mesquite::TRIANGLE, &vNodes[0]);
   }

    if (optmesh == NULL) return 1;

    Vector3D normal(0, 0, 1);
    Vector3D point(0, 0, 0);
    PlanarDomain mesh_plane(normal, point);

    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio inverse_mean_ratio(err);
    // sets the objective function template
    LPtoPTemplate obj_func(&inverse_mean_ratio, 2, err);
    // creates the optimization procedures
    FeasibleNewton f_newton(&obj_func);
    //performs optimization globally
    f_newton.use_global_patch();
    // creates a termination criterion and
    // add it to the optimization procedure
    // outer loop: default behavior: 1 iteration
    // inner loop: stop if gradient norm < eps
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_gradient_L2_norm(1e-4);
    f_newton.set_inner_termination_criterion(&tc_inner);

    // creates a quality assessor
    QualityAssessor m_ratio_qa(&inverse_mean_ratio);

    // creates an instruction queue
    InstructionQueue queue;
    queue.add_quality_assessor(&m_ratio_qa, err);
    queue.set_master_quality_improver(&f_newton, err);
    queue.add_quality_assessor(&m_ratio_qa, err);

    // do optimization of the mesh_set
    queue.run_instructions(optmesh, &mesh_plane, err);

    if (err)
    {
        std::cout << err << std::endl;
        return 2;
    }

    mesh->setCoordsArray(vCoords);
    delete optmesh;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Jaal::MeshOptimization::untangle(Mesh *mesh)
{
    MsqError err;
    unsigned long int numnodes = mesh->getSize(0);
    unsigned long int numfaces = mesh->getSize(2);

    mesh->search_boundary();

    vector<int> vfixed(numnodes);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        assert(vertex->getID() == i);
        vfixed[i] = 0;
        if (vertex->isBoundary())
        {
            vfixed[i] = 1;
        }
    }

    vector<unsigned long int> vNodes;
    vNodes.reserve(4 * numfaces);

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        for (int j = 0; j < face->getSize(0); j++)
        {
            Vertex *v = face->getNodeAt(j);
            vNodes.push_back(v->getID());
        }
    }

    vector<double> vCoords = mesh->getCoordsArray();

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
    ConjugateGradient pass1(&obj_func, err);
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
    pass1.set_outer_termination_criterion(&sc1);
    pass2.set_inner_termination_criterion(&sc3);

    // adds 1 pass of pass1 to mesh_set1
    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;
    //queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
    //queue1.add_preconditioner(pass1,err);MSQ_CHKERR(err);
    //queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
    queue1.set_master_quality_improver(&pass1, err);
    if (err) return 1;
    queue1.add_quality_assessor(&stop_qa2, err);
    if (err) return 1;

    // launches optimization on mesh_set1
    queue1.run_instructions(optmesh, &msq_geom, err);
    if (err) return 1;

    mesh->setCoordsArray(vCoords);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////


using namespace Mesquite;
using Mesquite::MsqError;

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
    FeasibleNewton* pass1 = new FeasibleNewton(obj_func, true);
    pass1->use_global_patch();
    if (err) return 1;

    QualityAssessor stop_qa(mean_ratio);

    // **************Set stopping criterion****************
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_vertex_movement(OF_value);
    if (err) return 1;
    TerminationCriterion tc_outer;
    tc_outer.add_iteration_limit(1);
    pass1->set_inner_termination_criterion(&tc_inner);
    pass1->set_outer_termination_criterion(&tc_outer);

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // adds 1 pass of pass1 to mesh_set1
    queue1.set_master_quality_improver(pass1, err);
    if (err) return 1;

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // launches optimization on mesh_set
    queue1.run_instructions(mesh, err);
    cout << " Error " << err << endl;
    if (err) return 1;

    MeshWriter::write_vtk(mesh, "feasible-newton-result.vtk", err);
    if (err) return 1;
    cout << "Wrote \"feasible-newton-result.vtk\"" << endl;

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
    SmartLaplacianSmoother* pass1 = new SmartLaplacianSmoother(obj_func);

    QualityAssessor stop_qa(mean_ratio);

    // **************Set stopping criterion****************
    TerminationCriterion tc_inner;
    tc_inner.add_absolute_vertex_movement(OF_value);
    TerminationCriterion tc_outer;
    tc_outer.add_iteration_limit(1);
    pass1->set_inner_termination_criterion(&tc_inner);
    pass1->set_outer_termination_criterion(&tc_outer);

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // adds 1 pass of pass1 to mesh_set
    queue1.set_master_quality_improver(pass1, err);
    if (err) return 1;

    queue1.add_quality_assessor(&stop_qa, err);
    if (err) return 1;

    // launches optimization on mesh_set
    queue1.run_instructions(mesh, err);
    if (err) return 1;

    MeshWriter::write_vtk(mesh, "smart-laplacian-result.vtk", err);
    if (err) return 1;
    cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;

    //print_timing_diagnostics( cout );
    return 0;
}

#endif

///////////////////////////////////////////////////////////////////////////////
