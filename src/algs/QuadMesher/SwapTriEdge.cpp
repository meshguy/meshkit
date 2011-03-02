#include "SwapTriEdge.hpp"

//#############################################################################

int SwapTriEdge::execute()
{
    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    num_edges_flipped = 0;

    size_t numfaces = mesh->getSize(2);
    while (1)
    {
        size_t ncount = 0;

        for (size_t i = 0; i < numfaces; i++)
        {
            Face *face = mesh->getFaceAt(i);
            face->setVisitMark(0);
        }

        for (size_t i = 0; i < numfaces; i++)
        {
            int err = atomicOp(mesh->getFaceAt(i));
            if (!err) ncount++;
        }

        if (ncount == 0) break;
        num_edges_flipped += ncount;
    }

    if (!relexist) mesh->clear_relations(0, 2);

    return num_edges_flipped;
}

//#############################################################################

int SwapTriEdge::atomicOp(const Face *face)
{
    for (int i = 0; i < 3; i++)
    {
        Vertex *v1 = face->getNodeAt((i + 1) % 3);
        Vertex *v2 = face->getNodeAt((i + 2) % 3);
        FlipEdge edge(v1, v2);
        if (is_edge_flip_allowed(edge))
        {
            commit(edge);
            return 0;
        }
    }
    return 1;
}

//#############################################################################

void SwapTriEdge::FlipEdge::process(Vertex *v1, Vertex *v2)
{
    assert(v1 != v2);

    this->setNodes(v1, v2);

    faces[0] = NULL;
    faces[1] = NULL;
    opposite_nodes[0] = NULL;
    opposite_nodes[1] = NULL;

    if (v1->isBoundary() && v2->isBoundary()) return;

    FaceSequence neighs = Mesh::getRelations112(v1, v2);

    assert(neighs.size() == 2);

    Vertex *ov1 = Face::opposite_node(neighs[0], v1, v2);
    Vertex *ov2 = Face::opposite_node(neighs[1], v1, v2);
    assert(ov1 && ov2);

    faces[0] = neighs[0];
    faces[1] = neighs[1];
    opposite_nodes[0] = ov1;
    opposite_nodes[1] = ov2;
}

//#############################################################################

bool SwapTriEdge::FlipEdge::isSharp(double featureAngle) const
{
    Vertex *v1 = getNodeAt(0);
    Vertex *v2 = getNodeAt(1);
    Vertex *ov1 = opposite_nodes[0];
    Vertex *ov2 = opposite_nodes[1];

    Vec3D A, B;
    A = Face::normal(ov1, v1, v2);
    B = Face::normal(ov2, v2, v1);

    double angle = Math::getVectorAngle(A, B, ANGLE_IN_DEGREES);
    if (angle > featureAngle && fabs(180 - angle) > 1.0E-06) return 1;

    return 0;
}

//#############################################################################

bool SwapTriEdge::FlipEdge::isConcave() const
{
    Vertex *v1 = getNodeAt(0);
    Vertex *v2 = getNodeAt(1);
    Vertex *ov1 = opposite_nodes[0];
    Vertex *ov2 = opposite_nodes[1];

    bool convex = Face::is_convex_quad(v1->getXYZCoords(),
                                       ov1->getXYZCoords(),
                                       v2->getXYZCoords(),
                                       ov2->getXYZCoords());
    if (!convex) return 1;

    return 0;
}


//*****************************************************************************

bool SwapTriEdge::is_edge_flip_allowed(const FlipEdge &edge) const
{
    Vertex *v1 = edge.getNodeAt(0);
    Vertex *v2 = edge.getNodeAt(1);

    if (v1->isBoundary() && v2->isBoundary()) return 0;

    if (edge.faces[0]->isVisited() || edge.faces[1]->isVisited()) return 0;

    Vertex *ov1 = edge.opposite_nodes[0];
    Vertex *ov2 = edge.opposite_nodes[1];

    double len1 = Vertex::length2(v1, v2);
    double len2 = Vertex::length2(ov1, ov2);
    if (len1 <= len2) return 0;

    if (edge.isConcave()) return 0;
    if (edge.isSharp(creaseAngle)) return 0;

    return 1;
}

//#############################################################################

int SwapTriEdge::commit(const FlipEdge &edge)
{
    Face *t1 = edge.faces[0];
    Face *t2 = edge.faces[1];
    Vertex *v1 = edge.getNodeAt(0);
    Vertex *v2 = edge.getNodeAt(1);
    Vertex *ov1 = edge.opposite_nodes[0];
    Vertex *ov2 = edge.opposite_nodes[1];

    v1->removeRelation2(t2);
    v2->removeRelation2(t1);
    ov1->addRelation2(t2);
    ov2->addRelation2(t1);

    NodeSequence vconn(3);
    vconn[0] = ov2;
    vconn[1] = ov1;
    vconn[2] = v1;
    t1->setNodes(vconn);

    vconn[0] = v2;
    vconn[1] = ov1;
    vconn[2] = ov2;
    t2->setNodes(vconn);

    t1->setVisitMark(1);
    t2->setVisitMark(1);

    return 0;
}

//#############################################################################

int VertexDegreeReduction::getVertexDegree(const Vertex *vertex) const
{
    assert(vertex);
    NodeSequence vneighs = vertex->getRelations0();
    return vneighs.size();
}

//#############################################################################

bool VertexDegreeReduction::is_edge_flip_allowed(const FlipEdge &flipedge) const
{
    if (flipedge.isBoundary()) return 0;
    if (flipedge.isConcave()) return 0;

    Vertex *v1 = flipedge.getNodeAt(0);
    Vertex *v2 = flipedge.getNodeAt(1);
    Vertex *ov1 = flipedge.opposite_nodes[0];
    Vertex *ov2 = flipedge.opposite_nodes[1];

    int d1 = getVertexDegree(v1);
    int d2 = getVertexDegree(v2);
    int d3 = getVertexDegree(ov1);
    int d4 = getVertexDegree(ov2);

    int relaxation_index = d1 + d2 - d3 - d4;

    if (relaxation_index < 3) return 0;
    return 1;
}

//#############################################################################

int VertexDegreeReduction::atomicOp(const Vertex *apexVertex)
{
    /*
        SharedEdgeMesh coveredges = facemesh->get_cover_edges(apexVertex);
        MeshUtil::make_chain(coveredges);

        SharedNodeSet nodeset = coveredges->getNodes();

        size_t nsuccess = 0;
        size_t numNeighs = nodeset->getSize();
        int istart = 0;
        if( numNeighs%2) istart = 1;

        //
        //1. Greedy algorithm don't work properly in many cases, therefore
        //   identify all the flippable edges first and later commit all of them
        //
        toCommit.clear();
        std::set<SharedVertex> visited;
        for( int i = istart; i < numNeighs; i=i+2) {
             SharedVertex dstVertex = nodeset->at(i%numNeighs);
             if( visited.find(dstVertex) == visited.end() ) {
                 FlipEdge flipedge(apexVertex,dstVertex);
                 if( is_edge_flip_allowed(flipedge) )  {
                     toCommit.push_back( flipedge);
                     visited.insert( flipedge.opposite_nodes[0] );
                     visited.insert( flipedge.opposite_nodes[1] );
                 }
             }

        }

        nsuccess = toCommit.size();
        for( int i = 0; i < nsuccess; i++)
             commit( toCommit[i] );

        if( nsuccess ) atomicOp(apexVertex);
        return nsuccess;
     */
}
///////////////////////////////////////////////////////////////////////////

int VertexDegreeReduction::execute()
{
    int relexist2 = mesh->build_relations(0, 2);
    int relexist0 = mesh->build_relations(0, 2);

    /*
      SharedNodeSet nodeset = facemesh->getNodes();
      std::sort( nodeset->container.begin(), nodeset->container.end(),
                 HighVertexDegreeCompare() );

     */
    size_t numNodes = mesh->getSize(0);

    int nsuccess = 0;
    for (size_t i = 0; i < numNodes; i++)
    {
        nsuccess += atomicOp(mesh->getNodeAt(i));
    }

    if (!relexist2) mesh->clear_relations(0, 2);
    if (!relexist0) mesh->clear_relations(0, 0);
    return nsuccess;

}

#ifdef CSV
///////////////////////////////////////////////////////////////////////////

void UnitTest::test_vertex_degree_reduction()
{
    int N = 100;
    double theta = 2.0 * M_PI / (double) N;
    Point3D p3d;
    SharedNodeSet nodeset = NodeSet::newObject();
    for (int i = 0; i < N; i++)
    {
        p3d[0] = cos(i * theta);
        p3d[1] = sin(i * theta);
        p3d[2] = 0.0;
        SharedVertex v = Vertex::newObject();
        v->setCoords(p3d);
        nodeset->add(v);
    }
    SharedEdgeMesh edgemesh = MeshUtil::create_close_contour(nodeset);

    boost::tuple<SharedFaceMesh, SharedVertex> product;
    product = MeshUtil::convex_star_triangulation(edgemesh);

    SharedFaceMesh facemesh = product.get < 0 > ();

    nodeset = facemesh->getNodes();
    nodeset->getRenumbered();

    VTKMeshExporter vtk;
    vtk.saveAs(facemesh, "vreduce_before.vtk");

    VertexDegreeReduction vd(facemesh);
    vd.execute();
    vtk.saveAs(facemesh, "vreduce_after.vtk");
}
///////////////////////////////////////////////////////////////////////////
#endif
