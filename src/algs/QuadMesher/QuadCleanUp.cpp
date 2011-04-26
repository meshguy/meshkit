#include "QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::automatic()
{
    int err, stage = 1;

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(100);
//
//  ****************************************************************************
// If the mesh is triangular first convert it into Quad mesh ....
//  ****************************************************************************
//
    if (mesh->isHomogeneous () == 3) {
        Tri2Quads t2quad;
        Mesh *quadmesh = t2quad.getQuadMesh( mesh, 1);
        delete mesh; // deletes the triangle mesh...
        mesh = quadmesh;
    }
    exit(0);

// Throughout the cleaning process, euler characteristic should remain same
    int euler1, euler0 = mesh->getEulerCharacteristic(); // Invariant

// Total area of the domain must be same ...
    double area0 = mesh->getSurfaceArea();

//  Ensure that there are irregular nodes in the mesh. if not, you are lucky and done.
    NodeSequence irreg_nodes = mesh->get_irregular_nodes( 4 );
    if( irreg_nodes.empty() ) {
        cout << "Great: There are no irregular nodes in the mesh" << endl;
        mopt.shape_optimize( mesh );
        return 0;
    }

//  ****************************************************************************
//  Triangle to Quad Transformation starts from here ...Input preparation..
//  ****************************************************************************

    assert( mesh->isHomogeneous () == 4);
    cout << " Input Mesh :    " << endl;
    cout << "      # Nodes : " << mesh->getSize(0) << endl;
    cout << "      # Faces : " << mesh->getSize(2) << endl;
    cout << "      # Components : " << mesh->getNumOfComponents() << endl;
    cout << "      # of irregular nodes " << irreg_nodes.size() << endl;
    cout << "      # of concave faces   " << mesh->count_concave_faces()    << endl;
    cout << "      Euler Characteristics : " << euler0 << endl;
    cout << "      Surface Area : " << area0 << endl;

//  ***************************************************************************
//  Mesh connectivity must be consistent:  Condition: Strict
//  ***************************************************************************

    if( !mesh->is_consistently_oriented() )
        mesh->make_consistently_oriented();

//  Initial mesh may have different connectivity. Condition Strict
    size_t ninvert  =  mesh->count_inverted_faces();
    size_t numfaces =  mesh->getSize(2);
    if( ninvert > 0.5*numfaces ) mesh->reverse();

//  ****************************************************************************
//  Initial mesh may have  doublets, remove them, they are troublesome: Condition Strict
//  ****************************************************************************

#ifdef DEBUG
    Jaal::set_doublet_tag(mesh);
    mesh->saveAs("stage0.dat");
#endif

    err = remove_interior_doublets();
    if( err ) {
        cout << "Fatal Error: There are interior doublets in the mesh " << endl;
        cout << "Check the mesh in : dbg.dat" << endl;
        mesh->saveAs( "dbg.dat");
        exit(0);
        return stage;
    }

    cout << "Info:  Stage "  << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces()    << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert( euler0 == euler1 );

//  ****************************************************************************
//  Check the boundary nodes connectivity, and ensure that all elements adjacent
//  to them are convex. Singlet must be called after doublets removal ...
//  ****************************************************************************
    err = remove_boundary_singlets();

#ifdef DEBUG
    Jaal::set_no_tags(mesh);
    mesh->saveAs("stage1.dat");
#endif

    swap_concave_faces();

#ifdef DEBUG
    Jaal::set_no_tags(mesh);
    mesh->saveAs("stage2.dat");
#endif

    if( mesh->count_concave_faces() )
        lapsmooth.convexify();

#ifdef DEBUG
    Jaal::set_no_tags(mesh);
    mesh->saveAs("stage3.dat");
#endif

    if( !mesh->count_concave_faces() )
        mopt.shape_optimize( mesh );

    if( mesh->count_concave_faces() )
        lapsmooth.convexify();

    if( mesh->count_concave_faces() )
        swap_concave_faces();

#ifdef DEBUG
    Jaal::set_no_tags(mesh);
    mesh->saveAs("stage4.dat");
#endif

    cout << "Info:  Stage "  << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces()    << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert( euler0 == euler1 );


//  ***************************************************************************
//  Tunnels may put constraints in the movements, so it is better to remove them,
//  Condition: Soft
//  ***************************************************************************

    err = remove_tunnels();
    cout << "Info:  Stage "  << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces()    << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert( euler0 == euler1 );

//  ***************************************************************************
//  Perform  vertex degree reduction with local edge swapping: Soft.
//  ***************************************************************************

    err = vertex_degree_reduction();
    cout << "Info:  Stage "  << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces()    << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert( euler0 == euler1 );
    mesh->get_topological_statistics();

#ifdef DEBUG
    Jaal::set_no_tags(mesh);
    mesh->saveAs("stage5.dat");
#endif

//  ***************************************************************************
//  Perform some local operations: Condition: Soft. Diamonds must be called after the
//  vertex deduction operation, otherewise, face-close operation might increase the
//  vertex degrees.
//  ***************************************************************************

    err = remove_diamonds();
    cout << "Info:  Stage "  << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces()    << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert( euler0 == euler1 );

#ifdef DEBUG
    Jaal::set_no_tags(mesh);
    mesh->saveAs("stage6.dat");
#endif

//  ***************************************************************************
//  Last attempt to make all elements convex ..
//  ***************************************************************************

    if( mesh->count_concave_faces() == 0)
        mopt.shape_optimize(mesh);

    cout << "Info:  Stage "  << stage++ << " :  Passed " << endl;
    cout << "# of irregular nodes " << mesh->count_irregular_nodes(4) << endl;
    cout << "# of concave faces   " << mesh->count_concave_faces()    << endl;
    euler1 = mesh->getEulerCharacteristic();
    assert( euler0 == euler1 );

//  ***************************************************************************
//  Cluster all irregular nodes away from the boundaries.
//  ***************************************************************************

// ****************************************************************************
// Perform Global remeshing to elimate 3 and 5 degree nodes ..
// ****************************************************************************
    err = remesh_defective_patches();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_regular_node_tag(Mesh *mesh)
{
    int relexist = mesh->build_relations(0, 2);

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        int nsize = vertex->getRelations2().size();
        if (nsize >= 6) vertex->setTag(0);
        if (nsize == 4) vertex->setTag(1);
        if (nsize <  4) vertex->setTag(2);
        if (nsize == 5) vertex->setTag(3);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::reduce_boundary_vertex_degree( Vertex *vertex )
{
  if (!vertex->isBoundary()) return 1;

  NodeSequence &vneighs = vertex->getRelations0();
  int vdegree  = vneighs.size();

  if ( vdegree <=  (vertex->get_ideal_vertex_degree( Face::QUADRILATERAL)+1) ) return 2;

  sort(vneighs.begin(), vneighs.end(), HighVertexDegreeCompare());

  for (int k = 0; k < vdegree; k++) {
       const NodeSequence &wneighs = vneighs[k]->getRelations0();
       if (!vneighs[k]->isBoundary() && wneighs.size() > 3) {
           SwapQuadEdge edge(mesh, vertex, vneighs[k]);
           int err = edge.apply_bound_rule();
           if (!err) return 0;
       }
   }
   return 3;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::boundary_vertex_degree_reduction_once()
{
    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();

    // Precondition : The mesh must be doublet free ...
    vector<Doublet> doublets = search_interior_doublets();
    assert(doublets.empty());

    size_t numnodes = mesh->getSize(0);

    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        int err =  reduce_boundary_vertex_degree( vertex );
        if( !err) ncount++;
    }

    if (ncount)
        cout << "Info: number of boundary edges swapped " << ncount << endl;

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    // Post-Condition : The mesh must be doublet free ...
    doublets = search_interior_doublets();
    assert(doublets.empty());

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////
int
QuadCleanUp::reduce_internal_vertex_degree( Vertex *vertex )
{
  if (vertex->isBoundary()) return 1;

  NodeSequence &vneighs = vertex->getRelations0();
  int nSize = vneighs.size();

  if( nSize <= 4) return 2;

  sort(vneighs.begin(), vneighs.end(), HighVertexDegreeCompare());

  for (int k = 0; k < nSize; k++) {
       const NodeSequence &wneighs = vneighs[k]->getRelations0();
       if (wneighs.size() > 3) {
           SwapQuadEdge edge(mesh, vertex, vneighs[k]);
           int err = edge.apply_reduce_degree_rule();
           if (!err) return 0;
       }
  }

  return 3;
}

///////////////////////////////////////////////////////////////////////////////

int 
QuadCleanUp::internal_vertex_degree_reduction_once()
{
    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();

    vector<Doublet> doublets = search_interior_doublets();
    assert(doublets.empty());

    size_t numnodes = mesh->getSize(0);

    NodeSequence nodes = mesh->getNodes();
    sort(nodes.begin(), nodes.end(), HighVertexDegreeCompare());

    size_t ncount = 0;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = nodes[i];
        int err = reduce_internal_vertex_degree( vertex );
        if( !err ) ncount++;
   }

    if (ncount)
      cout << "Info: number of internal edges swapped " << ncount << endl;

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    doublets = search_interior_doublets();
    assert(doublets.empty());

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::vertex_degree_reduction()
{

    int ncount1 = 0;
    int ncount2 = 0;
    while (1) {
        ncount1 = 0;
        ncount2 = 0;
        ncount1 = boundary_vertex_degree_reduction_once();
        ncount2 = internal_vertex_degree_reduction_once();
        int err = lapsmooth->execute();
        if (err || (ncount1 + ncount2) == 0) break;
    }

    if( mesh->count_concave_faces() == 0)
        mopt.shape_optimize(mesh);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

FaceSequence
QuadCleanUp::search_flat_quads()
{
    size_t numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    FaceSequence flatQ, neighs;

    int edgefaces[4];
    for (size_t iface = 0; iface < numfaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
        assert(face);
        if (face->getSize(0) == 4) {
            int boundary = 1;
            for (int j = 0; j < 4; j++) {
                Vertex *vb = face->getNodeAt(j);
                if (!vb->isBoundary()) {
                    boundary = 0;
                    break;
                }
            }

            if (boundary) {
                int bound_edges = 0;
                for (int j = 0; j < 4; j++) {
                    Vertex *v0 = face->getNodeAt((j + 0) % 4);
                    Vertex *v1 = face->getNodeAt((j + 1) % 4);
                    neighs = Mesh::getRelations112(v0, v1);
                    if (neighs.size() == 1)
                        bound_edges++;
                    edgefaces[j] = neighs.size();
                }

                if (bound_edges == 3) {
                    /*
                     Point3D v1 = make_vector( neighs[0], node );
                     Point3D v2 = make_vector( neighs[1], node );
                     double  angle = getVectorAngle(v1,v2);
                     if( angle > cutOffAngle )
                     degree2nodes.push_back(node);
                     flatQ.push_back(face);
                     */
                }

            }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    cout << "Number of flat Quads " << flatQ.size() << endl;
    return flatQ;
}

///////////////////////////////////////////////////////////////////////////////

FaceSequence
QuadCleanUp::search_restricted_faces()
{
    size_t numfaces = mesh->getSize(2);

    mesh->search_boundary();

    FaceSequence restricted_faces;

    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        face->setTag(0);
        Vertex *v0 = face->getNodeAt(0);
        Vertex *v1 = face->getNodeAt(1);
        Vertex *v2 = face->getNodeAt(2);
        Vertex *v3 = face->getNodeAt(3);
        if ((v0->isBoundary() && v2->isBoundary()) || (v1->isBoundary() && v3->isBoundary())) {
            restricted_faces.push_back(face);
            face->setTag(1);
        }
    }

    cout << "Info: Number of restricted faces detected : " << restricted_faces.size() << endl;

    return restricted_faces;
}

////////////////////////////////////////////////////////////////////////////////

int
RestrictedEdge::build()
{
    assert(connect[0] != NULL);
    assert(connect[1] != NULL);

    Vertex *resnode = connect[0];
    Vertex *bndnode = connect[1];

    assert(!resnode->isBoundary());
    assert(bndnode->isBoundary());

    FaceSequence vneighs = bndnode->getRelations2();
    if (vneighs.size() != 2) return 2;

    adjFaces[0] = vneighs[0];
    adjFaces[1] = vneighs[1];

    Face *face0 = vneighs[0];
    if (face0->isRemoved()) return 3;

    int pos0 = face0->getPosOf(bndnode);
    assert(pos0 >= 0);
    Vertex *v0 = face0->getNodeAt((pos0 + 2) % 4);

    Face *face1 = vneighs[1];
    if (face1->isRemoved()) return 4;
    int pos1 = face1->getPosOf(bndnode);
    assert(pos1 >= 0);
    Vertex *v2 = face1->getNodeAt((pos1 + 2) % 4);

    double area_before = face0->getArea() + face1->getArea();

    // One new node will be inserted ::
    Vertex *vnew = Vertex::newObject();
    Point3D xyz = Vertex::mid_point(v0, v2);
    vnew->setXYZCoords(xyz);
    newNodes.resize(1);
    newNodes[0] = vnew;

    // Three new faces will be created ...
    newFaces.resize(3);

    Face *newface;
    NodeSequence connect(4);

    // Face : 1
    connect[0] = resnode;
    connect[1] = v0;
    connect[2] = vnew;
    connect[3] = v2;

    newface = Face::newObject();
    newface->setNodes(connect);
    newFaces[0] = newface;

    // Face : 2
    Vertex *v3 = NULL;
    if (face0->getNodeAt((pos0 + 1) % 4) == resnode)
        v3 = face0->getNodeAt((pos0 + 3) % 4);

    if (face0->getNodeAt((pos0 + 3) % 4) == resnode)
        v3 = face0->getNodeAt((pos0 + 1) % 4);

    assert(v3);
    connect[0] = vnew;
    connect[1] = v0;
    connect[2] = v3;
    connect[3] = bndnode;

    newface = Face::newObject();
    newface->setNodes(connect);
    newFaces[1] = newface;

    // Face : 3
    Vertex *v4 = NULL;
    if (face1->getNodeAt((pos1 + 1) % 4) == resnode)
        v4 = face1->getNodeAt((pos1 + 3) % 4);

    if (face1->getNodeAt((pos1 + 3) % 4) == resnode)
        v4 = face1->getNodeAt((pos1 + 1) % 4);

    assert(v4);
    connect[0] = vnew;
    connect[1] = bndnode;
    connect[2] = v4;
    connect[3] = v2;

    newface = Face::newObject();
    newface->setNodes(connect);
    newFaces[2] = newface;

    double area_after = 0.0;
    for (int i = 0; i < 3; i++)
        area_after += newFaces[i]->getArea();

    if (fabs(area_after - area_before) > 1.0E-10) {
        delete newFaces[0];
        delete newFaces[1];
        delete newFaces[2];
        newFaces.clear();
        return 4;
    }

    for (int i = 0; i < 3; i++) {
        if (!newFaces[i]->isConvex()) {
            delete newFaces[0];
            delete newFaces[1];
            delete newFaces[2];
            newFaces.clear();
            return 5;
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::cleanup_internal_boundary_face()
{
    int relexist = mesh->build_relations(0, 2);

    size_t numfaces = mesh->getSize(2);

    FaceSequence boundfaces;

    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        if (face->hasBoundaryNode()) {
            int nfnodes = face->getSize(0);
            Vertex *v0 = NULL;
            Vertex *v1 = NULL;
            for (int j = 0; j < nfnodes; j++) {
                Vertex *v = face->getNodeAt(j);
                if (v->isBoundary()) {
                    v0 = face->getNodeAt((j + 1) % nfnodes);
                    v1 = face->getNodeAt((j + 3) % nfnodes);
                    if (!v0->isBoundary() && !v1->isBoundary()) {
                        FaceClose closeface(mesh, face, v0, v1);
                        break;
                    }
                }
            }
        }
    }

    mesh->prune();

    if (!relexist)
        mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////

int QuadCleanUp::refine_degree3_faces()
{
    int relexist2 = mesh->build_relations(0, 2);

    size_t numnodes = mesh->getSize(0);

    Vertex *vertex;
    FaceSequence vfaces;
    for (size_t i = 0; i < numnodes; i++) {
        vertex = mesh->getNodeAt(i);
        vfaces = vertex->getRelations2();
        int nsize = vfaces.size();
        if ( nsize == 3) {
            for (int j = 0; j < nsize; j++) {
                vertex = vfaces[j]->getNodeAt(0);
                if (vertex->getRelations2().size() > 4) continue;

                vertex = vfaces[j]->getNodeAt(1);
                if (vertex->getRelations2().size() > 4) continue;

                vertex = vfaces[j]->getNodeAt(2);
                if (vertex->getRelations2().size() > 4) continue;

                vertex = vfaces[j]->getNodeAt(3);
                if (vertex->getRelations2().size() > 4) continue;

                mesh->refine_quad15(vfaces[j]);
            }
        }
    }
    mesh->prune();
    mesh->collect_garbage();

    if (!relexist2)
        mesh->clear_relations(0, 2);

    return 0;
}

///////////////////////////////////////////////////////////////////////

void QuadCleanUp::report()
{
    if (mesh == NULL) return;
    cout << "Info: Reporting mesh information " << endl;

}

///////////////////////////////////////////////////////////////////////
void QuadCleanUp :: build_irregular_nodes_set()
{
    // Get all the irregular nodes from the mesh ( only the internal ones );
    irregular_nodes_set.clear();
    NodeSequence nset = mesh->get_irregular_nodes(4 );
    for( size_t i = 0; i < nset.size(); i++)
        irregular_nodes_set.insert( nset[i] );

}
///////////////////////////////////////////////////////////////////////

int QuadCleanUp :: remove_tunnels()
{
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

#ifdef REMOVE_LATER

int QuadCleanUp::has_interior_nodes_degree_345()
{
    int relexist = mesh->build_relations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    size_t numnodes = mesh->getSize(0);

    vector<int> degree(numnodes);
    FaceSequence neighs;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        neighs = v->getRelations2();
        if( neighs.size() < 3 || neighs.size() > 5 ) return 0;
    }

    return 1;
}

int
QuadCleanUp::free_restricted_nodes_once()
{
    vector<Vertex*> vrestrict = search_restricted_nodes ();

    if (vrestrict.size () == 0) return 0;

    int relexist0 = mesh->build_relations (0, 0);
    int relexist2 = mesh->build_relations (0, 2);

    vector<Vertex*> vneighs;
    vector<RestrictedEdge> redges;

    for (size_t i = 0; i < vrestrict.size (); i++) {
        Vertex *vertex = vrestrict[i];
        vneighs = vertex->getRelations0 ();
        for (size_t j = 0; j < vneighs.size (); j++) {
            if (vneighs[j]->isBoundary ()) {
                RestrictedEdge edge (mesh, vertex, vneighs[j]);
                int err = edge.build ();
                if (!err) {
                    redges.push_back (edge);
                }
            }
        }
    }

    if (redges.size ())
        cout << "Info:  # of valid restricted edges : " << redges.size () << endl;

    int ncount = 0;
    for (size_t i = 0; i < redges.size (); i++) {
        int err = redges[i].commit ();
        if (!err) ncount++;
    }

    if (ncount) {
        mesh->prune ();
        mesh->enumerate (0);
        mesh->enumerate (2);
        cout << "Info: # of restricted edges committed : " << ncount << endl;
    }

    if (!relexist0)
        mesh->clear_relations (0, 0);

    if (!relexist2)
        mesh->clear_relations (0, 2);

    lapsmooth->execute ();

    return ncount;
    return 1;
}

////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::free_restricted_nodes()
{
    int ncount;
    while (1) {
        ncount = free_restricted_nodes_once();
        if (ncount == 0) break;
    }

    if (!mesh->isConsistentlyOriented())
        mesh->makeConsistentlyOriented();
}

void
QuadCleanUp::cleanup_boundary(double cutOffAngle)
{
    int rel0exist = mesh->build_relations(0, 0);
    int rel2exist = mesh->build_relations(0, 2);

    mesh->search_boundary();
    mesh->setFeatureLength();

    Point3D pf, pb, pn;
    double flen, len, t;

    size_t numnodes = mesh->getSize(0);

    NodeSequence bound_neighs, free_neighs;
    FaceSequence vfaces;

    NodeSequence updated;

    for (int iter = 0; iter < 10; iter++) {
        updated.clear();

        for (size_t i = 0; i < numnodes; i++) {
            Vertex *bndnode = mesh->getNodeAt(i);
            if (bndnode->isBoundary()) {
                bound_neighs = bndnode->getRelations0();
                for (size_t j = 0; j < bound_neighs.size(); j++) {
                    Vertex *freenode = bound_neighs[j];
                    if (!freenode->isBoundary()) {
                        free_neighs = freenode->getRelations0();
                        int ncount = 0;
                        for (size_t k = 0; k < free_neighs.size(); k++)
                            if (free_neighs[k]->isBoundary()) ncount++;
                        if (ncount == 1) {
                            flen = bndnode->getFeatureLength();
                            len = Vertex::length(bndnode, freenode);
                            if (len > 2.0 * flen) {
                                pf = freenode->getXYZCoords();
                                pb = bndnode->getXYZCoords();
                                t = 0.75;
                                pn[0] = t * pf[0] + (1 - t) * pb[0];
                                pn[1] = t * pf[1] + (1 - t) * pb[1];
                                pn[2] = t * pf[2] + (1 - t) * pb[2];
                                freenode->setXYZCoords(pn);
                                vfaces = freenode->getRelations2();
                                bool pass = 1;
                                for (size_t iface = 0; iface < vfaces.size(); iface++) {
                                    if (!vfaces[iface]->isConvex()) {
                                        pass = 0.0;
                                        break;
                                    }
                                }
                                if (!pass) {
                                    freenode->setXYZCoords(pf);
                                } else
                                    updated.push_back(freenode);
                            }
                        }
                    }
                }
            }
        }

        if (updated.empty()) break;

        for (size_t i = 0; i < updated.size(); i++)
            updated[i]->setConstrainedMark(1);

        lapsmooth->execute();

        for (size_t i = 0; i < updated.size(); i++)
            updated[i]->setConstrainedMark(0);

    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);

    return;



    ///////////////////////////////////////////////////////////////
    // First try to handle flat node...
    ///////////////////////////////////////////////////////////////

    NodeSequence degree2nodes;


    Vertex* node, *onode;
    for (size_t i = 0; i < numnodes; i++) {
        node = mesh->getNodeAt(i);
        if (node->isBoundary()) {
            neighs = node->getRelations0();
            if (neighs.size() == 2) {
                if (neighs[0]->isBoundary() && neighs[1]->isBoundary()) {
                    Point3D v1 = make_vector(neighs[0], node);
                    Point3D v2 = make_vector(neighs[1], node);
                    double angle = Math::getVectorAngle(v1, v2);
                    if (angle > cutOffAngle)
                        degree2nodes.push_back(node);
                }
            }
        }
    }


    relexist = mesh->build_relations(0, 2);

    Face *boundface;
    FaceSequence faceneighs;
    for (size_t i = 0; i < degree2nodes.size(); i++) {
        node = degree2nodes[i];
        faceneighs = node->getRelations2();
        if (faceneighs.size() == 1) {
            boundface = faceneighs[0];
            if (boundface->getSize(0) == 4) {
                int j = boundface->getPosOf(node);
                onode = boundface->getNodeAt((j + 2) % 4);
                if (!onode->isBoundary())
                    insert_doublet(boundface, node, onode);
            }
        }
    }

    mesh->prune();
    mesh->enumerate(0);
    mesh->enumerate(2);

    if (!relexist)
        mesh->clear_relations(0, 2);

    mesh->setBoundaryKnown(0);
}

#endif

///////////////////////////////////////////////////////////////////////

/*
int QuadCleanUp::apply_shift_node3_rule(Vertex *vertex)
{
    int layerID = vertex->getLayerID();

    FaceSequence vfaces = vertex->getRelations2();
    if (vfaces.size() != 3) return 1;

    Face *face = NULL;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        if (vfaces[i]->getLayerID() == layerID)
        {
            face = vfaces[i];
            break;
        }
    }
    if (face == NULL) return 1;

    int pos = face->getPosOf(vertex);
    assert(pos >= 0);

    Vertex *opp_vertex = face->getNodeAt((pos + 2) % 4);

    vfaces = opp_vertex->getRelations2();
    if (vfaces.size() == 4)
        return refine_3444_pattern(face, pos);

    if (vfaces.size() == 5)
        return refine_3454_pattern(face, pos);

    if (vfaces.size() > 5)
    {
        while (1)
        {
            int err = apply_advance_front_excess_rule(opp_vertex);
            if (err) break;
        }
        vfaces = opp_vertex->getRelations2();
        if (vfaces.size() < 6)
            return apply_shift_node3_rule(vertex);
    }

    return 2;
}

///////////////////////////////////////////////////////////////////////

int QuadCleanUp::refine_3454_pattern(const Face *face, int pos)
{
    Point3D xyz;
    NodeSequence localnodes(13);

    localnodes[0] = face->getNodeAt((pos + 0) % 4);
    localnodes[1] = face->getNodeAt((pos + 1) % 4);
    localnodes[2] = face->getNodeAt((pos + 2) % 4);
    localnodes[3] = face->getNodeAt((pos + 3) % 4);

    int layerid = localnodes[0]->getLayerID();

    if (localnodes[1]->getLayerID() != layerid) return 1;
    if (localnodes[2]->getLayerID() <= layerid) return 1;
    if (localnodes[3]->getLayerID() != layerid) return 1;

    FaceSequence vfaces;
    vfaces = localnodes[0]->getRelations2();
    if (vfaces.size() != 3) return 1;

    vfaces = localnodes[1]->getRelations2();
    if (vfaces.size() != 4) return 1;

    vfaces = localnodes[2]->getRelations2();
    if (vfaces.size() != 5) return 1;

    vfaces = localnodes[3]->getRelations2();
    if (vfaces.size() != 4) return 1;

    Face *neigh1 = NULL;
    Face *neigh2 = NULL;

    FaceSequence adjFaces;

    adjFaces = Mesh::getRelations112(localnodes[1], localnodes[2]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[4] = Face::opposite_node(neigh1, localnodes[2]);
    localnodes[5] = Face::opposite_node(neigh1, localnodes[1]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[5]);
    localnodes[11] = Vertex::newObject();
    localnodes[11]->setXYZCoords(xyz);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[5]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == neigh1) neigh2 = adjFaces[1];
    if (adjFaces[1] == neigh1) neigh2 = adjFaces[0];
    localnodes[8] = Face::opposite_node(neigh2, localnodes[2]);
    localnodes[9] = Face::opposite_node(neigh2, localnodes[5]);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[3]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[6] = Face::opposite_node(neigh1, localnodes[3]);
    localnodes[7] = Face::opposite_node(neigh1, localnodes[2]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[6]);
    localnodes[12] = Vertex::newObject();
    localnodes[12]->setXYZCoords(xyz);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[6]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == neigh1) neigh2 = adjFaces[1];
    if (adjFaces[1] == neigh1) neigh2 = adjFaces[0];
    localnodes[10] = Face::opposite_node(neigh2, localnodes[2]);

    Face * newfaces[7];
    NodeSequence connect(4);

    connect[0] = localnodes[0];
    connect[1] = localnodes[1];
    connect[2] = localnodes[11];
    connect[3] = localnodes[2];
    newfaces[0] = Face::newObject();
    newfaces[0]->setNodes(connect);

    connect[0] = localnodes[1];
    connect[1] = localnodes[4];
    connect[2] = localnodes[5];
    connect[3] = localnodes[11];
    newfaces[1] = Face::newObject();
    newfaces[1]->setNodes(connect);

    connect[0] = localnodes[11];
    connect[1] = localnodes[5];
    connect[2] = localnodes[8];
    connect[3] = localnodes[9];
    newfaces[2] = Face::newObject();
    newfaces[2]->setNodes(connect);

    connect[0] = localnodes[0];
    connect[1] = localnodes[2];
    connect[2] = localnodes[12];
    connect[3] = localnodes[3];
    newfaces[3] = Face::newObject();
    newfaces[3]->setNodes(connect);

    connect[0] = localnodes[3];
    connect[1] = localnodes[12];
    connect[2] = localnodes[6];
    connect[3] = localnodes[7];
    newfaces[4] = Face::newObject();
    newfaces[4]->setNodes(connect);

    connect[0] = localnodes[6];
    connect[1] = localnodes[12];
    connect[2] = localnodes[9];
    connect[3] = localnodes[10];
    newfaces[5] = Face::newObject();
    newfaces[5]->setNodes(connect);

    connect[0] = localnodes[2];
    connect[1] = localnodes[11];
    connect[2] = localnodes[9];
    connect[3] = localnodes[12];
    newfaces[6] = Face::newObject();
    newfaces[6]->setNodes(connect);

    // Do backup of Coordinates. Only 11 nodes to be backed up.
    Point3D backupCoords[11];
    for (int i = 0; i < 11; i++)
        backupCoords[i] = localnodes[i]->getXYZCoords();

    Face * backupFaces[5];
    vfaces = localnodes[2]->getRelations2();
    for (int i = 0; i < 5; i++)
    {
        backupFaces[i] = vfaces[i];
        mesh->remove(vfaces[i]); // Goes to garbage but not deallocated.
    }

    // Update new nodes and faces ...
    mesh->addNode(localnodes[11]);
    mesh->addNode(localnodes[12]);
    for (int i = 0; i < 7; i++)
        mesh->addFace(newfaces[i]);

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(localnodes);

    set<Face*> faces_to_check;
    for (int i = 0; i < 13; i++)
    {
        FaceSequence vfaces = localnodes[i]->getRelations2();
        for (int j = 0; j < vfaces.size(); j++)
            faces_to_check.insert(vfaces[j]);
    }

    bool pass = 1;
    set<Face*>::const_iterator siter;
    for (siter = faces_to_check.begin(); siter != faces_to_check.end(); ++siter)
    {
        Face *f = *siter;
        if (f->invertedAt() >= 0)
        {
            pass = 0;
            break;
        }
    }

    if (!pass)
    {
        for (int i = 0; i < 7; i++)
            mesh->remove(newfaces[i]);

        for (int i = 0; i < 5; i++)
            mesh->addFace(backupFaces[i]); // Reactivated now ...

        for (int i = 0; i < 11; i++)
            localnodes[i]->setXYZCoords(backupCoords[i]);

        mesh->remove(localnodes[11]);
        mesh->remove(localnodes[12]);
        return 1;
    }

    // Update front levels ...
    localnodes[11]->setLayerID(layerid + 1);
    localnodes[12]->setLayerID(layerid + 1);

    newfaces[0]->setLayerID(layerid);
    newfaces[1]->setLayerID(layerid);
    newfaces[2]->setLayerID(layerid + 1);

    newfaces[3]->setLayerID(layerid);
    newfaces[4]->setLayerID(layerid);
    newfaces[5]->setLayerID(layerid + 1);

    newfaces[6]->setLayerID(layerid + 1);

    return 0;
}

///////////////////////////////////////////////////////////////////////

int QuadCleanUp::refine_3444_pattern(const Face *face, int pos)
{
    Point3D xyz;
    NodeSequence localnodes(12);

    localnodes[0] = face->getNodeAt((pos + 0) % 4);
    localnodes[1] = face->getNodeAt((pos + 1) % 4);
    localnodes[2] = face->getNodeAt((pos + 2) % 4);
    localnodes[3] = face->getNodeAt((pos + 3) % 4);

    int layerid = localnodes[0]->getLayerID();

    if (localnodes[1]->getLayerID() != layerid) return 1;
    if (localnodes[2]->getLayerID() <= layerid) return 1;
    if (localnodes[3]->getLayerID() != layerid) return 1;

    FaceSequence vfaces;
    vfaces = localnodes[0]->getRelations2();
    if (vfaces.size() != 3) return 1;

    vfaces = localnodes[1]->getRelations2();
    if (vfaces.size() != 4) return 1;

    vfaces = localnodes[2]->getRelations2();
    if (vfaces.size() != 4) return 1;

    vfaces = localnodes[3]->getRelations2();
    if (vfaces.size() != 4) return 1;

    Face *neigh1 = NULL;
    Face *neigh2 = NULL;

    FaceSequence adjFaces;

    adjFaces = Mesh::getRelations112(localnodes[1], localnodes[2]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[4] = Face::opposite_node(neigh1, localnodes[2]);
    localnodes[5] = Face::opposite_node(neigh1, localnodes[1]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[5]);
    localnodes[9] = Vertex::newObject();
    localnodes[9]->setXYZCoords(xyz);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[5]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == neigh1) neigh2 = adjFaces[1];
    if (adjFaces[1] == neigh1) neigh2 = adjFaces[0];
    localnodes[8] = Face::opposite_node(neigh2, localnodes[2]);

    adjFaces = Mesh::getRelations112(localnodes[2], localnodes[3]);
    if (adjFaces.size() != 2) return 1;
    if (adjFaces[0] == face) neigh1 = adjFaces[1];
    if (adjFaces[1] == face) neigh1 = adjFaces[0];
    localnodes[6] = Face::opposite_node(neigh1, localnodes[3]);
    localnodes[7] = Face::opposite_node(neigh1, localnodes[2]);
    xyz = Vertex::mid_point(localnodes[2], localnodes[6]);
    localnodes[11] = Vertex::newObject();
    localnodes[11]->setXYZCoords(xyz);

    xyz = Vertex::mid_point(localnodes[2], localnodes[8]);
    localnodes[10] = Vertex::newObject();
    localnodes[10]->setXYZCoords(xyz);

    Face * newfaces[7];
    NodeSequence connect(4);

    connect[0] = localnodes[0];
    connect[1] = localnodes[1];
    connect[2] = localnodes[9];
    connect[3] = localnodes[2];
    newfaces[0] = Face::newObject();
    newfaces[0]->setNodes(connect);

    connect[0] = localnodes[1];
    connect[1] = localnodes[4];
    connect[2] = localnodes[5];
    connect[3] = localnodes[9];
    newfaces[1] = Face::newObject();
    newfaces[1]->setNodes(connect);

    connect[0] = localnodes[9];
    connect[1] = localnodes[5];
    connect[2] = localnodes[8];
    connect[3] = localnodes[10];
    newfaces[2] = Face::newObject();
    newfaces[2]->setNodes(connect);

    connect[0] = localnodes[0];
    connect[1] = localnodes[2];
    connect[2] = localnodes[11];
    connect[3] = localnodes[3];
    newfaces[3] = Face::newObject();
    newfaces[3]->setNodes(connect);

    connect[0] = localnodes[3];
    connect[1] = localnodes[11];
    connect[2] = localnodes[6];
    connect[3] = localnodes[7];
    newfaces[4] = Face::newObject();
    newfaces[4]->setNodes(connect);

    connect[0] = localnodes[6];
    connect[1] = localnodes[11];
    connect[2] = localnodes[10];
    connect[3] = localnodes[8];
    newfaces[5] = Face::newObject();
    newfaces[5]->setNodes(connect);

    connect[0] = localnodes[2];
    connect[1] = localnodes[9];
    connect[2] = localnodes[10];
    connect[3] = localnodes[11];
    newfaces[6] = Face::newObject();
    newfaces[6]->setNodes(connect);

    Point3D backupCoords[9];
    for (int i = 0; i < 9; i++)
        backupCoords[i] = localnodes[i]->getXYZCoords();

    Face * backupFaces[4];
    vfaces = localnodes[2]->getRelations2();
    for (int i = 0; i < 4; i++)
    {
        backupFaces[i] = vfaces[i];
        mesh->remove(vfaces[i]); // Send to garbage, but don't delete ..
    }

    // Update the data structures ...
    mesh->addNode(localnodes[9]);
    mesh->addNode(localnodes[10]);
    mesh->addNode(localnodes[11]);
    for (int i = 0; i < 7; i++)
        mesh->addFace(newfaces[i]);

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(localnodes);

    FaceSet faces_to_check;
    for (int i = 0; i < 12; i++)
    {
        FaceSequence vfaces = localnodes[i]->getRelations2();
        for (int j = 0; j < vfaces.size(); j++)
            faces_to_check.insert(vfaces[j]);
    }

    bool pass = 1;
    set<Face*>::const_iterator siter;
    for (siter = faces_to_check.begin(); siter != faces_to_check.end(); ++siter)
    {
        Face *f = *siter;
        if (f->invertedAt() >= 0)
        {
            pass = 0;
            break;
        }
    }

    if (!pass)
    {
        for (int i = 0; i < 7; i++)
            mesh->remove(newfaces[i]);

        for (int i = 0; i < 4; i++)
            mesh->addFace(backupFaces[i]);

        for (int i = 0; i < 9; i++)
            localnodes[i]->setXYZCoords(backupCoords[i]);

        mesh->remove(localnodes[9]);
        mesh->remove(localnodes[10]);
        mesh->remove(localnodes[11]);

        return 1;
    }

    // Update front levels ...
    localnodes[9]->setLayerID(layerid + 1);
    localnodes[10]->setLayerID(layerid + 2);
    localnodes[11]->setLayerID(layerid + 1);

    newfaces[0]->setLayerID(layerid);
    newfaces[1]->setLayerID(layerid);
    newfaces[2]->setLayerID(layerid + 1);

    newfaces[3]->setLayerID(layerid);
    newfaces[4]->setLayerID(layerid);
    newfaces[5]->setLayerID(layerid + 1);

    newfaces[6]->setLayerID(layerid + 1);

    return 0;
}
*/

