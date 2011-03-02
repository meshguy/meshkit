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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        int nsize = vertex->getRelations2().size();
        if (nsize >= 6) vertex->setTag(0);
        if (nsize == 4) vertex->setTag(1);
        if (nsize < 4) vertex->setTag(2);
        if (nsize == 5) vertex->setTag(3);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);
}

///////////////////////////////////////////////////////////////////////////////

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
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
        {
           NodeSequence vneighs = vertex->getRelations0();
           int vdegree  = vneighs.size();
           if ( vdegree > vertex->get_ideal_vertex_degree() )
            {
                for (int k = 0; k < vdegree; k++)
                {
                    NodeSequence wneighs = vneighs[k]->getRelations0();
                    if (!vneighs[k]->isBoundary() && wneighs.size() > 3)
                    {
                        SwapQuadEdge edge(mesh, vertex, vneighs[k]);
                        int err = edge.apply_bound_rule();
                        if (!err)
                        {
                           ncount++;
                           break;
                        }
                    }
               }
            }
        }
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

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = nodes[i];
        if (!vertex->isBoundary())
        {
            NodeSequence vneighs = vertex->getRelations0();
            if (vneighs.size() > 4)
            {
                sort(vneighs.begin(), vneighs.end(), HighVertexDegreeCompare());
                for (int k = 0; k < vneighs.size(); k++)
                {
                    NodeSequence wneighs = vneighs[k]->getRelations0();
                    if (wneighs.size() > 3)
                    {
                        SwapQuadEdge edge(mesh, vertex, vneighs[k]);
                        int err = edge.apply_reduce_degree_rule();
                        if (!err)
                        {
                            ncount++;
                            break;
                        }
                    }
                }
            }

        }
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
    while (1)
    {
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
    //  Public Function ...
    size_t numfaces = mesh->getSize(2);

    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    FaceSequence flatQ, neighs;

    int edgefaces[4];
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        assert(face);
        if (face->getSize(0) == 4)
        {
            int boundary = 1;
            for (int j = 0; j < 4; j++)
            {
                Vertex *vb = face->getNodeAt(j);
                if (!vb->isBoundary())
                {
                    boundary = 0;
                    break;
                }
            }

            if (boundary)
            {
                int bound_edges = 0;
                for (int j = 0; j < 4; j++)
                {
                    Vertex *v0 = face->getNodeAt((j + 0) % 4);
                    Vertex *v1 = face->getNodeAt((j + 1) % 4);
                    neighs = Mesh::getRelations112(v0, v1);
                    if (neighs.size() == 1)
                        bound_edges++;
                    edgefaces[j] = neighs.size();
                }

                if (bound_edges == 3)
                {
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

    int ncount;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(0);
        Vertex *v0 = face->getNodeAt(0);
        Vertex *v1 = face->getNodeAt(1);
        Vertex *v2 = face->getNodeAt(2);
        Vertex *v3 = face->getNodeAt(3);
        if ((v0->isBoundary() && v2->isBoundary()) || (v1->isBoundary() && v3->isBoundary()))
        {
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

    if (fabs(area_after - area_before) > 1.0E-10)
    {
        delete newFaces[0];
        delete newFaces[1];
        delete newFaces[2];
        newFaces.clear();
        return 4;
    }

    for (int i = 0; i < 3; i++)
    {
        if (!newFaces[i]->isConvex())
        {
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

    /*
int
QuadCleanUp::free_restricted_nodes_once()
{
      vector<Vertex*> vrestrict = search_restricted_nodes ();

      if (vrestrict.size () == 0) return 0;

      int relexist0 = mesh->build_relations (0, 0);
      int relexist2 = mesh->build_relations (0, 2);

      vector<Vertex*> vneighs;
      vector<RestrictedEdge> redges;

      for (size_t i = 0; i < vrestrict.size (); i++)
        {
          Vertex *vertex = vrestrict[i];
          vneighs = vertex->getRelations0 ();
          for (size_t j = 0; j < vneighs.size (); j++)
            {
              if (vneighs[j]->isBoundary ())
                {
                  RestrictedEdge edge (mesh, vertex, vneighs[j]);
                  int err = edge.build ();
                  if (!err)
                    {
                      redges.push_back (edge);
                    }
                }
            }
        }

      if (redges.size ())
        cout << "Info:  # of valid restricted edges : " << redges.size () << endl;

      int ncount = 0;
      for (size_t i = 0; i < redges.size (); i++)
        {
          int err = redges[i].commit ();
          if (!err) ncount++;
        }

      if (ncount)
        {
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
    while (1)
    {
        ncount = free_restricted_nodes_once();
        if (ncount == 0) break;
    }

    if (!mesh->isConsistentlyOriented())
        mesh->makeConsistentlyOriented();
}
     */

////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::cleanup_internal_boundary_face()
{
    int relexist = mesh->build_relations(0, 2);

    size_t numfaces = mesh->getSize(2);

    FaceSequence boundfaces;

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if (face->hasBoundaryNode())
        {
            int nfnodes = face->getSize(0);
            Vertex *v0 = NULL;
            Vertex *v1 = NULL;
            for (int j = 0; j < nfnodes; j++)
            {
                Vertex *v = face->getNodeAt(j);
                if (v->isBoundary())
                {
                    v0 = face->getNodeAt((j + 1) % nfnodes);
                    v1 = face->getNodeAt((j + 3) % nfnodes);
                    if (!v0->isBoundary() && !v1->isBoundary())
                    {
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

#ifdef REMOVE_LATER

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

    for (int iter = 0; iter < 10; iter++)
    {
        updated.clear();

        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *bndnode = mesh->getNodeAt(i);
            if (bndnode->isBoundary())
            {
                bound_neighs = bndnode->getRelations0();
                for (size_t j = 0; j < bound_neighs.size(); j++)
                {
                    Vertex *freenode = bound_neighs[j];
                    if (!freenode->isBoundary())
                    {
                        free_neighs = freenode->getRelations0();
                        int ncount = 0;
                        for (size_t k = 0; k < free_neighs.size(); k++)
                            if (free_neighs[k]->isBoundary()) ncount++;
                        if (ncount == 1)
                        {
                            flen = bndnode->getFeatureLength();
                            len = Vertex::length(bndnode, freenode);
                            if (len > 2.0 * flen)
                            {
                                pf = freenode->getXYZCoords();
                                pb = bndnode->getXYZCoords();
                                t = 0.75;
                                pn[0] = t * pf[0] + (1 - t) * pb[0];
                                pn[1] = t * pf[1] + (1 - t) * pb[1];
                                pn[2] = t * pf[2] + (1 - t) * pb[2];
                                freenode->setXYZCoords(pn);
                                vfaces = freenode->getRelations2();
                                bool pass = 1;
                                for (size_t iface = 0; iface < vfaces.size(); iface++)
                                {
                                    if (!vfaces[iface]->isConvex())
                                    {
                                        pass = 0.0;
                                        break;
                                    }
                                }
                                if (!pass)
                                {
                                    freenode->setXYZCoords(pf);
                                }
                                else
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
    for (size_t i = 0; i < numnodes; i++)
    {
        node = mesh->getNodeAt(i);
        if (node->isBoundary())
        {
            neighs = node->getRelations0();
            if (neighs.size() == 2)
            {
                if (neighs[0]->isBoundary() && neighs[1]->isBoundary())
                {
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
    for (size_t i = 0; i < degree2nodes.size(); i++)
    {
        node = degree2nodes[i];
        faceneighs = node->getRelations2();
        if (faceneighs.size() == 1)
        {
            boundface = faceneighs[0];
            if (boundface->getSize(0) == 4)
            {
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

///////////////////////////////////////////////////////////////////////

int QuadCleanUp::refine_degree3_faces()
{
    int relexist2 = mesh->build_relations(0, 2);

    size_t numnodes = mesh->getSize(0);

    Vertex *vertex;
    FaceSequence vfaces;
    for (size_t i = 0; i < numnodes; i++)
    {
        vertex = mesh->getNodeAt(i);
        vfaces = vertex->getRelations2();
        if (vfaces.size() == 3)
        {
            for (int j = 0; j < vfaces.size(); j++)
            {
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

/*

int
SwapEdge::get_boundary_nodes_chain()
{
    vector<Edge> bndedges;
    bndedges.reserve(6);

    Edge sharededge(connect[0], connect[1]);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Vertex *vf0 = adjFaces[i]->getNodeAt((j + 0) % 4);
            Vertex *vf1 = adjFaces[i]->getNodeAt((j + 1) % 4);
            Edge edge(vf0, vf1);
            if (!edge.isSameAs(sharededge))
                bndedges.push_back(edge);
        }
    }

    Mesh::make_chain(bndedges);
    bound_nodes = Mesh::chain_nodes(bndedges);
    assert(bound_nodes.size() == 6);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void SwapEdge::update_front()
{
    int layerid = connect[0]->getLayerID();
    connect[1]->setLayerID(layerid + 1);
    adjFaces[0]->setLayerID(layerid);
    adjFaces[1]->setLayerID(layerid);
}

////////////////////////////////////////////////////////////////////////////////

void
SwapEdge::backup()
{
    bkp_data.diagonalConnect[0] = connect[0];
    bkp_data.diagonalConnect[1] = connect[1];
    bkp_data.face1Connect = adjFaces[0]->getNodes();
    bkp_data.face2Connect = adjFaces[1]->getNodes();
    for (int i = 0; i < 6; i++)
    {
        Vertex *v = bound_nodes[i];
        bkp_data.pCoords[v] = v->getXYZCoords();
    }
}
////////////////////////////////////////////////////////////////////////////////

void
SwapEdge::rollback()
{
    connect[0]->removeRelation0(connect[1]);
    connect[1]->removeRelation0(connect[0]);

    // Remove the vertex-face relationship.
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->removeRelation2(adjFaces[0]);
    }

    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->removeRelation2(adjFaces[1]);
    }

    // Change the connectivity of the two quads.
    adjFaces[0]->setNodes(bkp_data.face1Connect);
    adjFaces[1]->setNodes(bkp_data.face2Connect);
    connect[0] = bkp_data.diagonalConnect[0];
    connect[1] = bkp_data.diagonalConnect[1];

    // Update the vertex-vertex relationship.
    connect[0]->addRelation0(connect[1]);
    connect[1]->addRelation0(connect[0]);

    // Update the vertex-face relationships.
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->addRelation2(adjFaces[0]);
    }

    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->addRelation2(adjFaces[1]);
    }

    for (int i = 0; i < 6; i++)
    {
        Vertex *v = bound_nodes[i];
        v->setXYZCoords(bkp_data.pCoords[v]);
    }

    update_front();
}
////////////////////////////////////////////////////////////////////////////////

int
SwapEdge::make_new_diagonal_at(int pos, bool bound_check)
{
    //
    // Given closednodes ( exactly six nodes) in an order ( clockwise or
    // counter-clockwise) and create two quads having common diagonal between
    // (pos, pos+3).
    //
    assert(pos >= 0 && pos < 6);

    if (bound_check)
    {
        if (bound_nodes[(pos + 0) % 6]->isBoundary()) return 1;
        if (bound_nodes[(pos + 3) % 6]->isBoundary()) return 1;
    }

    // First perform full backup of the data structures.
    backup();

    // Remove the vertex-vertex relationship.
    connect[0]->removeRelation0(connect[1]);
    connect[1]->removeRelation0(connect[0]);

    // Remove the vertex-face relationship.
    assert(adjFaces[0]);
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->removeRelation2(adjFaces[0]);
    }

    assert(adjFaces[1]);
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->removeRelation2(adjFaces[1]);
    }

    // Change the connectivity of the two quads.
    vector<Vertex*> qConnect(4);
    qConnect[0] = bound_nodes[(pos + 0) % 6];
    qConnect[1] = bound_nodes[(pos + 1) % 6];
    qConnect[2] = bound_nodes[(pos + 2) % 6];
    qConnect[3] = bound_nodes[(pos + 3) % 6];
    connect[0] = qConnect[0];
    adjFaces[0]->setNodes(qConnect);

    qConnect[0] = bound_nodes[(pos + 3) % 6];
    qConnect[1] = bound_nodes[(pos + 4) % 6];
    qConnect[2] = bound_nodes[(pos + 5) % 6];
    qConnect[3] = bound_nodes[(pos + 6) % 6];
    connect[1] = qConnect[0];
    adjFaces[1]->setNodes(qConnect);

    // Update the vertex-vertex relationship.
    connect[0]->addRelation0(connect[1]);
    connect[1]->addRelation0(connect[0]);

    // Update the vertex-face relationships.
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[0]->getNodeAt(i);
        v->addRelation2(adjFaces[0]);
    }

    for (int i = 0; i < 4; i++)
    {
        Vertex *v = adjFaces[1]->getNodeAt(i);
        v->addRelation2(adjFaces[1]);
    }

    // The change may produce concave or tangled mesh. Do some local
    // Laplacian smoothing at the six nodes only. Hopefully, the elements
    // will be acceptable. If not, we need to rollback this operation.

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(bound_nodes);

    // Look at the neighbors, if some elements is inverted. Here we create
    // a set of neighbours, to avoid duplicate checking.
    set<Face*> neighSet;
    for (int i = 0; i < 6; i++)
    {
        Vertex *v = bound_nodes[i];
        if (!v->isBoundary())
        {
            vector<Face*> vfaces = v->getRelations2();
            for (int j = 0; j < vfaces.size(); j++)
                neighSet.insert(vfaces[j]);
        }
    }
    assert(!neighSet.empty());

    set<Face*>::const_iterator it;
    for (it = neighSet.begin(); it != neighSet.end(); ++it)
    {
        Face *face = *it;
        if (face->invertedAt() >= 0)
        {
            rollback();
            return 1;
        }
    }
    update_front();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
QuadEdge::commit()
{
    if (newFaces.empty()) return 1;

    for (int i = 0; i < newFaces.size(); i++)
        if (newFaces[i] == NULL) return 2;

    if (adjFaces[0]->isRemoved() || adjFaces[1]->isRemoved())
    {
        for (int i = 0; i < newNodes.size(); i++)
            delete newNodes[i];
        newNodes.clear();

        for (int i = 0; i < newFaces.size(); i++)
            delete newFaces[i];
        newFaces.clear();
        return 3;
    }

    Vertex *vertex;
    for (int i = 0; i < 4; i++)
    {
        vertex = adjFaces[0]->getNodeAt(i);
        if (vertex->isVisited()) return 4;
        vertex = adjFaces[1]->getNodeAt(i);
        if (vertex->isVisited()) return 4;
    }

    adjFaces[0]->setRemoveMark(1);
    adjFaces[1]->setRemoveMark(1);

    assert(mesh);

    for (int i = 0; i < newNodes.size(); i++)
        mesh->addNode(newNodes[i]);

    for (int i = 0; i < newFaces.size(); i++)
        mesh->addFace(newFaces[i]);

    for (int i = 0; i < 4; i++)
    {
        vertex = adjFaces[0]->getNodeAt(i);
        vertex->setVisitMark(1);
        vertex = adjFaces[1]->getNodeAt(i);
        vertex->setVisitMark(1);
    }

    return 0;
}
*/

/////////////////////////////////////////////////////////////////////

/*
int
SwapEdge::apply_reduce_degree_rule()
{
    if (build_boundary() != 0) return 1;

    vector<Vertex*> neighs;

    neighs = connect[0]->getRelations0();
    int d1 = neighs.size();
    if (d1 < 4) return 1;

    neighs = connect[1]->getRelations0();
    int d2 = neighs.size();
    if (d2 < 4) return 1;

    int start_pos = this->getPosOf(connect[0]);
    assert(start_pos >= 0);

    for (int i = 0; i < 2; i++)
    {
        int pos = (start_pos + i + 1) % 6;
        Vertex *v0 = bound_nodes[ pos ];
        if (v0->isBoundary()) continue;
        neighs = v0->getRelations0();
        int d3 = neighs.size();

        Vertex *v3 = bound_nodes[(pos + 3) % 6];
        if (v3->isBoundary()) continue;
        neighs = v3->getRelations0();
        int d4 = neighs.size();

        if (SwapEdge::is_topologically_valid_swap(d1, d2, d3, d4))
        {
            int err = make_new_diagonal_at(pos);
            if (!err) return 0;
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_concave_rule()
{
    if (build_boundary() != 0) return 1;

    int start_pos = getPosOf(connect[0]);
    for (int i = 0; i < 2; i++)
    {
        int err = make_new_diagonal_at((start_pos + 1 + i) % 6);
        if (!err) return 0;
    }
    return 1;
}

/////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_bound_rule()
{
    if (!connect[0]->isBoundary()) return 1;
    if (connect[1]->isBoundary()) return 1;

    if (build_boundary() != 0) return 2;

    for (int i = 0; i < 3; i++)
    {
        Vertex *v0 = bound_nodes[(i + 0) % 6];
        Vertex *v3 = bound_nodes[(i + 3) % 6];
        if ((!v0->isBoundary()) && (!v3->isBoundary()))
        {
            int err = make_new_diagonal_at(i, 0);
            if (!err) return 0;
        }
    }
    return 3;
}

/////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_advance_front_rule()
{
    if (build_boundary() != 0) return 1;

    int layerid = connect[0]->getLayerID();
    for (int i = 0; i < 3; i++)
    {
        Vertex *v0 = bound_nodes[(i + 0) % 6];
        Vertex *v3 = bound_nodes[(i + 3) % 6];
        if ((v0->getLayerID() > layerid) && (v3->getLayerID() > layerid))
        {
            int err = make_new_diagonal_at(i);
            if (!err) return 0;
        }
    }
    return 3;
}

/////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_singlet_rule(Vertex *singlet)
{
    assert(connect[0]->isBoundary());

    if (connect[1]->isBoundary()) return 1;

    if (build_boundary() != 0) return 1;
    ////////////////////////////////////////////////////////////////////////////
    // Objective :: Swap quads common diagonal such that the resulting diagonal
    //              contains the Singlet node.
    ////////////////////////////////////////////////////////////////////////////
    assert(QuadCleanUp::isSinglet(singlet));

    // First Vertex must be boundary

    if (adjFaces[0] == NULL || adjFaces[1] == NULL) return 2;

    assert(QuadCleanUp::hasSinglet(adjFaces[0]));

    // Make sure that other face doesn't have a singlet
    if (QuadCleanUp::hasSinglet(adjFaces[1])) return 3;

    //Find the position of singlet in the contour ...
    int pos = this->getPosOf(singlet);
    assert(pos >= 0);

    // Create new quads whose common diagonal must contain the singlet.
    return make_new_diagonal_at(pos, 0);
}

/////////////////////////////////////////////////////////////////////////////

int
SwapEdge::apply_deficient_rule(Vertex *vertex)
{
    assert(firstFace);
    if (vertex->isBoundary()) return 1;

    int d1 = connect[0]->getRelations2().size();
    int d2 = connect[1]->getRelations2().size();

    // Note 1: Just near the deficient node, there must be atleast 5 nodes in
    // order to swap the edge,otherwise, things may not converge, on the
    // same level, the deficient node will keep changing the location.
    //
    // Note 2: The node connect[1] is opposite to the deficient node and
    // most probably on the next level. If we keep moving the deficient
    // node inside the domain, then it is fine to have its degree to 4, so
    // that the deficient node will be generated there.
    vector<Face*> faces;
    if (d1 == 4 && d2 > 3)
    {
        int layerid = connect[0]->getLayerID();
        faces = vertex->getRelations2();
        if (faces.size() == 3)
        {
            Face *f0 = firstFace;
            faces = Mesh::getRelations112(connect[0], connect[1]);
            assert(faces.size() == 2);
            Face *f1 = NULL;
            if (faces[0] == f0) f1 = faces[1];
            if (faces[1] == f0) f1 = faces[0];
            assert(f1);
            int pos = f1->queryPosOf(connect[0]);
            Vertex *vopp = f1->getNodeAt((pos + 2) % 4);
            if (connect[1]->getLayerID() == layerid + 1 && vopp->getLayerID() == layerid + 1)
            {
                SwapEdge edge(mesh, connect[1], vopp, f1);
                int err = edge.apply_deficient_rule(connect[0]);
                if (err) return 1;
            }
        }
    }

    d1 = connect[0]->getRelations2().size();
    d2 = connect[1]->getRelations2().size();

    if (d1 < 5 && d2 < 3) return 2;

    // Having these conditions set, now build the structure.
    if (build_boundary() != 0)
    {
        return 3;
    }

    //Find the position of triplet in the contour ...
    int pos = this->getPosOf(vertex);
    assert(pos >= 0);

    // Create new quads whose common diagonal must contain the singlet.
    return make_new_diagonal_at(pos);
}
*/

/////////////////////////////////////////////////////////////////////////////

/*
int
QuadCleanUp::swap_concave_faces()
{
    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
    }

    int ncount = 0;
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        int pos = face->invertedAt();
        if (pos >= 0)
        {
            for (int i = 0; i < 2; i++)
            {
                Vertex *v0 = face->getNodeAt((pos + 1 + i) % 4);
                Vertex *v1 = face->getNodeAt((pos + 2 + i) % 4);
                SwapEdge edge(mesh, v0, v1);
                int err = edge.apply_concave_rule();
                if (!err)
                {
                    ncount++;
                    break;
                }
            }
        }
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: # of swapped edges " << ncount << endl;
    }

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    return ncount;
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
int
QuadCleanUp::apply_advance_front_singlet_rule(Vertex *singlet)
{
    int err;
    if (singlet->getLayerID() != 0) return 1;
    if (!QuadCleanUp::isSinglet(singlet)) return 1;

    vector<Face*> vfaces = singlet->getRelations2();
    assert(vfaces.size() == 1);

    Face *face = vfaces[0];
    int pos = face->queryPosOf(singlet);

    Vertex *v0 = face->getNodeAt((pos + 0) % 4);
    Vertex *v1 = face->getNodeAt((pos + 1) % 4);
    Vertex *v2 = face->getNodeAt((pos + 2) % 4);
    Vertex *v3 = face->getNodeAt((pos + 3) % 4);

    SwapEdge edge1(mesh, v1, v2);
    err = edge1.apply_singlet_rule(singlet);
    if (!err) return 0;

    SwapEdge edge2(mesh, v3, v2);
    err = edge2.apply_singlet_rule(singlet);
    if (!err) return 0;

    return 1;
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
int
QuadCleanUp::apply_advance_front_triplet_rule(Vertex *vertex)
{
    // ****************************************************************************
    // Implementation ideas:
    //
    // This is the case of 5-434 transition, Where 5 degree node is opposite
    // to 3 degree node and with first swap, it becomes 4-335 or 4-534 and with
    // the second swap, it becomes 3-444 and therefore, the 3 degree node moves
    // inside the domain.
    //
    // Having degree 5 is a somewhat strict condition, because degree four is
    // sufficient, but it will create doublet, and as per our design goals, we
    // do not want to introduce any new doublet in the mesh. In future, we may
    // change the conditions, if that is useful.
    //
    // ****************************************************************************
    if (vertex->isBoundary()) return 1;

    vector<Vertex*> vnodes = vertex->getRelations0();
    if (vnodes.size() != 3) return 1;

    int layerid = vertex->getLayerID();
    vector<Face*> vfaces = vertex->getRelations2();
    for (int j = 0; j < 3; j++)
    {
        Face *face = vfaces[j];
        if (face->getLayerID() == layerid)
        {
            int pos = face->queryPosOf(vertex);
            Vertex *v1 = face->getNodeAt((pos + 1) % 4);
            Vertex *v2 = face->getNodeAt((pos + 2) % 4);
            Vertex *v3 = face->getNodeAt((pos + 3) % 4);
            // If v1 and v3 are at the lower level, they are already processed.
            if (v1->getLayerID() < layerid || v3->getLayerID() < layerid) return 1;
            int d1 = v1->getRelations2().size();
            int d2 = v2->getRelations2().size();
            int d3 = v3->getRelations2().size();
            if (d3 > d1)
                std::swap(v1, v3);

            //  The opposite node must have degree at least five.
            if (max(d1, d3) <= 4 && d2 < 5) return 1;

            //  The opposite node must be at higher level than at (v0, v1, v3)
            if (v2->getLayerID() > layerid)
            {
                SwapEdge edge1(mesh, v1, v2, face);
                int err1 = edge1.apply_deficient_rule(vertex);
                if (!err1)
                {
                    return 0;
                }

                SwapEdge edge2(mesh, v3, v2, face);
                int err2 = edge2.apply_deficient_rule(vertex);
                if (!err2 == 0)
                {
                    return 0;
                }
            }
        }
    }

    return 1;

}
*/

///////////////////////////////////////////////////////////////////////////////

/*
int
QuadCleanUp::apply_advance_front_bridge_rule(Vertex *v0, Vertex *v1)
{
    int layerid = v0->getLayerID();
    if (v1->getLayerID() != layerid) return 1;

    // Although it is checked once again, but it is essential.
    vector<Face*> adjFaces;
    adjFaces = v0->getRelations2();
    if (adjFaces.size() != 3) return 2;

    adjFaces = v1->getRelations2();
    if (adjFaces.size() != 3) return 3;

    // Our assumption is that all the edges are simple.
    adjFaces = Mesh::getRelations112(v0, v1);
    if (adjFaces.size() != 2) return 4;

    if (adjFaces[0]->getLayerID() == adjFaces[1]->getLayerID()) return 5;

    // Check for the internal face.
    Face *internal_face = NULL;
    internal_face = adjFaces[0]->getLayerID() > adjFaces[1]->getLayerID() ? adjFaces[0] : adjFaces[1];

    // Swap the opposite edge:
    Vertex *v2, *v3;
    Face::opposite_nodes(internal_face, v0, v1, v2, v3);

    // Try swapping at the first node.
    int err;
    SwapEdge edge1(mesh, v2, v3, internal_face);
    err = edge1.apply_deficient_rule(v0);
    if (!err) return 0;

    // Try swapping at the second node.
    SwapEdge edge2(mesh, v2, v3, internal_face);
    err = edge2.apply_deficient_rule(v1);
    if (!err) return 0;

    return 1;

}
*/

///////////////////////////////////////////////////////////////////////////////

/*
int
QuadCleanUp::apply_advance_front_excess_rule(Vertex *vertex)
{
    if (vertex->isBoundary()) return 1;

    int layerid = vertex->getLayerID();
    vector<Vertex*> vneighs = vertex->getRelations0();
    int degree = vneighs.size();

    if (degree < 5) return 1;

    int ncount = 0;
    for (int k = 0; k < degree; k++)
        if (vneighs[k]->getLayerID() > layerid) ncount++;

    if (ncount < 2) return 1;

    for (int k = 0; k < degree; k++)
    {
        if (vneighs[k]->getLayerID() > layerid)
        {
            SwapEdge edge(mesh, vertex, vneighs[k]);
            int err = edge.apply_advance_front_rule();
            if (!err) return 0;
        }
    }

    return 1;
}
*/

///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
int
QuadCleanUp::advance_front_edges_swap_once(int layerid)
{
    int err;
    if (mesh->getAdjTable(0, 0))
        mesh->clear_relations(0, 0);

    if (mesh->getAdjTable(0, 2))
        mesh->clear_relations(0, 2);

    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    //
    // The input mesh should be doublet free and the output mesh will
    // be doublet free.
    //
    vector<Doublet> doublets = search_interior_doublets();
    assert(doublets.empty());

    // If the boundary is unknown ...
    if (layerid == 0)
        mesh->search_boundary();

    // We need atleast two layer to work with.
    size_t numnodes = mesh->getSize(0);

    // Check how many iiregular nodes in the present layer...
    size_t num_irregular_nodes = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitMark(0);
        if (vertex->getLayerID() == layerid)
        {
            if (vertex->getRelations2().size() != 4)
                num_irregular_nodes++;
        }
    }

    size_t ncount = 0;

    //
    // If the layer is boundary, then highest priority must be given to remove
    // the singlets. Many singlets can be removed by swapping the edges. There are
    // some cases, where "Swapping" may not remove singlets. There are two such
    // scenarios, which are handled by calling "remove_boundary_singlet" method.
    //

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
        {
            err = apply_advance_front_singlet_rule(vertex);
            if (!err) ncount++;
        }
    }

    // Second priority is given to bridges in the mesh. Removal of bridges
    // create a three degree node, which is the next case.

    if (ncount == 0)
    {
        /*
           vector<Edge33> bridges = search_edges33();
           for (size_t i = 0; i < bridges.size(); i++){
                Vertex *v0 = bridges[i].getNodeAt(0);
                Vertex *v1 = bridges[i].getNodeAt(1);
                if( v0->getLayerID() == layerid  && v1->getLayerID() == layerid ) {
                    err = apply_advance_front_bridge_rule(v0, v1);
                    if( !err ) ncount++;
                    v0->setTag(3);
                    v1->setTag(3);
                    mesh->saveAs("dbg.dat");
                    exit(0);
                }
           }
         */
    }

    // Third Priority is given to degree three nodes. if the adjacent nodes are
    // regular, then two swaps are required. In the first step, the degree of the
    // adjacent node is increased by swapping from outer layer, and then swapping
    // is done to make 3 degree node to regular node.

    if (ncount == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);

            if (vertex->getLayerID() == layerid)
            {
                int err = apply_advance_front_triplet_rule(vertex);
                if (!err) ncount++;
            }
        }
    }

    // Final case, this is the most intuitive and general; swapping is done to
    // reduce the vertex degree.
    if (ncount == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);
            int lid = vertex->getLayerID();
            if (lid >= layerid)
            {
                while (1)
                {
                    int err = apply_advance_front_excess_rule(vertex);
                    if (!err)
                        ncount++;
                    else
                        break;
                }
            }
        }
    }

#ifdef DEBUG
    if (ncount)
        cout << "Info: Layer " << layerid << " number of layer edges swapped " << ncount << endl;
#endif

    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    mesh->collect_garbage();

    doublets = search_interior_doublets();
    assert(doublets.empty());

    return ncount;
}
#endif

////////////////////////////////////////////////////////////////////////////////

/*
void
QuadCleanUp::advancing_front_edges_swap()
{
    // ****************************************************************************
    // Implementation ideas:
    // In general, irregular nodes are scattered around in the domain ( even after
    // initial cleanup operations. With the advance front swapping, we try to
    // (1) clean the mesh starting from the boundary (2) cluster irregular nodes
    // deep inside the mesh ( far from the boundary ) with the hope that clustering
    // will provide more opportunities for cleanup operations ( for example, enable
    // diamonds, bridges etc).
    //
    // There are two side-effects of this procedure.
    // 1)  It may create very high valance nodes in the last layers, which we may
    //     not be able to clean. ( this is not a stopper, as we can call vertex
    //     reduction modules, which will again scatter the irregular nodes.
    // 2)  The position of irregular nodes may depends on the density of mesh.
    //
    //
    // ***************************************************************************

    size_t ncount = 0;
    int nfronts;

    nfronts = mesh->setWavefront(0);
    nfronts = mesh->setWavefront(2);

    int istart = 1;
    int iend = nfronts;

    cout << "Temporay Changed " << endl;
    iend = 1;

    for (int ilayer = istart; ilayer < iend; ilayer++)
    {
        // All nodes below a current level never move or modified, therefore
        // make them inactive for smooting operations.
        //
        size_t numnodes = mesh->getSize(0);
        for (size_t j = 0; j < numnodes; j++)
        {
            Vertex *v = mesh->getNodeAt(j);
            if (v->getLayerID() < ilayer)
                v->setConstrainedMark(1);
            else
                v->setConstrainedMark(0);
        }

        int ncount1 = advance_front_edges_swap_once(ilayer);
        //          lapsmooth->execute ();
        ncount += ncount1;
    }

}
*/
////////////////////////////////////////////////////////////////////
/*
vector<Doublet>
QuadCleanUp::search_interior_doublets()
{
    //
    // An interior doublet is a vertex, which is shared by two
    // face neighbours. They are undesirables in the quadmesh.
    // as it would mean the angle is 180 between some adjacent
    // edges...

    // This module finds doublets in the interiour mesh...

    size_t numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    vDoublets.clear();

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (isDoublet(vertex))
        {
            Doublet newdoublet(mesh, vertex);
            if (newdoublet.isSafe())
            {
                newdoublet.makeShield();
                vDoublets.push_back(newdoublet);
            }
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    if (vDoublets.size())
    {
        cout << "Info: Number of interior doublets Detected : "
                << vDoublets.size() << endl;
    }

    return vDoublets;
}
*/

//////////////////////////////////////////////////////////////////////////

/*
vector<Singlet>
QuadCleanUp::search_boundary_singlets()
{
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
    mesh->setFeatureAngles();

    assert(mesh->getAdjTable(0, 2));

    for (size_t i = 0; i < mesh->getSize(2); i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }
   

    size_t numnodes = mesh->getSize(0);
    vSinglets.clear();
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        if (isSinglet(v))
        {
            if (v->getFeatureAngle() >= mesh->getMaxFeatureAngle())
            {
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
*/

///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
int
QuadCleanUp::remove_boundary_singlets_once()
{
    int relexist = mesh->build_relations(0, 2);

    if (vSinglets.empty())
        search_boundary_singlets();

    // Simplest case: just remove it by refinement;
    int ncount = 0;
    for (int i = 0; i < vSinglets.size(); i++)
    {
        int err = vSinglets[i].remove();
        if (!err) ncount++;
    }

    /* Reconsider the following code later.

    // Necessary for type2 changes.
    size_t numfaces = mesh->getSize (2);
    for (size_t i = 0; i < numfaces; i++)
      {
        Face *face = mesh->getFaceAt (i);
        face->setVisitMark (0);
      }

    for (int i = 0; i < vSinglets.size (); i++)
      vSinglets[i].update_type1 ();

    for (int i = 0; i < vSinglets.size (); i++)
      vSinglets[i].update_type2 ();

      for (int i = 0; i < vSinglets.size (); i++)
        vSinglets[i].update_type3 ();

      size_t numnodes = mesh->getSize (0);
      for (size_t i = 0; i < numnodes; i++)
        {
          Vertex *vertex = mesh->getNodeAt (i);
          vertex->setVisitMark (0);
        }

      for (int i = 0; i < vSinglets.size (); i++)
        {
          int err = vSinglets[i].commit ();
          if (!err) ncount++;
        }
     */

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        mesh->collect_garbage();
        mesh->setBoundaryKnown(0);

    }

    vSinglets.clear();

    if (!relexist)
        mesh->clear_relations(0, 2);


    return ncount;
}
#endif
//////////////////////////////////////////////////////////////////////////

/*
int Singlet::remove()
{
    vector<Face*> vfaces = vertex->getRelations2();
    if (vfaces.size() > 1) return 1;

    return mesh->refine_quad15(vfaces[0]);
}
*/
//////////////////////////////////////////////////////////////////////////

/*
int
QuadCleanUp::remove_boundary_singlets()
{
    mesh->setFeatureAngles();

    int ncount = 0;
    while (1)
    {
        size_t nremoved = remove_boundary_singlets_once();
        if (nremoved == 0) break;
        if (!nremoved) return ncount;
        ncount += nremoved;
    }

    if (!mesh->isConsistentlyOriented())
        mesh->makeConsistentlyOriented();

    return 0;
}
*/

////////////////////////////////////////////////////////////////////

/*
int FaceClose::remove()
{
    return 1;
}
*/

////////////////////////////////////////////////////////////////////

/*
bool
FaceClose::isSafe() const
{
    if (face->isRemoved()) return 0;

    if (vertex0->isRemoved()) return 0;
    if (vertex2->isRemoved()) return 0;

    if (vertex0->isBoundary()) return 0;
    if (vertex2->isBoundary()) return 0;

    int vpos = face->queryPosOf(vertex0);
    assert(vpos >= 0);

    if (face->getNodeAt((vpos + 2) % 4) != vertex2)
    {
        cout << "Warning: Face-open requires opposite vertices " << endl;
        cout << "Debug  : Face is : " << face->getNodeAt(0)->getID() << " "
                << face->getNodeAt(1)->getID() << " "
                << face->getNodeAt(2)->getID() << " "
                << face->getNodeAt(3)->getID() << endl;
        cout << "Opposite ends are " << vertex0->getID() << " " << vertex2->getID() << endl;
        return 1;
    }

    assert(face->hasNode(vertex0));
    assert(face->hasNode(vertex2));

    vector<Face*> neighs;
    neighs = vertex0->getRelations2();
    for (size_t j = 0; j < neighs.size(); j++)
    {
        if (neighs[j] != face)
        {
            int val0 = neighs[j]->hasNode(vertex0);
            int val1 = neighs[j]->hasNode(vertex2);
            if (val0 + val1 == 2) return 0;
        }
    }

    neighs = vertex2->getRelations2();
    for (size_t j = 0; j < neighs.size(); j++)
    {
        if (neighs[j] != face)
        {
            int val0 = neighs[j]->hasNode(vertex0);
            int val1 = neighs[j]->hasNode(vertex2);
            if (val0 + val1 == 2) return 0;
        }
    }

    return 1;
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
int
FaceClose::build()
{
    if (!isSafe()) return 1;

    assert(face->getSize(0) == 4);

    Vertex *v0 = vertex0;
    Vertex *v2 = vertex2;

    FaceSequence neighs = face->getRelations202();
    NodeSet nodeset;

    for (size_t i = 0; i < neighs.size(); i++)
    {
        Face *f = neighs[i];
        if (f->isRemoved()) return 1;
        for (int j = 0; j < 4; j++)
        {
            Vertex *v = neighs[i]->getNodeAt(j);
            if (v->isRemoved()) return 1;
            nodeset.insert(v);
        }
    }

    NodeSequence localnodes;
    localnodes.reserve(nodeset.size());

    FaceSet faces_to_check;
    map<Vertex*, Point3D> backupCoords;
    set<Vertex*> ::const_iterator it;
    for (it = nodeset.begin(); it != nodeset.end(); ++it)
    {
        Vertex *v = *it;
        backupCoords[v] = v->getXYZCoords();
        localnodes.push_back(v);
        FaceSequence vneighs = v->getRelations2();
        for (int i = 0; i < vneighs.size(); i++)
        {
            Face *f = vneighs[i];
            if (f->isRemoved()) return 1;
            faces_to_check.insert(f);
        }
    }
    faces_to_check.erase(face);

    // Add a new vertex in the mesh...
    Point3D p3d = Vertex::mid_point(v0, v2);

    // Both nodes will merge at the same point...
    vertex0->setXYZCoords(p3d);
    vertex2->setXYZCoords(p3d);

    // Don't move vertex0 and vertex2, because they will be merged.
    int fixnodes[] = {0, 0};
    if (!vertex0->isConstrained())
    {
        vertex0->setConstrainedMark(1);
        fixnodes[0] = 1;
    }

    if (!vertex2->isConstrained())
    {
        vertex2->setConstrainedMark(1);
        fixnodes[1] = 1;
    }

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(localnodes);

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

    // Roll back the positions
    if (!pass)
    {
        if (fixnodes[0]) vertex0->setConstrainedMark(0);
        if (fixnodes[1]) vertex2->setConstrainedMark(0);
        map<Vertex*, Point3D>::const_iterator miter;
        for (miter = backupCoords.begin(); miter != backupCoords.end(); ++miter)
        {
            Vertex *v = miter->first;
            v->setXYZCoords(miter->second);
        }
        return 1;
    }

    replacedNode = Vertex::newObject();
    replacedNode->setXYZCoords(p3d);

    return 0;
}
*/

////////////////////////////////////////////////////////////////////////////////

int QuadCleanUp :: remove_tunnels()
{
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

/*
Vertex*
QuadCleanUp::insert_doublet(Face *face, Vertex *v0, Vertex *v2)
{

    //  Public Function ...
    Point3D p3d = Vertex::mid_point(v0, v2);

    Vertex *doublet = Vertex::newObject();
    doublet->setXYZCoords(p3d);

    Vertex *o1 = NULL, *o2 = NULL;

    vector<Vertex*> connect = face->getNodes();
    for (size_t i = 0; i < connect.size(); i++)
    {
        if (connect[i] == v0)
        {
            o1 = connect[(i + 1) % 4];
            o2 = connect[(i + 3) % 4];
            break;
        }
    }

    mesh->addNode(doublet);

    connect[0] = doublet;
    connect[1] = v0;
    connect[2] = o1;
    connect[3] = v2;

    Face *newquad1 = Face::newObject();
    newquad1->setNodes(connect);
    mesh->addFace(newquad1);

    face->setRemoveMark(1);

    connect[0] = doublet;
    connect[1] = v2;
    connect[2] = o2;
    connect[3] = v0;

    Face *newquad2 = Face::newObject();
    newquad2->setNodes(connect);
    mesh->addFace(newquad2);

    if (mesh->getAdjTable(0, 2))
    {
        v0->removeRelation2(face);
        o1->removeRelation2(face);
        v2->removeRelation2(face);
        o2->removeRelation2(face);

        v0->addRelation2(newquad1);
        v0->addRelation2(newquad2);

        v2->addRelation2(newquad1);
        v2->addRelation2(newquad2);

        doublet->addRelation2(newquad1);
        doublet->addRelation2(newquad2);

        o1->addRelation2(newquad1);
        o2->addRelation2(newquad2);
    }

    return doublet;
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face)
{

    //  Public Function ...
    if (face->getSize(0) != 4)
        return NULL;
    Vertex *v0 = face->getNodeAt(0);
    Vertex *v2 = face->getNodeAt(2);
    return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_boundary_doublet(Face *face)
{

    //  Public Function ...
    if (!face->isBoundary())
        return NULL;

    if (face->getSize(0) != 4)
        return NULL;

    int ncount = 0;
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = face->getNodeAt(i);
        ncount += v->isBoundary();
    }

    if (ncount != 3)
        return NULL;

    Vertex *v0 = NULL, *v2 = NULL;
    for (int i = 0; i < 4; i++)
    {
        Vertex *v = face->getNodeAt(i);
        if (!v->isBoundary())
        {
            v0 = face->getNodeAt((i + 0) % 4);
            v2 = face->getNodeAt((i + 2) % 4);
        }
    }

    return insert_doublet(face, v0, v2);
}

////////////////////////////////////////////////////////////////////////////////
*/

