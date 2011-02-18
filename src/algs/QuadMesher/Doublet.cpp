#include <meshkit/QuadCleanUp.hpp>

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_doublet_tag(Mesh *mesh)
{
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (QuadCleanUp::isDoublet(vertex))
        {
            vertex->setTag(0);
        }
        else
            vertex->setTag(1);
    }

    if (!relexist2)
        mesh->clear_relations(0, 2);
}

//////////////////////////////////////////////////////////////////////

bool
Doublet::isSafe() const
{
    assert(vertex);

    if (vertex->isRemoved()) return 0;

    FaceSequence apexfaces = vertex->getRelations2();
    for (size_t i = 0; i < apexfaces.size(); i++)
    {
        Face *face = apexfaces[i];
        if (face->isRemoved()) return 0;
        if (face->isVisited()) return 0;
    }
    return 1;
}

////////////////////////////////////////////////////////////////////

void
Doublet::makeShield()
{
    FaceSequence apexfaces = vertex->getRelations2();
    for (size_t i = 0; i < apexfaces.size(); i++)
    {
        Face *face = apexfaces[i];
        face->setVisitMark(1);
    }
}

///////////////////////////////////////////////////////////////////////////////

vector<Doublet>
QuadCleanUp::search_interior_doublets()
{
    //
    ///////////////////////////////////////////////////////////////////////////
    // An interior doublet is a vertex, which is shared by two face neighbours. 
    // They are undesirables in the quadmesh as it would mean the angle is 180 
    // between some adjacent edges...
    // 
    ///////////////////////////////////////////////////////////////////////////

    int relexist2 = mesh->build_relations(0, 2);

    size_t numnodes = mesh->getSize(0);

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

    if (!relexist2)
        mesh->clear_relations(0, 2);

    if (vDoublets.size())
    {
        cout << "Info: Number of interior doublets Detected : "
                << vDoublets.size() << endl;
    }

    return vDoublets;
}

//////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face, Vertex *v0, Vertex *v2)
{
    //Create new vertex at the center of (v0,v2)
    Point3D p3d = Vertex::mid_point(v0, v2);

    Vertex *doublet = Vertex::newObject();
    doublet->setXYZCoords(p3d);

    int pos = face->getPosOf( v0 );
    if( pos < 0) return NULL;

    assert( v2 == face->getNodeAt( (pos+2)%4 ) );
 
    Vertex *o1 = face->getNodeAt( (pos+1)%4 );
    Vertex *o2 = face->getNodeAt( (pos+3)%4 );
   
    //  Creating a doublet in the mesh changes:
    //  (1)  insert new node
    //  (2)  one old face is removed
    //  (3)  two new faces inserted.
    //

    NodeSequence connect(4);

    mesh->addNode(doublet);
    connect[0] = doublet;
    connect[1] = v0;
    connect[2] = o1;
    connect[3] = v2;

    Face *newquad1 = Face::newObject();
    newquad1->setNodes(connect);
    mesh->addFace(newquad1);

    connect[0] = doublet;
    connect[1] = v2;
    connect[2] = o2;
    connect[3] = v0;

    Face *newquad2 = Face::newObject();
    newquad2->setNodes(connect);
    mesh->addFace(newquad2);

    mesh->remove( face );

    return doublet;
}

////////////////////////////////////////////////////////////////////////////////

Vertex*
QuadCleanUp::insert_doublet(Face *face)
{
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
    ////////////////////////////////////////////////////////////////////////////
    //
    //                                  X  ( Internal Vertex)
    //                               .      .
    //	                           .           .				
    //                           .               .
    //                         .                   .
    //               ********X...........X...........X************************
    //        Boundary                 Singlet                Boundary 
    //
    //   There is one quad on the boundary with three nodes on the boundary and
    //   one internal nodes. In order to remove the singlet node, we artifically
    //   create one doublet between the internal node and the singlet node.
    ////////////////////////////////////////////////////////////////////////////
  
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

int
Doublet::remove()
{
    if (vertex->isBoundary())
        return 1;

    FaceSequence neighs = vertex->getRelations2();
    assert(neighs.size() > 0);

    if (neighs.size() != 2)
        return 1;

    assert(neighs[0] != neighs[1]);

    Vertex *d1 = NULL, *d2 = NULL, *o1 = NULL, *o2 = NULL;

    NodeSequence connect = neighs[0]->getNodes();
    if (connect.size() != 4)
        return 1;

    if (neighs[0]->isRemoved())
        return 2;
    if (neighs[1]->isRemoved())
        return 2;

    for (size_t i = 0; i < connect.size(); i++)
    {
        if (connect[i] == vertex)
        {
            d1 = connect[(i + 1) % 4];
            o1 = connect[(i + 2) % 4];
            d2 = connect[(i + 3) % 4];
            break;
        }
    }

    connect = neighs[1]->getNodes();
    if (connect.size() != 4)
        return 1;

    for (size_t i = 0; i < connect.size(); i++)
    {
        if (connect[i] == vertex)
        {
            o2 = connect[(i + 2) % 4];
            break;
        }
    }

    assert(d1);
    assert(d2);
    assert(o1);
    assert(o2);
    assert(o1 != o2);

    d1->removeRelation2(neighs[0]);
    d1->removeRelation2(neighs[1]);

    d2->removeRelation2(neighs[0]);
    d2->removeRelation2(neighs[1]);

    o1->removeRelation2(neighs[0]);
    o2->removeRelation2(neighs[1]);

    connect[0] = d1;
    connect[1] = o1;
    connect[2] = d2;
    connect[3] = o2;

    Face *newquad = Face::newObject();
    newquad->setNodes(connect);
    mesh->addFace(newquad);

    d1->addRelation2(newquad);
    d2->addRelation2(newquad);

    o1->addRelation2(newquad);
    o2->addRelation2(newquad);

    neighs[0]->setStatus( MeshEntity::REMOVE);
    neighs[1]->setStatus( MeshEntity::REMOVE);
    vertex->setStatus( MeshEntity::REMOVE);

    return 0;
}

///////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_doublets_once()
{
    int relexist = mesh->build_relations(0, 2);

    search_interior_doublets();

    int ncount = 0;
    for (size_t i = 0; i < vDoublets.size(); i++)
    {
        int err = vDoublets[i].remove();
        if (!err)
            ncount++;
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: Number of doublets committed : " << ncount << endl;
        mesh->collect_garbage();
    }

    mesh->setBoundaryKnown(0);

    if (!relexist)
        mesh->clear_relations(0, 2);

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_interior_doublets()
{
    if (!mesh->isPruned()) mesh->prune();

    mesh->search_boundary();

    // It is possible that removal of doublet may create singlets on the boundary
    // so, it is better to call doublets first and then call to singlet removal next.

    while (1)
    {
        int ncount = remove_doublets_once();
        if (ncount == 0) break;
    }

    if (vDoublets.empty()) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

