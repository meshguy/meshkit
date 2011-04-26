#include "QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_diamond_tag(Mesh *mesh)
{
    QuadCleanUp qClean(mesh);
    vector<Diamond> diamonds = qClean.search_diamonds();

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setTag(1);
    }

    set<Face*>::const_iterator it;
    for (size_t i = 0; i < diamonds.size(); i++)
    {
        Face *face = diamonds[i].face;
        face->setTag(0);
    }
}

////////////////////////////////////////////////////////////////////
bool
FaceClose::isSafe() const
{
    if (face->isRemoved()) return 0;

    if (vertex0->isRemoved()) return 0;
    if (vertex2->isRemoved()) return 0;

    if (vertex0->isBoundary()) return 0;
    if (vertex2->isBoundary()) return 0;

    int vpos = face->getPosOf(vertex0);
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

    FaceSequence neighs;
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

///////////////////////////////////////////////////////////////////////////////

int
FaceClose::build()
{
    if (!isSafe()) return 1;

    assert(face->getSize(0) == 4);

    Vertex *v0 = vertex0;
    Vertex *v2 = vertex2;

    FaceSequence neighs = face->getRelations202();
    NodeSet nodeset;

    size_t nSize = neighs.size();
    for (size_t i = 0; i < nSize; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Vertex *v = neighs[i]->getNodeAt(j);
            nodeset.insert(v);
        }
    }

    NodeSequence localnodes;
#ifdef SEQUENCE_IS_VECTOR
    localnodes.reserve(nodeset.size());
#endif

    FaceSet faces_to_check;
    map<Vertex*, Point3D> backupCoords;
    set<Vertex*> ::const_iterator it;
    for (it = nodeset.begin(); it != nodeset.end(); ++it)
    {
        Vertex *v = *it;
        backupCoords[v] = v->getXYZCoords();
        localnodes.push_back(v);
        FaceSequence vneighs = v->getRelations2();
        nSize = vneighs.size();
        for (size_t i = 0; i < nSize; i++)
            faces_to_check.insert( vneighs[i] );
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

////////////////////////////////////////////////////////////////////////////////

int
FaceClose::commit()
{
    if (replacedNode == NULL) return 1;

    const FaceSequence &vr0 = vertex0->getRelations2();
    const FaceSequence &vr2 = vertex2->getRelations2();

    mesh->addNode(replacedNode);
    
    NodeSequence fnodes;

    size_t nSize = vr0.size();
    for (size_t i = 0; i < nSize; i++)
    {
        if (vr0[i] != face)
        {
            fnodes = vr0[i]->getNodes();
            int nv = fnodes.size();
            for( int j = 0; j < nv; j++) 
               if( fnodes[j] == vertex0) fnodes[j] = replacedNode;
            mesh->deactivate(vr0[i]);  // remove all existing relationships.
            vr0[i]->setNodes( fnodes );
            mesh->reactivate(vr0[i]);  // rebuilds new relationships.
        }
    }

    nSize = vr2.size();
    for (size_t i = 0; i < vr2.size(); i++)
    {
        if (vr2[i] != face)
        {
            fnodes = vr2[i]->getNodes();
            int nv = fnodes.size();
            for( int j = 0; j < nv; j++) 
               if( fnodes[j] == vertex2) fnodes[j] = replacedNode;
            mesh->deactivate(vr2[i]); // remove all existing relatioships
            vr2[i]->setNodes( fnodes );
            mesh->reactivate(vr2[i]); // rebuilds new relationships..
        }
    }

    // Two nodes and face go away from the mesh..
    mesh->remove(face);
    mesh->remove(vertex0);
    mesh->remove(vertex2);

    replacedNode = NULL; // So that destructor can delete if this is not used.

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool
QuadCleanUp::isDiamond(Face *face, int &pos, int type)
{
    pos = -1;
    const Vertex *v0 = face->getNodeAt(0);
    const Vertex *v1 = face->getNodeAt(1);
    const Vertex *v2 = face->getNodeAt(2);
    const Vertex *v3 = face->getNodeAt(3);

    FaceSequence neighs;

    neighs = v0->getRelations2();
    int d0 = neighs.size();

    neighs = v1->getRelations2();
    int d1 = neighs.size();

    neighs = v2->getRelations2();
    int d2 = neighs.size();

    neighs = v3->getRelations2();
    int d3 = neighs.size();

    // Boundary Cases ...
    if (v0->isBoundary() || v2->isBoundary())
    {
        if( d0 <= v0->get_ideal_vertex_degree( Face::QUADRILATERAL ) ) return 0;
        if( d2 <= v2->get_ideal_vertex_degree( Face::QUADRILATERAL ) ) return 0;
        if (!v1->isBoundary() && !v3->isBoundary())
        {
            if (d1 == 3 && d3 == 3)
            {
                pos = 1;
                return 1;
            }
        }
    }

    if (v1->isBoundary() || v3->isBoundary())
    {
        if( d1 <= v1->get_ideal_vertex_degree( Face::QUADRILATERAL ) ) return 0;
        if( d3 <= v3->get_ideal_vertex_degree( Face::QUADRILATERAL ) ) return 0;

        if (!v0->isBoundary() && !v2->isBoundary())
        {
            if ((d0 == 3 && d2 == 3))
            {
                pos = 0;
                return 1;
            }
        }
    }

    if (v0->isBoundary()) return 0;
    if (v1->isBoundary()) return 0;
    if (v2->isBoundary()) return 0;
    if (v3->isBoundary()) return 0;

    if ((d0 == 3 && d2 == 3))
    {
        if( d1 < 5 && d3 < 5)  return 0;
        pos = 0;
        return 1;
    }

    if (d1 == 3 && d3 == 3)
    {
        if( d0 < 5 && d2 < 5)  return 0;
        pos = 1;
        return 1;
    }

/*
    if (d1 == 5 && d3 == 5)
    {
        if ((d0 == 3 && d2 == 4) || (d0 == 4 && d2 == 3))
        {
            pos = 0;
            return 1;
        }
    }

    if (d0 == 5 && d2 == 5)
    {
        if ((d1 == 3 && d3 == 4) || (d1 == 4 && d3 == 3))
        {
            pos = 1;
            return 1;
        }
    }
*/

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

/*
int Diamond::remove()
{
    assert(face);
    assert(mesh);
    assert(vertex0);
    assert(vertex2);

    if (face->isRemoved()) return 1;
    if (vertex0->isRemoved()) return 1;
    if (vertex2->isRemoved()) return 1;

    vector<Face*> neighs;
    neighs = vertex0->getRelations2();
    if (neighs.size() != 3) return 2;

    neighs = vertex2->getRelations2();
    if (neighs.size() != 3) return 2;

    FaceClose op(mesh, face, vertex0, vertex2);
    int err = op.remove();

    return err;
}
*/



////////////////////////////////////////////////////////////////////

vector<Diamond>
QuadCleanUp::search_diamonds(int type)
{
    if (!mesh->isPruned()) mesh->prune();

    size_t numfaces = mesh->getSize(2);

    mesh->search_boundary();

    int relexist = mesh->build_relations(0, 2);

    int pos;

    assert(mesh->isBoundaryKnown());

    vDiamonds.clear();
    for (size_t iface = 0; iface < numfaces; iface++)
    {
        Face *face = mesh->getFaceAt(iface);
        if (isDiamond(face, pos))
        {
            Diamond diamond(mesh, face, pos);
            vDiamonds.push_back(diamond);
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    if (vDiamonds.size())
        cout << "Info: number of diamonds detected " << vDiamonds.size() << endl;

    return vDiamonds;
}

////////////////////////////////////////////////////////////////////

/*
int
Diamond::isSafe()
{
    vector<Face*> vr0 = vertex0->getRelations2();
    for (size_t i = 0; i < vr0.size(); i++)
        if (vr0[i]->isVisited()) return 0;

    vector<Face*> vr2 = vertex2->getRelations2();
    for (size_t i = 0; i < vr2.size(); i++)
        if (vr2[i]->isVisited()) return 0;

    return 1;

}

/////////////////////////////////////////////////////////////////////////

int
Diamond::makeShield()
{
    vector<Face*> vr0 = vertex0->getRelations2();
    for (size_t i = 0; i < vr0.size(); i++)
        vr0[i]->setVisitMark(1);

    vector<Face*> vr2 = vertex2->getRelations2();
    for (size_t i = 0; i < vr2.size(); i++)
        vr2[i]->setVisitMark(1);
    return 0;
}
*/

////////////////////////////////////////////////////////////////////////////////

int
Diamond::build()
{
    faceclose = new FaceClose(mesh, face, vertex0, vertex2);
    return faceclose->build();
}

///////////////////////////////////////////////////////////////////////////////

int
Diamond::commit()
{
    if (faceclose)
    {
        int err = faceclose->commit();
        delete faceclose;
        faceclose = NULL;
        return err;
    }

    return 1;
}
///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds_once()
{
    int rel0exist = mesh->build_relations(0, 0); // Need for laplace smoothing
    int rel2exist = mesh->build_relations(0, 2);

    // Essential step because only the faces are updated only.
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    int err, pos, ncount = 0;
    size_t faceid = 0;
    while(1) {
        if( faceid >= mesh->getSize(2) ) break;
        Face *face = mesh->getFaceAt(faceid++);
        if (isDiamond(face, pos))
        {
            Diamond diamond(mesh, face, pos);
            err = diamond.build();
            if (!err)
            {
                err = diamond.commit();
                if (!err) ncount++;
            }
        }
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: number of diamonds removed : " << ncount << endl;
        mesh->collect_garbage();
    }

    if (!rel0exist) mesh->clear_relations(0, 0);
    if (!rel2exist) mesh->clear_relations(0, 2);

    vDiamonds.clear();

    return ncount;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_diamonds_in_layer(int layerid)
{
    int rel0exist = mesh->build_relations(0, 0); // Need for laplace smoothing
    int rel2exist = mesh->build_relations(0, 2);

    // Essential step because only the faces are updated only.
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }

    int err, pos, ncount = 0;

    size_t faceid = 0;
    while (1)
    {
        if (faceid >= mesh->getSize(2)) break;
        Face *face = mesh->getFaceAt(faceid);
        if (face->getLayerID() == layerid)
        {
            if (isDiamond(face, pos))
            {
                Diamond diamond(mesh, face, pos);
                err = diamond.build();
                if (!err)
                {
                    err = diamond.commit();
                    if (!err) ncount++;
                }
            }
        }
        faceid++;
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: number of diamonds removed : " << ncount << endl;
        mesh->collect_garbage();
    }

    if (!rel0exist)
        mesh->clear_relations(0, 0);

    if (!rel2exist)
        mesh->clear_relations(0, 2);

    vDiamonds.clear();

    return ncount;
}

////////////////////////////////////////////////////////////////////

int 
QuadCleanUp::remove_diamonds()
{
    if (!mesh->isPruned()) mesh->prune();

    mesh->search_boundary();

    vector<Doublet> doublets;
    doublets = search_interior_doublets();
    size_t num_doublets_start = doublets.size();
    if (num_doublets_start)
    {
        cout << "Info: Mesh has interior doublets. Removing them " << endl;
        remove_interior_doublets();
        doublets = search_interior_doublets();
        num_doublets_start = doublets.size();
    }

    while (1)
    {
        int ncount = remove_diamonds_once();
        if (ncount == 0) break;
    }

    doublets = search_interior_doublets();
    size_t num_doublets_end = doublets.size();
    if (num_doublets_end > num_doublets_start)
    {
        cout << "Warning: diamonds removal created doublet in the mesh: Now removing them " << endl;
        remove_interior_doublets();
    }
  
    if( mesh->count_concave_faces() == 0) 
        mopt.shape_optimize(mesh);

    return 0;
}

////////////////////////////////////////////////////////////////////

