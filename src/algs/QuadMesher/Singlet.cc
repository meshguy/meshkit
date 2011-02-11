#include "meshkit/QuadCleanUp.h"
#include <sstream>
#include <assert.h>

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_singlet_tag(Mesh *mesh)
{
    size_t numnodes = mesh->getSize(0);

    int relexist = mesh->build_relations(0, 2);

    if (!mesh->isBoundaryKnown());
    mesh->search_boundary();

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (QuadCleanUp::isSinglet(vertex))
        {
            vertex->setTag(0);
        }
        else
            vertex->setTag(1);
    }

    if (!relexist)
        mesh->clear_relations(0, 2);
}

////////////////////////////////////////////////////////////////////

vector<Singlet>
QuadCleanUp::search_boundary_singlets()
{
    vSinglets.clear();
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

    assert(mesh->getAdjTable(0, 2));

    for (size_t i = 0; i < mesh->getSize(2); i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitMark(0);
    }
   
    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        if (isSinglet(v))
        {
           Singlet newsinglet(mesh, v);
           vSinglets.push_back(newsinglet);
        }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    if (vSinglets.size())
        cout << "Info: Number of Singlets detected " << vSinglets.size() << endl;

    return vSinglets;
}

///////////////////////////////////////////////////////////////////////////////

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

/*

    // Necessary for type2 changes.
    size_t numfaces = mesh->getSize (2);
    for (size_t i = 0; i < numfaces; i++)
      {
        Face *face = mesh->getFaceAt (i);
        face->setVisitMark (0);
      }


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
        mesh->setBoundaryKnown(0);
    }

    vSinglets.clear();

    if (!relexist)
        mesh->clear_relations(0, 2);

    return ncount;
}
//////////////////////////////////////////////////////////////////////////
int Singlet::remove()
{
  int err = update_type1();
  if( err == 0) return 0;

  return remove_by_refinement();
  return 1;
}

//////////////////////////////////////////////////////////////////////////

int Singlet::remove_by_refinement()
{
    //
    // It must be the last resort, if everything fails. This method is in
    // general bad, because it create 4 new elements and can create many
    // irrgular nodes ( maximum 8 nodes ) 
    //
    FaceSequence vfaces = vertex->getRelations2();
    if (vfaces.size() > 1) return 1;

    return mesh->refine_quad15(vfaces[0]);
}
//////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_boundary_singlets()
{
    int ncount = 0;
    while (1)
    {
        size_t nremoved = remove_boundary_singlets_once();
        if (nremoved == 0) break;
        ncount += nremoved;
    }

    if (!mesh->is_consistently_oriented())
        mesh->make_consistently_oriented();

    mesh->collect_garbage();

    if( mesh->count_concave_faces() == 0) 
        mopt.shape_optimize(mesh);

    return 0;
}

////////////////////////////////////////////////////////////////////

int
Singlet::update_type1()
{
    // If any neighbour to the singlet has excess degree then perhaps
    // by Swapping an edge, singlet can be removed.
    // 
    //                             X
    //                           *   *
    //  Can this be swapped ?  *       *  Can this be swapped ?
    //                       *           *
    //                     X ***** X******X  ( Ihave extra degree ) 
    //
    /////////////////////////////////////////////////////////////////////

    if (!active) return 1;

    // Find out the face containing the singlet first.
    FaceSequence vfaces = vertex->getRelations2();

    if (vfaces.size() != 1) return 1;

    Face *f0 = vfaces[0];
    int pos = f0->getPosOf(vertex);

    // Try the first edge
    Vertex *nextvertex0, *nextvertex1;

    nextvertex0 = f0->getNodeAt((pos + 1) % 4);
    nextvertex1 = f0->getNodeAt((pos + 2) % 4);
    assert(nextvertex0->isBoundary());
    vfaces = nextvertex0->getRelations2();
    if (vfaces.size() > 2)
    {
        type = 1;
        SwapQuadEdge edge(mesh, nextvertex0, nextvertex1, f0);
        int err = edge.apply_singlet_rule(vertex);
        if (!err)
        {
            active = 0;
            return 0;
        }
    }
    
    // If failed, try the other edge ...
    nextvertex0 = f0->getNodeAt((pos + 3) % 4);
    nextvertex1 = f0->getNodeAt((pos + 2) % 4);
    assert(nextvertex0->isBoundary());
    vfaces = nextvertex0->getRelations2();

    if (vfaces.size() > 2)
    {
        type = 1;
        SwapQuadEdge edge(mesh, nextvertex0, nextvertex1, f0);
        int err = edge.apply_singlet_rule(vertex);
        if (!err)
        {
            active = 0;
            return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
Singlet::update_type2()
{
    Break();
/*
    if (!active) return 1;

    if (!vertex->isBoundary()) return 1;

    vector<Face*> vfaces = vertex->getRelations2();

    if (vfaces.size() != 1) return 1;

    Face *f0 = vfaces[0];
    Face *f1 = NULL;

    Vertex * v[6];
    v[0] = vertex;
    v[1] = NULL;
    v[2] = NULL;
    v[3] = NULL;
    v[4] = NULL;
    v[5] = NULL;

    int pos = f0->getPosOf(vertex);
    Vertex *v1 = f0->getNodeAt((pos + 1) % 4);
    Vertex *v2 = f0->getNodeAt((pos + 2) % 4);
    Vertex *v3 = f0->getNodeAt((pos + 3) % 4);

    vector<Face*> neighs;

    bool found = 0;

    if (!found)
    {
        vfaces = v1->getRelations2();
        if (vfaces.size() == 2)
        {
            neighs = Mesh::getRelations112(v1, v2);
            assert(neighs.size() > 0);
            if (neighs.size() == 2)
            {
                if (neighs[0] == f0) f1 = neighs[1];
                if (neighs[1] == f0) f1 = neighs[0];
                if (!QuadCleanUp::hasSinglet(f1))
                {
                    type = 2;
                    pos = f1->getPosOf(v1);

                    v[1] = v1;
                    v[2] = v2;
                    v[3] = v3;
                    v[5] = f1->getNodeAt((pos + 2) % 4);

                    if (f1->getNodeAt((pos + 1) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 3) % 4);
                        found = 1;
                    }

                    if (f1->getNodeAt((pos + 3) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 1) % 4);
                        found = 1;
                    }
                }
            } // End checking first vertex
        }
    }

    if (!found)
    {
        vfaces = v3->getRelations2();
        if (vfaces.size() == 2)
        {
            neighs = Mesh::getRelations112(v2, v3);
            assert(neighs.size() > 0);
            if (neighs.size() == 2)
            {
                if (neighs[0] == f0) f1 = neighs[1];
                if (neighs[1] == f0) f1 = neighs[0];
                if (!QuadCleanUp::hasSinglet(f1))
                {
                    type = 2;
                    int pos = f1->getPosOf(v3);
                    v[1] = v3;
                    v[2] = v2;
                    v[3] = v1;
                    v[5] = f1->getNodeAt((pos + 2) % 4);

                    if (f1->getNodeAt((pos + 1) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 3) % 4);
                        found = 1;
                    }

                    if (f1->getNodeAt((pos + 3) % 4) == v2)
                    {
                        v[4] = f1->getNodeAt((pos + 1) % 4);
                        found = 1;
                    }
                }
            } // End checking second vertex
        }
    }

    if (!found) return 1;

    if (f0->isVisited() || f1->isVisited()) return 2;

    assert(v[0]->isBoundary());
    assert(v[1]->isBoundary());
    assert(v[3]->isBoundary());

    Point3D xyz;

    // Two new nodes will be created ...
    Vertex *doublet1 = Vertex::newObject();
    xyz = Vertex::mid_point(v[0], v[2]);
    doublet1->setXYZCoords(xyz);
    mesh->addNode(doublet1);

    Vertex *doublet2 = Vertex::newObject();
    xyz = Vertex::mid_point(v[1], v[5]);
    doublet2->setXYZCoords(xyz);
    mesh->addNode(doublet2);

    // Two old faces will be updated and two new faces will be created ...

    Face *newface;
    vector<Vertex*> connect(4);

    connect[0] = v[3];
    connect[1] = v[0];
    connect[2] = doublet1;
    connect[3] = v[2];
    f0->setNodes(connect);
    f0->setVisitMark(1);

    connect[0] = v[1];
    connect[1] = v[4];
    connect[2] = v[5];
    connect[3] = doublet2;
    f1->setNodes(connect);
    f1->setVisitMark(1);

    connect[0] = v[0];
    connect[1] = v[1];
    connect[2] = doublet2;
    connect[3] = doublet1;

    newface = Face::newObject();
    newface->setNodes(connect);
    mesh->addFace(newface);

    connect[0] = doublet1;
    connect[1] = doublet2;
    connect[2] = v[5];
    connect[3] = v[2];

    newface = Face::newObject();
    newface->setNodes(connect);
    mesh->addFace(newface);
*/

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
Singlet::update_type3()
{
    Break();

    if (!active) return 1;

    if (!vertex->isBoundary()) return 1;

    FaceSequence vfaces = vertex->getRelations2();

    if (vfaces.size() != 1) return 1;

    Face *f0 = vfaces[0];

    int pos = f0->getPosOf(vertex);

    Vertex *v0 = f0->getNodeAt((pos + 0) % 4);
    Vertex *v2 = f0->getNodeAt((pos + 2) % 4);

    if (v2->isBoundary()) return 1;

    Point3D p0 = v0->getXYZCoords();
    Point3D p2 = v2->getXYZCoords();

    vfaces = v2->getRelations2();
    double area0 = 0.0;
    for (size_t i = 0; i < vfaces.size(); i++)
        area0 += vfaces[i]->getArea();

    v2->setXYZCoords(p0);

    bool pass = 1;
    double area1 = 0.0;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        if (vfaces[i] != f0)
        {
            if (!vfaces[i]->isConvex())
            {
                pass = 0;
                break;
            }
            area1 += vfaces[i]->getArea();
        }
    }

    if (pass && fabs(area1 - area0) > 1.0E-06) pass = 0;

    if (!pass)
    {
        v2->setXYZCoords(p2);
        return 2;
    }

    Vertex *vnew = Vertex::newObject();
    vnew->setXYZCoords(p0);
    vnew->setBoundaryMark(vertex->getBoundaryMark());
    newNodes.resize(1);
    newNodes[0] = vnew;

#ifdef SEQUENCE_IS_VECTOR
    newFaces.reserve(vfaces.size() - 1);
#endif

    NodeSequence connect;
    for (size_t i = 0; i < vfaces.size(); i++)
    {
        if (vfaces[i] != f0)
        {
            connect = vfaces[i]->getNodes();
            for (int j = 0; j < 4; j++)
            {
                if (connect[j] == v2)
                    connect[j] = vnew;
            }
            Face *newface = Face::newObject();
            newface->setNodes(connect);
            newFaces.push_back(newface);
        }
    }
    oldFaces = vfaces;

    oldNodes.resize(1);
    oldNodes[0] = vertex;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
Singlet::clear()
{

    for (size_t j = 0; j < newNodes.size(); j++)
        delete newNodes[j];
    newNodes.clear();
    for (size_t j = 0; j < newFaces.size(); j++)
        delete newFaces[j];
    newFaces.clear();

}

///////////////////////////////////////////////////////////////////////////////

int
Singlet::commit()
{
    for (size_t i = 0; i < oldFaces.size(); i++)
    {
        if (oldFaces[i]->isRemoved())
        {
            clear();
            return 1;
        }
    }

    if (!active) return 2;

    for (size_t i = 0; i < oldNodes.size(); i++)
        oldNodes[i]->setStatus( MeshEntity::REMOVE);
    oldNodes.clear();

    for (size_t i = 0; i < oldFaces.size(); i++)
        oldFaces[i]->setStatus(MeshEntity::REMOVE);
    oldFaces.clear();

    for (size_t i = 0; i < newNodes.size(); i++)
        mesh->addNode(newNodes[i]);
    newNodes.clear();

    for (size_t i = 0; i < newFaces.size(); i++)
        mesh->addFace(newFaces[i]);
    newFaces.clear();

    return 0;

}

///////////////////////////////////////////////////////////////////////////////

