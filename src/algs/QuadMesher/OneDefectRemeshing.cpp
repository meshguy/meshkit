#include "QuadCleanUp.hpp"

using namespace Jaal;

ObjectPool<Vertex> OneDefectPatch :: nodePool;
ObjectPool<Face>   OneDefectPatch :: facePool;

///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::setTags()
{
    size_t nSize;

    nSize = mesh->getSize(0);
    for (size_t i = 0; i < nSize; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        v->setTag(0);
    }

    nSize = mesh->getSize(2);
    for (size_t i = 0; i < nSize; i++)
    {
        Face *f = mesh->getFaceAt(i);
        f->setTag(0);
    }

    nSize = nodepath.size();
    for (size_t i = 0; i < nSize; i++)
        nodepath[i]->setTag(1);

    nSize = bound_nodes.size();
    for (size_t i = 0; i <  nSize; i++)
        bound_nodes[i]->setTag(3);

    set<Vertex*>::const_iterator vit;
    for (vit = corners.begin(); vit != corners.end(); ++vit)
        (*vit)->setTag(2);


    set<Face*>::const_iterator it;
    for (it = faces.begin(); it != faces.end(); ++it)
        (*it)->setTag(1);
}

///////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::isSafe()
{
    FaceSet::const_iterator it;
    for (it = faces.begin(); it != faces.end(); ++it)
    {
        if ((*it)->isRemoved()) return 0;
    }

    // Check if the previous updates modified the splitting node degree.
    if (quad_splitting_node)
    {
        const FaceSequence &vfaces = quad_splitting_node->getRelations2();
        int nSize = vfaces.size();
        if ( nSize != quad_splitting_node_degree) return 0;
    }

    return 1;
}
///////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::getPosOf(const Vertex *v)
{
    size_t nSize = bound_nodes.size();
    for (size_t i = 0; i <  nSize; i++)
        if (bound_nodes[i] == v) return i;

//  setTags();
    mesh->saveAs("dbg.dat");
    cout << "Error: Vertex not found on the boundary " << endl;
    exit(0);

    return -1;
}

///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::start_boundary_loop_from(Vertex *vmid)
{
    assert(corners.find(vmid) != corners.end());

    NodeSequence::iterator middle;
    middle = find(bound_nodes.begin(), bound_nodes.end(), vmid);
    assert(middle != bound_nodes.end());

    std::rotate(bound_nodes.begin(), middle, bound_nodes.end());
    assert(bound_nodes[0] == vmid);

    set_boundary_segments();
}
///////////////////////////////////////////////////////////////////////////////

NodeSequence OneDefectPatch::get_bound_nodes(const Vertex *src, const Vertex *dst)
{
    int start_pos = getPosOf(src);
    int end_pos = getPosOf(dst);
    int nsize = bound_nodes.size();

    assert(nsize > 1);

    if (end_pos == 0) end_pos = nsize;
    assert(end_pos > start_pos);

    NodeSequence seq(end_pos - start_pos + 1);
    int index = 0;
    for (int i = start_pos; i <= end_pos; i++)
        seq[index++] = bound_nodes[i % nsize];

    return seq;
}

///////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::has_irregular_node_on_first_segment() const
{
    int start_pos = cornerPos[0];
    int end_pos = start_pos + segSize[0] - 1;

    FaceSequence vfaces;
    for (int i = 1; i < end_pos - 1; i++)
    {
        vfaces = bound_nodes[i]->getRelations2();
        if (vfaces.size() < 4) return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::is_simply_connected()
{
    vector<Edge> edges;
    edges.reserve(4 * faces.size());

    FaceSet::const_iterator it;
    for (it = faces.begin(); it != faces.end(); ++it)
    {
        Face *face = *it;
        for (int i = 0; i < 4; i++)
        {
            Vertex *v0 = face->getNodeAt((i + 0) % 4);
            Vertex *v1 = face->getNodeAt((i + 1) % 4);
            Edge edge(v0, v1);
            int found = 0;
            for (size_t j = 0; j < edges.size(); j++)
            {
                if (edges[j].isSameAs(edge))
                {
                    found = 1;
                    break;
                }
            }
            if (!found) edges.push_back(edge);
        }
    }

    size_t V = inner_nodes.size() + bound_nodes.size();
    size_t E = edges.size();
    size_t F = faces.size();

    if (V - E + F == 1) return 1;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::reorient_4_sided_loop()
{
    //Always remeshable, nothing has to be done...
    if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3])) return 0;

    //////////////////////////////////////////////////////////////////////////
    // Defination:  A four sided convex loop has four segments.
    // Objectives:  A four sided convex loop must be orietned such that
    //   1.  First segment must be smaller than the third one, because
    //       we need to create a triangle patch based at segment#3.  
    //
    //   2.  If there are two choices, then the side having irregular
    //       node must be given higher priority. ( Here irregulaty means
    //       that vertex valency < 4 ).
    //
    //   3.  A side having less number of nodes on the first segment than
    //       the third is given preference.
    // 
    // Pre-Conditions  :  A loop must be oriented ( CW or CCW ).
    //
    // Date: 17th Nov. 2010.
    //////////////////////////////////////////////////////////////////////////

    Vertex *start_corner = NULL;

    if (segSize[0] == segSize[2])
    {
        if (min(segSize[1], segSize[3]) == 2) return 1;
        //  Either Segment 2 or 3 must be starting node
        if (segSize[1] < segSize[3])
            start_corner = bound_nodes[ cornerPos[1] ];
        else
            start_corner = bound_nodes[ cornerPos[3] ];
        start_boundary_loop_from(start_corner);
    }

    if (min(segSize[0], segSize[2]) == 2) return 1;

    if (segSize[2] < segSize[0])
    {
        start_corner = bound_nodes[ cornerPos[2] ];
        start_boundary_loop_from(start_corner);
    }

    // By this stage, the loop must be reoriented correctly.
    assert(segSize[0] < segSize[2]);

    // Great, we found one irregular node on the first boundary...
    if (has_irregular_node_on_first_segment()) return 1;

    // If the segment 2 and 3 have same size, Alas, nothing can be done. 
    if (segSize[1] == segSize[3]) return 1;

    if (min(segSize[1], segSize[3]) == 2) return 1;

    if (segSize[3] < segSize[1])
    {
        start_corner = bound_nodes[ cornerPos[3] ];
        start_boundary_loop_from(start_corner);
    }
    else
    {
        start_corner = bound_nodes[ cornerPos[1] ];
        start_boundary_loop_from(start_corner);
    }

    //
    // Note that we didn't check for irregular node here. So if this segment
    // has at least one irregular node, then we are lucky. Otherwise decision
    // to remesh it done based wthere remeshing will result in the reduction 
    // of irregular nodes in patch.
    //
    assert(segSize[0] < segSize[2]);
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh_3_sided_patch()
{
    Vertex *c0, *c1, *c2;

    c0 = bound_nodes[ cornerPos[0] ];
    c1 = bound_nodes[ cornerPos[1] ];
    c2 = bound_nodes[ cornerPos[2] ];

    NodeSequence anodes, bnodes, cnodes;
    anodes = get_bound_nodes(c0, c1);
    bnodes = get_bound_nodes(c1, c2);
    cnodes = get_bound_nodes(c2, c0);

    int err;
    err = remesh_tri_loop(mesh, anodes, bnodes, cnodes, NULL, newnodes, newfaces);
    return err;
}
//////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh_4_sided_patch()
{
    int err;
    size_t nSize;
    Vertex *c0, *c1, *c2, *c3, *ta, *tb, *tc;

    // Collect landmark nodes on the patch..
    c0 = bound_nodes[ cornerPos[0] ];
    c1 = bound_nodes[ cornerPos[1] ];
    c2 = bound_nodes[ cornerPos[2] ];
    c3 = bound_nodes[ cornerPos[3] ];
    ta = quad_splitting_node;

    if (ta == NULL)
    {
        if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3]))
        {
            NodeSequence anodes = get_bound_nodes(c0, c1);
            NodeSequence bnodes = get_bound_nodes(c1, c2);
            NodeSequence cnodes = get_bound_nodes(c2, c3);
            NodeSequence dnodes = get_bound_nodes(c3, c0);
            err = remesh_quad_loop(mesh, anodes, bnodes, cnodes, dnodes, newnodes, newfaces, 0);
            return err;
        }
    }

    assert(ta != NULL);

    NodeSequence a1nodes, a2nodes, b1nodes, c0nodes, c1nodes, c2nodes,
            abnodes, canodes, bcnodes, d1nodes;

    a1nodes = get_bound_nodes(c0, ta);
    a2nodes = get_bound_nodes(ta, c1);

    tb = bound_nodes[cornerPos[2] + a2nodes.size() - 1 ];
    tc = bound_nodes[cornerPos[3] - a1nodes.size() + 1];
    c0nodes = get_bound_nodes(tb, tc);

    abnodes = linear_interpolation(ta, tb, segSize[1]);

    nSize = abnodes.size();
    for (size_t  i = 1; i < nSize - 1; i++)
    {
        mesh->addNode(abnodes[i]);
        newnodes.push_back(abnodes[i]);
    }

    canodes = linear_interpolation(tc, ta, segSize[3]);

    nSize = canodes.size();
    for (size_t i = 1; i < nSize - 1; i++)
    {
        mesh->addNode(canodes[i]);
        newnodes.push_back(canodes[i]);
    }

    bcnodes = get_bound_nodes(tb, tc);
    b1nodes = get_bound_nodes(c1, c2);
    c2nodes = get_bound_nodes(c2, tb);
    c1nodes = get_bound_nodes(tc, c3);
    d1nodes = get_bound_nodes(c3, c0);

    NodeSequence nnodes;
    FaceSequence nfaces;

    err = remesh_tri_loop(mesh, abnodes, bcnodes, canodes, NULL, nnodes, nfaces, 0);
    if (!err)
    {
        nSize = nnodes.size();
        for (size_t i = 0; i < nSize; i++)
            newnodes.push_back(nnodes[i]);

        nSize = nfaces.size();
        for (size_t i = 0; i < nSize; i++)
            newfaces.push_back(nfaces[i]);
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a1nodes, canodes, c1nodes, d1nodes, nnodes, nfaces, 0);
        if (!err)
        {
            nSize = nnodes.size();
            for (size_t i = 0; i < nSize; i++)
                newnodes.push_back(nnodes[i]);

            nSize = nfaces.size();
            for (size_t i = 0; i < nSize; i++)
                newfaces.push_back(nfaces[i]);
        }
    }

    if (!err)
    {
        err = remesh_quad_loop(mesh, a2nodes, b1nodes, c2nodes, abnodes, nnodes, nfaces, 0);
        if (!err)
        {
            nSize = nnodes.size();
            for (size_t i = 0; i < nSize; i++)
                newnodes.push_back(nnodes[i]);

            nSize = nfaces.size();
            for (size_t i = 0; i < nSize; i++)
                newfaces.push_back(nfaces[i]);
        }
    }

    return err;

    return 0;
}

//////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh_5_sided_patch()
{

    Vertex *c0, *c1, *c2, *c3, *c4;
    c0 = bound_nodes[ cornerPos[0] ];
    c1 = bound_nodes[ cornerPos[1] ];
    c2 = bound_nodes[ cornerPos[2] ];
    c3 = bound_nodes[ cornerPos[3] ];
    c4 = bound_nodes[ cornerPos[4] ];

    NodeSequence anodes, bnodes, cnodes, dnodes, enodes;
    anodes = get_bound_nodes(c0, c1);
    bnodes = get_bound_nodes(c1, c2);
    cnodes = get_bound_nodes(c2, c3);
    dnodes = get_bound_nodes(c3, c4);
    enodes = get_bound_nodes(c4, c0);

    int err;
    err = remesh_penta_loop(mesh, anodes, bnodes, cnodes, dnodes, enodes, NULL, newnodes, newfaces);

    return err;
}

//////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::is_quad_breakable_at(const Vertex *ta)
{

    // Landmark nodes on the boundary. There are seven such nodes.
    Vertex *c0, *c1, *c2, *c3, *tb, *tc;
    c0 = bound_nodes[ cornerPos[0] ];
    c1 = bound_nodes[ cornerPos[1] ];
    c2 = bound_nodes[ cornerPos[2] ];
    c3 = bound_nodes[ cornerPos[3] ];

    NodeSequence a1nodes, a2nodes, c0nodes;
    int segments[3], partSegments[6];

    a1nodes = get_bound_nodes(c0, ta);
    a2nodes = get_bound_nodes(ta, c1);
    tb = bound_nodes[cornerPos[2] + a2nodes.size() - 1 ];
    tc = bound_nodes[cornerPos[3] - a1nodes.size() + 1];
    c0nodes = get_bound_nodes(tb, tc);
    segments[0] = c0nodes.size() - 1;
    segments[1] = segSize[3] - 1;
    segments[2] = segSize[2] - 1;

    if (Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments))
        return 0;

    return 1;
}
//////////////////////////////////////////////////////////////////////////////

bool OneDefectPatch::is_4_sided_convex_loop_quad_meshable()
{
    assert(corners.size() == 4);

    // It is always quadrangulable.
    if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3])) return 1;

    int nr0 = count_irregular_nodes(0);
    int nr1 = count_irregular_nodes(1);

    //
    // If there are only two irregular nodes in the patch ( inner + boundary)
    // then it may be quadrangable, but it will not reduce number of irregular
    // nodes ( one irregular node and one irregular node on the boundary 
    // will be created. Since there is no reduction in #of irregular nodes,
    // it is not profitable, and we just return it.
    //
    // If there are more than two irrgular nodes, then we can ensure that
    // that there will at most two irregular nodes in the patch ( one inside
    // and probably one on the boundary.
    //

    if (nr0 + nr1 == 2) return 0;

    // The loop is reoriented so that the first edge is smaller than the
    // opposite edge ....
    if (reorient_4_sided_loop()) return 0;

    // We need at least 3 points on the first seqment ...
    if (segSize[0] < 3) return 0;

    // Landmark nodes on the boundary. There are seven such nodes.
    Vertex *c0, *c1, *c2, *c3, *ta, *tb, *tc;
    c0 = bound_nodes[ cornerPos[0] ];
    c1 = bound_nodes[ cornerPos[1] ];
    c2 = bound_nodes[ cornerPos[2] ];
    c3 = bound_nodes[ cornerPos[3] ];

    // First let us prioritieze vertices on the first segment that will
    // be tried for splitting the region. First of all, if #of irregular
    // nodes inside the region is 2, then we must have one irregular
    // node on the boundary. But even if the #irregular nodes inside
    // the region > 2, we still prefer irregular node in the boundary
    // as it will reduce irregularity by one.

    // Both end nodes on the segment are not included in splitting.

    // Method: Start node(1,...n-1) and check it could lead to acceptable
    // splitting. A proper splitting means that the triangle sandwitched
    // between two quads regions is "quad-meshable" with one defect.
    //
    // Future Work: It is possible that the triangle region is not "one
    // defect" meshable, but with more than one defect. We will extend
    // this case later.
    //

    NodeSequence trynodes;
    FaceSequence vfaces;
    for (int i = 1; i < segSize[0] - 1; i++)
    {
        vfaces = bound_nodes[i]->getRelations2();
        if (vfaces.size() == 3)
            trynodes.push_back(bound_nodes[i]);
    }

    if (trynodes.empty() && nr0 == 2)
    {
        // Patch not remeshable because there are no irregular node on the boundary
        return 0;
    }

    if (nr0 > 2)
    {
        for (int i = 1; i < segSize[0] - 1; i++)
        {
            vfaces = bound_nodes[i]->getRelations2();
            if (vfaces.size() == 4)
                trynodes.push_back(bound_nodes[i]);
        }
    }

    NodeSequence a1nodes, a2nodes, c0nodes;
    int segments[3], partSegments[6];

    int nSize = trynodes.size();
    for (int i = 0; i < nSize; i++)
    {
        ta = trynodes[i];
        a1nodes = get_bound_nodes(c0, ta);
        a2nodes = get_bound_nodes(ta, c1);

        tb = bound_nodes[cornerPos[2] + a2nodes.size() - 1 ];
        tc = bound_nodes[cornerPos[3] - a1nodes.size() + 1];

        c0nodes = get_bound_nodes(tb, tc);
        segments[0] = c0nodes.size() - 1;
        segments[1] = segSize[3] - 1;
        segments[2] = segSize[2] - 1;

        if (Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments))
        {
            FaceSequence vfaces = ta->getRelations2();
            quad_splitting_node = ta;
            quad_splitting_node_degree = vfaces.size();
            return 1;
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::init_blob()
{
    inner_nodes.clear();
    faces.clear();

    FaceSequence vfaces;
    size_t nSize;
    if (apex)
    {
        inner_nodes.insert(apex);
        vfaces = apex->getRelations2();
        nSize  = vfaces.size();
        for (size_t i = 0; i < nSize; i++)
            faces.insert(vfaces[i]);
    }
    else
    {
        nSize = nodepath.size();
        for (size_t i = 0; i < nSize; i++)
        {
            Vertex *v = nodepath[i];
            inner_nodes.insert(v);
            vfaces = v->getRelations2();
            for (size_t i = 0; i < vfaces.size(); i++)
                faces.insert(vfaces[i]);
        }
    }

    // Not needed, used for debugging...
    seednodes.clear();
    set<Vertex*>::const_iterator nit;
    for (nit = inner_nodes.begin(); nit != inner_nodes.end(); ++nit)
        seednodes.push_back(*nit);

    seedfaces.clear();
    set<Face*>::const_iterator fit;
    for (fit = faces.begin(); fit != faces.end(); ++fit)
        seedfaces.push_back(*fit);

    create_boundary();

    nSize = bound_nodes.size();
    for (size_t i = 0; i < nSize; i++)
    {
        vfaces = bound_nodes[i]->getRelations2();
        if (vfaces.size() != 4)
        {
            expand_blob(bound_nodes[i]);
            break;
        }
    }
}

////////////////////////////////////////////////////////////////////

size_t OneDefectPatch::count_irregular_nodes(int where)
{
    set<Vertex*>::const_iterator it;
    assert(mesh->getAdjTable(0, 2));

    size_t ncount = 0;

    if (where == 0)
    {
        for (it = inner_nodes.begin(); it != inner_nodes.end(); ++it)
        {
            Vertex *v = *it;
            if (v->getRelations2().size() != 4) ncount++;
        }
    }
    else
    {
        for (size_t i = 0; i < bound_nodes.size(); i++)
        {
            Vertex *v = bound_nodes[i];
            if (!v->isBoundary() && v->getRelations2().size() == 3) ncount++;
        }
    }

    return ncount;
}
////////////////////////////////////////////////////////////////////

int
OneDefectPatch::create_boundary()
{
    corners.clear();
    boundary.clear();
    Vertex *vertex;

    assert(!faces.empty());

    // We need to rebuild relations locally to identfy corners and boundary.
    FaceSet::const_iterator fiter;
    std::map<Vertex*, FaceSet> relations02;

    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        Face *face = *fiter;
        for (int j = 0; j < face->getSize(0); j++)
        {
            vertex = face->getNodeAt(j);
            relations02[vertex].insert(face);
        }
    }

    // A boundary edge must have exactly one face neighbor...
    FaceSequence faceneighs;
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        Face *face = *fiter;
        for (int j = 0; j < 4; j++)
        {
            Vertex *v0 = face->getNodeAt((j + 0) % 4);
            Vertex *v1 = face->getNodeAt((j + 1) % 4);
            faceneighs.clear();
            assert(relations02[v0].size() > 0);
            assert(relations02[v1].size() > 0);
            set_intersection(relations02[v0].begin(), relations02[v0].end(),
                             relations02[v1].begin(), relations02[v1].end(),
                             back_inserter(faceneighs));
            if (faceneighs.size() == 1)
            {
                Edge newedge(v0, v1);
                boundary.push_back(newedge);
            }
        }
    }

    // Sequence the chain and start from one of the corner...
    int err = Mesh::make_chain(boundary);
    if (err) return 2;

    //
    // Identify corners in the mesh.
    // Should we check only one vertex per edge ?
    //
    set<Face*> neighs;

    size_t nSize = boundary.size();
    for (size_t k = 0; k < nSize; k++)
    {
        vertex = boundary[k].getNodeAt(0);
        neighs = relations02[vertex];
        if (neighs.size() == 1) corners.insert(vertex);

        vertex = boundary[k].getNodeAt(1);
        neighs = relations02[vertex];
        if (neighs.size() == 1) corners.insert(vertex);
    }

    // We don't handle these cases...
    //   if (corners.size() < 3 || corners.size() > 5) return 2;

    // Start the chain from one of the corners.
    err = Mesh::rotate_chain(boundary, *corners.begin());
    if (err) return 3;

    // Collect the sequence of boundary nodes...,
    bound_nodes.resize( nSize );
    for (size_t k = 0; k < nSize; k++)
        bound_nodes[k] = boundary[k].getNodeAt(0); // Only the first node.

    //
    // Collect the inner nodes of the patch. These nodes will be deleted, if
    // the remesh operation is successful...
    //
    inner_nodes.clear();
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        Face *face = *fiter;
        for (int j = 0; j < 4; j++)
        {
            Vertex *v = face->getNodeAt(j);
            if (find(bound_nodes.begin(), bound_nodes.end(), v) == bound_nodes.end())
                inner_nodes.insert(v);
        }
    }

    // Split the boundary loop into segments.
    // (i.e. End of the segments are the corners identified earlier )
    set_boundary_segments();

    return 0;
}

////////////////////////////////////////////////////////////////////

void OneDefectPatch::set_boundary_segments()
{
    size_t nSize = corners.size();
    // Although this stage will not come in this algorithm...
    if ( nSize == 0) return;

    cornerPos.resize( nSize + 1);

    set<Vertex*>::const_iterator it;
    int index = 0;
    for (it = corners.begin(); it != corners.end(); ++it)
    {
        cornerPos[index++] = getPosOf(*it);
    }

    cornerPos[corners.size()] = bound_nodes.size();

    sort(cornerPos.begin(), cornerPos.end());

    segSize.resize( nSize );

    for (size_t i = 0; i < nSize; i++)
        segSize[i] = cornerPos[(i + 1)] - cornerPos[i] + 1;
}

////////////////////////////////////////////////////////////////////

int OneDefectPatch::get_topological_outer_angle(Vertex *vertex)
{
    FaceSequence vfaces = vertex->getRelations2();

    // How many faces are outside the regions.
    int ncount = 0;
    for (size_t i = 0; i < vfaces.size(); i++)
        if (faces.find(vfaces[i]) == faces.end()) ncount++;

    return ncount - 2;
}
////////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::expand_blob(Vertex *vertex)
{
    FaceSequence vfaces = vertex->getRelations2();
    for (size_t i = 0; i < vfaces.size(); i++)
        faces.insert(vfaces[i]);
}
////////////////////////////////////////////////////////////////////////////////

int
OneDefectPatch::build_remeshable_boundary()
{
    assert(mesh->getAdjTable(0, 2));

    /////////////////////////////////////////////////////////////////////////////
    // A valid boundary patch has the following propeties:
    //
    // 1) The loop is closed.
    // 2) The region is simply connected. The Euler characteristic must be
    //    one, which mean there are no holes in the region.
    // 3) The loop has 3, 4 or 5 corner nodes.
    // 4) The loop form almost convex boundary.
    //
    // If these conditions are met, we may try to remesh the region with
    // quadrilateral elements. Sometimes it may not be possible to remesh
    // a reason, or the resulting elements have unacceptable quality that
    // we would like to avoid.
    //
    // A resulting mesh MUST decrease the irregular nodes count.
    //////////////////////////////////////////////////////////////////////////
    init_blob();

    if (count_irregular_nodes(0) < 2) return 1;

    int err;
    while (1)
    {
        err = create_boundary();
        if (err) return 2;

        // Should not cover the entire original mesh.
        if (faces.size() == mesh->getSize(2)) return 3;

        // Expand the concave blob.
        int topo_convex_region = 1;
        size_t nSize = bound_nodes.size();
        for (size_t i = 0; i < nSize; i++)
        {
            if (!bound_nodes[i]->isBoundary())
            {
                int topo_angle = get_topological_outer_angle(bound_nodes[i]);
                if (topo_angle < 0)
                {
                    expand_blob(bound_nodes[i]);
                    topo_convex_region = 0;
                }
            }
        }

        if (topo_convex_region)
        {

            int nsides = corners.size();

            // We restrict ourselves with 3-4-5 sides patches.
            if (nsides > 5) return 6;

            int segments[5], partSegments[10];

            int meshable = 0;
            switch (nsides)
            {
            case 3:
                segments[0] = segSize[0] - 1;
                segments[1] = segSize[1] - 1;
                segments[2] = segSize[2] - 1;
                meshable = Face::is_3_sided_convex_loop_quad_meshable(segments, partSegments);
                break;
            case 4:
                meshable = is_4_sided_convex_loop_quad_meshable();
                break;

            case 5:
                segments[0] = segSize[0] - 1;
                segments[1] = segSize[1] - 1;
                segments[2] = segSize[2] - 1;
                segments[3] = segSize[3] - 1;
                segments[4] = segSize[4] - 1;
                meshable = Face::is_5_sided_convex_loop_quad_meshable(segments, partSegments);
                break;
            }

            if (meshable) return 0;

            // Bad method for the time being.
            for (size_t i = 0; i < bound_nodes.size(); i++)
                expand_blob(bound_nodes[i]);
        }
    }

    // You should never come at this stage...
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::backup()
{
    // Node Coordinates ...
    set<Vertex*>::const_iterator it;

    backupCoords.reserve(inner_nodes.size());

    backupCoords.clear();
    for (it = inner_nodes.begin(); it != inner_nodes.end(); ++it)
    {
        Vertex *v = *it;
        assert(!v->isRemoved());
        backupCoords.push_back(v->getXYZCoords());
    }

    // Face Connectivity information ...
    backupConnect.clear();
    set<Face*>::const_iterator fiter;

#ifdef SEQUENCE_IS_VECTOR
    backupConnect.reserve(4 * faces.size());
#endif

    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        Face *face = *fiter;
        assert(!face->isRemoved());
        backupConnect.push_back(face->getNodeAt(0));
        backupConnect.push_back(face->getNodeAt(1));
        backupConnect.push_back(face->getNodeAt(2));
        backupConnect.push_back(face->getNodeAt(3));
    }
}

////////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::rollback()
{
    assert(!inner_nodes.empty());

    set<Vertex*>::const_iterator it;
    int index = 0;
    for (it = inner_nodes.begin(); it != inner_nodes.end(); ++it)
    {
        Vertex *v = *it;
        v->setXYZCoords(backupCoords[index++]);
        mesh->reactivate(v);
        assert(!v->isRemoved());
        assert(mesh->contains(v));
    }

    assert(!faces.empty());

    set<Face*>::const_iterator fiter;
    NodeSequence connect(4);
    index = 0;
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
    {
        connect[0] = backupConnect[index++];
        connect[1] = backupConnect[index++];
        connect[2] = backupConnect[index++];
        connect[3] = backupConnect[index++];
        for (int i = 0; i < 4; i++)
            assert(!connect[i]->isRemoved());
        Face *face = *fiter;
        face->setNodes(connect);
        mesh->reactivate(face);
        assert(!face->isRemoved());
        assert(mesh->contains(face));
    }

}

////////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::pre_remesh()
{
    //
    // Before remeshing deactivate all the internal faces and nodes so that all
    // the relations are cleared. They are not removed because if something goes
    // wrong, the data structures are restored using rollback operations.
    //
    FaceSet::const_iterator fiter;
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
        mesh->deactivate(*fiter);

    irregular_nodes_removed.clear();
    NodeSet::const_iterator niter;
    for (niter = inner_nodes.begin(); niter != inner_nodes.end(); ++niter)
    {
        if (!QuadCleanUp::isRegular(*niter)) irregular_nodes_removed.push_back(*niter);
        mesh->deactivate(*niter);
    }

    // At this stage, we will have an empty patch ...
}
////////////////////////////////////////////////////////////////////////////////

void OneDefectPatch::post_remesh()
{
    FaceSet::const_iterator fiter;
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
        mesh->remove(*fiter);

    NodeSet::const_iterator niter;
    for (niter = inner_nodes.begin(); niter != inner_nodes.end(); ++niter)
        mesh->remove(*niter);

    faces.clear();
    for (size_t i = 0; i < newfaces.size(); i++)
        faces.insert(newfaces[i]);
    newfaces.clear();

    inner_nodes.clear();
    for (size_t i = 0; i < newnodes.size(); i++)
        inner_nodes.insert(newnodes[i]);
    newnodes.clear();

    int ncount = 0;
    new_defective_node = NULL;
    for (niter = inner_nodes.begin(); niter != inner_nodes.end(); ++niter)
    {
        Vertex *v = *niter;
        if (!QuadCleanUp::isRegular(v))
        {
            ncount++;
            new_defective_node = v;
        }
    }
    assert(ncount < 2);

    assert(mesh->is_consistently_oriented());

    // Perform some local smoothing and ensure that elements are valid.
    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(100);

    NodeSet nset0 = inner_nodes;
    size_t nSize = bound_nodes.size();
    for (size_t i = 0; i < nSize; i++)
        nset0.insert(bound_nodes[i]);

    NodeSet nset1;

    double minq;
    vector<double> quality(faces.size());
    NodeSequence localnodes;

    for (niter = nset0.begin(); niter != nset0.end(); ++niter)
    {
        nset1.insert(*niter);
        NodeSequence vneighs = (*niter)->getRelations0();
        for (size_t j = 0; j < vneighs.size(); j++)
            nset1.insert(vneighs[j]);
    }

    localnodes.clear();
    for (niter = nset1.begin(); niter != nset1.end(); ++niter)
        localnodes.push_back(*niter);

    lapsmooth.localized_at(localnodes);

    int index = 0;
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter)
        quality[index++] = (*fiter)->getAspectRatio();
    minq = *min_element(quality.begin(), quality.end());

    cout << "Min element quality of the patch " << minq << endl;

    if (minq < 0.2) 
       lapsmooth.execute();
}
////////////////////////////////////////////////////////////////////////////////

int OneDefectPatch::remesh()
{
    // Check if any node or element is modifed by previous operations.
    if (!isSafe()) return 1;

    if (!is_simply_connected()) return 2;

    newnodes.clear();
    newfaces.clear();

    if (corners.size() > 5) cout << "Warning: Corners more than 5 " << endl;

    // backup all nodes and faces within the patch.
    backup();

    int nirregular0 = count_irregular_nodes(0) + count_irregular_nodes(1);

    // Do some pre-processing before remeshing the patch
    pre_remesh();

    int err;
    switch (corners.size())
    {
    case 3:
        err = remesh_3_sided_patch();
        break;
    case 4:
        err = remesh_4_sided_patch();
        break;
    case 5:
        err = remesh_5_sided_patch();
        break;
    }

    if (err)
    {
        for (size_t i = 0; i < newfaces.size(); i++)
            mesh->remove(newfaces[i]);
        newfaces.clear();

        for (size_t i = 0; i < newnodes.size(); i++)
            mesh->remove(newnodes[i]);
        newnodes.clear();

        rollback();
        return 1;
    }

    // Do some post-processing of the patch
    post_remesh();

    int nirregular1 = count_irregular_nodes(0) + count_irregular_nodes(1);

    // Our method must stricly decrease the number of irregular nodes,
    // otherwise it will be go into infinite loop.

    assert(nirregular1 < nirregular0);

    return err;
}

////////////////////////////////////////////////////////////////////////////////

OneDefectPatch *QuadCleanUp::build_one_defect_patch(Vertex *vertex)
{
    OneDefectPatch *patch = NULL;

    if( vertex == NULL ) {
        size_t numNodes = mesh->getSize(0);
        for( size_t i = 0; i < numNodes; i++) {
             size_t id = Math::random_value( size_t(0), numNodes-1);
             cout << " Search " << i << " From " << numNodes << endl;
             Vertex *v = mesh->getNodeAt(id);
             patch =  build_one_defect_patch( v );
             if( patch) return patch;
        }
        cout << "Sorry: Couldn't find any defective remeshable patch" << endl;
        return NULL;
    }

    if (isRegular(vertex)) return NULL;

    MeshFilter *filter = new FirstIrregularNode();
    DijkstraShortestPath djk(mesh, 1);

    NodeSequence nodepath = djk.getPath(vertex, filter);
    delete filter;

    if (isRegular(nodepath.front())) return NULL;
    if (isRegular(nodepath.back())) return NULL;

    patch = new OneDefectPatch(mesh, nodepath);
    int err = patch->build_remeshable_boundary();

    if (err)
    {
        delete patch;
        return NULL;
    }

    return patch;
}

////////////////////////////////////////////////////////////////////////////////
/*
vector<OneDefectPatch>
QuadCleanUp::search_one_defect_patches()
{
    vDefectPatches.clear();

    int relexist = mesh->build_relations(0, 2);
    mesh->search_boundary();

    size_t numnodes = mesh->getSize(0);

    region_search_method = 1;

    assert(mesh->getAdjTable(0, 2));

    int err;

    if (region_search_method == 0)
    {
        vector<QTrack> qpath = Jaal::generate_quad_irregular_graph(mesh);
        for (int i = 0; i < qpath.size(); i++)
        {
           OneDefectPatch patch(mesh, qpath[i].sequence);
           err = patch.build_remeshable_boundary();
           if (!err) vDefectPatches.push_back(patch);
        }
    }

    if (region_search_method == 1)
    {
        MeshFilter *filter =  new FirstIrregularNode();
        DijkstraShortestPath djk(mesh, 1);

        if (irregular_nodes_set.empty()) build_irregular_nodes_set();
       
        NodeSet ::const_iterator it;
        for (it = irregular_nodes_set.begin(); it != irregular_nodes_set.end(); ++it)
        {
            Vertex *vertex1 = *it;
            if( !isRegular( vertex1 )  ) 
            {
                NodeSequence nodeseq = mesh->get_breadth_first_ordered_nodes(vertex1, filter);
                if( !nodeseq.empty()) {
                     Vertex *vertex2 = nodeseq.back();
                     assert( !isRegular( vertex2 ) );
                     NodeSequence nodepath = djk.getPath(vertex1, vertex2);
                     OneDefectPatch patch(mesh, nodepath);
                     err = patch.build_remeshable_boundary();
                     if (!err)
                     {
                         vDefectPatches.push_back(patch); // Only one pickup.
                         break;
                     }
                  }
                }
            }
    }

    if (!relexist)
        mesh->clear_relations(0, 2);

    return vDefectPatches;
}
 */

/////////////////////////////////////////////////////////////////////////////

int QuadCleanUp::remesh_defective_patches()
{
    int relexist0 = mesh->build_relations(0, 0);
    int relexist2 = mesh->build_relations(0, 2);

    mesh->search_boundary();

    assert(mesh->getAdjTable(0, 2));

    cout << "# of Irregular nodes before cleanup "
            << mesh->count_irregular_nodes(4) << endl;

    LaplaceSmoothing lapsmooth(mesh);
    LaplaceWeight *lapweight = NULL;
    lapsmooth.setMethod(0);
    lapweight = new LaplaceAreaWeight();
    lapsmooth.setWeight(lapweight);
    lapsmooth.setNumIterations(5); // Don't do excess smoothing at every stage.

    region_search_method = 1;
    lapsmooth.setNumIterations(5); // Don't do excess smoothing at every stage.

    OneDefectPatch *patch;
    while (1)
    {
        build_irregular_nodes_set();

        size_t ncount = 0;
        while (1)
        {
            if (irregular_nodes_set.empty()) break;
            Vertex *vertex = *irregular_nodes_set.begin();
            irregular_nodes_set.erase(vertex);
            patch = build_one_defect_patch(vertex);
            if (patch)
            {
                int err = patch->remesh();
                if (!err)
                {
                    NodeSequence vrem = patch->get_irregular_nodes_removed();
                    size_t nSize = vrem.size();
                    for (size_t j = 0; j < nSize; j++)
                        irregular_nodes_set.erase(vrem[j]);
                    Vertex *vtx = patch->get_new_defective_node();
                    if (vtx) irregular_nodes_set.insert(vtx);
                    ncount++;
                    lapsmooth.execute();
                }
                delete patch;
            }
        }
        if (ncount == 0) break;
    }
    mopt.shape_optimize(mesh);


    if (!relexist0)
        mesh->clear_relations(0, 0);

    if (!relexist2)
        mesh->clear_relations(0, 2);

    delete lapweight;

    return 0;
}
///////////////////////////////////////////////////////////////////////

