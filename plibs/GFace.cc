#include "GFace.h"
#include <sstream>
#include <fstream>
#include <limits>
#include <set>

void iGeom_getEntAdj_Ext(iGeom_Instance instance,
                         iBase_EntityHandle entity_handle,
                         int to_dimension,
                         std::vector<iBase_EntityHandle> &entities,
                         int *err)
{
    SimpleArray<iBase_EntityHandle> edges;
    iGeom_getEntAdj(instance, entity_handle, iBase_EDGE, ARRAY_INOUT(edges), err);

    // First Edge is always correct. It is the reference.
    set<iBase_EntityHandle> edgeset;
    for (int i = 1; i < edges.size(); i++) edgeset.insert(edges[i]);

    iBase_EntityHandle start_vertex, next_vertex;
    SimpleArray<iBase_EntityHandle> edgenodes;
    iGeom_getEntAdj(instance, edges[0], iBase_VERTEX, ARRAY_INOUT(edgenodes), err);

    start_vertex = edgenodes[0];
    if (edgenodes.size() > 1) next_vertex = edgenodes[1];

    entities.clear();
    entities.push_back(start_vertex);

    iBase_EntityHandle curr_edge;
    for (int i = 1; i < edges.size(); i++)
    {
        entities.push_back(next_vertex);

        BOOST_FOREACH(curr_edge, edgeset)
        {
            iGeom_getEntAdj(instance, curr_edge, iBase_VERTEX, ARRAY_INOUT(edgenodes), err);
            if (edgenodes[0] == next_vertex)
            {
                edges[i] = curr_edge;
                next_vertex = edgenodes[1];
                edgeset.erase(curr_edge);
                break;
            }
            if (edgenodes[1] == next_vertex)
            {
                edges[i] = curr_edge;
                next_vertex = edgenodes[0];
                edgeset.erase(curr_edge);
                break;
            }
        }
    }
    assert(next_vertex == start_vertex);
    assert(entities.size() == edges.size());

    if (to_dimension == iBase_VERTEX) return;

    entities.clear();
    for (int i = 0; i < edges.size(); i++)
        entities.push_back(edges[i]);

    if (to_dimension == iBase_EDGE) return;

    cout << "Warning: Not yet implemented " << endl;
    exit(0);
}

///////////////////////////////////////////////////////////////////////////////

int GFace::get_hex_face_number(const SimpleArray<iBase_EntityHandle> &cellnodes,
                               const SimpleArray<iBase_EntityHandle> &facenodes,
                               int &side_no, int &orientation)
{
    static vector< vector<iBase_EntityHandle> > cellfaces(6);
    for (int i = 0; i < 6; i++) cellfaces[i].resize(4);

    cellfaces[0][0] = cellnodes[0];
    cellfaces[0][1] = cellnodes[3];
    cellfaces[0][2] = cellnodes[7];
    cellfaces[0][3] = cellnodes[4];

    cellfaces[1][0] = cellnodes[1];
    cellfaces[1][1] = cellnodes[2];
    cellfaces[1][2] = cellnodes[6];
    cellfaces[1][3] = cellnodes[5];

    cellfaces[2][0] = cellnodes[0];
    cellfaces[2][1] = cellnodes[1];
    cellfaces[2][2] = cellnodes[5];
    cellfaces[2][3] = cellnodes[4];

    cellfaces[3][0] = cellnodes[3];
    cellfaces[3][1] = cellnodes[2];
    cellfaces[3][2] = cellnodes[6];
    cellfaces[3][3] = cellnodes[7];

    cellfaces[4][0] = cellnodes[0];
    cellfaces[4][1] = cellnodes[1];
    cellfaces[4][2] = cellnodes[2];
    cellfaces[4][3] = cellnodes[3];

    cellfaces[5][0] = cellnodes[4];
    cellfaces[5][1] = cellnodes[5];
    cellfaces[5][2] = cellnodes[6];
    cellfaces[5][3] = cellnodes[7];

    side_no = -1;
    vector<iBase_EntityHandle>::const_iterator istart, iend;
    for (int i = 0; i < 6; i++)
    {
        istart = cellfaces[i].begin();
        iend = cellfaces[i].end();
        int found = 1;
        for (int j = 0; j < 4; j++)
        {
            if (std::find(istart, iend, facenodes[j]) == iend)
            {
                found = 0;
                break;
            }
        }
        if (found) side_no = i;
    }

    assert(side_no >= 0 && side_no < 6);

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////

GFace::GFace(iGeom_Instance g, iBase_EntityHandle h,
             iMesh_Instance m, iRel_Instance a, iRel_RelationHandle r)
{
    geometry = g;
    gFaceHandle = h;
    kdtree = NULL;

    int err;
    iGeom_getEntUVRange(geometry, gFaceHandle, &umin, &vmin, &umax, &vmax, &err);
    assert(!err);

    iGeom_getEntBoundBox(geometry, gFaceHandle, &xmin, &ymin, &zmin,
                         &xmax, &ymax, &zmax, &err);

    periodic[0] = 0;
    periodic[1] = 0;
    int in_u, in_v;
    iGeom_isEntPeriodic(geometry, gFaceHandle, &in_u, &in_v, &err);
    assert(!err);
    if (in_u) periodic[0] = 1;
    if (in_v) periodic[1] = 1;

    xlength = fabs(xmax - xmin);
    ylength = fabs(ymax - ymin);
    zlength = fabs(zmax - zmin);
    maxlength = max(xlength, max(ylength, zlength));

    uvCoords = NULL;
    numNeighs = 1;
    annIdx.resize(numNeighs);
    anndist.resize(numNeighs);
    build_kdtree();

    mesh = m;
    assoc = a;
    rel = r;
}

///////////////////////////////////////////////////////////////////////////////

GFace::~GFace()
{
    delete_kdtree();
}

///////////////////////////////////////////////////////////////////////////////

bool GFace::hasSeam() const
{
    int err, sense_out;

    SimpleArray<iBase_EntityHandle> edgeHandles;
    iGeom_getEntAdj(geometry, gFaceHandle, iBase_EDGE, ARRAY_INOUT(edgeHandles), &err);
    assert(!err);

    for (int i = 0; i < edgeHandles.size(); i++)
    {
        iGeom_getEgFcSense(geometry, edgeHandles[i], gFaceHandle, &sense_out, &err);
        if (sense_out == 0) return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void GFace::build_kdtree()
{
    delete_kdtree(); // If Allocated earlier

    int err;

    int N = 1000;
    double x, y, z;

    double du = (umax - umin) / (double) N;
    double dv = (vmax - vmin) / (double) N;

    int numnodes = (N + 1)*(N + 1);
    kdnodes = annAllocPts(numnodes, 3);

    uvCoords = new double[2 * numnodes];

    int index = 0;
    for (int j = 0; j < N + 1; j++)
    {
        double v = vmin + j*dv;
        if (v > vmax) v = vmax;
        for (int i = 0; i < N + 1; i++)
        {
            double u = umin + i*du;
            if (u > umax) u = umax;
            iGeom_getEntUVtoXYZ(geometry, gFaceHandle, u, v, &x, &y, &z, &err);
            kdnodes[index][0] = x;
            kdnodes[index][1] = y;
            kdnodes[index][2] = z;
            uvCoords[2 * index + 0] = u;
            uvCoords[2 * index + 1] = v;
            index++;
        }
    }

    kdtree = new ANNkd_tree(kdnodes, numnodes, 3);
    assert(kdtree);
}

///////////////////////////////////////////////////////////////////////////////

void GFace::delete_kdtree()
{
    if (kdtree == NULL) return;

    annDeallocPts(kdnodes);
    delete [] uvCoords;
    delete kdtree;
    kdtree = NULL;
    annClose();
}

///////////////////////////////////////////////////////////////////////////////

Point2D GFace::getUVRange(int i) const
{
    assert(kdtree);
    static Point2D p2d;
    if (i == 0)
    {
        p2d[0] = umin;
        p2d[1] = umax;
        return p2d;
    }

    p2d[0] = vmin;
    p2d[1] = vmax;
    return p2d;

}

///////////////////////////////////////////////////////////////////////////////

double GFace::getGeodesicLength(const Point2D &u0, const Point2D &uN) const
{
    int nsub = 2;

    int err;
    double u, v, du, dv;

    double x0, y0, z0;
    double x1, y1, z1;
    double dx, dy, dz, dl, arclen, maxlen;

    maxlen = 0.0;
    for (int ilevel = 0; ilevel < 10; ilevel++)
    {
        du = (uN[0] - u0[0]) / (double) nsub;
        dv = (uN[1] - u0[1]) / (double) nsub;

        iGeom_getEntUVtoXYZ(geometry, gFaceHandle, u0[0], u0[1], &x0, &y0, &z0, &err);
        assert(!err);

        arclen = 0.0;
        for (int i = 0; i < nsub; i++)
        {
            u = u0[0] + (i + 1) * du;
            v = u0[1] + (i + 1) * dv;
            iGeom_getEntUVtoXYZ(geometry, gFaceHandle, u, v, &x1, &y1, &z1, &err);
            assert(!err);
            dx = x1 - x0;
            dy = y1 - y0;
            dz = z1 - z0;
            dl = sqrt(dx * dx + dy * dy + dz * dz);
            arclen += dl;
            x0 = x1;
            y0 = y1;
            z0 = z1;
        }
        if (fabs(arclen - maxlen) < 1.0E-06) break;
        maxlen = arclen;
        nsub *= 2;
    }
    return arclen;
}

///////////////////////////////////////////////////////////////////////////////

Vec3D GFace::getNormal(const Point2D &uv) const
{
    assert(kdtree);

    int err;
    double nx, ny, nz;
    iGeom_getEntNrmlUV(geometry, gFaceHandle, uv[0], uv[1],
                       &nx, &ny, &nz, &err);
    assert(!err);

    double mag = 1.0 / sqrt(nx * nx + ny * ny + nz * nz);

    Vec3D avec;
    avec[0] = nx*mag;
    avec[1] = ny*mag;
    avec[2] = nz*mag;

    return avec;
}
///////////////////////////////////////////////////////////////////////////////

std::pair<Vec3D, Vec3D> GFace::getFirstDer(const Point2D &uv) const
{
    int err;
    static SimpleArray<double> uDeriv, vDeriv;
    iGeom_getEnt1stDrvt(geometry, gFaceHandle, uv[0], uv[1],
                        ARRAY_INOUT(uDeriv), ARRAY_INOUT(vDeriv),
                        &err);
    assert(!err);

    Vec3D du;
    du[0] = uDeriv[0];
    du[1] = uDeriv[1];
    du[2] = uDeriv[2];

    Vec3D dv;
    dv[0] = vDeriv[0];
    dv[1] = vDeriv[1];
    dv[2] = vDeriv[2];

    return std::make_pair(du, dv);
}

///////////////////////////////////////////////////////////////////////////////

Point3D GFace::getXYZCoords(const Point2D &uv) const
{
    assert(kdtree);

    static Point3D p3d;
    int err;

    double x, y, z;
    iGeom_getEntUVtoXYZ(geometry, gFaceHandle, uv[0], uv[1], &x, &y, &z, &err);
    assert(!err);

    p3d[0] = x;
    p3d[1] = y;
    p3d[2] = z;

    return p3d;
}

///////////////////////////////////////////////////////////////////////////////

Point2D GFace::getUVCoords(const Point3D &xyz) const
{
    assert(kdtree);

    Point2D uv;
    double tol = 1.0E-06;
    int err;
    double x, y, z, u, v, dx, dy, dz, derr;

    double xq = xyz[0];
    double yq = xyz[1];
    double zq = xyz[2];

    if (xq < xmin || xq > xmax)
        cout << "Warning: Query point outside X Range " << endl;

    if (yq < ymin || yq > ymax)
        cout << "Warning: Query point outside Y Range " << endl;

    if (zq < zmin || zq > zmax)
        cout << "Warning: Query point outside Z Range " << endl;

    iGeom_getEntXYZtoUV(geometry, gFaceHandle, xq, yq, zq, &u, &v, &err);
    assert(!err);

    iGeom_getEntUVtoXYZ(geometry, gFaceHandle, u, v, &x, &y, &z, &err);
    assert(!err);

    dx = fabs(x - xq);
    dy = fabs(y - yq);
    dz = fabs(z - zq);
    derr = dx * dx + dy * dy + dz*dz;

    if (derr < tol * tol)
    {
        uv[0] = u;
        uv[1] = v;
        return uv;
    }

    double queryPoint[3], eps = 0.0;
    double xon, yon, zon;
    double dist1, dist2;

    queryPoint[0] = xq;
    queryPoint[1] = yq;
    queryPoint[2] = zq;

    kdtree->annkSearch(queryPoint, numNeighs, &annIdx[0], &anndist[0], eps);

    int index = annIdx[0];
    dist1 = anndist[0];

    x = kdnodes[index][0];
    y = kdnodes[index][1];
    z = kdnodes[index][2];
    u = uvCoords[2 * index + 0];
    v = uvCoords[2 * index + 1];

    double uguess = u;
    double vguess = v;
    
    iGeom_getEntXYZtoUVHint( geometry, gFaceHandle, x, y, z, &uguess, &vguess, &err);
    iGeom_getEntUVtoXYZ(geometry, gFaceHandle, uguess, vguess, &xon, &yon, &zon, &err);
    dx = queryPoint[0] - xon;
    dy = queryPoint[1] - yon;
    dz = queryPoint[2] - zon;
    dist2 =  dx*dx + dy*dy + dz*dz;

    if( dist2 < dist1 ) {
        u = uguess;
        v = vguess;
    }

    uv[0] = u;
    uv[1] = v;
    return uv;
}

///////////////////////////////////////////////////////////////////////////////

Point2D GFace::getUVCoords(const Point3D &xyz, const Point2D &nearto) const
{
    assert(kdtree);

    Point2D uv;
    double tol = 1.0E-06;
    int err;
    double x, y, z, u, v, dx, dy, dz, derr;

    double xq = xyz[0];
    double yq = xyz[1];
    double zq = xyz[2];

    if (xq < xmin || xq > xmax)
        cout << "Warning: Query point outside X Range " << endl;

    if (yq < ymin || yq > ymax)
        cout << "Warning: Query point outside Y Range " << endl;

    if (zq < zmin || zq > zmax)
        cout << "Warning: Query point outside Z Range " << endl;

    u = nearto[0];
    v = nearto[1];
    iGeom_getEntXYZtoUVHint(geometry, gFaceHandle, xq, yq, zq, &u, &v, &err);
    assert(!err);

    iGeom_getEntUVtoXYZ(geometry, gFaceHandle, u, v, &x, &y, &z, &err);
    assert(!err);

    dx = fabs(x - xq);
    dy = fabs(y - yq);
    dz = fabs(z - zq);
    derr = dx * dx + dy * dy + dz*dz;

    if (derr < tol * tol)
    {
        uv[0] = u;
        uv[1] = v;
        return uv;
    }

    assert(1);

    double queryPoint[3], eps = 0.0;
    double xon, yon, zon;
    double dist1, dist2;

    queryPoint[0] = xq;
    queryPoint[1] = yq;
    queryPoint[2] = zq;

    kdtree->annkSearch(queryPoint, numNeighs, &annIdx[0], &anndist[0], eps);

    int index = annIdx[0];
    dist1 = anndist[0];

    x = kdnodes[index][0];
    y = kdnodes[index][1];
    z = kdnodes[index][2];
    u = uvCoords[2 * index + 0];
    v = uvCoords[2 * index + 1];

    uv[0] = u;
    uv[1] = v;
    return uv;
}

///////////////////////////////////////////////////////////////////////////////

int GFace::getUVCoords(const Point2D &uvstart, double dist, Point2D &uvguess) const
{
    Point2D uv;
    Point3D p0 = getXYZCoords(uvstart);
    Point3D p1 = getXYZCoords(uvguess);

    double dx, dy, dz, dl;
    double u = uvguess[0];
    double v = uvguess[1];

    int ntrials = 0;
    while (1)
    {
        dx = p1[0] - p0[0];
        dy = p1[1] - p0[1];
        dz = p1[2] - p0[2];
        dl = sqrt(dx * dx + dy * dy + dz * dz);

        if (fabs(dl - dist) < 1.0E-15) break;

        u = uvstart[0] + (u - uvstart[0]) * dist / dl;
        v = uvstart[1] + (v - uvstart[1]) * dist / dl;
        uv[0] = u;
        uv[1] = v;

        p1 = getXYZCoords(uv);

        if (ntrials++ == 100)
        {
            cout << " Warning: UV  Search failed" << endl;
            break;
        }
    }

    uvguess[0] = u;
    uvguess[1] = v;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

Point3D GFace::getClosestPoint(const Point3D &qPoint) const
{
    int err;
    Point3D p3d;
    double p[3], eps = 0.0;
    double xon, yon, zon;
    double xon1, yon1, zon1;
    double dx, dy, dz, dist1, dist2;

    dist1 = std::numeric_limits<double>::max();

    if (!kdtree)
    {
        xon = qPoint[0];
        yon = qPoint[1];
        zon = qPoint[2];
    }

    iGeom_getEntClosestPt(geometry, gFaceHandle, qPoint[0], qPoint[1], qPoint[2],
                          &xon, &yon, &zon, &err);

    if (!err)
    {
        dx = qPoint[0] - xon;
        dy = qPoint[1] - yon;
        dz = qPoint[2] - zon;
        dist1 = dx * dx + dy * dy + dz*dz;
        p3d[0] = xon;
        p3d[1] = yon;
        p3d[2] = zon;
        return p3d;
    }

    cout << " Searching from k-d tree " << endl;

    if (kdtree)
    {
        p[0] = qPoint[0];
        p[1] = qPoint[1];
        p[2] = qPoint[2];
        kdtree->annkSearch(p, numNeighs, &annIdx[0], &anndist[0], eps);
        dist2 = anndist[0];

        int index = annIdx[0];
        double u = uvCoords[2 * index + 0];
        double v = uvCoords[2 * index + 1];

        iGeom_getEntUVtoXYZ(geometry, gFaceHandle, u, v, &xon1, &yon1, &zon1, &err);
        dx = qPoint[0] - xon1;
        dy = qPoint[1] - yon1;
        dz = qPoint[2] - zon1;
        dist2 = dx * dx + dy * dy + dz*dz;

        if (dist2 < dist1)
        {
            //          cout << " With iGeom " << dist1 << " With KD Tree : " << dist2 << endl;
            xon = xon1;
            yon = yon1;
            zon = zon1;
        }
    }

    p3d[0] = xon;
    p3d[1] = yon;
    p3d[2] = zon;
    return p3d;
}

///////////////////////////////////////////////////////////////////////////////

void GFace::projectEdgeHigherOrderNodes(const vector<double> &gnodes)
{
    int err;
    gllnodes = gnodes;

    if (mesh == 0)
    {
        cout << " Warning: No mesh is present in geometric face : Higher order nodes projection not done " << endl;
        return;
    }

    if (assoc == 0)
    {
        cout << " Warning: No assoc is present in geometric face : Higher order nodes projection not done " << endl;
        return;
    }

    if (rel == 0)
    {
        cout << " Warning: No relation is present in geometric face : Higher order nodes projection not done " << endl;
        return;
    }

    const char *tag1 = "HO_POINTS";
    iMesh_getTagHandle(mesh, tag1, &horder_tag, &err, strlen(tag1));
    assert(!err);

    int nsize = gllnodes.size();
    arclength_ratio.resize(nsize);
    for (int i = 0; i < nsize; i++)
        arclength_ratio[i] = 0.5 * (1.0 + gllnodes[i]);

    iBase_EntitySetHandle meshSet;
    iRel_getEntSetAssociation(assoc, rel, gFaceHandle, 0, &meshSet, &err);
    assert(!err);

    iBase_TagHandle dim_tag;
    const char *tag2 = "GEOM_DIMENSION";
    iMesh_getTagHandle(mesh, tag2, &dim_tag, &err, strlen(tag2));
    assert(!err);

    int geom_dim;
    iMesh_getEntSetIntData(mesh, meshSet, dim_tag, &geom_dim, &err);
    assert(!err);

    if (geom_dim == 2)
    {
        SimpleArray<iBase_EntityHandle> mEdges;
        iMesh_getEntities(mesh, meshSet, iBase_EDGE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mEdges), &err);
        for (int i = 0; i < mEdges.size(); i++) projectEdgeHigherOrderNodes(mEdges[i]);
    }
}

///////////////////////////////////////////////////////////////////////////////

void GFace::projectEdgeHigherOrderNodes(iBase_EntityHandle mEdgeHandle)
{
    int err;
    SimpleArray<iBase_EntityHandle> edgeNodes;

    iMesh_getEntAdj(mesh, mEdgeHandle, iBase_VERTEX, ARRAY_INOUT(edgeNodes), &err);
    //
    // end nodes parametric coordinates can be ambiguous in the case of periodic
    // surfaces, therefore, estimate UV near to the ends first.
    // In this case, at least one vertex of the edge will be on the surface.
    // i.e. both the nodes cann't be on the boundary curves. If that were the
    // case, ambiguities can occur.
    //

    Point3D p0, pN, pnear, pon;
    Point2D uv, uv0, uvN, uvnear;

    iMesh_getVtxCoord(mesh, edgeNodes[0], &p0[0], &p0[1], &p0[2], &err);
    iMesh_getVtxCoord(mesh, edgeNodes[1], &pN[0], &pN[1], &pN[2], &err);

    double u, dl, arclen, arcdist, unear = 0.98;

    pnear[0] = (1 - unear) * p0[0] + unear * pN[0];
    pnear[1] = (1 - unear) * p0[1] + unear * pN[1];
    pnear[2] = (1 - unear) * p0[2] + unear * pN[2];
    uvnear = getUVCoords(pnear);
    uv0 = getUVCoords(p0, uvnear);

    pnear[0] = (1 - unear) * pN[0] + unear * p0[0];
    pnear[1] = (1 - unear) * pN[1] + unear * p0[1];
    pnear[2] = (1 - unear) * pN[2] + unear * p0[2];
    uvnear = getUVCoords(pnear);
    uvN = getUVCoords(pN, uvnear);

    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;
    iMesh_getData(mesh, mEdgeHandle, horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
    assert(!err);

    HO_Points *hopoints = (HO_Points *) tag_val;
    iBase_EntityHandle *nodesOnEdge = hopoints->nodeHandles;
    int numHPoints = hopoints->nx;
    int nhalf = numHPoints / 2;

    arclen = getGeodesicLength(uv0, uvN);

    // From the start->nhalf
    arcdist = 0.0;
    uvnear[0] = uv0[0];
    uvnear[1] = uv0[1];

    for (int k = 1; k < nhalf; k++)
    {
        iBase_EntityHandle currvertex = nodesOnEdge[k];
        dl = arclength_ratio[k] * arclen - arcdist;
        u = gllnodes[k];
        uv[0] = 0.5 * (1 - u) * uv0[0] + 0.5 * (1 + u) * uvN[0];
        uv[1] = 0.5 * (1 - u) * uv0[1] + 0.5 * (1 + u) * uvN[1];
//      getUVCoords(uvnear, dl, uv);
        pon = getXYZCoords(uv);
        iMesh_setVtxCoord(mesh, currvertex, pon[0], pon[1], pon[2], &err);
        uvnear[0] = uv[0];
        uvnear[1] = uv[1];
        arcdist += dl;
    }

    // From the end->nhalf
    arcdist = 0.0;
    uvnear[0] = uvN[0];
    uvnear[1] = uvN[1];
    for (int k = 0; k < nhalf; k++)
    {
        iBase_EntityHandle currvertex = nodesOnEdge[numHPoints - 1 - k];
        dl = arclength_ratio[k] * arclen - arcdist;
        u = gllnodes[k];
        uv[0] = 0.5 * (1 - u) * uvN[0] + 0.5 * (1 + u) * uv0[0]; // careful, u0, uN swapped
        uv[1] = 0.5 * (1 - u) * uvN[1] + 0.5 * (1 + u) * uv0[1]; // careful, u0, uN swapped
//      getUVCoords(uvnear, dl, uv);
        pon = getXYZCoords(uv);
        iMesh_setVtxCoord(mesh, currvertex, pon[0], pon[1], pon[2], &err);
        uvnear[0] = uv[0];
        uvnear[1] = uv[1];
        arcdist += dl;
    }

    if (numHPoints % 2)
    {
        int midpos = nhalf;
        iBase_EntityHandle currvertex = nodesOnEdge[midpos];
        u = gllnodes[midpos + 1];
        uv[0] = 0.5 * (1 - u) * uv0[0] + 0.5 * (1 + u) * uvN[0];
        uv[1] = 0.5 * (1 - u) * uv0[1] + 0.5 * (1 + u) * uvN[1];
        pon = getXYZCoords(uv);
        iMesh_setVtxCoord(mesh, currvertex, pon[0], pon[1], pon[2], &err);
    }
}

////////////////////////////////////////////////////////////////////////////////

void GFace::projectFaceHigherOrderNodes(const vector<double> &gnodes)
{
    int err;
    gllnodes = gnodes;

    if (mesh == 0)
    {
        cout << " Warning: No mesh is present in geometric face : Higher order nodes projection not done " << endl;
        return;
    }

    if (assoc == 0)
    {
        cout << " Warning: No assoc is present in geometric face : Higher order nodes projection not done " << endl;
        return;
    }

    if (rel == 0)
    {
        cout << " Warning: No relation is present in geometric face : Higher order nodes projection not done " << endl;
        return;
    }

    const char *tag1 = "HO_POINTS";
    iMesh_getTagHandle(mesh, tag1, &horder_tag, &err, strlen(tag1));
    assert(!err);

    int nsize = gllnodes.size();
    arclength_ratio.resize(nsize);
    for (int i = 0; i < nsize; i++)
        arclength_ratio[i] = 0.5 * (1.0 + gllnodes[i]);

    iBase_EntitySetHandle meshSet;
    iRel_getEntSetAssociation(assoc, rel, gFaceHandle, 0, &meshSet, &err);
    assert(!err);

    iBase_TagHandle dim_tag;
    const char *tag2 = "GEOM_DIMENSION";
    iMesh_getTagHandle(mesh, tag2, &dim_tag, &err, strlen(tag2));
    assert(!err);

    int geom_dim;
    iMesh_getEntSetIntData(mesh, meshSet, dim_tag, &geom_dim, &err);
    assert(!err);

    if (geom_dim == 2)
    {
        SimpleArray<iBase_EntityHandle> mFaces;
        iMesh_getEntities(mesh, meshSet, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mFaces), &err);
        for (int i = 0; i < mFaces.size(); i++) projectFaceHigherOrderNodes(mFaces[i]);
    }
}

////////////////////////////////////////////////////////////////////////////////

void GFace::projectFaceHigherOrderNodes(iBase_EntityHandle mFaceHandle)
{
    int err;

    int offset, nx, ny, numHPoints;

    SimpleArray<iBase_EntityHandle> faceNodes;
    iMesh_getEntAdj(mesh, mFaceHandle, iBase_VERTEX, ARRAY_INOUT(faceNodes), &err);

    char *tag_val = NULL;
    int tag_val_allocated, tag_val_size;
    iMesh_getData(mesh, mFaceHandle, horder_tag, &tag_val, &tag_val_allocated, &tag_val_size, &err);
    assert(!err);
    HO_Points *hopoints = (HO_Points *) tag_val;
    nx = hopoints->nx;
    ny = hopoints->ny;
    numHPoints = nx*ny;

    iBase_EntityHandle *nodeHandles = hopoints->nodeHandles;

    vector<double> u(numHPoints);
    vector<double> v(numHPoints);

    iBase_EntityHandle currvertex;

    Point3D p3d, p0, p1, corners[4], pnear;
    Point2D uv, uvnear;
    double  t = 0.95;
    iMesh_getVtxCoord(mesh, faceNodes[0], &corners[0][0], &corners[0][1], &corners[0][2], &err);
    iMesh_getVtxCoord(mesh, faceNodes[1], &corners[1][0], &corners[1][1], &corners[1][2], &err);
    iMesh_getVtxCoord(mesh, faceNodes[2], &corners[2][0], &corners[2][1], &corners[2][2], &err);
    iMesh_getVtxCoord(mesh, faceNodes[3], &corners[3][0], &corners[3][1], &corners[3][2], &err);

    p0     = linear_interpolation01( corners[0], corners[1], t );
    p1     = linear_interpolation01( corners[0], corners[3], t );
    pnear  = linear_interpolation01( p0, p1, 0.50);
    uvnear = getUVCoords(pnear); 

    offset = 0;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = getUVCoords(p3d, uvnear);
    u[offset] = uv[0];
    v[offset] = uv[1];

    p0     = linear_interpolation01( corners[1], corners[0], t );
    p1     = linear_interpolation01( corners[1], corners[2], t );
    pnear  = linear_interpolation01( p0, p1, 0.50);
    uvnear = getUVCoords(pnear); 

    offset = nx - 1;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = getUVCoords(p3d, uvnear);
    u[offset] = uv[0];
    v[offset] = uv[1];

    p0     = linear_interpolation01( corners[3], corners[2], t );
    p1     = linear_interpolation01( corners[3], corners[0], t );
    pnear  = linear_interpolation01( p0, p1, 0.50);
    uvnear = getUVCoords(pnear); 

    offset = (ny - 1) * nx;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = getUVCoords(p3d, uvnear);
    u[offset] = uv[0];
    v[offset] = uv[1];

    p0     = linear_interpolation01( corners[2], corners[1], t );
    p1     = linear_interpolation01( corners[2], corners[3], t );
    pnear  = linear_interpolation01( p0, p1, 0.50);
    uvnear = getUVCoords(pnear); 

    offset = nx * ny - 1;
    currvertex = nodeHandles[offset];
    iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
    uv = getUVCoords(p3d, uvnear);
    u[offset] = uv[0];
    v[offset] = uv[1];

    uvnear[0] = u[0];
    uvnear[1] = v[0];
    for (int i = 1; i < nx - 1; i++)
    {
        offset = i;
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = getUVCoords(p3d, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
        uvnear    = uv;
    }

    uvnear[0] = u[nx*(ny-1)];
    uvnear[1] = v[nx*(ny-1)];
    for (int i = 1; i < nx - 1; i++)
    {
        offset = i + (ny - 1) * nx;
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = getUVCoords(p3d, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
        uvnear    = uv;
    }

    uvnear[0] = u[0];
    uvnear[1] = v[0];
    for (int j = 1; j < ny - 1; j++)
    {
        offset = j*nx;
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = getUVCoords(p3d, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
    }


    uvnear[0] = u[(nx-1)];
    uvnear[1] = v[(nx-1)];
    for (int j = 1; j < ny - 1; j++)
    {
        offset = j * nx + (nx - 1);
        currvertex = nodeHandles[offset];
        iMesh_getVtxCoord(mesh, currvertex, &p3d[0], &p3d[1], &p3d[2], &err);
        uv = getUVCoords(p3d, uvnear);
        u[offset] = uv[0];
        v[offset] = uv[1];
    }

    TFIMap::blend_from_edges(u, gllnodes, gllnodes);
    TFIMap::blend_from_edges(v, gllnodes, gllnodes);

    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            offset = j * nx + i;
            uv[0] = u[offset];
            uv[1] = v[offset];
            Point3D pon = getXYZCoords(uv);
            iMesh_setVtxCoord(mesh, nodeHandles[offset], pon[0], pon[1], pon[2], &err);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////

void GFace::saveAs(const string &filename) const
{
    ofstream ofile(filename.c_str(), ios::out);
    if (ofile.fail()) return;

    int nx = 51;
    int ny = 51;

    double du = (umax - umin) / (double) (nx - 1);
    double dv = (vmax - vmin) / (double) (ny - 1);

    double x, y, z;
    ofile << "#Nodes " << nx * ny << endl;
    int err, index = 0;
    for (int j = 0; j < ny; j++)
    {
        double v = vmin + j*dv;
        for (int i = 0; i < nx; i++)
        {
            double u = umin + i*du;
            iGeom_getEntUVtoXYZ(geometry, gFaceHandle, u, v, &x, &y, &z, &err);
            ofile << index++ << " " << x << " " << y << " " << z << endl;
        }
    }
}

