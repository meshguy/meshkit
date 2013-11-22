#include "ITAP_NetGen_EdgeMesh.h"

void ITAP_NetGen_EdgeMesh :: discretize_close_edge(iBase_EntityHandle edgeHandle, int numEdges)
{
    int err, status;
    double x, y, z, u, v, du;

    int numNodes = numEdges;

    double umin, umax;
    iGeom_getEntURange( geometry, edgeHandle, &umin, &umax, &err);

    du = (umax - umin) / (double) numEdges;

    vector<int> nodes( numNodes );
    vector<double> param( numNodes + 1);

    Point<3> p3d;

    param[0] = umin;
    for (int i = 1; i < numNodes; i++)
    {
        u = umin + i*du;
        iGeom_getEntUtoXYZ(geometry, edgeHandle, u, &x, &y, &z, &err);
        p3d(0) = x;
        p3d(1) = y;
        p3d(2) = z;
        MeshPoint mp( p3d );
        nodes[i] = mesh.AddPoint(mp);
        param[i] = u;
        global_meshedge_id++;
    }
    param[numNodes] = umax;

    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj(geometry, edgeHandle, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

    int vertexID;
    iGeom_getIntData( geometry, gNodes[0], id_tag, &vertexID, &err);
    nodes[0] = vertexID;

    int edgeID;
    iGeom_getIntData( geometry, edgeHandle, id_tag, &edgeID, &err);

    SimpleArray<iBase_EntityHandle> gFaces;
    iGeom_getEntAdj(geometry, edgeHandle, iBase_FACE, ARRAY_INOUT(gFaces), &err);

    int faceID;
    iGeom_getIntData( geometry, gFaces[0], id_tag, &faceID, &err);

    vector<double> faceU( numNodes+1 );
    vector<double> faceV( numNodes+1 );

    SearchUV  searchUV;
    for( int i = 0; i < numNodes+1; i++) {
        iGeom_getEntUtoXYZ(geometry, edgeHandle, param[i], &x, &y, &z, &err);
        searchUV.getUV( geometry, gFaces[0], x, y, z, u, v);
        faceU[i] = u;
        faceV[i] = v;
    }

    Segment seg;
    for (int i = 0; i < numEdges; i++)
    {
        seg.p1 = nodes[i];
        seg.p2 = nodes[(i+1)%numNodes];

        seg.edgenr = global_meshedge_id++;
        seg.si = faceID;

        seg.epgeominfo[0].dist = param[i];
        seg.epgeominfo[0].edgenr = edgeID;
        seg.epgeominfo[0].u = faceU[i];
        seg.epgeominfo[0].v = faceV[i];

        seg.epgeominfo[1].dist = param[i+1];
        seg.epgeominfo[1].edgenr = edgeID;
        seg.epgeominfo[1].u = faceU[i+1];
        seg.epgeominfo[1].v = faceV[i+1];

        mesh.AddSegment(seg);
     }
}

////////////////////////////////////////////////////////////////////////////////

void ITAP_NetGen_EdgeMesh :: discretize_open_edge( iBase_EntityHandle edgeHandle, int numEdges )
{
    int err, status, numNodes;
    double x, y, z, u, v, du;

    double umin, umax;
    iGeom_getEntURange( geometry, edgeHandle, &umin, &umax, &err);
    du = (umax - umin) / (double) numEdges;

    numNodes = numEdges + 1;
    
    Point<3> p3d;

    vector<int>  nodes( numNodes );
    vector<double> param( numNodes );

    param[0] = umin;
    for (int i = 1; i < numNodes - 1; i++)
    {
        u = umin + i*du;
        iGeom_getEntUtoXYZ(geometry, edgeHandle, u, &x, &y, &z, &err);
        p3d(0) = x;
        p3d(1) = y;
        p3d(2) = z;
        MeshPoint mp( p3d );
        nodes[i] = mesh.AddPoint(mp);
        param[i] = u;
        global_meshedge_id++;
    }
    param[numNodes-1] = umax;

    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj(geometry, edgeHandle, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

    int vertexID;
    iGeom_getIntData( geometry, gNodes[0], id_tag, &vertexID, &err);
    nodes[0] = vertexID;

    iGeom_getIntData( geometry, gNodes[1], id_tag, &vertexID, &err);
    nodes[numNodes-1] = vertexID;

    int edgeID;
    iGeom_getIntData( geometry, edgeHandle, id_tag, &edgeID, &err);

    SimpleArray<iBase_EntityHandle> gFaces;
    iGeom_getEntAdj(geometry, edgeHandle, iBase_FACE, ARRAY_INOUT(gFaces), &err);

    int faceID;
    iGeom_getIntData( geometry, gFaces[0], id_tag, &faceID, &err);

    vector<double> faceU( numNodes+1 );
    vector<double> faceV( numNodes+1 );

    SearchUV  searchUV;
    for( int i = 0; i < numNodes+1; i++) {
        iGeom_getEntUtoXYZ(geometry, edgeHandle, param[i], &x, &y, &z, &err);
        searchUV.getUV( geometry, gFaces[0], x, y, z, u, v);
        faceU[i] = u;
        faceV[i] = v;
    }

    Segment seg;
    for (int i = 0; i < numEdges; i++)
    {
        seg.p1 = nodes[i];
        seg.p2 = nodes[(i+1)%numNodes];

        seg.edgenr = global_meshedge_id++;
        seg.si = faceID;

        seg.epgeominfo[0].dist = param[i];
        seg.epgeominfo[0].edgenr = edgeID;
        seg.epgeominfo[0].u = faceU[i];
        seg.epgeominfo[0].v = faceV[i];

        seg.epgeominfo[1].dist = param[i+1];
        seg.epgeominfo[1].edgenr = edgeID;
        seg.epgeominfo[1].u = faceU[i+1];
        seg.epgeominfo[1].v = faceV[i+1];

        mesh.AddSegment(seg);
     }

}

////////////////////////////////////////////////////////////////////////////////
void ITAP_NetGen_EdgeMesh :: discretize(iBase_EntityHandle edgeHandle, int numEdges)
{
    int err;

    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj(geometry, edgeHandle, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

    switch (gNodes.size())
    {
    case 1:
        discretize_close_edge(edgeHandle, numEdges);
        break;
    case 2:
        discretize_open_edge(edgeHandle, numEdges);
        break;
    default:
        cout << "Fatal Error: Invalid Geometric edge " << endl;
        exit(0);
    }
}

/////////////////////////////////////////////////////////////////////////////

/*
void DivideEdge( iGeom_Instance &geometry, iBase_EntityHandle &edgeHandle, 
                 vector<MeshPoint> & ps, vector<double> &params, Mesh & mesh)
{
    int err;
    double s0, s1;
    double maxh = mparam.maxh;
    int nsubedges = 1;
    double pnt[3], oldpnt[3];
    double svalue[DIVIDE_EDGE_SECTIONS];

    iGeom_getEntURange( geometry, edgeHandle, &s0, &s1, &err);

    double hvalue[DIVIDE_EDGE_SECTIONS + 1];
    hvalue[0] = 0;

    pnt = c->Value(s0);

    double olddist = 0;
    double dist = 0;

    for (int i = 1; i <= DIVIDE_EDGE_SECTIONS; i++)
    {
        oldpnt = pnt;
        pnt = c->Value(s0 + (i / double(DIVIDE_EDGE_SECTIONS))*(s1 - s0));
        hvalue[i] = hvalue[i - 1] +
                1.0 / mesh.GetH(Point3d(pnt.X(), pnt.Y(), pnt.Z())) *
                pnt.Distance(oldpnt);

        olddist = dist;
        dist = pnt.Distance(oldpnt);
    }

    nsubedges = max(1, int(floor(hvalue[DIVIDE_EDGE_SECTIONS] + 0.5)));

    ps.SetSize(nsubedges - 1);
    params.SetSize(nsubedges + 1);

    int i = 1;
    int i1 = 0;
    do
    {
        if (hvalue[i1] / hvalue[DIVIDE_EDGE_SECTIONS] * nsubedges >= i)
        {
            params[i] = s0 + (i1 / double(DIVIDE_EDGE_SECTIONS))*(s1 - s0);
            pnt = c->Value(params[i]);
            ps[i - 1] = MeshPoint(Point3d(pnt.X(), pnt.Y(), pnt.Z()));
            i++;
        }

        i1++;
        if (i1 > DIVIDE_EDGE_SECTIONS)
        {
            nsubedges = i;
            ps.SetSize(nsubedges - 1);
            params.SetSize(nsubedges + 1);
            cout << "divide edge: local h too small" << endl;
        }

    }
    while (i < nsubedges);

    params[0] = s0;
    params[nsubedges] = s1;

    if (params[nsubedges] <= params[nsubedges - 1])
    {
        cout << "CORRECTED" << endl;
        ps.SetSize(nsubedges - 2);
        params.SetSize(nsubedges);
        params[nsubedges] = s1;
    }
}
*/

/////////////////////////////////////////////////////////////////////////////

void ITAP_NetGen_EdgeMesh :: execute()
{
    int err;

    iBase_EntitySetHandle rootSet;
    iGeom_getRootSet(geometry, &rootSet, &err);

    const char *tag = "GLOBAL_ID";
    int namelen = strlen(tag);
    iGeom_getTagHandle(geometry, tag, &id_tag, &err, namelen);

    int numNodes, numEdges, numFaces, numCells;

    // First use Geometric Nodes.
    SimpleArray<iBase_EntityHandle> geomNodes;
    iGeom_getEntities( geometry, rootSet, iBase_VERTEX, ARRAY_INOUT( geomNodes), &err);

    int index, vertexID;
    Point<3> p3d;
    double x, y, z;
    for( int i = 0; i < geomNodes.size(); i++) 
    {
       iGeom_getVtxCoord( geometry, geomNodes[i], &x, &y, &z, &err);
       iGeom_getIntData( geometry, geomNodes[i], id_tag, &vertexID, &err);
       p3d(0) = x;
       p3d(1) = y;
       p3d(2) = z;
       MeshPoint mp( p3d );
       index = mesh.AddPoint(mp);
       assert( index == vertexID );
    }
    
    SimpleArray<iBase_EntityHandle> geomFaces;
    iGeom_getEntities( geometry, rootSet, iBase_FACE, ARRAY_INOUT(geomFaces), &err);

    numFaces = geomFaces.size();

    vector<int> face2solid[2];
    face2solid[0].resize(numFaces);
    face2solid[1].resize(numFaces);
    for (int i = 0; i < numFaces; i++)
    {
        face2solid[0][i] = 0;
        face2solid[1][i] = 0;
    } 

    SimpleArray<iBase_EntityHandle> faceCells;
    iGeom_getEntities( geometry, rootSet, iBase_REGION, ARRAY_INOUT( faceCells ), &err);

    int edgeID, faceID, cellID;
    for( int i = 0; i < numFaces; i++) {
         iGeom_getEntAdj(geometry, geomFaces[i], iBase_REGION, ARRAY_INOUT(faceCells), &err);
         iGeom_getIntData( geometry, geomFaces[i], id_tag, &faceID, &err);
         for( int j  = 0; j < faceCells.size(); j++) {
              iGeom_getIntData( geometry, faceCells[j], id_tag, &cellID, &err);
              if (face2solid[0][faceID - 1] == 0)
                   face2solid[0][faceID - 1] = cellID - 1;
              else
                   face2solid[1][faceID - 1] = cellID - 1;
         }
    }

    for( int i = 0; i < numFaces; i++) 
    {
         int cell0 = face2solid[0][i];
         int cell1 = face2solid[1][i];
         mesh.AddFaceDescriptor(FaceDescriptor(i, cell0, cell1, 0));

    }

    SimpleArray<iBase_EntityHandle> geomEdges;
    iGeom_getEntities( geometry, rootSet, iBase_EDGE, ARRAY_INOUT( geomEdges), &err);

    numEdges = geomEdges.size();

    vector<double> elen(numEdges);
    SimpleArray<double> measure;

    for (int i = 0; i < numEdges; i++) 
    {
        iGeom_measure(geometry, &geomEdges[i], 1, ARRAY_INOUT(measure), &err);
        elen[i] = measure[0];
    }

    vector<double> tmpelen;
    tmpelen = elen;
    sort( tmpelen.begin(), tmpelen.end() );
    double mean_edge_length = tmpelen[numEdges/2];

    cout << " Minimum Edge Length : " << *min_element( elen.begin(), elen.end() ) << endl;
    cout << " Maximum Edge Length : " << *max_element( elen.begin(), elen.end() ) << endl;
    cout << " Mean Length         : " <<  mean_edge_length << endl;

    double spacing = mean_edge_length /(double) 10.0;

    int numSegments;

    for( int i  = 0; i < numEdges; i++) 
    {
         if( elen[i] < 0.5*mean_edge_length)
             numSegments = 1;
        else
             numSegments  = elen[i]/ spacing;

        discretize(geomEdges[i], numSegments);
    }

   mesh.CalcSurfacesOfNode();
}

