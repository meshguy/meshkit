#include "SurfMesher.h"
#include "SearchUV.h"
#include <assert.h>
#include <fstream>
#include <limits>
#include <math.h>
#include <sstream>
#include <algorithm>

#include <map>
using namespace std;

#include <iostream>
using namespace std;

double get_edge_length(const Vertex *v0, const Vertex *v1)
{
    double dx = v0->xyzCoords[0] - v1->xyzCoords[0];
    double dy = v0->xyzCoords[1] - v1->xyzCoords[1];
    double dz = v0->xyzCoords[2] - v1->xyzCoords[2];

    return sqrt(dx * dx + dy * dy + dz * dz);
}

int laplacian_uv_smoothing(iGeom_Instance &geom, iBase_EntityHandle &gFace,
                           SurfaceMesh &mesh)
{
    int err;
    int maxIterations = 1000;
    double t = 0.25;
    std::map<Vertex*, vector<Vertex*> > vmap;

    vector<Vertex*> neighs;
    for (int i = 0; i < mesh.faces.size(); i++)
    {
        Face face = mesh.faces[i];
        for (int j = 0; j < 3; j++)
        {
            Vertex *v1 = face.getVertex((j + 1) % 3);
            Vertex *v2 = face.getVertex((j + 2) % 3);

            if (!v1->onBoundary)
            {
                neighs = vmap[v1];
                if (find(neighs.begin(), neighs.end(), v2) == neighs.end())
                    vmap[v1].push_back(v2);
            }

            if (!v2->onBoundary)
            {
                neighs = vmap[v2];
                if (find(neighs.begin(), neighs.end(), v1) == neighs.end())
                    vmap[v2].push_back(v1);
            }
        }
    }

    vector<Vertex*> movableNodes;
    std::map<Vertex*, vector<Vertex*> > ::const_iterator it;
    for (it = vmap.begin(); it != vmap.end(); ++it)
    {
        movableNodes.push_back(it->first);
    }

    int numNeighs, numNodes = movableNodes.size();

    vector<double> newCoords;
    newCoords.resize(2 * numNodes);

    double usum, vsum, wght, sum_wght;
    double dx, dy, dz, du, dv, dist;
    double x, y, z, u, v;

    int iter = 0;
    while (1)
    {
        for (int i = 0; i < numNodes; i++)
        {
            Vertex *apex = movableNodes[i];
            vector<Vertex*> &neighs = vmap[apex];
            numNeighs = neighs.size();
            assert(numNeighs);

            usum = 0.0;
            vsum = 0.0;
            sum_wght = 0.0;
            for (int j = 0; j < numNeighs; j++)
            {
                dx = apex->xyzCoords[0] - neighs[j]->xyzCoords[0];
                dy = apex->xyzCoords[1] - neighs[j]->xyzCoords[1];
                dz = apex->xyzCoords[2] - neighs[j]->xyzCoords[2];
                wght = sqrt(dx * dx + dy * dy + dz * dz);
                usum += wght * neighs[j]->uvCoords[0];
                vsum += wght * neighs[j]->uvCoords[1];
                sum_wght += wght;
            }
            assert(sum_wght > 0.0);

            usum /= sum_wght;
            vsum /= sum_wght;

            newCoords[2 * i + 0] = usum;
            newCoords[2 * i + 1] = vsum;
        }

        double maxerror = 0.0;
        for (int i = 0; i < numNodes; i++)
        {
            Vertex *apex = movableNodes[i];
            du = newCoords[2 * i + 0] - apex->uvCoords[0];
            dv = newCoords[2 * i + 1] - apex->uvCoords[1];
            dist = sqrt(du * du + dv * dv);
            maxerror = max(dist, maxerror);
        }
        cout << " Max error : " << maxerror << endl;

        if (maxerror < 1.0E-10) return 0;

        for (int i = 0; i < numNodes; i++)
        {
            Vertex *apex = movableNodes[i];
            apex->uvCoords[0] = (1.0 - t) * apex->uvCoords[0] + t * newCoords[2 * i + 0];
            apex->uvCoords[1] = (1.0 - t) * apex->uvCoords[1] + t * newCoords[2 * i + 1];
            u = apex->uvCoords[0];
            v = apex->uvCoords[1];
            iGeom_getEntUVtoXYZ(geom, gFace, u, v, &x, &y, &z, &err);
            apex->xyzCoords[0] = x;
            apex->xyzCoords[1] = y;
            apex->xyzCoords[2] = z;
        }

        if (iter++ == maxIterations) return 1;
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////

int SurfaceMesher :: discretize(GeomMesh &geomesh, iBase_EntityHandle gFace)
{
    int status, err;
    iMesh_Instance mesh = geomesh.mesh;
    iGeom_Instance geom = geomesh.geom;
    iRel_Instance assoc = geomesh.assoc;
    iRel_RelationHandle relation = geomesh.relation02;


    const char *tag1 = "GLOBAL_ID";
    int namelen = strlen(tag1);
    iGeom_getTagHandle(geom, tag1, &geom_id_tag, &err, namelen);

    double umin, umax, vmin, vmax;
    iGeom_getEntUVRange(geom, gFace, &umin, &vmin, &umax, &vmax, &err);

    SimpleArray<iBase_EntityHandle> gEdges, mEdges, mNodes, edgeNodes;
    iGeom_getEntAdj(geom, gFace, iBase_EDGE, ARRAY_INOUT(gEdges), &err);

    int numGeomEdges = gEdges.size();

    vector<Edge> boundSegments;
    vector<iBase_EntityHandle> nodeHandles;
    vector<iBase_EntitySetHandle> edgeMeshSets;
    std::map<iBase_EntityHandle, Vertex*> uvVertex;

    edgeMeshSets.resize(numGeomEdges);

    int id = 0;
    int edgeID;
    double x, y, z, u, uc, v;
    double xg, yg, zg;
    Edge newedge;

    SearchUV searchUV;
    for (int i = 0; i < numGeomEdges; i++)
    {
        iRel_getEntSetAssociation(assoc, relation, gEdges[i], 0, &edgeMeshSets[i], 
                                  &err);
        iGeom_getIntData(geom, gEdges[i], geom_id_tag, &edgeID, &err);

        mEdges.clear();
        iMesh_getEntities( mesh, edgeMeshSets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES, 
                           ARRAY_INOUT(mEdges), &err);
        for (int j = 0; j < mEdges.size(); j++)
        {
            edgeNodes.clear();
            iMesh_getEntAdj( mesh, mEdges[j], iBase_VERTEX, ARRAY_INOUT(edgeNodes), 
                             &err);
            for (int k = 0; k < edgeNodes.size(); k++)
            {
                if (uvVertex.find(edgeNodes[k]) == uvVertex.end())
                {
                    iMesh_getVtxCoord(mesh, edgeNodes[k], &x, &y, &z, &err);
                    searchUV.getUV(geom, gFace, x, y, z, u, v);
                    Vertex *vertex = new Vertex;
                    vertex->id = id;
                    vertex->onBoundary   = 1;
                    vertex->uvCoords[0]  = u;
                    vertex->uvCoords[1]  = v;
                    vertex->xyzCoords[0] = x;
                    vertex->xyzCoords[1] = y;
                    vertex->xyzCoords[2] = z;
                    uvVertex[edgeNodes[k]] = vertex;
                    nodeHandles.push_back(edgeNodes[k]);
                    id++;
                }
            }
            newedge.connect[0] = uvVertex[ edgeNodes[0] ];
            newedge.connect[1] = uvVertex[ edgeNodes[1] ];
            newedge.geomEdgeID = edgeID;
            boundSegments.push_back(newedge);
        }
    }

    int numBoundNodes = uvVertex.size();

    SurfaceMesh surfmesh = getMesh(boundSegments);
    //exit(0);

    iBase_EntityHandle newHandle;

    vector<Vertex*> interior = getInsertedNodes();

    for (int i = 0; i < interior.size(); i++)
    {
        Vertex *vertex = interior[i];
        double u = vertex->uvCoords[0];
        double v = vertex->uvCoords[1];
        iGeom_getEntUVtoXYZ(geom, gFace, u, v, &x, &y, &z, &err);
        assert(!err);
        vertex->xyzCoords[0] = x;
        vertex->xyzCoords[1] = y;
        vertex->xyzCoords[2] = z;
        iMesh_createVtx(mesh, x, y, z, &newHandle, &err);
        assert(!err);
        nodeHandles.push_back(newHandle);
    }

    vector<Face> surfElements = surfmesh.faces;
    int numFaces = surfElements.size();

    assert(numFaces);

    vector<iBase_EntityHandle> fConnect, faceHandles;
    faceHandles.resize(numFaces);

    for (int i = 0; i < numFaces; i++)
    {
        Face face = surfElements[i];
        int numNodes = face.getNumNodes();
        fConnect.resize(numNodes);
        for (int j = 0; j < numNodes; j++)
        {
            Vertex *vertex = face.getVertex(j);
            fConnect[j] = nodeHandles[vertex->id];
        }
        switch (numNodes)
        {
        case 3:
            iMesh_createEnt( mesh, iMesh_TRIANGLE, &fConnect[0], numNodes, 
                             &newHandle, &status, &err);
            break;
        case 4:
            iMesh_createEnt(mesh, iMesh_QUADRILATERAL, &fConnect[0], numNodes, 
                            &newHandle, &status, &err);
            break;
        default:
            iMesh_createEnt(mesh, iMesh_POLYGON, &fConnect[0], numNodes, &newHandle,
                            &status, &err);
            break;
        }
        faceHandles[i] = newHandle;
    }

    iBase_EntitySetHandle entitySet;
    iMesh_createEntSet(mesh, 1, &entitySet, &err);
    iMesh_addEntArrToSet(mesh, &faceHandles[0], numFaces, entitySet, &err);

    int faceID;
    iGeom_getIntData(geom, gFace, geom_id_tag, &faceID, &err);

    iBase_TagHandle mesh_id_tag;
    iMesh_getTagHandle(mesh, tag1, &mesh_id_tag, &err, namelen);
    iMesh_setEntSetIntData(mesh, entitySet, mesh_id_tag, faceID, &err);

    iBase_TagHandle dim_tag;
    const char *tag2 = "GEOM_DIMENSION";
    namelen = strlen(tag2);
    iMesh_getTagHandle(mesh, tag2, &dim_tag, &err, namelen);

    int dim = 2;
    iMesh_setEntSetIntData(mesh, entitySet, dim_tag, dim, &err);

    iRel_setEntSetAssociation(assoc, relation, gFace, entitySet, &err);

    for (int i = 0; i < gEdges.size(); i++)
        iMesh_addPrntChld(mesh, entitySet, edgeMeshSets[i], &err);
}

////////////////////////////////////////////////////////////////////////////////

int SurfaceMesher :: saveAs(const vector<Edge> &segments, const string &filename)
{
    std::map<int, Vertex*> vmap;

    for (int i = 0; i < segments.size(); i++)
    {
        Vertex *v1 = segments[i].connect[0];
        Vertex *v2 = segments[i].connect[1];
        vmap[v1->id] = v1;
        vmap[v2->id] = v2;
    }
    ofstream ofile(filename.c_str(), ios::out);

    ofile << "OFF " << endl;
    ofile << vmap.size() << " 0 " << segments.size() << endl;

    for (int i = 0; i < vmap.size(); i++)
    {
        Vertex *v1 = vmap[i];
        ofile << v1->uvCoords[0] << " "
                << v1->uvCoords[1] << " 0 " << endl;
    }

    for (int i = 0; i < segments.size(); i++)
    {
        Vertex *v1 = segments[i].connect[0];
        Vertex *v2 = segments[i].connect[1];
        ofile << v1->id << " " << v2->id << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////

int SurfaceMesher :: saveAs(SurfaceMesh &uvmesh, const string &filename)
{
    ofstream ofile(filename.c_str(), ios::out);

    ofile << "OFF " << endl;
    ofile << uvmesh.nodes.size() << " " << uvmesh.faces.size() << " 0 " << endl;

    for (int i = 0; i < uvmesh.nodes.size(); i++)
    {
        Vertex *vertex = uvmesh.nodes[i];
        ofile << vertex->xyzCoords[0] << " "
                << vertex->xyzCoords[1] << " "
                << vertex->xyzCoords[2] << endl;
    }
    /*
        for (int i = 0; i < uvmesh.nodes.size(); i++)
        {
            Vertex *vertex = uvmesh.nodes[i];
            ofile << vertex->uvCoords[0] << " "
                  << vertex->uvCoords[1] << "  0 " << endl;
        }
     */

    for (int i = 0; i < uvmesh.faces.size(); i++)
    {
        Face face = uvmesh.faces[i];
        ofile << " 3 " << face.connect[0]->id << " "
                << face.connect[1]->id << " "
                << face.connect[2]->id << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////
#ifdef TRIANGLE
SurfaceMesher* TriangleSurfaceMesher::newObject()
{
    SurfaceMesher *mesher = new TriangleSurfaceMesher;
    return mesher;
}
#endif

////////////////////////////////////////////////////////////////////////////////
#ifdef NETGEN
SurfaceMesher* NetGenSurfaceMesher::newObject()
{
    SurfaceMesher *mesher = new NetGenSurfaceMesher;
    return mesher;
}
#endif

////////////////////////////////////////////////////////////////////////////////

SurfaceMesher* TetGenSurfaceMesher::newObject()
{
    SurfaceMesher *mesher = new TetGenSurfaceMesher;
    return mesher;
}

////////////////////////////////////////////////////////////////////////////////

SurfaceMesher* SurfaceMesherFactory::getProduct(const string &product)
{
#ifdef NETGEN
    if (product == "NetGen") return NetGenSurfaceMesher::newObject();
#endif

    if (product == "TetGen") return TetGenSurfaceMesher::newObject();

#ifdef TRIANGLE
    if (product == "Triangle") return TriangleSurfaceMesher::newObject();
#endif

    return NULL;
}

////////////////////////////////////////////////////////////////////////////////
#ifdef NETGEN
SurfaceMesh NetGenSurfaceMesher::getMesh(const vector<Edge> &edges) const
{
    int  err;
    Mesh  mesh;

    int faceID;
    iGeom_getIntData(geometry, faceHandle, geom_id_tag, &faceID, &err);
    assert( !err );

    SimpleArray<iBase_EntityHandle> faceCells;
    iGeom_getEntAdj(geometry, faceHandle, iBase_REGION, ARRAY_INOUT(faceCells), &err);
    assert( !err );
  
    assert( faceCells.size() > 0 && faceCells.size() <= 2);

    int region0, region1;
    if( faceCells.size() == 1) {
        iGeom_getIntData(geometry, faceCells[0], geom_id_tag, &region0, &err);
        region1 = -1;
    }

    if( faceCells.size() == 2) {
        iGeom_getIntData(geometry, faceCells[0], geom_id_tag, &region0, &err);
        iGeom_getIntData(geometry, faceCells[1], geom_id_tag, &region1, &err);
    }
 
    mesh.AddFaceDescriptor( FaceDescriptor( faceID, region0, region1, 0));

    SurfaceMesh surfmesh;

    insertedNodes.clear();

    size_t numEdges = edges.size();
    if (numEdges < 3) return surfmesh;

    std::map<int, Vertex*> vmap;

    for (int i = 0; i < numEdges; i++)
    {
        Vertex *v1 = edges[i].connect[0];
        Vertex *v2 = edges[i].connect[1];

        if (vmap.find(v1->id) == vmap.end()) vmap[v1->id] = v1;
        if (vmap.find(v2->id) == vmap.end()) vmap[v2->id] = v2;
    }

    size_t numNodes = vmap.size();

    for (int i = 0; i < numNodes; i++)
    {
        if (vmap.find(i) == vmap.end())
        {
            cout << "Error: Node IDs are not continuous: " << endl;
            return surfmesh;
        }
    }

    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    iGeom_getEntBoundBox( geometry, faceHandle, &xmin, &ymin, &zmin, 
                          &xmax, &ymax, &zmax, &err);

    Point<3> minCorner, maxCorner;
    minCorner(0) = xmin;
    minCorner(1) = ymin;
    minCorner(2) = zmin;

    maxCorner(0) = xmax;
    maxCorner(1) = ymax;
    maxCorner(2) = zmax;

    double dx = fabs(xmax-xmin);
    double dy = fabs(ymax-ymin);
    double dz = fabs(zmax-zmin);
    double maxlen = std::max( dx, std::max(dy,dz) );

    Box<3> bbox(minCorner,maxCorner);

//  int projecttype = PLANESPACE;    
    int projecttype = PARAMETERSPACE;

    Meshing2 *ngmesher;
    ngmesher = new ITAP_NetGen_SurfaceMesh(geometry, faceHandle, bbox, projecttype);
    assert(ngmesher);

    vector<Vertex*> boundNodes( numNodes );

    Point<3> apoint;
    for (int i = 0; i < numNodes; i++)
    {
        Vertex *vertex = vmap[i];
        boundNodes[i] = vertex;
        apoint(0) = vertex->xyzCoords[0];
        apoint(1) = vertex->xyzCoords[1];
        apoint(2) = vertex->xyzCoords[2];
        MeshPoint mp(apoint);
        mesh.AddPoint(mp);
        ngmesher->AddPoint(mp, i+1);
    }

    Segment seg;
    for (int i = 0; i < numEdges; i++)
    {
        Vertex *v1 = edges[i].connect[0];
        Vertex *v2 = edges[i].connect[1];
        seg[0] = v1->id+1;
        seg[1] = v2->id+1;
        seg.si = faceID;

        seg.epgeominfo[0].edgenr = edges[i].geomEdgeID;
        seg.epgeominfo[0].u = v1->uvCoords[0];
        seg.epgeominfo[0].v = v1->uvCoords[1];

        seg.epgeominfo[1].edgenr = edges[i].geomEdgeID;
        seg.epgeominfo[1].u = v2->uvCoords[0];
        seg.epgeominfo[1].v = v2->uvCoords[1];

        mesh.AddSegment(seg);
    }

    PointGeomInfo gi0, gi1;
    for (int i = 0; i < numEdges; i++)
    {
        Vertex *v1 = edges[i].connect[0];
        Vertex *v2 = edges[i].connect[1];

        gi0.trignum = faceID;
        gi0.u = v1->uvCoords[0];
        gi0.v = v1->uvCoords[1];

        gi1.trignum = faceID;
        gi1.u = v2->uvCoords[0];
        gi1.v = v2->uvCoords[1];

        int p1 = v1->id+1;
        int p2 = v2->id+1;
        ngmesher->AddBoundaryElement( p1, p2, gi0, gi1 );
    }

    SimpleArray<double> faceArea;
    iGeom_measure( geometry, &faceHandle, 1, ARRAY_INOUT(faceArea), &err);

    ngmesher->SetMaxArea( faceArea[0] );
   
//  double maxh = mparam.maxh;
    double maxh = maxlen;

    MESHING2_RESULT res = ngmesher->GenerateMesh(mesh, maxh, faceID);

    /*
    vector<Face> surfTriangles;
        int nstart = in.numberofpoints;
        numNodes = out.numberofpoints;

        double uv[2];
        index = nstart;
        for (int i = nstart; i < numNodes; i++)
        {
            uv[0] = out.pointlist[2 * i + 0];
            uv[1] = out.pointlist[2 * i + 1];
            Vertex *vertex = new Vertex;
            vertex->uvCoords[0] = uv[0];
            vertex->uvCoords[1] = uv[1];
            vertex->id = i;
            vmap[i] = vertex;
            insertedNodes.push_back(vertex);
        }

        int numFaces = out.numberoftriangles;
        assert(numFaces);
        surfTriangles.resize(numFaces);

        vector<Vertex*> connect(3);
        for (int i = 0; i < numFaces; i++)
        {
            int n1 = out.trianglelist[3 * i + 0];
            connect[0] = vmap[n1];
            int n2 = out.trianglelist[3 * i + 1];
            connect[1] = vmap[n2];
            int n3 = out.trianglelist[3 * i + 2];
            connect[2] = vmap[n3];
            surfTriangles[i].connect = connect;
        }

        if (in.holelist) free(in.holelist);
        if (out.pointlist) free(out.pointlist);
        if (out.trianglelist) free(out.trianglelist);

        vector<Vertex*> nodes;
        if (numNodes) nodes.resize(numNodes);
        for (int i = 0; i < numNodes; i++) nodes[i] = vmap[i];

        surfmesh.faces = surfTriangles;
        surfmesh.nodes = nodes;
     */

    return surfmesh;
}
#endif

////////////////////////////////////////////////////////////////////////////////
#ifdef TRIANGLE
SurfaceMesh TriangleSurfaceMesher::getMesh(const vector<Edge> &edges) const
{
    SurfaceMesh surfmesh;

    insertedNodes.clear();

    vector<Face> surfTriangles;
    struct triangulateio in, out; // triangle sw data structure


    size_t numEdges = edges.size();
    if (numEdges < 3) return surfmesh;

    std::map<int, Vertex*> vmap;
    for (int i = 0; i < numEdges; i++)
    {
        Vertex *v1 = edges[i].connect[0];
        Vertex *v2 = edges[i].connect[1];
        if (vmap.find(v1->id) == vmap.end()) vmap[v1->id] = v1;
        if (vmap.find(v2->id) == vmap.end()) vmap[v2->id] = v2;
    }

    size_t numNodes = vmap.size();

    for (int i = 0; i < numNodes; i++)
    {
        if (vmap.find(i) == vmap.end())
        {
            cout << "Error: Node IDs are not continuous: " << endl;
            return surfmesh;
        }
    }

    vector<double> uvCoords( 2*numNodes);
    vector<Vertex*> boundNodes( numNodes);

    size_t index = 0;
    for (int i = 0; i < numNodes; i++)
    {
        Vertex *vertex = vmap[i];
        boundNodes[i] = vertex;
        uvCoords[index++] = vertex->uvCoords[0];
        uvCoords[index++] = vertex->uvCoords[1];
    }

    vector<int> segments( 2*numEdges);

    index = 0;
    double sum = 0.0;
    for (int i = 0; i < numEdges; i++)
    {
        Vertex *v1 = edges[i].connect[0];
        Vertex *v2 = edges[i].connect[1];
        segments[index++] = v1->id;
        segments[index++] = v2->id;
        double du = v1->uvCoords[0] - v2->uvCoords[0];
        double dv = v1->uvCoords[1] - v2->uvCoords[1];
        sum += sqrt(du * du + dv * dv);
    }
    double avg_length = sum / (double) numEdges;
    double area = 2.0 * avg_length*avg_length;

    // The following stuff is for "Triangle" software ...
    in.numberofpoints = uvCoords.size() / 2;
    in.pointlist = &uvCoords[0];
    in.numberofsegments = segments.size() / 2;
    in.segmentlist = &segments[0];
    in.numberofholes = 0;
    in.holelist = NULL;

    in.numberofpointattributes = 0;
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    in.segmentmarkerlist = NULL;
    in.triangleattributelist = NULL;
    in.numberofregions = 0;
    in.regionlist = NULL;

    // These needs to declared, I wasted full one day to find out this 
    // restriction

    out.pointlist = NULL;
    out.pointmarkerlist = NULL;
    out.trianglelist = NULL;
    out.neighborlist = NULL;
    out.segmentlist = NULL;
    out.segmentmarkerlist = NULL;
    out.edgelist = NULL;
    out.edgemarkerlist = NULL;

    cout << "Triangulating UV Space " << endl;

    std::ostringstream cmd;
    //  cmd << "q10.0BCPpzYQ" << "a" << area;
    cmd << "BCPpzYQ";

    char triSwitches[100];
    strcpy(triSwitches, cmd.str().c_str());

    triangulate(triSwitches, &in, &out, (struct triangulateio *) NULL);

    int nstart = in.numberofpoints;
    numNodes = out.numberofpoints;

    double uv[2];
    index = nstart;
    for (int i = nstart; i < numNodes; i++)
    {
        uv[0] = out.pointlist[2 * i + 0];
        uv[1] = out.pointlist[2 * i + 1];
        Vertex *vertex = new Vertex;
        vertex->uvCoords[0] = uv[0];
        vertex->uvCoords[1] = uv[1];
        vertex->id = i;
        vmap[i] = vertex;
        insertedNodes.push_back(vertex);
    }

    int numFaces = out.numberoftriangles;
    assert(numFaces);
    surfTriangles.resize(numFaces);

    vector<Vertex*> connect(3);
    for (int i = 0; i < numFaces; i++)
    {
        int n1 = out.trianglelist[3 * i + 0];
        connect[0] = vmap[n1];
        int n2 = out.trianglelist[3 * i + 1];
        connect[1] = vmap[n2];
        int n3 = out.trianglelist[3 * i + 2];
        connect[2] = vmap[n3];
        surfTriangles[i].connect = connect;
    }

    if (in.holelist) free(in.holelist);
    if (out.pointlist) free(out.pointlist);
    if (out.trianglelist) free(out.trianglelist);

    vector<Vertex*> nodes;
    if (numNodes) nodes.resize(numNodes);
    for (int i = 0; i < numNodes; i++) nodes[i] = vmap[i];

    surfmesh.faces = surfTriangles;
    surfmesh.nodes = nodes;

    return surfmesh;
}
#endif
///////////////////////////////////////////////////////////////////////////////

FlipEdge EdgeFlipping::build_edge(Vertex *v0, Vertex *v1)
{
    assert(v0 != NULL && v1 != NULL);
    assert(v0 != v1);

    FlipEdge flipedge;

    flipedge.connect[0] = v0;
    flipedge.connect[1] = v1;

    std::set<Face*> setA = vFaceMap[v0];
    std::set<Face*> setB = vFaceMap[v1];

    vector<Face*> setC;
    set_intersection(setA.begin(), setA.end(),
                     setB.begin(), setB.end(),
                     inserter(setC, setC.begin()));

    flipedge.faces[0] = NULL;
    flipedge.faces[1] = NULL;

    int numNeighs = setC.size();

    switch (numNeighs)
    {
    case 1:
        flipedge.faces[0] = setC[0];
        flipedge.faces[1] = NULL;
        break;
    case 2:
        flipedge.faces[0] = setC[0];
        flipedge.faces[1] = setC[1];
        break;
    default:
        cout << " Fatal Error: Invalid number of edge neighbours " << numNeighs << endl;
        exit(0);
    }

    if (flipedge.faces[0])
        flipedge.oppositeNodes[0] = opposite_vertex(flipedge.faces[0], v0, v1);

    if (flipedge.faces[1])
        flipedge.oppositeNodes[1] = opposite_vertex(flipedge.faces[1], v0, v1);

    return flipedge;
}

///////////////////////////////////////////////////////////////////////////////

void EdgeFlipping::build_relations()
{
    int numFaces = mesh.faces.size();
    faces.resize(numFaces);

    for (int i = 0; i < numFaces; i++)
    {
        faces[i] = new Face;
        faces[i]->connect = mesh.faces[i].connect;
    }

    vFaceMap.clear();

    for (int i = 0; i < numFaces; i++)
    {
        Face *face = faces[i];
        assert(face->getNumNodes() == 3);
        for (int j = 0; j < face->getNumNodes(); j++)
        {
            Vertex *vertex = face->getVertex(j);
            vFaceMap[vertex].insert(face);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

Vertex* EdgeFlipping::opposite_vertex(Face *face, Vertex *v0, Vertex *v1)
{
    assert(face->getNumNodes() == 3);

    for (int i = 0; i < 3; i++)
    {
        Vertex *av = face->getVertex((i + 1) % 3);
        Vertex *bv = face->getVertex((i + 2) % 3);
        if (av == v0 && bv == v1) return face->getVertex(i);
        if (av == v1 && bv == v0) return face->getVertex(i);
    }
    cout << "Fatal Error: Invalid vertex search " << endl;
    exit(0);
    return NULL;
}

///////////////////////////////////////////////////////////////////////////////

bool FlipEdge::isFlipAllowed() const
{
    if (faces[0] == NULL) return 0;
    if (faces[1] == NULL) return 0;

    if (oppositeNodes[0] == NULL) return 0;
    if (oppositeNodes[1] == NULL) return 0;

    double dist0 = get_edge_length(connect[0], connect[1]);
    double dist1 = get_edge_length(oppositeNodes[0], oppositeNodes[1]);
    if (dist0 > dist1) return 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int EdgeFlipping::commit(const FlipEdge &edge)
{
    Face *t1 = edge.faces[0];
    Face *t2 = edge.faces[1];
    Vertex *v1 = edge.connect[0];
    Vertex *v2 = edge.connect[1];
    Vertex *ov1 = edge.oppositeNodes[0];
    Vertex *ov2 = edge.oppositeNodes[1];

    vFaceMap[v1].erase(t2);
    vFaceMap[v2].erase(t1);

    vFaceMap[ov1].insert(t2);
    vFaceMap[ov2].insert(t1);

    vector<Vertex*> vconn(3);
    vconn[0] = ov2;
    vconn[1] = ov1;
    vconn[2] = v1;
    t1->connect = vconn;

    vconn[0] = v2;
    vconn[1] = ov1;
    vconn[2] = ov2;
    t2->connect = vconn;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void EdgeFlipping::saveAs(const string &filename)
{
    ofstream ofile(filename.c_str(), ios::out);

    ofile << "OFF " << endl;
    ofile << mesh.nodes.size() << " " << mesh.faces.size() << " 0 " << endl;

    for (int i = 0; i < mesh.nodes.size(); i++)
    {
        Vertex *vertex = mesh.nodes[i];
        ofile << vertex->uvCoords[0] << " "
                << vertex->uvCoords[1] << "  0 " << endl;
    }

    for (int i = 0; i < faces.size(); i++)
    {
        Face *face = faces[i];
        ofile << " 3 " << face->connect[0]->id << " "
                << face->connect[1]->id << " "
                << face->connect[2]->id << endl;
    }
}

void EdgeFlipping::execute()
{
    build_relations();

    FlipEdge flipedge;
    int numFaces = faces.size();

    int progress, index = 0;
    while (1)
    {
        progress = 0;
        for (int i = 0; i < numFaces; i++)
        {
            int nConn = faces[i]->getNumNodes();
            assert(nConn == 3);
            for (int j = 0; j < nConn; j++)
            {
                Vertex *v0 = faces[i]->getVertex((j + 1) % 3);
                Vertex *v1 = faces[i]->getVertex((j + 2) % 3);
                flipedge = build_edge(v0, v1);
                if (flipedge.isFlipAllowed())
                {
                    cout << "Commiting : " << endl;
                    commit(flipedge);
                    /*
                                        ostringstream oss;
                                        oss << "flipedge" << index << ".off";
                                        saveAs( oss.str() );
                                        index++;
                     */
                    progress = 0;
                    break;
                }
            }
        }
        if (!progress) break;
    }
}
