#include "VolMesher.h"

///////////////////////////////////////////////////////////////////////////////

#ifdef NETGEN
VolumeMesher* NetGenVolumeMesher::newObject()
{
    VolumeMesher *mesher = new NetGenVolumeMesher;
    return mesher;
}
#endif

///////////////////////////////////////////////////////////////////////////////

VolumeMesher* TetGenVolumeMesher::newObject()
{
    VolumeMesher *mesher = new TetGenVolumeMesher;
    return mesher;
}

///////////////////////////////////////////////////////////////////////////////

VolumeMesher* VolumeMesherFactory::getProduct(const string &product)
{
#ifdef NETGEN
    if (product == "NetGen") return NetGenVolumeMesher::newObject();
#endif

    if (product == "TetGen") return TetGenVolumeMesher::newObject();

    return NULL;
}

///////////////////////////////////////////////////////////////////////////////

int VolumeMesher :: discretize(GeomMesh &geomesh, iBase_EntityHandle gCell)
{
    int status, err;

    iMesh_Instance mesh = geomesh.mesh;
    iGeom_Instance geom = geomesh.geom;
    iRel_Instance assoc = geomesh.assoc;
    iRel_RelationHandle relation02 = geomesh.relation02;

    VolumeMesher *volmesher = geomesh.volMesher;
    assert(volmesher);

    SimpleArray<iBase_EntityHandle> gFaces, mFaces, mModes, faceNodes;
    iGeom_getEntAdj(geom, gCell, iBase_FACE, ARRAY_INOUT(gFaces), &err);

    int numGeomFaces = gFaces.size();

    std::map<iBase_EntityHandle, Vertex*> uvVertex;
    double x, y, z;
    Face newface;
    int id = 0;
    vector<Face> boundFaces;

    vector<iBase_EntityHandle> nodeHandles;
    vector<iBase_EntitySetHandle> faceMeshSets(numGeomFaces);
    for (int i = 0; i < numGeomFaces; i++)
    {
        iRel_getEntSetAssociation( assoc, relation02, gFaces[i], 0, 
                                   &faceMeshSets[i], &err);

        mFaces.clear();
        iMesh_getEntities(mesh, faceMeshSets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mFaces), &err);

        for (int j = 0; j < mFaces.size(); j++)
        {
            faceNodes.clear();
            iMesh_getEntAdj( mesh, mFaces[j], iBase_VERTEX, ARRAY_INOUT(faceNodes), 
                             &err);
            for (int k = 0; k < faceNodes.size(); k++)
            {
                if (uvVertex.find(faceNodes[k]) == uvVertex.end())
                {
                    iMesh_getVtxCoord(mesh, faceNodes[k], &x, &y, &z, &err);
                    Vertex *vertex = new Vertex;
                    vertex->id = id++;
                    uvVertex[faceNodes[k]] = vertex;
                    nodeHandles.push_back(faceNodes[k]);
                }
            }

            for (int k = 0; k < faceNodes.size(); k++)
                newface.connect[k] = uvVertex[ faceNodes[k] ];
            boundFaces.push_back(newface);
        }
    }

    int numBoundNodes = uvVertex.size();

    VolumeMesh volmesh = getMesh(boundFaces);

    vector<Vertex*> interior = volmesher->getInsertedNodes();

    iBase_EntityHandle newHandle;

    for (int i = 0; i < interior.size(); i++)
    {
        Vertex *vertex = interior[i];
        x = vertex->xyzCoords[0];
        y = vertex->xyzCoords[0];
        z = vertex->xyzCoords[0];
        iMesh_createVtx(mesh, x, y, z, &newHandle, &err);
        nodeHandles.push_back(newHandle);
    }

    vector<Cell> volElements = volmesh.cells;
    int numCells = volElements.size();

    vector<iBase_EntityHandle> vConnect, cellHandles;
    cellHandles.resize(numCells);

    for (int i = 0; i < numCells; i++)
    {
        Cell cell = volElements[i];
        int numNodes = cell.getNumNodes();
        vConnect.resize(numNodes);
        for (int j = 0; j < numNodes; j++)
        {
            Vertex *vertex = cell.getVertex(j);
            vConnect[j] = nodeHandles[vertex->id];
        }
        switch (numNodes)
        {
        case 4:
            iMesh_createEnt( mesh, iMesh_TETRAHEDRON, &vConnect[0], numNodes, 
                             &newHandle, &status, &err);
            break;
        case 8:
            iMesh_createEnt( mesh, iMesh_HEXAHEDRON, &vConnect[0], numNodes, 
                             &newHandle, &status, &err);
            break;
        }
        cellHandles[i] = newHandle;
    }

    iBase_EntitySetHandle entitySet;
    iMesh_createEntSet(mesh, 1, &entitySet, &err);
    iMesh_addEntArrToSet(mesh, &cellHandles[0], numCells, entitySet, &err);

    iBase_TagHandle geom_id_tag, mesh_id_tag, dim_tag;
    const char *tag1 = "GLOBAL_ID";
    int namelen = strlen(tag1);
    iGeom_getTagHandle(geom, tag1, &geom_id_tag, &err, namelen);

    int cellID;
    iGeom_getIntData(geom, gCell, geom_id_tag, &cellID, &err);

    iMesh_getTagHandle(mesh, tag1, &mesh_id_tag, &err, namelen);
    iMesh_setEntSetIntData(mesh, entitySet, mesh_id_tag, cellID, &err);

    const char *tag2 = "GEOM_DIMENSION";
    namelen = strlen(tag2);
    iMesh_getTagHandle(mesh, tag2, &dim_tag, &err, namelen);
    assert(!err);

    int dim = 3;
    iMesh_setEntSetIntData(mesh, entitySet, dim_tag, dim, &err);

    iRel_setEntSetAssociation(assoc, relation02, gCell, entitySet, &err);

    for (int i = 0; i < gFaces.size(); i++)
        iMesh_addPrntChld(mesh, entitySet, faceMeshSets[i], &err);
}

////////////////////////////////////////////////////////////////////////////////

#ifdef NETGEN
VolumeMesh NetGenVolumeMesher::getMesh(const vector<Face> &boundFaces) const
{
    double xyzCoords[3];
    int numNodes, numFaces, numCells, numBoundNodes, numBoundFaces;
    vector<int> connect;
    VolumeMesh volmesh;

    numBoundFaces = boundFaces.size();

    std::map<int, Vertex*> vmap;
    for (int i = 0; i < numBoundFaces; i++)
    {
        Face face = boundFaces[i];
        assert(face.getNumNodes() == 3);
        for (int j = 0; j < face.getNumNodes(); j++)
        {
            Vertex *v = face.connect[j];
            if (vmap.find(v->id) == vmap.end()) vmap[v->id] = v;
        }
    }
    numBoundNodes = vmap.size();

    for (int i = 0; i < numBoundNodes; i++)
    {
        if (vmap.find(i) == vmap.end())
        {
            cout << "Error: Node IDs are not continuous: " << endl;
            return volmesh;
        }
    }

    vector<Vertex*> nodes(numBoundNodes);
    for (int i = 0; i < numBoundNodes; i++)
    {
        Vertex *vertex = vmap[i];
        vertex->id += 1; // Counting starts from 1 in netgen
        nodes[i] = vertex;
    }
    vmap.clear();

    Ng_Mesh *ngmesh = Ng_NewMesh();

    for (int i = 0; i < numBoundNodes; i++)
    {
        Vertex *vertex = nodes[i];
        assert(vertex);
        Ng_AddPoint(ngmesh, vertex->xyzCoords);
    }

    connect.resize(3);
    for (int i = 0; i < numBoundFaces; i++)
    {
        connect[0] = boundFaces[i].connect[0]->id;
        connect[1] = boundFaces[i].connect[1]->id;
        connect[2] = boundFaces[i].connect[2]->id;
        Ng_AddSurfaceElement(ngmesh, NG_TRIG, &connect[0]);
    }

    Ng_Meshing_Parameters meshParams;
    Ng_GenerateVolumeMesh(ngmesh, &meshParams);

    numNodes = Ng_GetNP(ngmesh);
    numFaces = Ng_GetNSE(ngmesh);
    numCells = Ng_GetNE(ngmesh);

    assert(numFaces == numBoundFaces);

    nodes.resize(numNodes);
    insertedNodes.reserve(numNodes - numBoundNodes);

    for (int i = 0; i < numBoundNodes; i++)
        nodes[i]->id -= 1; // Since it was locally incremented by 1

    int index = 0;
    for (int i = numBoundNodes; i < numNodes; i++)
    {
        Ng_GetPoint(ngmesh, i + 1, xyzCoords);
        Vertex *vertex = new Vertex;
        vertex->id = i;

        vertex->xyzCoords[0] = xyzCoords[0];
        vertex->xyzCoords[1] = xyzCoords[1];
        vertex->xyzCoords[2] = xyzCoords[2];

        nodes[i] = vertex;
        insertedNodes.push_back(vertex);
    }

    volmesh.cells.resize(numCells);

    connect.resize(100);
    for (int i = 0; i < numCells; i++)
    {
        int vol_elem_type = Ng_GetVolumeElement(ngmesh, i + 1, &connect[0]);
        assert(vol_elem_type == NG_TET);
        volmesh.cells[i].connect.resize(4);
        volmesh.cells[i].connect[0] = nodes[ connect[0] ];
        volmesh.cells[i].connect[1] = nodes[ connect[1] ];
        volmesh.cells[i].connect[2] = nodes[ connect[2] ];
        volmesh.cells[i].connect[3] = nodes[ connect[3] ];
    }

    Ng_DeleteMesh(ngmesh);

    return volmesh;
}
#endif
///////////////////////////////////////////////////////////////////////////////

VolumeMesh TetGenVolumeMesher::getMesh(const vector<Face> &boundFaces) const
{
    double xyzCoords[3];
    int numNodes, numFaces, numCells, numBoundNodes, numBoundFaces;
    vector<int> connect;
    VolumeMesh volmesh;

    numBoundFaces = boundFaces.size();

    std::map<int, Vertex*> vmap;
    for (int i = 0; i < numBoundFaces; i++)
    {
        Face face = boundFaces[i];
        assert(face.getNumNodes() == 3);
        for (int j = 0; j < face.getNumNodes(); j++)
        {
            Vertex *v = face.connect[j];
            if (vmap.find(v->id) == vmap.end()) vmap[v->id] = v;
        }
    }
    numBoundNodes = vmap.size();

    for (int i = 0; i < numBoundNodes; i++)
    {
        if (vmap.find(i) == vmap.end())
        {
            cout << "Error: Node IDs are not continuous: " << endl;
            return volmesh;
        }
    }

    vector<Vertex*> nodes(numBoundNodes);
    for (int i = 0; i < numBoundNodes; i++)
    {
        Vertex *vertex = vmap[i];
        nodes[i] = vertex;
    }
    vmap.clear();

    tetgenio in, out;
    tetgenio::facet *f;
    tetgenio::polygon *p;
    int i;

    in.firstnumber = 0;

    in.numberofpoints = numBoundNodes;
    in.pointlist = new REAL[3 * numBoundNodes];
    for (int i = 0; i < numBoundNodes; i++)
    {
        in.pointlist[3 * i + 0] = nodes[i]->xyzCoords[0];
        in.pointlist[3 * i + 1] = nodes[i]->xyzCoords[1];
        in.pointlist[3 * i + 2] = nodes[i]->xyzCoords[2];
    }

    in.numberoffacets = numBoundFaces;
    in.facetlist = new tetgenio::facet[numBoundFaces];
    in.facetmarkerlist = NULL;

    for (int i = 0; i < numBoundFaces; i++)
    {
        f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[3];
        p->vertexlist[0] = boundFaces[i].connect[0]->id;
        p->vertexlist[1] = boundFaces[i].connect[1]->id;
        p->vertexlist[2] = boundFaces[i].connect[2]->id;
    }

    tetrahedralize("pq1.414Y", &in, &out);

    numNodes = out.numberofpoints;
    nodes.reserve(numNodes);

    int index = 0;
    for (int i = numBoundNodes; i < numNodes; i++)
    {
        Vertex *v = new Vertex;
        v->id = i;
        v->xyzCoords[0] = out.pointlist[3 * i + 0];
        v->xyzCoords[1] = out.pointlist[3 * i + 1];
        v->xyzCoords[2] = out.pointlist[3 * i + 2];
        nodes.push_back(v);
    }

    numCells = out.numberoftetrahedra;

    volmesh.cells.resize(numCells);
    for (int i = 0; i < numCells; i++)
    {
        volmesh.cells[i].connect.resize(4);
        volmesh.cells[i].connect[0] = nodes[ out.tetrahedronlist[4 * i + 0] ];
        volmesh.cells[i].connect[1] = nodes[ out.tetrahedronlist[4 * i + 1] ];
        volmesh.cells[i].connect[2] = nodes[ out.tetrahedronlist[4 * i + 2] ];
        volmesh.cells[i].connect[2] = nodes[ out.tetrahedronlist[4 * i + 3] ];
    }
}

///////////////////////////////////////////////////////////////////////////////

void example_volume_mesher(const string &filename)
{
    ifstream ifile(filename.c_str(), ios::in);
    if (ifile.fail()) return;

    string str;
    ifile >> str;
    assert(str == "OFF");

    int numNodes, numFaces, numEdges;

    ifile >> numNodes >> numFaces >> numEdges;

    assert(numNodes && numFaces);

    vector<Vertex*> nodes(numNodes);

    double x, y, z;
    for (int i = 0; i < numNodes; i++)
    {
        ifile >> x >> y >> z;
        Vertex *v = new Vertex;
        v->xyzCoords[0] = x;
        v->xyzCoords[1] = y;
        v->xyzCoords[2] = z;
        v->id = i;
        nodes[i] = v;
    }

    vector<Face> faces(numFaces);

    int nConn, n0, n1, n2;
    for (int i = 0; i < numFaces; i++)
    {
        ifile >> nConn >> n0 >> n1 >> n2;
        assert(nConn == 3);
        faces[i].connect.resize(3);
        faces[i].connect[0] = nodes[n0];
        faces[i].connect[1] = nodes[n1];
        faces[i].connect[2] = nodes[n2];
    }

    VolumeMesher  *volmesher;
    VolumeMesh     tetmesh;

    cout << "Volume meshing with TetGen : " << endl;
    volmesher = new  TetGenVolumeMesher;
    tetmesh = volmesher->getMesh(faces);
    cout << " Number of Tets : " << tetmesh.getSize(3) << endl;
    delete volmesher;

#ifdef NETGEN
    cout << "Volume meshing with NetGen : " << endl;
    volmesher = new  NetGenVolumeMesher;
    tetmesh = volmesher->getMesh(faces);
    cout << " Number of Tets : " << tetmesh.getSize(3) << endl;
    delete volmesher;
#endif

}
///////////////////////////////////////////////////////////////////////////////
