
iBase_EntityHandle
Mesh::get_MOAB_Handle(iMesh_Instance imesh, Vertex *vertex)
{
    int err;

    iBase_EntityHandle newHandle = vertex->get_MOAB_Handle();
    if (newHandle)
        return newHandle;

    Point3D p = vertex->getXYZCoords();
    iMesh_createVtx(imesh, p[0], p[1], p[2], &newHandle, &err);
    assert(!err);
    vertex->set_MOAB_Handle(newHandle);

    return newHandle;
}

///////////////////////////////////////////////////////////////////////////////

iBase_EntityHandle
Mesh::get_MOAB_Handle(iMesh_Instance imesh, Face *face)
{
    int status, err;

    // Not a good way. Shouldn't delete it if the connectivity is unchanged.

    iBase_EntityHandle newHandle = face->get_MOAB_Handle();
    if (newHandle)
        iMesh_deleteEnt(imesh, newHandle, &err);

    vector<iBase_EntityHandle> connect;

    int nnodes = face->getSize(0);
    connect.resize(nnodes);

    for (int j = 0; j < nnodes; j++)
    {
        Vertex *v = face->getNodeAt(j);
        connect[j] = get_MOAB_Handle(imesh, v);
    }

    switch (nnodes)
    {
    case 3:
        iMesh_createEnt(imesh, iMesh_TRIANGLE, &connect[0], nnodes, &newHandle,
                        &status, &err);
        assert(!err);
        break;
    case 4:
        iMesh_createEnt(imesh, iMesh_QUADRILATERAL, &connect[0], nnodes,
                        &newHandle, &status, &err);
        assert(!err);
        break;
    default:
        iMesh_createEnt(imesh, iMesh_POLYGON, &connect[0], nnodes, &newHandle,
                        &status, &err);
        assert(!err);
        break;
    }

    face->set_MOAB_Handle(newHandle);

    return newHandle;
}

///////////////////////////////////////////////////////////////////////////////

int
Mesh::toMOAB(iMesh_Instance &imesh, iBase_EntitySetHandle entitySet)
{
    int err;

    if (imesh == 0)
    {
        int optlen = 0;
        char *options = NULL;
        iMesh_newMesh(options, &imesh, &err, optlen);
        assert(!err);
    }

    search_boundary();

    const char *tagname = "fixed";
    int namelen = strlen(tagname);

    iBase_TagHandle idtag;

    iMesh_getTagHandle(imesh, tagname, &idtag, &err, namelen);

    if (err)
        iMesh_createTag(imesh, tagname, 1, iBase_INTEGER, &idtag, &err, namelen);

    assert(!err);

    iBase_EntityHandle newHandle;

    size_t numnodes = getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = getNodeAt(i);
        newHandle = get_MOAB_Handle(imesh, v);
        assert(newHandle);
        if (entitySet)
            iMesh_addEntToSet(imesh, newHandle, entitySet, &err);
        int bmark = v->getBoundaryMark();
        iMesh_setIntData(imesh, newHandle, idtag, bmark, &err);
        assert(!err);
    }

    size_t numfaces = getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *f = getFaceAt(i);
        newHandle = get_MOAB_Handle(imesh, f);
        if (entitySet)
            iMesh_addEntToSet(imesh, newHandle, entitySet, &err);
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

int
Mesh::fromMOAB(iMesh_Instance imesh, iBase_EntitySetHandle entitySet)
{
    int err;
    int numNodes, numFaces;

    if (entitySet == 0)
        iMesh_getRootSet(imesh, &entitySet, &err);

    iMesh_getNumOfType(imesh, entitySet, iBase_VERTEX, &numNodes, &err);
    assert(!err);

    if (numNodes == 0)
    {
        cout << "Warning: There are no nodes in iMesh " << endl;
        return 1;
    }

    reserve(numNodes, 0);

    SimpleArray<iBase_EntityHandle> nodeHandles;
    iMesh_getEntities(imesh, entitySet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                      ARRAY_INOUT(nodeHandles), &err);

    typedef map<iBase_EntityHandle, Vertex*> Moab2JaalNodeMap;
    Moab2JaalNodeMap moab2jaalNodes;

    Point3D p3d;
    double x, y, z;
    for (int i = 0; i < numNodes; i++)
    {
        Vertex *vtx = Vertex::newObject();
        moab2jaalNodes[nodeHandles[i]] = vtx;
        iMesh_getVtxCoord(imesh, nodeHandles[i], &x, &y, &z, &err);
        p3d[0] = x;
        p3d[1] = y;
        p3d[2] = z;
        vtx->setXYZCoords(p3d);
        vtx->setID(i);
        vtx->set_MOAB_Handle(nodeHandles[i]);
        addNode(vtx);
    }

    iMesh_getNumOfType(imesh, entitySet, iBase_FACE, &numFaces, &err);

    if (numNodes == 0)
    {
        cout << "Warning: There are no faces in iMesh " << endl;
        return 1;
    }

    reserve(numFaces, 2);

    NodeSequence connect(3);
    SimpleArray<iBase_EntityHandle> tfaceHandles, qfaceHandles, facenodes;

    iMesh_getEntities(imesh, entitySet, iBase_FACE, iMesh_TRIANGLE, ARRAY_INOUT(
                                                                                tfaceHandles), &err);

    size_t numTris = tfaceHandles.size();
    if (numTris)
    {
        connect.resize(3);
        for (size_t i = 0; i < numTris; i++)
        {
            iMesh_getEntAdj(imesh, tfaceHandles[i], iBase_VERTEX, ARRAY_INOUT(
                                                                              facenodes), &err);
            for (int j = 0; j < 3; j++)
                connect[j] = moab2jaalNodes[facenodes[j]];
            Face *face = new Face;
            face->setNodes(connect);
            face->set_MOAB_Handle(tfaceHandles[i]);
            addFace(face);
        }
        facenodes.clear();
    }

    iMesh_getEntities(imesh, entitySet, iBase_FACE, iMesh_QUADRILATERAL,
                      ARRAY_INOUT(qfaceHandles), &err);

    size_t numQuads = qfaceHandles.size();
    if (numQuads)
    {
        connect.resize(4);
        for (size_t i = 0; i < numQuads; i++)
        {
            iMesh_getEntAdj(imesh, qfaceHandles[i], iBase_VERTEX, ARRAY_INOUT(
                                                                              facenodes), &err);
            for (int j = 0; j < 4; j++)
                connect[j] = moab2jaalNodes[facenodes[j]];
            Face *face = new Face;
            face->setNodes(connect);
            face->set_MOAB_Handle(qfaceHandles[i]);
            addFace(face);
        }
    }
}

vector<int>
Jaal::getVertexFaceDegrees(iMesh_Instance &imesh)
{
    Mesh *jmesh = new Mesh;
    jmesh->fromMOAB(imesh);
    vector<int> quality = jmesh->get_topological_statistics();
    delete jmesh;
    return quality;
}

#endif

