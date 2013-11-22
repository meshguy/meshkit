#include "EdgeMesher.h"

iBase_EntitySetHandle
EdgeMesher :: discretize_close_edge(iBase_EntityHandle gEdge, int numEdges)
{
    int err, status;
    double x, y, z, umin, umax;

    iMesh_Instance mesh = geomesh.mesh;
    iGeom_Instance geom = geomesh.geom;
    iRel_Instance assoc = geomesh.assoc;
    iRel_RelationHandle relation00 = geomesh.relation00;

    iGeom_getEntURange(geom, gEdge, &umin, &umax, &err);

    int numNodes = numEdges;

    vector<iBase_EntityHandle> nodeHandles;
    nodeHandles.resize(numNodes);

    double u, du = (umax - umin) / (double) numEdges;

    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj(geom, gEdge, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

    iRel_getEntEntAssociation(assoc, relation00, gNodes[0], 0, &nodeHandles[0], &err);

    iBase_EntityHandle newHandle;
    for (int i = 1; i < numNodes; i++)
    {
        u = umin + i*du;
        iGeom_getEntUtoXYZ(geom, gEdge, u, &x, &y, &z, &err);
        iMesh_createVtx(mesh, x, y, z, &newHandle, &err);
        nodeHandles[i] = newHandle;
    }

    vector<iBase_EntityHandle> edgeHandles, eConnect(2);
    edgeHandles.resize(numEdges);

    for (int i = 0; i < numEdges; i++)
    {
        eConnect[0] = nodeHandles[i];
        eConnect[1] = nodeHandles[(i + 1) % numNodes];
        iMesh_createEnt( mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &newHandle, 
                         &status, &err);
        assert(!err);
        edgeHandles[i] = newHandle;
    }

    iBase_EntitySetHandle entitySet;
    iMesh_createEntSet(mesh, 1, &entitySet, &err);
    iMesh_addEntArrToSet(mesh, &edgeHandles[0], numEdges, entitySet, &err);

    return entitySet;
}

////////////////////////////////////////////////////////////////////////////////

iBase_EntitySetHandle
EdgeMesher :: discretize_open_edge(iBase_EntityHandle gEdge, int numEdges )
{
    int err, status, numNodes;
    iBase_EntityHandle newHandle;
    double x, y, z, u, du, umin, umax;

    iMesh_Instance mesh = geomesh.mesh;
    iGeom_Instance geom = geomesh.geom;
    iRel_Instance assoc = geomesh.assoc;
    iRel_RelationHandle relation00 = geomesh.relation00;

    iGeom_getEntURange(geom, gEdge, &umin, &umax, &err);
    //
    // Presently. Very simple discretization, but later replaced with curvature
    // based discretization.

    du = (umax - umin) / (double) numEdges;

    numNodes = numEdges + 1;
    vector<iBase_EntityHandle> nodeHandles;
    nodeHandles.resize(numNodes);

    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj(geom, gEdge, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

    iRel_getEntEntAssociation( assoc, relation00, gNodes[0], 0, &nodeHandles[0], 
                               &err);

    for (int i = 1; i < numNodes - 1; i++)
    {
        u = umin + i*du;
        iGeom_getEntUtoXYZ(geom, gEdge, u, &x, &y, &z, &err);
        iMesh_createVtx(mesh, x, y, z, &newHandle, &err);
        nodeHandles[i] = newHandle;
    }

    iRel_getEntEntAssociation( assoc, relation00, gNodes[1], 0, 
                               &nodeHandles[numNodes - 1], &err);

    vector<iBase_EntityHandle> edgeHandles, eConnect(2);
    edgeHandles.resize(numEdges);

    for (int i = 0; i < numEdges; i++)
    {
        eConnect[0] = nodeHandles[i];
        eConnect[1] = nodeHandles[i + 1];
        iMesh_createEnt( mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &newHandle, 
                         &status, &err);
        assert(!err);
        edgeHandles[i] = newHandle;
    }

    iBase_EntitySetHandle entitySet;
    iMesh_createEntSet(mesh, 1, &entitySet, &err);
    assert(!err);
    iMesh_addEntArrToSet(mesh, &edgeHandles[0], numEdges, entitySet, &err);
    assert(!err);

    return entitySet;
}

////////////////////////////////////////////////////////////////////////////////

int 
EdgeMesher :: discretize(GeomMesh &gmsh, iBase_EntityHandle gEdge, int numEdges)
{
    geomesh = gmsh;

    iMesh_Instance mesh = geomesh.mesh;
    iGeom_Instance geom = geomesh.geom;
    iRel_Instance assoc = geomesh.assoc;
    iRel_RelationHandle relation = geomesh.relation02;

    int err;
    SimpleArray<iBase_EntityHandle> gNodes;
    iGeom_getEntAdj(geom, gEdge, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

    iBase_EntitySetHandle entitySet;

    switch (gNodes.size())
    {
    case 1:
        entitySet = discretize_close_edge(gEdge, numEdges);
        break;
    case 2:
        entitySet = discretize_open_edge(gEdge, numEdges);
        break;
    default:
        cout << "Fatal Error: Invalid Geometric edge " << endl;
        exit(0);
    }

    iBase_TagHandle geom_id_tag, mesh_id_tag, dim_tag;

    const char *tag1 = "GLOBAL_ID";
    int namelen = strlen(tag1);
    iGeom_getTagHandle(geom, tag1, &geom_id_tag, &err, namelen);

    int edgeID;
    iGeom_getIntData(geom, gEdge, geom_id_tag, &edgeID, &err);

    iMesh_getTagHandle(mesh, tag1, &mesh_id_tag, &err, namelen);
    iMesh_setEntSetIntData(mesh, entitySet, mesh_id_tag, edgeID, &err);

    const char *tag2 = "GEOM_DIMENSION";
    namelen = strlen(tag2);
    iMesh_getTagHandle(mesh, tag2, &dim_tag, &err, namelen);
    assert(!err);

    int dim = 1;
    iMesh_setEntSetIntData(mesh, entitySet, dim_tag, dim, &err);

    iRel_setEntSetAssociation(assoc, relation, gEdge, entitySet, &err);
}

