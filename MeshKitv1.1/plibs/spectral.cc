#include "itaps.h"

///////////////////////////////////////////////////////////////////////////////

void getCentroid(iMesh_Instance &mesh, SimpleArray<iBase_EntityHandle> &nodeHandles,
                 double &xc, double &yc, double &zc)
{
    int err;
    xc = 0.0;
    yc = 0.0;
    zc = 0.0;

    int numNodes = nodeHandles.size();

    double x, y, z;
    for (int i = 0; i < numNodes; i++)
    {
        iMesh_getVtxCoord(mesh, nodeHandles[i], &x, &y, &z, &err);
        xc += x;
        yc += y;
        zc += z;
    }
    xc /= (double) numNodes;
    yc /= (double) numNodes;
    zc /= (double) numNodes;
}

///////////////////////////////////////////////////////////////////////////////

void spectral_edge(GeomMesh &geomesh, iBase_EntityHandle edgeHandle,
                   iBase_TagHandle horder_tag, 
                   vector<iBase_EntityHandle> &hOrderNodes, int projectOnGeometry)
{
    int err;
    iMesh_Instance mesh = geomesh.mesh;

    SimpleArray<iBase_EntityHandle> edgeNodes;
    iMesh_getEntAdj(mesh, edgeHandle, iBase_VERTEX, ARRAY_INOUT(edgeNodes), &err);

    double xc, yc, zc;
    getCentroid(mesh, edgeNodes, xc, yc, zc);

    iBase_EntityHandle newHandle;
    iMesh_createVtx(mesh, xc, yc, zc, &newHandle, &err);
    iMesh_setEHData(mesh, edgeHandle, horder_tag, newHandle, &err);

    hOrderNodes.resize(1);
    hOrderNodes[0] = newHandle;

}
///////////////////////////////////////////////////////////////////////////////

void spectral_face(GeomMesh &geomesh, iBase_EntityHandle faceHandle,
                   iBase_TagHandle horder_tag, 
                   vector<iBase_EntityHandle> &hOrderNodes, int projectOnGeometry)
{
    int err;
    iMesh_Instance mesh = geomesh.mesh;

    SimpleArray<iBase_EntityHandle> faceNodes;
    iMesh_getEntAdj(mesh, faceHandle, iBase_VERTEX, ARRAY_INOUT(faceNodes), &err);

    double xc, yc, zc;
    getCentroid(mesh, faceNodes, xc, yc, zc);

    iBase_EntityHandle newHandle;
    iMesh_createVtx(mesh, xc, yc, zc, &newHandle, &err);
    iMesh_setEHData(mesh, faceHandle, horder_tag, newHandle, &err);

    hOrderNodes.resize(1);
    hOrderNodes[0] = newHandle;

}
///////////////////////////////////////////////////////////////////////////////

void spectral_cell( GeomMesh &geomesh, iBase_EntityHandle cellHandle,
                    iBase_TagHandle horder_tag, 
                    vector<iBase_EntityHandle> &hOrderNodes)
{
    int err;
    iMesh_Instance mesh = geomesh.mesh;

    SimpleArray<iBase_EntityHandle> cellNodes;
    iMesh_getEntAdj(mesh, cellHandle, iBase_VERTEX, ARRAY_INOUT(cellNodes), &err);

    double xc, yc, zc;
    getCentroid(mesh, cellNodes, xc, yc, zc);

    iBase_EntityHandle newHandle;
    iMesh_createVtx(mesh, xc, yc, zc, &newHandle, &err);
    iMesh_setEHData(mesh, cellHandle, horder_tag, newHandle, &err);

    hOrderNodes.resize(1);
    hOrderNodes[0] = newHandle;

}
///////////////////////////////////////////////////////////////////////////////

void generate_spectral_elements(GeomMesh &geomesh, int numNodes)
{
    int err;
    iMesh_Instance mesh = geomesh.mesh;
    iGeom_Instance geom = geomesh.geom;
    iRel_Instance assoc = geomesh.assoc;
    iRel_RelationHandle relation02 = geomesh.relation02;
    iBase_EntitySetHandle meshRootSet = geomesh.meshRootSet;
    iBase_EntitySetHandle geomRootSet = geomesh.geomRootSet;

    iBase_EntitySetHandle entitySet;
    vector<iBase_EntityHandle> nodeHandles, hOrderNodes;
    SimpleArray<iBase_EntityHandle> gEdges, gFaces;
    SimpleArray<iBase_EntityHandle> mEdges, mFaces, mCells;

    //////////////////////////////////////////////////////////////////////////
    // Create Tag
    //////////////////////////////////////////////////////////////////////////

    const char *tagname = "HO_POINTS";
    int namelen = strlen(tagname);

    iBase_TagHandle horder_tag;
    iMesh_createTag( mesh, tagname, 1, iBase_ENTITY_HANDLE, &horder_tag, &err, 
                     namelen);

    //////////////////////////////////////////////////////////////////////////
    // Step-1::  place higher order nodes on the edges lies on geometric edges
    //////////////////////////////////////////////////////////////////////////

    iGeom_getEntities(geom, geomRootSet, iBase_EDGE, ARRAY_INOUT(gEdges), &err);
    vector<iBase_EntitySetHandle> edgeMeshSets;

    int numGeoEdges = gEdges.size();
    for (int i = 0; i < numGeoEdges; i++)
    {
        iRel_getEntSetAssociation( assoc, relation02, gEdges[i], 0, 
                                   &edgeMeshSets[i], &err);

        mEdges.clear();
        iMesh_getEntities( mesh, edgeMeshSets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES, 
                           ARRAY_INOUT(mEdges), &err);
        for (int j = 0; j < mEdges.size(); j++)
        {
            spectral_edge(geomesh, mEdges[j], horder_tag, hOrderNodes, 1);
            iMesh_createEntSet(mesh, 1, &entitySet, &err);
            iMesh_addEntArrToSet( mesh, &hOrderNodes[0], (int) hOrderNodes.size(), 
                                  entitySet, &err);
        }
    }
    ///////////////////////////////////////////////////////////////////////////
    // Step-2:  place higher order nodes on the edges lies on geometric faces.
    ///////////////////////////////////////////////////////////////////////////

    iGeom_getEntities(geom, geomRootSet, iBase_FACE, ARRAY_INOUT(gFaces), &err);
    vector<iBase_EntitySetHandle> faceMeshSets;

    int numGeoFaces = gFaces.size();
    for (int i = 0; i < numGeoFaces; i++)
    {
        iRel_getEntSetAssociation(assoc, relation02, gFaces[i], 0, &faceMeshSets[i], &err);

        mFaces.clear();
        iMesh_getEntities( mesh, faceMeshSets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES, 
                           ARRAY_INOUT(mFaces), &err);
        for (int j = 0; j < mFaces.size(); j++)
        {
            mEdges.clear();
            iMesh_getEntAdj(mesh, mFaces[j], iBase_EDGE, ARRAY_INOUT(mEdges), &err);

            for (int k = 0; k < mEdges.size(); k++)
            {
                spectral_edge(geomesh, mEdges[k], horder_tag, hOrderNodes, 1);
                iMesh_createEntSet(mesh, 1, &entitySet, &err);
                iMesh_addEntArrToSet(mesh, &hOrderNodes[0], (int) hOrderNodes.size(), entitySet, &err);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Step-3: place higher order nodes on interior edges.
    ////////////////////////////////////////////////////////////////////////////

    iMesh_getEntities(mesh, meshRootSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(mCells), &err);

    for (int i = 0; i < mCells.size(); i++)
    {
        mEdges.clear();
        iMesh_getEntAdj(mesh, mCells[i], iBase_EDGE, ARRAY_INOUT(mEdges), &err);

        for (int j = 0; j < mEdges.size(); j++)
        {
            spectral_edge(geomesh, mEdges[j], horder_tag, hOrderNodes, 0);
            iMesh_createEntSet(mesh, 1, &entitySet, &err);
            iMesh_addEntArrToSet( mesh, &hOrderNodes[0], (int) hOrderNodes.size(), 
                                  entitySet, &err);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step-4:  place higher order nodes on geometric faces.
    ///////////////////////////////////////////////////////////////////////////

    iGeom_getEntities(geom, geomRootSet, iBase_FACE, ARRAY_INOUT(gFaces), &err);

    numGeoFaces = gFaces.size();
    for (int i = 0; i < numGeoFaces; i++)
    {
        iRel_getEntSetAssociation(assoc, relation02, gFaces[i], 0, &faceMeshSets[i], &err);

        mFaces.clear();
        iMesh_getEntities( mesh, faceMeshSets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES, 
                           ARRAY_INOUT(mFaces), &err);
        for (int j = 0; j < mFaces.size(); j++)
        {
            spectral_face(geomesh, mFaces[j], horder_tag, hOrderNodes, 1);
            iMesh_createEntSet(mesh, 1, &entitySet, &err);
            iMesh_addEntArrToSet( mesh, &hOrderNodes[0], (int) hOrderNodes.size(), 
                                  entitySet, &err);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step-5:   place higher order nodes on interior faces.
    ///////////////////////////////////////////////////////////////////////////

    iMesh_getEntities( mesh, meshRootSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, 
                       ARRAY_INOUT(mCells), &err);

    for (int i = 0; i < mCells.size(); i++)
    {
        mEdges.clear();
        iMesh_getEntAdj(mesh, mCells[i], iBase_FACE, ARRAY_INOUT(mFaces), &err);

        for (int j = 0; j < mFaces.size(); j++)
        {
            spectral_face(geomesh, mFaces[j], horder_tag, hOrderNodes, 0);
            iMesh_createEntSet(mesh, 1, &entitySet, &err);
            iMesh_addEntArrToSet( mesh, &hOrderNodes[0], (int) hOrderNodes.size(), 
                                  entitySet, &err);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step-6:  rearrage nodes on each face.
    ///////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // Step-7:  place higher order nodes in the interior of cells.
    //////////////////////////////////////////////////////////////////////////
    //
    iMesh_getEntities( mesh, meshRootSet, iBase_REGION, iMesh_ALL_TOPOLOGIES, 
                       ARRAY_INOUT(mCells), &err);

    for (int i = 0; i < mCells.size(); i++)
    {
        spectral_cell(geomesh, mCells[i], horder_tag, hOrderNodes);
        iMesh_createEntSet(mesh, 1, &entitySet, &err);
        iMesh_addEntArrToSet( mesh, &hOrderNodes[0], (int) hOrderNodes.size(), 
                              entitySet, &err);
    }

    //
    // Step-8: Rearrange higher order nodes in each cells.
    //
}

///////////////////////////////////////////////////////////////////////////////

