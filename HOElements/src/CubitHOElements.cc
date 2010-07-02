#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <sstream>

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include "SimpleArray.hpp"

#include "SpectralElements.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int readSimpleMesh(const char *filename, iMesh_Instance &mesh)
{
    int err;

    int optlen = 0;
    char *options = NULL;
    iMesh_newMesh(options, &mesh, &err, optlen);
    assert(!err);

    iBase_TagHandle idtag;
    iMesh_getTagHandle(mesh, "GLOBAL_ID", &idtag, &err, strlen("GLOBAL_ID"));
    assert(!err);

    ifstream infile(filename, ios::in);

    if (infile.fail())
    {
        cout << "Error: cann't read file " << filename << endl;
        return 1;
    }

    string s;
    infile >> s;
    assert(s == "#Nodes");

    int numNodes;
    infile >> numNodes;

    vector<iBase_EntityHandle> nodeHandles;
    nodeHandles.resize(numNodes);

    int gid;
    double x, y, z;

    iBase_EntityHandle newhandle;
    int id, index = 0;
    for (int i = 0; i < numNodes; i++)
    {
        infile >> id >> x >> y >> z;
        iMesh_createVtx(mesh, x, y, z, &newhandle, &err);
        nodeHandles[index++] = newhandle;
    }

    infile >> s;
    assert(s == "#Hex");
    int numCells;
    infile >> numCells;

    vector<iBase_EntityHandle> cellHandles;
    cellHandles.resize(numCells);

    vector<iBase_EntityHandle> connect(8);
    int status, n0, n1, n2, n3, n4, n5, n6, n7;
    index = 0;
    for (int i = 0; i < numCells; i++)
    {
        infile >> n0 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7;
        connect[0] = nodeHandles[n0];
        connect[1] = nodeHandles[n1];
        connect[2] = nodeHandles[n3];
        connect[3] = nodeHandles[n2];
        connect[4] = nodeHandles[n4];
        connect[5] = nodeHandles[n5];
        connect[6] = nodeHandles[n7];
        connect[7] = nodeHandles[n6];
        iMesh_createEnt(mesh, iMesh_HEXAHEDRON, &connect[0], 8, &newhandle, &status, &err);
        cellHandles[index++] = newhandle;
    }
}

////////////////////////////////////////////////////////////////////////////////

void buildAssociations(iGeom_Instance &geom, iMesh_Instance &mesh,
                       iRel_Instance &assoc, iRel_RelationHandle &rel)
{

    int err, namelen;

    SimpleArray<iBase_EntitySetHandle> entitySets;
    iBase_EntitySetHandle geom_root_set, mesh_root_set;

    // Get the root sets of the geometry and mesh.
    iMesh_getRootSet(mesh, &mesh_root_set, &err);
    iGeom_getRootSet(geom, &geom_root_set, &err);

    iBase_TagHandle geom_id_tag, mesh_id_tag, geom_dim_tag;
    const char *tag1 = "GLOBAL_ID";
    namelen = strlen(tag1);
    iGeom_getTagHandle(geom, tag1, &geom_id_tag, &err, namelen);
    iMesh_getTagHandle(mesh, tag1, &mesh_id_tag, &err, namelen);

    const char *tag2 = "GEOM_DIMENSION";
    namelen = strlen(tag2);
    iMesh_getTagHandle(mesh, tag2, &geom_dim_tag, &err, namelen);
    assert(!err);

    iRel_newAssoc(0, &assoc, &err, 0);

    iRel_createAssociation(assoc,
                           geom, 0, iRel_IGEOM_IFACE,
                           mesh, 2, iRel_IMESH_IFACE, &rel, &err);
    assert(!err);

    // Get all the entitySet in the mesh
    iMesh_getEntSets(mesh, mesh_root_set, 0, ARRAY_INOUT(entitySets), &err);

    int ncount;
    iBase_EntityHandle gEntity;

    // Map all the geometric edges
    SimpleArray<iBase_EntityHandle> gEdges;
    iGeom_getEntities(geom, geom_root_set, iBase_EDGE, ARRAY_INOUT(gEdges), &err);
    assert(!err);

    int geom_id, geom_dim;
    std::map<int, iBase_EntityHandle> mapEdges;
    for (int i = 0; i < gEdges.size(); i++)
    {
        iGeom_getIntData(geom, gEdges[i], geom_id_tag, &geom_id, &err);
        mapEdges[geom_id] = gEdges[i];
    }

    // Map all the geometric faces ...
    SimpleArray<iBase_EntityHandle> gFaces;
    iGeom_getEntities(geom, geom_root_set, iBase_FACE, ARRAY_INOUT(gFaces), &err);
    assert(!err);

    std::map<int, iBase_EntityHandle> mapFaces;
    for (int i = 0; i < gFaces.size(); i++)
    {
        iGeom_getIntData(geom, gFaces[i], geom_id_tag, &geom_id, &err);
        mapFaces[geom_id] = gFaces[i];
    }

    // Map all the geometric cells ...
    SimpleArray<iBase_EntityHandle> gCells;
    iGeom_getEntities(geom, geom_root_set, iBase_REGION, ARRAY_INOUT(gCells), &err);
    assert(!err);

    std::map<int, iBase_EntityHandle> mapCells;
    for (int i = 0; i < gCells.size(); i++)
    {
        iGeom_getIntData(geom, gCells[i], geom_id_tag, &geom_id, &err);
        mapCells[geom_id] = gCells[i];
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Create Edge Assocations:
    ///////////////////////////////////////////////////////////////////////////////
    cout << " Building Edge Associations " << endl;

    int numEdges;
    iMesh_getNumOfType(mesh, mesh_root_set, iBase_EDGE, &numEdges, &err);

    SimpleArray<iBase_EntityHandle> mEdges;

    int numAssociations = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mEdges.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mEdges), &err);

        if (mEdges.size() && (mEdges.size() != numEdges))
        {
            ncount += mEdges.size();

            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_id_tag, &geom_id, &err);
            assert(!err);

            iMesh_getEntSetIntData(mesh, entitySets[i], geom_dim_tag, &geom_dim, &err);
            assert(!err);

            gEntity = 0;
            switch (geom_dim)
            {
            case 1:

                if (mapEdges.find(geom_id) != mapEdges.end())
                {
                    gEntity = mapEdges[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                    gEntity = mapFaces[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Face not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric Cell not found : " << geom_id << endl;
                    exit(0);
                }
                break;
            default:
                cout << "Error: Invalid geometric dimension " << geom_dim << endl;
                exit(0);
            }

            if (gEntity)
                iRel_setEntSetAssociation(assoc, rel, gEntity, entitySets[i], &err);
        }
    }

    if (numAssociations != mapEdges.size())
        cout << "Warning: There are more edge entitySet than geometric edges " << endl;

    //////////////////////////////////////////////////////////////////////////////
    // Face Association
    //////////////////////////////////////////////////////////////////////////////
    cout << " Building Face Associations " << endl;

    SimpleArray<iBase_EntityHandle> mFaces;

    int numFaces;
    iMesh_getNumOfType(mesh, mesh_root_set, iBase_FACE, &numFaces, &err);

    int mesh_faceid;

    numAssociations = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mFaces.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mFaces), &err);

        if (mFaces.size() && (mFaces.size() != numFaces))
        {
            ncount += mFaces.size();
            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_id_tag, &geom_id, &err);
            assert(!err);
            iMesh_getEntSetIntData(mesh, entitySets[i], geom_dim_tag, &geom_dim, &err);
            assert(!err);

            gEntity = 0;
            switch (geom_dim)
            {
            case 2:
                if (mapFaces.find(geom_id) != mapFaces.end())
                {
                    gEntity = mapFaces[geom_id];
                    numAssociations++;
                }
                else
                {
                    cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                    exit(0);
                }
                break;
            case 3:
                if (mapCells.find(geom_id) != mapCells.end())
                    gEntity = mapCells[geom_id];
                else
                {
                    cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                    exit(0);
                }
                break;
            }
            if (gEntity)
                iRel_setEntSetAssociation(assoc, rel, gEntity, entitySets[i], &err);
        }
    }

    if (numAssociations != mapFaces.size())
        cout << "Warning: There are more face entitySet than geometric faces " << endl;

    //////////////////////////////////////////////////////////////////////////////
    // Cell Association
    //////////////////////////////////////////////////////////////////////////////

    SimpleArray<iBase_EntityHandle> mCells;

    int mesh_cellid;

    int numCells;
    iMesh_getNumOfType(mesh, mesh_root_set, iBase_REGION, &numCells, &err);

    ncount = 0;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mCells.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_REGION, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mCells), &err);

        if (mCells.size() && (mCells.size() != numCells))
        {
            ncount += mCells.size();
            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_id_tag, &geom_id, &err);
            assert(!err);

            if (mapCells.find(geom_id) != mapCells.end())
            {
                if (mapCells.find(geom_id) != mapCells.end())
                {
                    gEntity = mapCells[geom_id];
                    iRel_setEntSetAssociation(assoc, rel, gEntity, entitySets[i], &err);
                }
            }
        }
    }

}
////////////////////////////////////////////////////////////////////////////////

int readGeometry(const char *filename, iGeom_Instance &geom)
{
    int err;
    iGeom_newGeom(NULL, &geom, &err, 0);

    iGeom_load(geom, filename, 0, &err, strlen(filename), 0);

    iBase_EntitySetHandle rootSet;
    iGeom_getRootSet(geom, &rootSet, &err);

    cout << "Model Contents " << endl;
    const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};

    int count;
    for (int i = 0; i <= 3; ++i)
    {
        iGeom_getNumOfType(geom, rootSet, i, &count, &err);
        std::cout << gtype[i] << count << std::endl;
    }

    iBase_TagHandle idtag;
    iGeom_getTagHandle(geom, "GLOBAL_ID", &idtag, &err, strlen("GLOBAL_ID"));
    assert(!err);

    int gid;
    set<int> idset;

    SimpleArray<iBase_EntityHandle> gEdges;
    iGeom_getEntities(geom, rootSet, iBase_EDGE, ARRAY_INOUT(gEdges), &err);

    int numEdges = gEdges.size();
    idset.clear();
    for (int i = 0; i < numEdges; i++)
    {
        iGeom_getIntData(geom, gEdges[i], idtag, &gid, &err);
        idset.insert(gid);
    }

    SimpleArray<iBase_EntityHandle> gFaces;
    iGeom_getEntities(geom, rootSet, iBase_FACE, ARRAY_INOUT(gFaces), &err);

    int numFaces = gFaces.size();
    idset.clear();
    for (int i = 0; i < numFaces; i++)
    {
        iGeom_getIntData(geom, gFaces[i], idtag, &gid, &err);
        idset.insert(gid);
        GFace currface(geom, gFaces[i]);
        ostringstream oss;
        oss << "gface" << i << ".dat";
        currface.saveAs(oss.str());
    }

    ofstream ofile( "modeledges.dat", ios::out);

    vector<Point3D> xyz;
    vector<double>  u;

    ofile << "#Nodes " << 100*gEdges.size() << endl;
    
    int index = 0;
    for (int i = 0; i < gEdges.size(); i++)
    {
        GEdge curredge(geom, gEdges[i]);
        curredge.uniform_u_discretization(100, xyz, u);
        for( int j = 0; j < xyz.size(); j++) 
             ofile << index++ << " " << xyz[j][0] << " " << xyz[j][1] << " " << xyz[j][2] << endl;
    }
  
}

///////////////////////////////////////////////////////////////////////////////

void saveCubitMesh(iGeom_Instance &geom, iMesh_Instance &mesh, const char *filename)
{

    ofstream ofile(filename, ios::out);
    if (ofile.fail()) return;

    int err;

    iBase_TagHandle mesh_idtag;
    iMesh_getTagHandle(mesh, "GLOBAL_ID", &mesh_idtag, &err, strlen("GLOBAL_ID"));
    assert(!err);

    iBase_TagHandle geom_idtag;
    iGeom_getTagHandle(geom, "GLOBAL_ID", &geom_idtag, &err, strlen("GLOBAL_ID"));
    assert(!err);

    const char *tag2 = "GEOM_DIMENSION";
    iBase_TagHandle geom_dim_tag;
    int namelen = strlen(tag2);
    iMesh_getTagHandle(mesh, tag2, &geom_dim_tag, &err, namelen);
    assert(!err);

    iBase_EntitySetHandle rootSet;
    iMesh_getRootSet(mesh, &rootSet, &err);

    SimpleArray<iBase_EntitySetHandle> entitySets;
    iMesh_getEntSets(mesh, rootSet, 0, ARRAY_INOUT(entitySets), &err);

    SimpleArray<iBase_EntityHandle> mNodes;
    iMesh_getEntities(mesh, rootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                      ARRAY_INOUT(mNodes), &err);

    int gid, gdim;
    double x, y, z;

    ofile << "#Nodes " << endl;
    for (int i = 0; i < mNodes.size(); i++)
    {
        iMesh_getVtxCoord(mesh, mNodes[i], &x, &y, &z, &err);
        iMesh_getIntData(mesh, mNodes[i], mesh_idtag, &gid, &err);
        ofile << gid << " " << fixed << x << " " << y << " " << z << endl;
    }
    ofile << endl;

    int numEdges;
    iMesh_getNumOfType(mesh, rootSet, iBase_EDGE, &numEdges, &err);

    int index = 0;

    int id1, id2;
    SimpleArray<iBase_EntityHandle> mEdges;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mEdges.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mEdges), &err);

        if (mEdges.size() && (mEdges.size() != numEdges))
        {
            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_idtag, &gid, &err);
            assert(!err);

            iMesh_getEntSetIntData(mesh, entitySets[i], geom_dim_tag, &gdim, &err);
            assert(!err);

            ofile << "#EntitySet EDGE : GEOM_ID " << gid << " GEOM_DIMENSION " << gdim << endl;

            SimpleArray<iBase_EntityHandle> edgeNodes;
            for (int i = 0; i < mEdges.size(); i++)
            {
                iMesh_getEntAdj(mesh, mEdges[i], iBase_VERTEX, ARRAY_INOUT(edgeNodes), &err);
                iMesh_getIntData(mesh, edgeNodes[0], mesh_idtag, &id1, &err);
                assert(!err);
                iMesh_getIntData(mesh, edgeNodes[1], mesh_idtag, &id2, &err);
                assert(!err);
                ofile << id1 << " " << id2 << endl;
            }
            ofile << endl;
        }
    }

    int numFaces;
    iMesh_getNumOfType(mesh, rootSet, iBase_FACE, &numFaces, &err);

    int id3, id4;
    SimpleArray<iBase_EntityHandle> mFaces;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mFaces.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mFaces), &err);

        if (mFaces.size() && (mFaces.size() != numFaces))
        {
            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_idtag, &gid, &err);
            assert(!err);

            iMesh_getEntSetIntData(mesh, entitySets[i], geom_dim_tag, &gdim, &err);
            assert(!err);

            ofile << "#EntitySet FACE : GEOM_ID " << gid << " GEOM_DIMENSION " << gdim << endl;

            SimpleArray<iBase_EntityHandle> faceNodes;
            for (int j = 0; j < mFaces.size(); j++)
            {
                faceNodes.clear();
                iMesh_getEntAdj(mesh, mFaces[j], iBase_VERTEX, ARRAY_INOUT(faceNodes), &err);
                for (int k = 0; k < faceNodes.size(); k++)
                {
                    iMesh_getIntData(mesh, faceNodes[k], mesh_idtag, &id1, &err);
                    ofile << id1 << " ";
                }
                ofile << endl;
            }
            ofile << endl;
        }
    }


    int numCells;
    iMesh_getNumOfType(mesh, rootSet, iBase_REGION, &numFaces, &err);

    int id5, id6, id7, id8;
    SimpleArray<iBase_EntityHandle> mCells;
    for (int i = 0; i < entitySets.size(); i++)
    {
        mCells.clear();
        iMesh_getEntities(mesh, entitySets[i], iBase_REGION, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(mCells), &err);

        if (mCells.size() && (mCells.size() != numCells))
        {
            iMesh_getEntSetIntData(mesh, entitySets[i], mesh_idtag, &gid, &err);
            assert(!err);

            ofile << "#EntitySet CELL : GEOM_ID " << gid << endl;

            SimpleArray<iBase_EntityHandle> cellNodes;
            for (int j = 0; j < mCells.size(); j++)
            {
                cellNodes.clear();
                iMesh_getEntAdj(mesh, mCells[j], iBase_VERTEX, ARRAY_INOUT(cellNodes), &err);
                for (int k = 0; k < cellNodes.size(); k++)
                {
                    iMesh_getIntData(mesh, cellNodes[k], mesh_idtag, &id1, &err);
                    ofile << id1 << " ";
                }
                ofile << endl;
            }
            ofile << endl;
        }
    }


}
///////////////////////////////////////////////////////////////////////////////

int readCubitMesh(const char *filename, iMesh_Instance &mesh)
{
    char *options = NULL;
    int optlen = 0;
    int err;

    iMesh_newMesh(options, &mesh, &err, optlen);

    SimpleArray<int> adjTable;
    iMesh_getAdjTable(mesh, ARRAY_INOUT(adjTable), &err);
    assert(!err);
    if (adjTable[5] == 0) adjTable[5] = 1;
    if (adjTable[10] == 0) adjTable[10] = 1;
    iMesh_setAdjTable(mesh, ARRAY_IN(adjTable), &err);
    if (err)
    {
        char descr[1000];
        int len = 1000;
        iMesh_getDescription(mesh, descr, &err, len);
        cout << descr << endl;
        exit(0);
    }

    iBase_EntitySetHandle rootSet;
    iMesh_getRootSet(mesh, &rootSet, &err);

    iMesh_load(mesh, rootSet, filename, options, &err, strlen(filename), optlen);

    if (err)
    {
        char desc[1024];
        iMesh_getDescription(mesh, desc, &err, 1024);
        cout << desc << endl;
        exit(0);
    }

    const char *outfile = "meshcubit.vtk";
    int namelen = strlen(outfile);
    iMesh_save(mesh, rootSet, outfile, NULL, &err, namelen, 0);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void usage()
{
    cout << "Executable : i:o:g:n:" << endl;
    cout << "   Where  i : input  file (must be cubit file *.cub )" << endl;
    cout << "          o : output file (must be *.dat) " << endl;
    cout << "          g : (0-1) Whether model has geometry : Default O " << endl;
    cout << "          n : Order of higher elements ( 3-32 ): Default 8 " << endl;
    cout << "Example hoelem -i mymodel.cub -o horder.dat  -g 1 -n 16 " << endl;
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        usage();
    }
    char *infile = 0,  *outfile = 0;

    int hasGeometry = 0, norder = 8;

    int c;
    while ((c = getopt(argc, argv, "i:o:g:n:")) != -1)
    {
        switch (c)
        {
        case 'i':
            infile = optarg;
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'g':
            hasGeometry = atoi(optarg);
            break;
        case 'n':
            norder = atoi(optarg);
            break;
        default:
            usage();
        }
    }

    if( infile == 0) {
        cout << " Error: No input file provided " << endl;
        return 1;
    }

    if( outfile == 0) {
        cout << " Error: No output file provided " << endl;
        return 1;
    }

    if( norder < 3 ) {
        cout << "Warning: Order less than 3, nothing needs to be done " << endl;
        return 2;
    }

    iGeom_Instance geom;

    if (hasGeometry) readGeometry(infile, geom);

    iMesh_Instance mesh;
    readCubitMesh(infile, mesh);

    iRel_Instance assoc;
    iRel_RelationHandle rel;

    SpectralElements spe;

    if (hasGeometry)
    {
        buildAssociations(geom, mesh, assoc, rel);
        spe.generate(mesh, norder, geom, assoc, rel);
    }
    else
        spe.generate(mesh, norder);

    spe.saveAs(outfile);

}

////////////////////////////////////////////////////////////////////////////////


