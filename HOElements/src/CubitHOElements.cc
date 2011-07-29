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

#include <meshkit/iGeom.hpp>
#include <meshkit/iMesh.hpp>
#include <meshkit/iRel.hpp>

#include "SpectralElements.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

iMesh* readSimpleMesh(const char *filename)
{
     ifstream infile(filename, ios::in);
     if (infile.fail()) {
          cout << "Error: cann't read file " << filename << endl;
          return NULL;
     }

     iMesh *mesh = new iMesh;

     int err;

     iBase_TagHandle idtag;
     err = mesh->getTagHandle("GLOBAL_ID", idtag);
     assert(!err);

     string s;
     infile >> s;
     assert(s == "#Nodes");

     int numNodes;
     infile >> numNodes;

     vector<iBase_EntityHandle> nodeHandles;
     nodeHandles.resize(numNodes);

     double x, y, z;

     iBase_EntityHandle newhandle;
     int id, index = 0;
     for (int i = 0; i < numNodes; i++) {
          infile >> id >> x >> y >> z;
          mesh->createVtx(x, y, z, newhandle);
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
     for (int i = 0; i < numCells; i++) {
          infile >> n0 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7;
          connect[0] = nodeHandles[n0];
          connect[1] = nodeHandles[n1];
          connect[2] = nodeHandles[n3];
          connect[3] = nodeHandles[n2];
          connect[4] = nodeHandles[n4];
          connect[5] = nodeHandles[n5];
          connect[6] = nodeHandles[n7];
          connect[7] = nodeHandles[n6];
          mesh->createEnt(iMesh_HEXAHEDRON, &connect[0], 8, newhandle);
          cellHandles[index++] = newhandle;
     }
     return mesh;
}

///////////////////////////////////////////////////////////////////////////////

iGeom *readGeometry(const char *filename )
{
     ifstream infile(filename, ios::in);
     if (infile.fail()) {
          cout << "Error: cann't read file " << filename << endl;
          return NULL;
     }

     int err;

     iGeom *geom = new iGeom();

     geom->load(filename);

     iBase_EntitySetHandle rootSet = geom->getRootSet();

     cout << "Model Contents " << endl;
//  const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};

     int count;

     geom->getNumOfType(rootSet, iBase_VERTEX, count);
     cout << "#Vertex  : " << count << endl;

     geom->getNumOfType(rootSet, iBase_EDGE, count);
     cout << "#Edges   : " << count << endl;

     geom->getNumOfType(rootSet, iBase_FACE, count);
     cout << "#Faces   : " << count << endl;

     geom->getNumOfType(rootSet, iBase_REGION, count);
     cout << "#Regions : "  << count << endl;

     iBase_TagHandle idtag;
     err = geom->getTagHandle("GLOBAL_ID", idtag );
     assert(!err);

     int gid;
     set<int> idset;

     vector<iBase_EntityHandle> gEdges;
     err = geom->getEntities(rootSet, iBase_EDGE, gEdges);

     int numEdges = gEdges.size();
     idset.clear();
     for (int i = 0; i < numEdges; i++) {
          err = geom->getIntData(gEdges[i], idtag, gid);
          assert( !err );
          idset.insert(gid);
     }

     vector<iBase_EntityHandle> gFaces;
     err = geom->getEntities(rootSet, iBase_FACE, gFaces);

     int numFaces = gFaces.size();
     idset.clear();
     for (int i = 0; i < numFaces; i++) {
          err = geom->getIntData(gFaces[i], idtag, gid);
          assert( !err );
          idset.insert(gid);
//      GFace currface(gFaces[i], geom);
//      ostringstream oss;
//      oss << "gface" << i << ".dat";
//      currface.saveAs(oss.str());
     }

     ofstream ofile( "modeledges.dat", ios::out);

     vector<Point3D> xyz;
     vector<double>  u;

     ofile << "#Nodes " << 100*gEdges.size() << endl;

     int index = 0;
     for (int i = 0; i < gEdges.size(); i++) {
          GEdge curredge(gEdges[i], geom);
          curredge.uniform_u_discretization(100, xyz, u);
          for( int j = 0; j < xyz.size(); j++)
               ofile << index++ << " " << xyz[j][0] << " " << xyz[j][1] << " " << xyz[j][2] << endl;
     }

     return geom;
     cout << " IDSET : " << idset.size() << endl;
     exit(0);
}

///////////////////////////////////////////////////////////////////////////////

void buildAssociations(iGeom *geom, iMesh *mesh, iRel *rel, iRel::PairHandle *&relPair)
{

     int err, namelen;

     vector<iBase_EntitySetHandle> entitySets;
     iBase_EntitySetHandle geom_root_set, mesh_root_set;

     // Get the root sets of the geometry and mesh.
     mesh_root_set = mesh->getRootSet();
     geom_root_set = geom->getRootSet();

     iBase_TagHandle geom_id_tag, mesh_id_tag, geom_dim_tag;
     const char *tag1 = "GLOBAL_ID";
     geom->getTagHandle(tag1, geom_id_tag);
     mesh->getTagHandle(tag1, mesh_id_tag);

     const char *tag2 = "GEOM_DIMENSION";
     mesh->getTagHandle(tag2, geom_dim_tag);
     assert(!err);

     err = rel->createPair(geom, iRel::ENTITY, iRel::IGEOM_IFACE, iRel::ACTIVE,
                           mesh, iRel::SET,    iRel::IMESH_IFACE, iRel::ACTIVE,
                           relPair);
     assert(!err);

     // Get all the entitySet in the mesh
     mesh->getEntSets(mesh_root_set, 0, entitySets);

     int ncount;
     iBase_EntityHandle gEntity;

     // Map all the geometric edges
     vector<iBase_EntityHandle> gEdges;
     geom->getEntities(geom_root_set, iBase_EDGE, gEdges);
     assert(!err);

     int geom_id, geom_dim;
     std::map<int, iBase_EntityHandle> mapEdges;
     for (int i = 0; i < gEdges.size(); i++) {
          geom->getIntData(gEdges[i], geom_id_tag, geom_id);
          mapEdges[geom_id] = gEdges[i];
     }

     // Map all the geometric faces ...
     vector<iBase_EntityHandle> gFaces;
     geom->getEntities(geom_root_set, iBase_FACE, gFaces);
     assert(!err);

     std::map<int, iBase_EntityHandle> mapFaces;
     for (int i = 0; i < gFaces.size(); i++) {
          geom->getIntData(gFaces[i], geom_id_tag, geom_id);
          mapFaces[geom_id] = gFaces[i];
     }

     // Map all the geometric cells ...
     vector<iBase_EntityHandle> gCells;
     geom->getEntities(geom_root_set, iBase_REGION, gCells);
     assert(!err);

     std::map<int, iBase_EntityHandle> mapCells;
     for (int i = 0; i < gCells.size(); i++) {
          geom->getIntData(gCells[i], geom_id_tag, geom_id);
          mapCells[geom_id] = gCells[i];
     }

     ///////////////////////////////////////////////////////////////////////////////
     // Create Edge Assocations:
     ///////////////////////////////////////////////////////////////////////////////
     cout << " Building Edge Associations " << endl;

     int numEdges;
     mesh->getNumOfType(mesh_root_set, iBase_EDGE, numEdges);

     vector<iBase_EntityHandle> mEdges;

     int numAssociations = 0;
     for (int i = 0; i < entitySets.size(); i++) {
          mEdges.clear();
          mesh->getEntities(entitySets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES, mEdges);

          if (mEdges.size() && (mEdges.size() != numEdges)) {
               ncount += mEdges.size();

               mesh->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
               assert(!err);

               mesh->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
               assert(!err);

               gEntity = 0;
               switch (geom_dim) {
               case 1:

                    if (mapEdges.find(geom_id) != mapEdges.end()) {
                         gEntity = mapEdges[geom_id];
                         numAssociations++;
                    } else {
                         cout << "Fatal Error: Geometric Edge not found : " << geom_id << endl;
                         exit(0);
                    }
                    break;
               case 2:
                    if (mapFaces.find(geom_id) != mapFaces.end())
                         gEntity = mapFaces[geom_id];
                    else {
                         cout << "Fatal Error: Geometric Face not found : " << geom_id << endl;
                         exit(0);
                    }
                    break;
               case 3:
                    if (mapCells.find(geom_id) != mapCells.end())
                         gEntity = mapCells[geom_id];
                    else {
                         cout << "Fatal Error: Geometric Cell not found : " << geom_id << endl;
                         exit(0);
                    }
                    break;
               default:
                    cout << "Error: Invalid geometric dimension " << geom_dim << endl;
                    exit(0);
               }

               if (gEntity)
                    err = relPair->setEntSetRelation(gEntity, entitySets[i]);
          }
     }

     if (numAssociations != mapEdges.size())
          cout << "Warning: There are more edge entitySet than geometric edges " << endl;

     //////////////////////////////////////////////////////////////////////////////
     // Face Association
     //////////////////////////////////////////////////////////////////////////////
     cout << " Building Face Associations " << endl;

     vector<iBase_EntityHandle> mFaces;

     int numFaces;
     mesh->getNumOfType(mesh_root_set, iBase_FACE, numFaces);

     int mesh_faceid;

     numAssociations = 0;
     for (int i = 0; i < entitySets.size(); i++) {
          mFaces.clear();
          mesh->getEntities(entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES, mFaces);

          if (mFaces.size() && (mFaces.size() != numFaces)) {
               ncount += mFaces.size();
               err = mesh->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
               assert(!err);
               err = mesh->getEntSetIntData(entitySets[i], geom_dim_tag, geom_dim);
               assert(!err);

               gEntity = 0;
               switch (geom_dim) {
               case 2:
                    if (mapFaces.find(geom_id) != mapFaces.end()) {
                         gEntity = mapFaces[geom_id];
                         numAssociations++;
                    } else {
                         cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                         exit(0);
                    }
                    break;
               case 3:
                    if (mapCells.find(geom_id) != mapCells.end())
                         gEntity = mapCells[geom_id];
                    else {
                         cout << "Fatal Error: Geometric face not found " << geom_id << endl;
                         exit(0);
                    }
                    break;
               }

               if (gEntity)
                    err = relPair->setEntSetRelation(gEntity, entitySets[i]);
          }
     }

     if (numAssociations != mapFaces.size())
          cout << "Warning: There are more face entitySet than geometric faces " << endl;

     //////////////////////////////////////////////////////////////////////////////
     // Cell Association
     //////////////////////////////////////////////////////////////////////////////

     vector<iBase_EntityHandle> mCells;

     int mesh_cellid;

     int numCells;
     mesh->getNumOfType(mesh_root_set, iBase_REGION, numCells);

     ncount = 0;
     for (int i = 0; i < entitySets.size(); i++) {
          mCells.clear();
          mesh->getEntities(entitySets[i], iBase_REGION, iMesh_ALL_TOPOLOGIES, mCells);

          if (mCells.size() && (mCells.size() != numCells)) {
               ncount += mCells.size();
               err = mesh->getEntSetIntData(entitySets[i], mesh_id_tag, geom_id);
               assert(!err);

               if (mapCells.find(geom_id) != mapCells.end()) {
                    if (mapCells.find(geom_id) != mapCells.end()) {
                         gEntity = mapCells[geom_id];
                         err = relPair->setEntSetRelation(gEntity, entitySets[i]);
                    }
               }
          }
     }

}
////////////////////////////////////////////////////////////////////////////////
void saveCubitMesh(iGeom *geom, iMesh *mesh, const char *filename)
{
     ofstream ofile(filename, ios::out);
     if (ofile.fail()) return;

     int err;

     iBase_TagHandle mesh_idtag;
     err = mesh->getTagHandle("GLOBAL_ID", mesh_idtag);
     assert(!err);

     iBase_TagHandle geom_idtag;
     err = geom->getTagHandle("GLOBAL_ID", geom_idtag );
     assert(!err);

     const char *tag2 = "GEOM_DIMENSION";
     iBase_TagHandle geom_dim_tag;
     err = mesh->getTagHandle(tag2, geom_dim_tag);
     assert(!err);

     iBase_EntitySetHandle rootSet = mesh->getRootSet();

     vector<iBase_EntitySetHandle> entitySets;
     mesh->getEntSets(rootSet, 0, entitySets);

     vector<iBase_EntityHandle> mNodes;
     mesh->getEntities(rootSet, iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                       mNodes);
     int gid, gdim;
     double x, y, z;

     ofile << "#Nodes " << endl;
     for (int i = 0; i < mNodes.size(); i++) {
          err = mesh->getVtxCoord(mNodes[i], x, y, z);
          err = mesh->getIntData( mNodes[i], mesh_idtag, gid);
          ofile << gid << " " << fixed << x << " " << y << " " << z << endl;
     }
     ofile << endl;

     int numEdges;
     mesh->getNumOfType(rootSet, iBase_EDGE, numEdges);

     int index = 0;

     int id1, id2;
     vector<iBase_EntityHandle> mEdges;
     for (int i = 0; i < entitySets.size(); i++) {
          mEdges.clear();
          mesh->getEntities(entitySets[i], iBase_EDGE, iMesh_ALL_TOPOLOGIES, mEdges);

          if (mEdges.size() && (mEdges.size() != numEdges)) {
               err = mesh->getEntSetIntData(entitySets[i], mesh_idtag, gid);
               assert(!err);

               err = mesh->getEntSetIntData(entitySets[i], geom_dim_tag, gdim);
               assert(!err);

               ofile << "#EntitySet EDGE : GEOM_ID " << gid << " GEOM_DIMENSION " << gdim << endl;

               vector<iBase_EntityHandle> edgeNodes;
               for (int i = 0; i < mEdges.size(); i++) {
                    mesh->getEntAdj(mEdges[i], iBase_VERTEX, edgeNodes);
                    err = mesh->getIntData(edgeNodes[0], mesh_idtag, id1);
                    assert(!err);
                    err = mesh->getIntData(edgeNodes[1], mesh_idtag, id2);
                    assert(!err);
                    ofile << id1 << " " << id2 << endl;
               }
               ofile << endl;
          }
     }


     int numFaces;
     mesh->getNumOfType(rootSet, iBase_FACE, numFaces);

     int id3, id4;
     vector<iBase_EntityHandle> mFaces;
     for (int i = 0; i < entitySets.size(); i++) {
          mFaces.clear();
          mesh->getEntities(entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES,mFaces);

          if (mFaces.size() && (mFaces.size() != numFaces)) {
               err = mesh->getEntSetIntData(entitySets[i], mesh_idtag, gid);
               assert(!err);

               err = mesh->getEntSetIntData(entitySets[i], geom_dim_tag, gdim);
               assert(!err);

               ofile << "#EntitySet FACE : GEOM_ID " << gid << " GEOM_DIMENSION " << gdim << endl;

               vector<iBase_EntityHandle> faceNodes;
               for (int j = 0; j < mFaces.size(); j++) {
                    faceNodes.clear();
                    mesh->getEntAdj(mFaces[j], iBase_VERTEX, faceNodes);
                    for (int k = 0; k < faceNodes.size(); k++) {
                         mesh->getIntData(faceNodes[k], mesh_idtag, id1);
                         ofile << id1 << " ";
                    }
                    ofile << endl;
               }
               ofile << endl;
          }
     }

     int numCells;
     mesh->getNumOfType(rootSet, iBase_REGION, numCells);

     int id5, id6, id7, id8;
     vector<iBase_EntityHandle> mCells;
     for (int i = 0; i < entitySets.size(); i++) {
          mCells.clear();
          mesh->getEntities(entitySets[i], iBase_REGION, iMesh_ALL_TOPOLOGIES, mCells);

          if (mCells.size() && (mCells.size() != numCells)) {
               mesh->getEntSetIntData(entitySets[i], mesh_idtag, gid);
               assert(!err);

               ofile << "#EntitySet CELL : GEOM_ID " << gid << endl;

               vector<iBase_EntityHandle> cellNodes;
               for (int j = 0; j < mCells.size(); j++) {
                    cellNodes.clear();
                    mesh->getEntAdj(mCells[j], iBase_VERTEX, cellNodes);
                    for (int k = 0; k < cellNodes.size(); k++) {
                         mesh->getIntData(cellNodes[k], mesh_idtag, id1);
                         ofile << id1 << " ";
                    }
                    ofile << endl;
               }
               ofile << endl;
          }
     }
}
///////////////////////////////////////////////////////////////////////////////

iMesh* readCubitMesh(const char *filename)
{
     int err;

     iMesh *mesh = new iMesh;

     vector<int> adjTable;
     /*
         mesh->getAdjTable(adjTable);

         assert(!err);
         if (adjTable[5] == 0) adjTable[5] = 1;
         if (adjTable[10] == 0) adjTable[10] = 1;

         iMesh_setAdjTable(mesh, ARRAY_IN(adjTable), &err);
         if (err)
         {
             char descr[1000];
             int len = 1000;
             iMesh_getDescription(mesh, descr, len);
             cout << descr << endl;
             exit(0);
         }

     */
     iBase_EntitySetHandle rootSet = mesh->getRootSet();
     mesh->load(rootSet, filename);

     const char *outfile = "meshcubit.vtk";
     mesh->save(rootSet, outfile);

     return mesh;
}

///////////////////////////////////////////////////////////////////////////////

void usage()
{
     cout << "Executable : i:o:g:n:" << endl;
     cout << "   Where  i : input  file (must be cubit file *.cub )" << endl;
     cout << "          o : output file (must be *.dat) " << endl;
     cout << "          g : (0-1) Whether model has geometry : Default 0 " << endl;
     cout << "          n : Order of higher elements ( 3-32 ): Default 8 " << endl;
     cout << "Example hoelem -i mymodel.cub -o horder.dat  -g 1 -n 16 " << endl;
     exit(0);
}
///////////////////////////////////////////////////////////////////////////////

#ifdef UNIT_TEST

int main(int argc, char **argv)
{
     if (argc == 1) {
          usage();
     }
     char *infile = 0,  *outfile = 0;

     int hasGeometry = 0, norder = 8;

     int c;
     while ((c = getopt(argc, argv, "i:o:g:n:")) != -1) {
          switch (c) {
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

     iGeom *geom = NULL;
     if (hasGeometry) geom = readGeometry(infile);

     /*
         iMesh *mesh = readCubitMesh(infile);

         iRel rel;
         iRel_PairHandle assoc;

         SpectralElements spe;

         if (hasGeometry)
         {
             buildAssociations(&geom, &mesh, &rel, assoc);
             spe.generate(&mesh, norder, &geom, &rel, assoc);
         }
         else
             spe.generate(&mesh, norder);

         spe.saveAs(outfile);
     */
}

////////////////////////////////////////////////////////////////////////////////

#endif


