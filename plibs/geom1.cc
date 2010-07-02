#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include "SimpleArray.hpp"

using namespace std;

struct GeomMesh
{
   GeomMesh() 
   {
       char  *options = NULL;
       int    optlen  = 0;
       int    err;
       iMesh_newMesh( options, &mesh, &err, optlen);  
       iMesh_getRootSet( mesh, &meshRootSet, &err);
       iRel_newAssoc( 0, &assoc, &err, 0);
       iRel_createAssociation( assoc, 
                          geom, 0, iRel_IGEOM_IFACE,
                          mesh, 2, iRel_IMESH_IFACE, &relation, &err); 
   }

   iMesh_Instance  mesh;
   iGeom_Instance  geom;
   iRel_Instance   assoc;
   iRel_RelationHandle relation;
   iBase_EntitySetHandle meshRootSet, geomRootSet;
};

///////////////////////////////////////////////////////////////////////////////

int readMesh(const string &filename, iMesh_Instance &mesh)
{

   // Read the mesh generated by OpenCascade OCCMeshGenerator and convert into
   // TetGen format which is supported by MOAB. 

   ofstream ofile;
   ifstream ifile( filename.c_str(), ios::in);
   if( ifile.fail() ) {
       cout << "Warning: Cann't open file " << filename << endl;
       return 1;
   }

   char  *options = NULL;
   int    optlen  = 0;
   int    err;
   iMesh_newMesh( options, &mesh, &err, optlen);  

   string str;
   ifile >> str;   assert( str == "mesh3d");

   ifile >> str;   assert( str == "dimension");

   int dimension;
   ifile >> dimension; assert( dimension == 3);
   
   ifile >> str;   assert( str == "geomtype");
   int geomtype;
   ifile >> geomtype;  assert( geomtype == 12);

   ifile >> str;
   char line[256];
   if( str[0] == '#') {
       ifile.getline( line, 256 );
   }

   ifile >> str; assert( str == "surfaceelementsuv");

   int numFaces;

   ifile >> numFaces;

   int surfnr, bcnr, domin, domout, np, vid;
   double u,v;

   ofile.open( "model.face", ios::out);
   ofile << numFaces <<  " 1 " <<  endl;

   map<int,int>  geomface;
   for( int i = 0; i < numFaces; i++) 
   {
        ifile >> surfnr >> bcnr >> domin >> domout >> np;
        surfnr = surfnr - 1;  // Problem in NetGen, but developers not ready to modify.
        ofile << i << " ";
        for( int j = 0; j < np; j++) {
            ifile >> vid;
            ofile << vid - 1 << " ";
        }
        ofile << surfnr << endl;
        geomface[i] = surfnr;
      
        for( int j = 0; j < np; j++) 
             ifile >> u >> v;
   }
   ofile.close();

   ifile >> str;
   if( str[0] = '#') {
       ifile.getline( line, 256 );
   }
   ifile >> str; assert( str == "volumeelements");
   
   int numCells;
   ifile >> numCells;

   ofile.open( "model.ele", ios::out);
   ofile << numCells <<  " 4  1 " <<  endl;

   int matnr;
   for( int i = 0; i < numCells; i++) {
        ifile >> matnr >> np;
        ofile << i << " ";
        for( int j = 0; j < np; j++)  {
             ifile >> vid;
             ofile << vid-1 << " ";
        }
        ofile << matnr << endl;
   }
   ofile.close();

   ifile >> str;
   if( str[0] = '#') {
       ifile.getline( line, 256 );
   }
   ifile >> str; assert( str == "edgesegmentsgi2");

   int numEdges;
   ifile >> numEdges;

   int surfid, dummy , p1, p2, trignum1, trignum2, ednr1, ednr2;
   double dist1, dist2;

   ofile.open( "model.edge", ios::out);
   ofile << numEdges <<  endl;

   for( int i = 0; i < numEdges; i++) 
   {
        ifile >> surfid >> dummy >> p1 >> p2 >> trignum1 >> trignum2 >> domin >> domout 
              >> ednr1  >> dist1 >> ednr2 >> dist2;
        ofile << i << " " << p1-1 << " " << p2-1 << endl;
   }
   ofile.close();

   ifile >> str;
   if( str[0] = '#') ifile.getline( line, 256 );

   ifile >> str; assert( str == "points");

   int numNodes;
   ifile >> numNodes;

   double  x, y, z;
   ofile.open( "model.node", ios::out);
   ofile << numNodes <<  " 3  0  0 " <<  endl;

   for( int i = 0; i <  numNodes; i++) 
   {
        ifile >> x >> y >> z;
        ofile << i << " " << x << " " << y << " " << z << endl;
   }
   ofile.close();

   iBase_EntitySetHandle rootSet;
   iMesh_getRootSet( mesh, &rootSet, &err);

   const char *infile = "model.ele";
   int  namelen = strlen( infile );

   iMesh_load( mesh, rootSet, infile , options, &err, namelen, optlen);

   if( err ) {
       char desc[1024];
       iMesh_getDescription( mesh, desc, &err, 1024 );
       cout << desc << endl;
       exit(0);
   }

  // Presently, the TetGen reader in MOAB created EntitySets, but doesn't
  // assign any Tag to the entitySet and it seems that the surfaceID ( attribute
  // on face in tetgen surface) are not stored. Therefore, we need to do
  // all these modification outside the reader.

  iBase_TagHandle  id_tag, dim_tag;
  const char *tag1 = "GLOBAL_ID";
  namelen  = strlen( tag1 );
  iMesh_getTagHandle(mesh, tag1,  &id_tag, &err, namelen );

  const char *tag2 = "GEOM_DIMENSION";
  namelen  = strlen( tag2 );
  iMesh_getTagHandle(mesh, tag2,  &dim_tag, &err, namelen ); assert( !err );

  SimpleArray<iBase_EntitySetHandle>  entitySets;
  iMesh_getEntSets( mesh, rootSet, 0, ARRAY_INOUT(entitySets), &err);

  //
  // Assign Geometric face to each entitySet in mesh. In order to do that, get any
  // face in an entitySet, and assign the geometric face of the mesh face to the
  // entitySet. 
  //
  // One problem( or feature ) of the iMesh_getEntSets is that it returns 
  // parent entitySet also. It must be excluded for the classifications, therefore,
  // I have checked "entitySize == numFaces" to determine which one is parent
  // set.
  //
  int gid, faceid;
  SimpleArray<iBase_EntityHandle> facehandles;

  iMesh_getNumOfType( mesh, rootSet, iBase_FACE, &numFaces, &err);

  int ncount = 0;
  for( int i = 0; i < entitySets.size(); i++) {
       facehandles.clear();
       iMesh_getEntities( mesh, entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(facehandles), &err);
       if( facehandles.size() && (facehandles.size() != numFaces )) 
       {
           ncount += facehandles.size();
           iMesh_getIntData( mesh, facehandles[0], id_tag, &gid, &err);
           faceid = geomface[gid]; // To which geometry face the mesh face belong.
           iMesh_setEntSetIntData( mesh, entitySets[i], id_tag,  faceid, &err);
           iMesh_setEntSetIntData( mesh, entitySets[i], dim_tag, 2, &err);
       }
  }
  assert( ncount == numFaces );

/*
  // Not necessary, only for visualization purpose.
  const char *outfile = "model.vtk";
  namelen = strlen( outfile );
  iMesh_save( mesh, rootSet, outfile , options, &err, namelen, optlen);
*/

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void buildAssociations( iGeom_Instance &geom, iMesh_Instance &mesh, iRel_Instance &assoc)
{
   int err, namelen;
   SimpleArray<iBase_EntitySetHandle>  entitySets;
   SimpleArray<iBase_EntityHandle>     facehandles;
   iBase_EntitySetHandle geom_root_set, mesh_root_set;

   // Get the root sets of the geometry and mesh.
   iMesh_getRootSet( mesh, &mesh_root_set, &err);
   iGeom_getRootSet( geom, &geom_root_set, &err);

   // Get the geometric faces and mesh faces from the root sets.
   SimpleArray<iBase_EntityHandle> geomfaces, meshfaces;
   iGeom_getEntities( geom, geom_root_set, iBase_FACE, ARRAY_INOUT(geomfaces), &err);
   iMesh_getEntities( mesh, mesh_root_set, iBase_FACE, iMesh_ALL_TOPOLOGIES, ARRAY_INOUT(meshfaces), &err);

   iBase_TagHandle  geom_id_tag, mesh_id_tag;
   const char *tag1 = "GLOBAL_ID";
   namelen  = strlen( tag1 );
   iGeom_getTagHandle( geom, tag1,  &geom_id_tag, &err, namelen );
   iMesh_getTagHandle( mesh, tag1,  &mesh_id_tag, &err, namelen );

   int mesh_faceid, geom_faceid;
   std::map<int, iBase_EntityHandle> mapfaces;

   for( int i = 0; i < geomfaces.size(); i++) {
        iGeom_getIntData( geom, geomfaces[i], geom_id_tag, &geom_faceid, &err);
        mapfaces[geom_faceid] = geomfaces[i];
   }

   // Each geometric faces associates with a set of mesh faces.
   iRel_RelationHandle rel1;
   iRel_createAssociation( assoc, 
                           geom, 0, iRel_IGEOM_IFACE, 
                           mesh, 2, iRel_IMESH_IFACE, &rel1, &err); assert( !err );

  iMesh_getEntSets( mesh, mesh_root_set, 0, ARRAY_INOUT(entitySets), &err);

  int numFaces; 
  iMesh_getNumOfType( mesh, mesh_root_set, iBase_FACE, &numFaces, &err);

  iBase_EntityHandle  geoentity;
  int ncount = 0;
  for( int i = 0; i < entitySets.size(); i++) 
  {
       facehandles.clear();
       iMesh_getEntities( mesh, entitySets[i], iBase_FACE, iMesh_ALL_TOPOLOGIES,
                          ARRAY_INOUT(facehandles), &err);

       if( facehandles.size()  && ( facehandles.size() != numFaces )) 
       {
           ncount += facehandles.size();

           iMesh_getEntSetIntData( mesh, entitySets[i], mesh_id_tag, &geom_faceid, &err); assert( !err );
           geoentity  = mapfaces[geom_faceid];
           iRel_setEntSetAssociation( assoc, rel1, geoentity, entitySets[i], &err );
       }
  }
}

////////////////////////////////////////////////////////////////////////////////

iBase_EntitySetHandle  discretize_close_edge( GeomMesh &geomesh, iBase_EntityHandle gEdge)
{
  int err;
  double x, y, z, u0, u1;

  iMesh_Instance mesh = geomesh.mesh;
  iGeom_Instance geom = geomesh.geom;

  SimpleArray<iBase_EntityHandle>  gNodes;
  iGeom_getEntAdj( geom, gEdge, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

  iGeom_getVtxCoord( geom, gNodes[0], &x, &y, &z, &err);
  iGeom_getEntXYZtoU( geom, gNodes[0], x, y, z, &u0, &err);

  iBase_EntitySetHandle  entitySet;
  iMesh_createEntSet( mesh, 1, &entitySet, &err);

  return entitySet;
}

////////////////////////////////////////////////////////////////////////////////

iBase_EntitySetHandle discretize_open_edge( GeomMesh &geomesh, iBase_EntityHandle gEdge)
{
  int err, status, numEdges = 10;
  iBase_EntityHandle newHandle;
  double x, y, z, u, du, umin, umax;

  iMesh_Instance mesh = geomesh.mesh;
  iGeom_Instance geom = geomesh.geom;

  iGeom_getEntURange( geom, gEdge, &umin, &umax, &err);

  du = (umax-umin)/(double)numEdges;

  vector<iBase_EntityHandle>  nodeHandles;
  nodeHandles.resize(numEdges+1);

  for( int i = 0; i < numEdges+1; i++) {
       u = umin + i*du;
       iGeom_getEntUtoXYZ( geom, gEdge, u, &x, &y, &z, &err);
       iMesh_createVtx(mesh, x, y, z, &newHandle, &err);
       nodeHandles[i] = newHandle;
  }
       
  vector<iBase_EntityHandle> edgeHandles, eConnect(2);
  edgeHandles.resize(numEdges);

  for( int i = 0; i < numEdges; i++) {
       eConnect[0] =  nodeHandles[i];
       eConnect[1] =  nodeHandles[i+1];
       iMesh_createEnt( mesh, iMesh_LINE_SEGMENT, &eConnect[0], 2, &newHandle, &status, &err);
       edgeHandles[i] = newHandle;
  }

  iBase_EntitySetHandle  entitySet;
  iMesh_createEntSet( mesh, 1, &entitySet, &err);
  iMesh_addEntArrToSet( mesh, &edgeHandles[0], numEdges, entitySet, &err);

  return entitySet;
}

////////////////////////////////////////////////////////////////////////////////

int discretize_geom_edge( GeomMesh &geomesh, iBase_EntityHandle gEdge)
{
  int  err;

  iMesh_Instance mesh = geomesh.mesh;
  iGeom_Instance geom = geomesh.geom;

  SimpleArray<iBase_EntityHandle>  gNodes;
  iGeom_getEntAdj( geom, gEdge, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

  iBase_EntitySetHandle  entitySet;

  switch( gNodes.size() )
  {
     case 1:
         entitySet = discretize_close_edge( geomesh, gEdge);
         break;
     case 2:
         entitySet = discretize_open_edge( geomesh, gEdge);
         break;
      default:
         cout << "Fatal Error: Invalid Geometric edge " << endl;
         exit(0);
   }

   iBase_TagHandle  geom_id_tag, mesh_id_tag, dim_tag;
   const char *tag1 = "GLOBAL_ID";
   int namelen  = strlen( tag1 );
   iGeom_getTagHandle( geom, tag1,  &geom_id_tag, &err, namelen );

   int edgeID;
   iGeom_getIntData( geom, gEdge, geom_id_tag, &edgeID, &err);

   iMesh_getTagHandle( mesh, tag1, &mesh_id_tag, &err, namelen );
   iMesh_setEntSetIntData( mesh, entitySet, mesh_id_tag, edgeID, &err);

   const char *tag2 = "GEOM_DIMENSION";
   namelen  = strlen( tag2 );
   iMesh_getTagHandle(mesh, tag2,  &dim_tag, &err, namelen ); assert( !err );

   int dim = 1;
   iMesh_setEntSetIntData( mesh, entitySet, dim_tag, dim, &err);

   iRel_Instance assoc = geomesh.assoc;
   iRel_RelationHandle relation = geomesh.relation;
   iRel_setEntSetAssociation( assoc, relation, gEdge, entitySet, &err);
}

////////////////////////////////////////////////////////////////////////////////

iBase_EntitySetHandle  
generate_surface_mesh(iGeom_Instance geom, vector<iBase_EntityHandle> &boundSegments )
{


}

////////////////////////////////////////////////////////////////////////////////

int discretize_geom_face( GeomMesh &geomesh, iBase_EntityHandle gFace)
{
  int  err;

  iMesh_Instance mesh = geomesh.mesh;
  iGeom_Instance geom = geomesh.geom;
  iRel_Instance assoc = geomesh.assoc;
  iRel_RelationHandle relation = geomesh.relation;

  SimpleArray<iBase_EntityHandle>  gEdges;
  iGeom_getEntAdj( geom, gFace, iBase_VERTEX, ARRAY_INOUT(gEdges), &err);

  vector<iBase_EntitySetHandle> edgeMeshSets;
  for(int i = 0; i < gEdges.size(); i++) 
      iRel_getEntSetAssociation(assoc, relation, gEdges[i], 0, &edgeMeshSets[i], &err);

  vector<iBase_EntityHandle> boundSegment;
  iBase_EntitySet entitySet = generate_surface_mesh(geom, boundSegments );

  iBase_TagHandle  geom_id_tag, mesh_id_tag, dim_tag;
  const char *tag1 = "GLOBAL_ID";
  int namelen  = strlen( tag1 );
  iGeom_getTagHandle( geom, tag1,  &geom_id_tag, &err, namelen );

  int faceID;
  iGeom_getIntData( geom, gFace, geom_id_tag, &edgeID, &err);

  iMesh_getTagHandle( mesh, tag1, &mesh_id_tag, &err, namelen );
  iMesh_setEntSetIntData( mesh, entitySet, mesh_id_tag, faceID, &err);

  const char *tag2 = "GEOM_DIMENSION";
  namelen  = strlen( tag2 );
  iMesh_getTagHandle(mesh, tag2,  &dim_tag, &err, namelen ); assert( !err );

  int dim = 2;
  iMesh_setEntSetIntData( mesh, entitySet, dim_tag, dim, &err);

  iRel_Instance assoc = geomesh.assoc;
  iRel_RelationHandle relation = geomesh.relation;
  iRel_setEntSetAssociation( assoc, relation, gFace, entitySet, &err);

  for(int i = 0; i < gEdges.size(); i++) 
      iMesh_addPrntChld( mesh, entitySet, edgeMeshSets[i], &err);
}

int discretize_geom_cell( GeomMesh &geomesh, iBase_EntityHandle gCell)
{
  int  err;

  iMesh_Instance mesh = geomesh.mesh;
  iGeom_Instance geom = geomesh.geom;
  iRel_Instance assoc = geomesh.assoc;
  iRel_RelationHandle relation = geomesh.relation;

  SimpleArray<iBase_EntityHandle>  gFaces;
  iGeom_getEntAdj( geom, gCell, iBase_FACE, ARRAY_INOUT(gFaces), &err);

  vector<iBase_EntitySetHandle> faceMeshSets;

  for(int i = 0; i < gFaces.size(); i++) 
      iRel_getEntSetAssociation(assoc, relation, gFaces[i], 0, &faceMeshSets[i], &err);

  vector<iBase_EntityHandle> boundFaces;
  iBase_EntitySet entitySet = generate_volume_mesh(geom, boundFaces);

  iBase_TagHandle  geom_id_tag, mesh_id_tag, dim_tag;
  const char *tag1 = "GLOBAL_ID";
  int namelen  = strlen( tag1 );
  iGeom_getTagHandle( geom, tag1,  &geom_id_tag, &err, namelen );

  int faceID;
  iGeom_getIntData( geom, gFace, geom_id_tag, &edgeID, &err);

  iMesh_getTagHandle( mesh, tag1, &mesh_id_tag, &err, namelen );
  iMesh_setEntSetIntData( mesh, entitySet, mesh_id_tag, faceID, &err);

  const char *tag2 = "GEOM_DIMENSION";
  namelen  = strlen( tag2 );
  iMesh_getTagHandle(mesh, tag2,  &dim_tag, &err, namelen ); assert( !err );

  int dim = 3;
  iMesh_setEntSetIntData( mesh, entitySet, dim_tag, dim, &err);

  iRel_Instance assoc = geomesh.assoc;
  iRel_RelationHandle relation = geomesh.relation;
  iRel_setEntSetAssociation( assoc, relation, gCell, entitySet, &err);

  for(int i = 0; i < gFaces.size(); i++) 
      iMesh_addPrntChld( mesh, entitySet, faceMeshSets[i], &err);
}


////////////////////////////////////////////////////////////////////////////////

int readGeometry( string &filename, iGeom_Instance &geom )
{
   string engine_opt = ";engine=OCC";

   int err;
   iGeom_newGeom( engine_opt.c_str(), &geom, &err, engine_opt.length() );
   iGeom_load( geom, &filename[0], 0, &err, filename.length(), 0 );

   iBase_EntitySetHandle rootSet;
   iGeom_getRootSet( geom, &rootSet, &err);

   cout << "Model Contents " << endl;
   const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};

   int count;
   for (int i = 0; i <= 3; ++i) {
    iGeom_getNumOfType( geom, rootSet, i, &count, &err );
    std::cout << gtype[i] << count << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////

iMesh_Instance Discretize_Geometric_Model( GeomMesh &geomesh)
{
  char  *options = NULL;
  int    optlen  = 0;
  int    err;

  iMesh_Instance mesh = geomesh.mesh;
  iGeom_Instance geom = geomesh.geom;

  iBase_EntitySetHandle rootSet = geomesh.meshRootSet;

  ///////////////////////////////////////////////////////////////////////////////
  // Step I: Collect all the corner nodes in Geometry:
  ///////////////////////////////////////////////////////////////////////////////
  SimpleArray<iBase_EntityHandle>  geoNodes;
  iGeom_getEntities( geom, rootSet, iBase_VERTEX, ARRAY_INOUT(geoNodes), &err);

  double x, y, z;
  iBase_EntityHandle newHandle;
  for( int i = 0; i < geoNodes.size(); i++) {
       iGeom_getVtxCoord(geom, geoNodes[i], &x, &y, &z, &err);
       iMesh_createVtx(mesh, x, y, z, &newHandle, &err);
 }

  ///////////////////////////////////////////////////////////////////////////////
  // Step II: Collect all the geometric edges and discretize them. 
  ///////////////////////////////////////////////////////////////////////////////
  SimpleArray<iBase_EntityHandle>  geoEdges;
  iGeom_getEntities( geom, rootSet, iBase_EDGE, ARRAY_INOUT(geoEdges), &err);

  for( int i = 0; i < geoEdges.size(); i++)  
       discretize_geom_edge(geomesh, geoEdges[i]);

 ///////////////////////////////////////////////////////////////////////////////
 //Step III: Collect all the geometric faces and discretize them.
 ///////////////////////////////////////////////////////////////////////////////
  SimpleArray<iBase_EntityHandle>  geoFaces;
  iGeom_getEntities( geom, rootSet, iBase_FACE, ARRAY_INOUT(geoFaces), &err);

  for( int i = 0; i < geoFaces.size(); i++)  
       discretize_geom_face( geomesh, geoFaces[i]);

 return mesh;
}

////////////////////////////////////////////////////////////////////////////////


int main( int argc, char **argv)
{
  int err;
/*
   int err, result;
   if( argc != 3 ) {
       cout << "Usage: executable geomfile (*.brep)  meshfile (*.vol) " << endl;
       return 1;
   }
*/

   iGeom_Instance geom;
   string geomfile =  argv[1];
   readGeometry( geomfile, geom );

/*
   GeomMesh geomesh;
   geomesh.geom = geom;

   iMesh_Instance mesh = Discretize_Geometric_Model( geomesh );

   const char *outfile = "model.vtk";
   int  namelen = strlen( outfile );

   const char *options = NULL;
   int  optlen = 0;

   iBase_EntitySetHandle rootSet;
   iMesh_getRootSet( mesh, &rootSet, &err);
   iMesh_save( mesh, rootSet, outfile , options, &err, namelen, optlen);
*/


   iMesh_Instance mesh;
   string meshfile = argv[1];

   readMesh( meshfile, mesh  );

   const char *outfile = "model.vtk";
   int  namelen = strlen( outfile );

   const char *options = NULL;
   int  optlen = 0;

   iBase_EntitySetHandle rootSet;
   iMesh_getRootSet( mesh, &rootSet, &err);
   iMesh_save( mesh, rootSet, outfile , options, &err, namelen, optlen);

   exit(0);


   iRel_Instance assoc;
   iRel_newAssoc( 0, &assoc, &result, 0);

   buildAssociations( geom, mesh, assoc);
}

////////////////////////////////////////////////////////////////////////////////


