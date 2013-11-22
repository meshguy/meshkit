#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>
#include <iomanip>
#include <limits>
#include <sstream>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <iGeom.h>
#include <iMesh.h>
#include <iRel.h>

#include "SimpleArray.hpp"

using namespace std;

#include <gmsh/Gmsh.h>
#include <gmsh/GModel.h>
#include <gmsh/meshGEdge.h>
#include <gmsh/meshGFace.h>

#include "ITAPVertex.h"
#include "ITAPEdge.h"
#include "ITAPFace.h"
#include "ITAPRegion.h"

map<MVertex*, int> edgeVertexMap;
int  vertexID = 0;

double  minchar_len = 0.02;

double maxErrorCurveEval = 0.0, maxErrorFaceEval = 0.0;

iBase_TagHandle  geom_id_tag;
//////////////////////////////////////////////////////////////////////////////

double edgelength( GEdge *edge)
{
   Range<double> urange = edge->parBounds(0);
   double l = edge->length( urange.low(), urange.high(), 10);
   return l;
}

void print_edge( GEdge *edge)
{
   Range<double> urange = edge->parBounds(0);

   cout << " Edge " << edge->tag()  
        << " End "  << edge->getBeginVertex()->tag() << " "
                    << edge->getEndVertex()->tag()  
        << " U "    << urange.low() << " " 
                    << urange.high() 
        << " Length  " << edgelength(edge) 
        << " MElements " << edge->getNumMeshElements() << endl;
}

int readGeometry(const string &filename, iGeom_Instance &geom)
{
    string engine_opt = ";engine=OCC";

    int err;
//  iGeom_newGeom(engine_opt.c_str(), &geom, &err, engine_opt.length());
    iGeom_newGeom( NULL, &geom, &err, 0);
    iGeom_load(geom, &filename[0], 0, &err, filename.length(), 0);

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

   const char *tag = "GLOBAL_ID";
   int namelen  = strlen( tag );
   iGeom_getTagHandle( geom, tag,  &geom_id_tag, &err, namelen );

}

//////////////////////////////////////////////////////////////////////////////
GModel *build_gmsh_model( string &filename )
{
  iGeom_Instance geom;
  readGeometry( filename, geom);

  int err;
  GModel *model = new GModel;

  iBase_EntitySetHandle rootSet;
  iGeom_getRootSet(geom, &rootSet, &err);

  /////////////////////////////////////////////////////////////////////////////
  //Geometric Nodes
  /////////////////////////////////////////////////////////////////////////////
  map<iBase_EntityHandle, GVertex*>  vertexMap;

  SimpleArray<iBase_EntityHandle> gNodes;
  iGeom_getEntities(geom, rootSet, iBase_VERTEX, ARRAY_INOUT(gNodes), &err);

  double dmax = 0.99*numeric_limits<double>::max();
  double xmin = dmax;
  double ymin = dmax;
  double zmin = dmax;

  double xmax = -1.0*dmax;
  double ymax = -1.0*dmax;
  double zmax = -1.0*dmax;

  int numNodes = gNodes.size();

  for( int i = 0; i < gNodes.size(); i++) 
  {
     GVertex *vertex = new ITAPVertex( model, geom, &gNodes[i] ); assert( vertex );
     vertexMap[gNodes[i]] = vertex;
     model->add( vertex );

     xmin = min( xmin, vertex->x() );
     ymin = min( ymin, vertex->y() );
     zmin = min( zmin, vertex->z() );

     xmax = max( xmax, vertex->x() );
     ymax = max( ymax, vertex->y() );
     zmax = max( zmax, vertex->z() );
  }

  double dx = fabs(xmax-xmin);
  double dy = fabs(ymax-ymin);
  double dz = fabs(zmax-zmin);
  double dlmax  = minchar_len*max( dx, max(dy,dz) );
  GmshSetOption("Mesh", "CharacteristicLengthMin", dlmax);

  /////////////////////////////////////////////////////////////////////////////
  // Add Geometric Edges
  /////////////////////////////////////////////////////////////////////////////
  SimpleArray<iBase_EntityHandle> gEdges, edgeNodes;
  iGeom_getEntities(geom, rootSet, iBase_EDGE, ARRAY_INOUT(gEdges), &err);

  int numEdges = gEdges.size();
  vector<double>  elen(numEdges);
 
  GVertex *v1, *v2;
  for( int i = 0; i < numEdges ; i++) 
  {
     edgeNodes.clear();
     iGeom_getEntAdj(geom, gEdges[i], iBase_VERTEX, ARRAY_INOUT(edgeNodes), &err);

     v1 = NULL;
     v2 = NULL;
     if( edgeNodes.size() == 1 ) {
         assert( vertexMap.find( edgeNodes[0] ) != vertexMap.end() );
         v1 = vertexMap[ edgeNodes[0] ];
         v2 = v1;
     }
 
     if( edgeNodes.size() == 2 ) {
         assert( vertexMap.find( edgeNodes[0] ) != vertexMap.end() );
         v1 = vertexMap[ edgeNodes[0] ];
         assert( vertexMap.find( edgeNodes[1] ) != vertexMap.end() );
         v2 = vertexMap[ edgeNodes[1] ];
     }
 
     ITAPEdge *e = new ITAPEdge( model, geom, &gEdges[i], v1, v2);
     model->add( e );
  }

  /////////////////////////////////////////////////////////////////////////////
  // Add Geometric Faces
  /////////////////////////////////////////////////////////////////////////////

  SimpleArray<iBase_EntityHandle> gFaces;
  iGeom_getEntities(geom, rootSet, iBase_FACE, ARRAY_INOUT(gFaces), &err);

  int numFaces = gFaces.size();
  for( int i = 0; i < numFaces; i++) 
  {
     ITAPFace *face = new ITAPFace(model, geom, &gFaces[i]); assert(face);
     model->add( face );
  }

  /////////////////////////////////////////////////////////////////////////////
  // Add Geometric Region
  /////////////////////////////////////////////////////////////////////////////
  SimpleArray<iBase_EntityHandle> gCells;
  iGeom_getEntities(geom, rootSet, iBase_REGION, ARRAY_INOUT(gCells), &err);

  int numCells = gCells.size();
  for( int i = 0; i < gCells.size(); i++) 
  {
     GRegion *r = new ITAPRegion( model, geom, &gCells[i]);
     model->add( r );
  }

  ////////////////////////////////////////////////////////////////////////////
  // Verify Model
  ////////////////////////////////////////////////////////////////////////////
  cout << " Verifying Model : " << endl;
  assert( model->getNumVertices() ==  numNodes );
  assert( model->getNumEdges()    ==  numEdges );
  assert( model->getNumFaces()    ==  numFaces );
  assert( model->getNumRegions()  ==  numCells );

  SimpleArray<iBase_EntityHandle> faceNodes, faceEdges;
  set<int> aSet, bSet;
  for( int i = 0; i < numFaces; i++) 
  {
     int faceID;
     iGeom_getIntData(geom, gFaces[i], geom_id_tag, &faceID, &err);
     GFace *face = model->getFaceByTag(faceID); assert( face );

     aSet.clear();
     faceNodes.clear();
     iGeom_getEntAdj(geom, gFaces[i], iBase_VERTEX, ARRAY_INOUT(faceNodes), &err);
     for( int j = 0; j < faceNodes.size(); j++)  {
          iGeom_getIntData(geom, faceNodes[j], geom_id_tag, &vertexID, &err);
          aSet.insert( vertexID );
     }

     bSet.clear();
     list<GVertex*> vertices = face->vertices();
     GVertex *vertex;
     BOOST_FOREACH( vertex, vertices) bSet.insert( vertex->tag() );

     if( aSet != bSet ) {
         cout << " Fatal Error: Vertex Set not matched : Face " << faceID << endl;
         set<int>::const_iterator sit;
         cout << " A Set : ";
         for( sit = aSet.begin(); sit != aSet.end(); ++sit) 
              cout << *sit << " ";
         cout << endl;

         cout << " B Set : ";
         for( sit = bSet.begin(); sit != bSet.end(); ++sit) 
              cout << *sit << " ";
         cout << endl;
         exit(0);
     }

     aSet.clear();
     faceEdges.clear();
     iGeom_getEntAdj(geom, gFaces[i], iBase_EDGE, ARRAY_INOUT(faceEdges), &err);
     int edgeID;
     for( int j = 0; j < faceEdges.size(); j++)  {
          iGeom_getIntData(geom, faceEdges[j], geom_id_tag, &edgeID, &err);
          aSet.insert( edgeID );
     }

     bSet.clear();
     list<GEdge*> fedge = face->edges();
     GEdge* edge;
     BOOST_FOREACH( edge, fedge) bSet.insert( edge->tag() );

     if( aSet != bSet ) {
         cout << " Fatal Error: Edge Set not matched : Face " << faceID << endl;
         set<int>::const_iterator sit;
         cout << " A Set : ";
         for( sit = aSet.begin(); sit != aSet.end(); ++sit) 
              cout << *sit << " ";
         cout << endl;

         cout << " B Set : ";
         for( sit = bSet.begin(); sit != bSet.end(); ++sit) 
              cout << *sit << " ";
         cout << endl;
         exit(0);
     }
   }

  ////////////////////////////////////////////////////////////////////////////
  // Edge Meshing ...
  ////////////////////////////////////////////////////////////////////////////
  cout << "Edge Meshing Starts .... " << endl;

  meshGEdge meshedge;
  for( int i = 0; i < numEdges; i++) 
  {
       int edgeID = i + 1;
       GEdge *edge = model->getEdgeByTag(edgeID); assert( edge );
       meshedge(edge);
//     save(edge);
  }
  model->writeVTK( "modeledges.vtk");
  cout << "Edge Meshing Completed .... " << endl;
 
  /////////////////////////////////////////////////////////////////////////////
  // Surface Meshing ...
  /////////////////////////////////////////////////////////////////////////////
  
  cout << "Surface Meshing Starts : " << endl;

//CSV: GMSH
  meshGFace meshface;
  for( int i = 0; i < numFaces; i++) 
  {
     int faceID = i+ 1;
     ITAPFace *face = (ITAPFace *)model->getFaceByTag(faceID); assert( face );
     save(face);
     face->create_kdtree();
     meshface(face);
     cout <<  "Face " << faceID << " Mesh Size : " << face->getNumMeshElements() << endl;
     face->delete_kdtree();
  }
  model->writeVTK( "modelfaces.vtk");
  cout << "Surface Meshing Completed : " << endl;

  exit(0);
  return model;

  /////////////////////////////////////////////////////////////////////////////
  // Volume Meshing ..
  /////////////////////////////////////////////////////////////////////////////

  return model;
}

///////////////////////////////////////////////////////////////////////////////

void save(GEdge *edge, iGeom_Instance &geom, iBase_EntityHandle edgeHandle)
{
  int edgeID = edge->tag();

  ostringstream oss;
  oss <<  "ModelMesh/meshedge" << edgeID << ".dat";

  string filename = oss.str();

  ofstream ofile( filename.c_str(), ios::out);

  int numMeshElements = edge->getNumMeshElements();

  set<MVertex*>  vset;
  for( int i = 0; i < numMeshElements; i++) {
       MElement *melem  = edge->getMeshElement(i);
       for( int j = 0; j < melem->getNumVertices(); j++) {
            MVertex *vertex = melem->getVertex(j);
            if( edgeVertexMap.find( vertex ) == edgeVertexMap.end() ) 
                edgeVertexMap[vertex] = vertexID++;
            vset.insert( vertex );
       }
  }

  ofile << "EdgeID  " << edgeID << endl;
  ofile << vset.size() <<  "  " << numMeshElements << endl;

  int err;

  double umin, umax;
  iGeom_getEntURange( geom, edgeHandle, &umin, &umax, &err);

  Range<double> ur = edge->parBounds(0);

  double x, y, z, u, u1, maxerror = 0.0;
  set<MVertex*>::const_iterator it;
  for( it = vset.begin(); it != vset.end(); ++it) 
  {
       MVertex *vertex = *it;
       int vid = edgeVertexMap[vertex];

       x = vertex->x();
       y = vertex->y();
       z = vertex->z();
       ofile <<  vid  << "  " << x << " " << y << " " << z << endl;
      
       reparamMeshVertexOnEdge( vertex, edge, u );

       iGeom_getEntUtoXYZ( geom, edgeHandle, u, &x, &y, &z, &err); assert( !err );

       double dx = fabs(vertex->x() - x );
       double dy = fabs(vertex->y() - y );
       double dz = fabs(vertex->z() - z );
       maxerror  = max( maxerror, max(dx, max(dy,dz)));
  }

  maxErrorCurveEval = max( maxerror, maxErrorCurveEval );

  for( int i = 0; i < numMeshElements; i++) {
       MElement *melem  = edge->getMeshElement(i);
       for( int j = 0; j < melem->getNumVertices(); j++) {
            MVertex *vertex = melem->getVertex(j);
            ofile << edgeVertexMap[vertex] << " ";
       }
       ofile << endl;
 }

}

////////////////////////////////////////////////////////////////////////////////

void save(GFace *face, iGeom_Instance &geom, iBase_EntityHandle faceHandle)
{
  int err;
  int faceID = face->tag();

  ostringstream oss;
  oss <<  "ModelMesh/meshface" << faceID << ".dat";

  string filename = oss.str();

  ofstream ofile( filename.c_str(), ios::out);

  Range<double> urange = face->parBounds(0);
  Range<double> vrange = face->parBounds(1);

  double eps = 1.0E-06;
  double umin, umax, vmin, vmax;
  iGeom_getEntUVRange( geom, faceHandle, &umin, &vmin, &umax, &vmax, &err);

  int numFaces  = face->getNumMeshElements();

  map<MVertex*,int> vmap;
  for( int i = 0; i < numFaces; i++) 
  {
       MElement *melem  = face->getMeshElement(i);
       for( int j = 0; j < melem->getNumVertices(); j++) {
            MVertex *vertex = melem->getVertex(j);
            if( vmap.find( vertex ) == vmap.end() ) {
                if( edgeVertexMap.find(vertex) != edgeVertexMap.end() ) 
                    vmap[vertex] = edgeVertexMap[vertex];
                else
                    vmap[vertex] = vertexID++;
            }
       }
  }

  int numNodes  = vmap.size();
  double u,v, u1, v1;

  ofile << "FaceID  " << faceID << endl;
  ofile << numNodes <<  "  " <<  numFaces<< endl;

  map<MVertex*,int>::const_iterator it1;

  map<int,MVertex*> idmap;
  for( it1 = vmap.begin(); it1 != vmap.end(); ++it1) 
  {
       MVertex *vertex = it1->first;
       int  id = it1->second;
       idmap[id] = vertex;
  }

  SPoint2   uv;

  double x, y, z, maxerror = 0.0;
  map<int,MVertex*> ::const_iterator it2;
  for( it2 = idmap.begin(); it2 != idmap.end(); ++it2) 
  {
       int  id = it2->first;
       MVertex *vertex = it2->second;
       x = vertex->x();
       y = vertex->y();
       z = vertex->z();
       ofile <<  id  << "  " << x << " " << y << " " << z << endl;

       reparamMeshVertexOnFace( vertex, face, uv );
       iGeom_getEntUVtoXYZ( geom, faceHandle, uv[0], uv[1], &x, &y, &z, &err);
       assert( !err );
       double dx = fabs(vertex->x() - x );
       double dy = fabs(vertex->y() - y );
       double dz = fabs(vertex->z() - z );
       maxerror  = max( maxerror, max(dx, max(dy,dz)));
  }

  maxErrorFaceEval = max( maxerror, maxErrorFaceEval );

  for( int i = 0; i < numFaces; i++) 
  {
       MElement *melem  = face->getMeshElement(i);
       for( int j = 0; j < melem->getNumVertices(); j++) {
            MVertex *vertex = melem->getVertex(j);
            ofile << vmap[vertex] << " ";
       }
       ofile << endl;
 }

}

//////////////////////////////////////////////////////////////////////////////

void build_occ_model( string filename, int whichone )
{
  int err;

  GModel *model  = new GModel;
 
  model->readOCCBREP( filename ); 

  assert( model->getNumVertices() );
  assert( model->getNumEdges() );
  assert( model->getNumFaces() );
  assert( model->getNumRegions() );

  double dmax = 0.99*numeric_limits<double>::max();
  double xmin = dmax;
  double ymin = dmax;
  double zmin = dmax;

  double xmax = -1.0*dmax;
  double ymax = -1.0*dmax;
  double zmax = -1.0*dmax;

  for( int i = 0; i < model->getNumVertices() ; i++) {
       GVertex *vertex = (GVertex *)model->getVertexByTag(i+1);
       xmin = min( xmin, vertex->x() );
       ymin = min( ymin, vertex->y() );
       zmin = min( zmin, vertex->z() );

       xmax = max( xmax, vertex->x() );
       ymax = max( ymax, vertex->y() );
       zmax = max( zmax, vertex->z() );
  }

  double dx = fabs(xmax-xmin);
  double dy = fabs(ymax-ymin);
  double dz = fabs(zmax-zmin);
  double dlmax  = minchar_len*max( dx, max(dy,dz) );

  GmshSetOption("Mesh", "CharacteristicLengthMin", dlmax);

  SVector3 deriv;

  int edgeID, ncount = 0;
  meshGEdge meshedge;
  for( int i = 0; i < model->getNumEdges() ; i++) {
       edgeID = i+ 1;
       GEdge *edge = model->getEdgeByTag(edgeID); assert( edge );
       meshedge(edge);
//     save(edge);
  }
  model->writeVTK( "modeledges.vtk");

//CSV  OCC
  int faceID;
  meshGFace meshface;
  for( int i = 92; i < model->getNumFaces() ; i++) {
       faceID = i + 1;
       GFace *face = model->getFaceByTag( faceID ); assert( face );
//     save(face);
       meshface(face);
       cout <<  "Face : " << faceID << " Mesh Elements " << face->getNumMeshElements() << endl;
       break;
  }

  model->writeVTK( "modelfaces.vtk");
  cout << "Surface Meshing Completed : " << endl;

  exit(0);
}


////////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv)
{
  int err;
  if( argc != 3 ) {
      cout << " Usage: executable geomfile  modeltype (0,1)" << endl;
      return 1;
  }


  int model_type = atoi( argv[2] );

  GmshInitialize();
  GmshSetOption("Mesh", "Algorithm", 5);

  string fname = string( argv[1] );

  if( model_type == 0) 
      build_gmsh_model( fname  );
  else
     build_occ_model( fname, model_type );

  GmshFinalize();
}

