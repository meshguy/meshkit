#include "Mesh.h"
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void update_vertex_position( Vertex *vertex, double &localerror, double &maxerror, double lambda = 0.5)
{
   localerror = 0.0;
   if( vertex->isConstrained() ) return;

   vector<Vertex*> neighs = vertex->getRelations0();

   Point3D newpos, xyz;
   newpos[0] = 0.0;
   newpos[1] = 0.0;
   newpos[2] = 0.0;
   for( size_t i = 0; i < neighs.size(); i++) 
   {
        xyz = neighs[i]->getXYZCoords();
        newpos[0] += xyz[0];
        newpos[1] += xyz[1];
        newpos[2] += xyz[2];
   }

   newpos[0] /= (double)neighs.size();
   newpos[1] /= (double)neighs.size();
   newpos[2] /= (double)neighs.size();

   xyz = vertex->getXYZCoords(); // Old position

   localerror =  max( localerror, fabs(newpos[0] - xyz[0] ));
   localerror =  max( localerror, fabs(newpos[1] - xyz[1] ));
   localerror =  max( localerror, fabs(newpos[2] - xyz[2] ));
   maxerror   =  max( maxerror, localerror);

   newpos[0] = (1.0-lambda)*xyz[0] + lambda*newpos[0];
   newpos[1] = (1.0-lambda)*xyz[1] + lambda*newpos[1];
   newpos[2] = (1.0-lambda)*xyz[2] + lambda*newpos[2];

   vertex->setXYZCoords(newpos);
}

///////////////////////////////////////////////////////////////////////////////

void Jaal::laplacian_smoothing( Mesh *mesh, int numIters, int verbose)
{
  if( verbose ) 
      cout << " Laplacian Smoothing ... " << endl;

  mesh->search_boundary();

  int relexist = mesh->build_relations(0,0);

  size_t  numnodes = mesh->getSize(0);
 
  double localerror, maxerror = 0.0;
  for( int iter = 0; iter < numIters; iter++) {
       maxerror = 0.0;
       for( int i = 0; i < numnodes; i++) 
            update_vertex_position( mesh->getNodeAt(i) , localerror, maxerror);
  }

  if( verbose ) 
      cout << "Info: Laplacian Maximum Error : " << maxerror << endl;

  if( !relexist) mesh->clear_relations(0,0);
}

///////////////////////////////////////////////////////////////////////////////

void Jaal::laplacian_smoothing( Mesh *mesh, vector<Vertex*> &vertexQ, int numIters, int verbose)
{

  /////////////////////////////////////////////////////////////////////////////
  // This is a local smoothing procedure. User specifies nodes which needs to be 
  // smoothed. So the smmoothing starts from the specified nodes and neighbouring
  // nodes are added if the error between the old and new position is more than
  // "tolerance" limit. This procedure is effective when the mesh has been globally
  // smoothed before calling this function, otherwise, it may behave like a global
  // smoothner...
  //
  // The idea of "local smoothing" is not new and numerous papers talk about it.
  // But I haven't seen its implementation anywhere before, so I am not citing
  // any reference for it.
  //
  // Also, if you are calling this function many times, for good optimization,
  // call mesh->build_relations(0,0) before calling this function, otherwise
  // for each invocation of this function, this function will be called,
  // relationships will be build and destroyed.
  //
  /////////////////////////////////////////////////////////////////////////////
  
  if( verbose ) 
      cout << " Laplacian Smoothing ... " << endl;

  mesh->search_boundary();

  int relexist = mesh->build_relations(0,0);

  size_t numnodes = mesh->getSize(0);

  for( size_t i = 0; i < numnodes; i++) {
       Vertex *vertex = mesh->getNodeAt(i);
       vertex->setVisitMark(0);
  }

  for( size_t i = 0; i < vertexQ.size(); i++)
       vertexQ[i]->setVisitMark(1);

  for( size_t i = 0; i < vertexQ.size(); i++) {
       Vertex *vertex = vertexQ[i];
       vertex->setVisitMark(1);
       if( vertex->isConstrained() ) {
           vector<Vertex*> neighs = vertex->getRelations0();
	   for( int j = 0; j < neighs.size(); j++) {
	        if( neighs[j]->isVisited() == 0) {
		    vertexQ.push_back( neighs[j] );
                    neighs[j]->setVisitMark(1);
                }
           }
       }
  }

  double localerror, maxerror = 0.0;
  for( int iter = 0; iter < numIters; iter++) {
       maxerror = 0.0;
       for( size_t i = 0; i < vertexQ.size(); i++) {
            Vertex *vertex = vertexQ[i];
            update_vertex_position( vertex, localerror, maxerror);
	    if( localerror > 1.0E-05) {
                vector<Vertex*> neighs = vertex->getRelations0();
	        for( int j = 0; j < neighs.size(); j++) {
	            if( neighs[j]->isVisited() == 0) {
		        vertexQ.push_back( neighs[j] );
                        neighs[j]->setVisitMark(1);
                    }
                }
             }
       }
  }

  if( verbose ) 
      cout << "Info: Laplacian Maximum Error : " << maxerror << endl;

  if( !relexist) mesh->clear_relations(0,0);
}

///////////////////////////////////////////////////////////////////////////////

void Jaal::laplacian_smoothing(iMesh_Instance imesh, int numIters)
{
      Mesh *jmesh = new Mesh;
      jmesh->fromMOAB(imesh);
      Jaal::laplacian_smoothing( jmesh, numIters);
      delete jmesh;
}
///////////////////////////////////////////////////////////////////////////////

