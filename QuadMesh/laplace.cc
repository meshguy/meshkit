#include "Mesh.h"
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void update_vertex_position( Vertex *vertex, double &maxerror, double lambda = 0.5)
{
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

   maxerror =  max( maxerror, fabs(newpos[0] - xyz[0] ));
   maxerror =  max( maxerror, fabs(newpos[1] - xyz[1] ));
   maxerror =  max( maxerror, fabs(newpos[2] - xyz[2] ));

   newpos[0] = (1.0-lambda)*xyz[0] + lambda*newpos[0];
   newpos[1] = (1.0-lambda)*xyz[1] + lambda*newpos[1];
   newpos[2] = (1.0-lambda)*xyz[2] + lambda*newpos[2];

   vertex->setXYZCoords(newpos);
}

///////////////////////////////////////////////////////////////////////////////

void Jaal::laplacian_smoothing( Mesh *mesh, int numIters)
{
  cout << " Laplacian Smoothing ... " << endl;
  mesh->search_boundary();

  int relexist = mesh->build_relations(0,0);

  int numnodes = mesh->getSize(0);
 
  double maxerror = 0.0;
  for( int iter = 0; iter < numIters; iter++) {
       maxerror = 0.0;
       for( int i = 0; i < numnodes; i++) 
            update_vertex_position( mesh->getNode(i) , maxerror);
  }

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

