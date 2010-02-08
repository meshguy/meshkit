#include "Mesh.h"

using namespace Jaal;

//##############################################################################
void skip_comments(FILE *f)
{
    // Skip comments in an ASCII file (lines beginning with #)
    int c;
    bool in_comment = false;
    while (1) {
         c = fgetc(f);
         if (c == EOF) return;
         if (in_comment) {
              if (c == '\n') in_comment = false;
         } else if (c == '#') {
            in_comment = true;
         } else if (!isspace(c)) {
            break;
         }
    }
    ungetc(c, f);
}
//##############################################################################

Mesh* Jaal::readOffData( const string &fname)
{
  ifstream infile( fname.c_str(), ios::in);
  if( infile.fail() ) {
      cout << "Warning: Cann't open file " << fname << endl;
      return NULL;
  }

  Mesh *mesh = new Mesh;
 
  //  The codelet is borrowed from TriMesh Software
  vector<int> facevtx;
  double  x, y, z;

  Vertex* vertex;
  vector<Vertex*> vnodes, connect(3);
  string str;

  infile >> str;
  assert( str == "OFF");

  int  numNodes, numFaces, numEdges;
  infile >> numNodes >> numFaces >> numEdges;

  Point3D p3d;
  for( int i = 0; i < numNodes; i++) {
       infile >> x >> y >> z;
       p3d[0] = x;
       p3d[1] = y;
       p3d[2] = z;
       vertex = Vertex::newObject();
       vertex->setXYZCoords(p3d);
       mesh->addNode( vertex );
  } 

  for( int i = 0; i < numFaces; i++) 
  {
       infile >> numNodes;

       facevtx.resize(numNodes);
       connect.resize(numNodes);
       for( int j = 0; j < numNodes; j++) 
            infile >> facevtx[j];

       switch( numNodes )
       {
          case 3:
              connect[0] = mesh->getNode(facevtx[0]);
              connect[1] = mesh->getNode(facevtx[1]);
              connect[2] = mesh->getNode(facevtx[2]);
              break;
          case 4:
              connect[0] = mesh->getNode(facevtx[0]);
              connect[1] = mesh->getNode(facevtx[1]);
              connect[2] = mesh->getNode(facevtx[2]);
              connect[3] = mesh->getNode(facevtx[3]);
              break;
          default:
              for( int j = 0; j < numNodes; j++) 
                   connect[j] = mesh->getNode(facevtx[j]);
              break; 
       }
       Face *face = new Face;
       face->setConnection( connect );
       mesh->addFace(face);
   }  

   return mesh;
}

//##############################################################################
