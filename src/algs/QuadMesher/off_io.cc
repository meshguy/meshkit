#include <meshkit/Mesh.h>

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

int Mesh ::read_off_format_data(const string &fname)
{
  ifstream infile( fname.c_str(), ios::in);
  if( infile.fail() ) {
      cout << "Warning: Cann't open file " << fname << endl;
      return 1;
  }

  //  The codelet is borrowed from TriMesh Software
  vector<int> facevtx;
  double  x, y, z;

  Vertex* vertex;
  NodeSequence vnodes, connect(3);
  string str;

  infile >> str;
  assert( str == "OFF");

  size_t  numNodes, numFaces, numEdges;
  infile >> numNodes >> numFaces >> numEdges;

  Point3D p3d;
  for( int i = 0; i < numNodes; i++) {
       infile >> x >> y >> z;
       p3d[0] = x;
       p3d[1] = y;
       p3d[2] = z;
       vertex = Vertex::newObject();
       vertex->setXYZCoords(p3d);
       addNode( vertex );
  } 

  for( size_t i = 0; i < numFaces; i++) 
  {
       infile >> numNodes;

       facevtx.resize(numNodes);
       connect.resize(numNodes);
       for( int j = 0; j < numNodes; j++) 
            infile >> facevtx[j];

       switch( numNodes )
       {
          case 3:
              connect[0] = getNodeAt(facevtx[0]);
              connect[1] = getNodeAt(facevtx[1]);
              connect[2] = getNodeAt(facevtx[2]);
              break;
          case 4:
              connect[0] = getNodeAt(facevtx[0]);
              connect[1] = getNodeAt(facevtx[1]);
              connect[2] = getNodeAt(facevtx[2]);
              connect[3] = getNodeAt(facevtx[3]);
              break;
          default:
              for( int j = 0; j < numNodes; j++) 
                   connect[j] = getNodeAt(facevtx[j]);
              break; 
       }
       Face *face = new Face;
       face->setNodes( connect );
       addFace(face);
   }  
   cout << "Reading Off file complete " << endl;
   return 0;
}

//##############################################################################
