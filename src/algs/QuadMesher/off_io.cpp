#include "Mesh.hpp"

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

int MeshImporter ::off_file(const string &fname)
{
  if( mesh == NULL ) mesh = new Mesh;

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
  for( size_t i = 0; i < numNodes; i++) {
       infile >> x >> y >> z;
       p3d[0] = x;
       p3d[1] = y;
       p3d[2] = z;
       vertex = Vertex::newObject();
       vertex->setXYZCoords(p3d);
       mesh->addNode( vertex );
  } 

  for( size_t i = 0; i < numFaces; i++) 
  {
       infile >> numNodes;

       facevtx.resize(numNodes);
       connect.resize(numNodes);
       for( size_t j = 0; j < numNodes; j++) 
            infile >> facevtx[j];

       switch( numNodes )
       {
          case 3:
              connect[0] = mesh->getNodeAt(facevtx[0]);
              connect[1] = mesh->getNodeAt(facevtx[1]);
              connect[2] = mesh->getNodeAt(facevtx[2]);
              break;
          case 4:
              connect[0] = mesh->getNodeAt(facevtx[0]);
              connect[1] = mesh->getNodeAt(facevtx[1]);
              connect[2] = mesh->getNodeAt(facevtx[2]);
              connect[3] = mesh->getNodeAt(facevtx[3]);
              break;
          default:
              for( size_t j = 0; j < numNodes; j++) 
                   connect[j] = mesh->getNodeAt(facevtx[j]);
              break; 
       }
       Face *face = new Face;
       face->setNodes( connect );
       mesh->addFace(face);
   }  

   return 0;
}

//##############################################################################

int 
MeshExporter ::off_file(Mesh *mesh, const string &s)
{
    if (!mesh->isPruned())
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
    }
    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) 
        return 1;

    if (!mesh->is_consistently_oriented())
    {
        cout << "Warning: Mesh is not conistently oriented " << endl;
    }

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    size_t nn = numnodes;
    ofile << "OFF" << endl;

    ofile << nn << " " << numfaces << " 0  " << endl;

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        const Point3D &p3d = v->getXYZCoords();
        ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

    NodeSequence oldConnect, newConnect;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if (face->getSize(0) == 4)
        {
            oldConnect = face->getNodes();
            Face::quad_tessalate(oldConnect, newConnect); // Because of OpenGL
        }
        else
            newConnect = face->getNodes();

        int nnodes = newConnect.size();
        ofile << nnodes << " ";
        for (int j = 0; j < nnodes; j++)
        {
            size_t vid = newConnect[j]->getID();
            if (vid >= numnodes)
            {
                assert(!face->isRemoved());
                assert(!newConnect[j]->isRemoved());
                cout << face->getID() << endl;
                cout << face->isRemoved() << endl;
                cout << "Vertex indexing out of range " << vid << endl;
                exit(0);
            }
            ofile << vid << " ";
        }
        ofile << endl;
    }
    return 0;
}
