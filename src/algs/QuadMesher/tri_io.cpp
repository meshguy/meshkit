#include <meshkit/Mesh.hpp>

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void MeshImporter ::readTriNodes( const string &fname)
{
  string filename = fname + ".node";

  ifstream infile( filename.c_str(), ios::in);
  if( infile.fail() )  {
      cout << "Warning: cann't open node file " << filename << endl;
      return;
  }

  cout << "Reading node file " << filename << endl;

  int id, numnodes, ndim, numattrib, boundflag, bmark;

  infile >> numnodes >> ndim >> numattrib >> boundflag;

  assert( ndim == 2 );
  assert( numattrib == 0);

  mesh->reserve(numnodes, 0);

  double x, y, z = 0.0;
  Point3D xyz;
  for( int i = 0; i < numnodes; i++)  {
       infile >> id >> x >> y ;
       if( ndim == 3) infile >> z;
       if( boundflag ) infile >> bmark;
       global2local[id] = i;
       xyz[0] = x;
       xyz[1] = y;
       xyz[2] = z;
       Vertex *v = Vertex::newObject();
       v->setID(i);
       v->setXYZCoords(xyz);
       mesh->addNode( v );
  }
}
///////////////////////////////////////////////////////////////////////////////
void MeshImporter::readTriEdges( const string &fname)
{
/*
  string filename = fname + ".edge";
  ifstream infile( filename.c_str(), ios::in);
  if( infile.fail() ) {
      cout << "Warning: cann't open node file " << filename << endl;
      return;
  }

  cout << "Reading edge file " << filename << endl;

  int id, numedges, boundflag;
  int n0, n1;

  infile >> numedges >> boundflag;
  edges.resize( numedges );

  for( int i = 0; i < numedges; i++) {
       infile >> id >> n0 >> n1;
       if( boundflag ) infile >> edges[i].gid;
       n0 = global2local[n0];
       n1 = global2local[n1];
       edges[i].connect[0] = n0;
       edges[i].connect[1] = n1;
  }
  */
}

///////////////////////////////////////////////////////////////////////////////

void MeshImporter ::readTriFaces( const string &fname)
{
  string filename = fname + ".ele";
  ifstream infile( filename.c_str(), ios::in);
  if( infile.fail() ) {
      cout << "Warning: cann't open node file " << filename << endl;
      return;
  }

  int id, numfaces, numelemnodes, boundflag, bmark;
  int n0, n1, n2, n3, facetype;
  
  infile >> numfaces >> numelemnodes >> boundflag;

  mesh->reserve( numfaces, 2);

  vector<PNode> connect(3);
  Face *newface;
  
  if( numelemnodes == 3 ) {
      facetype = 3;
      connect.resize(3);
      for( int i = 0; i < numfaces; i++) {
           infile >> id >> n0 >> n1 >> n2;
           if( boundflag )  infile >> bmark;
           n0 = global2local[n0];
           n1 = global2local[n1];
           n2 = global2local[n2];

           connect[0] = mesh->getNodeAt(n0);
           connect[1] = mesh->getNodeAt(n1);
           connect[2] = mesh->getNodeAt(n2);
	   newface = new Face;
	   newface->setNodes(connect);
	   mesh->addFace( newface );
       }
  }

  if( numelemnodes == 4 ) {
      facetype = 4;
      connect.resize(4);
      for( int i = 0; i < numfaces; i++)  {
           infile >> id >> n0 >> n1 >> n2 >> n3;
           if( boundflag )  infile >> bmark;
           n0 = global2local[n0];
           n1 = global2local[n1];
           n2 = global2local[n2];
           n3 = global2local[n3];
           connect[0] = mesh->getNodeAt(n0);
           connect[1] = mesh->getNodeAt(n1);
           connect[2] = mesh->getNodeAt(n2);
           connect[3] = mesh->getNodeAt(n3);
	   newface = new Face;
	   newface->setNodes(connect);
	   mesh->addFace( newface );
       }
  }
  
}

///////////////////////////////////////////////////////////////////////////////
int MeshImporter ::triangle_file( const string &fname)
{
   readTriNodes( fname );
   readTriEdges( fname );
   readTriFaces( fname );
   mesh->enumerate(2);
   global2local.clear();
   return 0;
}

///////////////////////////////////////////////////////////////////////////////
