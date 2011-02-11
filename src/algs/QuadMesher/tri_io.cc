#include "Mesh.h"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void Mesh::readNodes( const string &fname)
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

  nodes.resize( numnodes );

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
       nodes[i] = Vertex::newObject();
       nodes[i]->setID(i);
       nodes[i]->setXYZCoords(xyz);
  }

  double xmin, xmax, ymin, ymax;
  xyz = nodes[0]->getXYZCoords();
  xmin = xyz[0];
  xmax = xyz[0];
  ymin = xyz[1];
  ymax = xyz[1];

  for( int i = 0; i < numnodes; i++)  {
        xyz = nodes[i]->getXYZCoords();
       xmin = min( xmin, xyz[0] );
       xmax = max( xmax, xyz[0] );
       ymin = min( ymin, xyz[1] );
       ymax = max( ymax, xyz[1] );
  }

  double xlen  = fabs(xmax-xmin);
  double ylen  = fabs(ymax-ymin);
  double scale = max( xlen, ylen );

  for( int i = 0; i < numnodes; i++)  
  {
       xyz = nodes[i]->getXYZCoords();
       xyz[0] /= scale;
       xyz[1] /= scale;
  }
}
///////////////////////////////////////////////////////////////////////////////
void Mesh::readEdges( const string &fname)
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

void Mesh::readFaces( const string &fname)
{
  string filename = fname + ".ele";
  ifstream infile( filename.c_str(), ios::in);
  if( infile.fail() ) {
      cout << "Warning: cann't open node file " << filename << endl;
      return;
  }

  cout << "Reading face file " << filename << endl;

  int id, numfaces, numelemnodes, boundflag, bmark;
  int n0, n1, n2, n3, facetype;
  
  infile >> numfaces >> numelemnodes >> boundflag;

  faces.resize( numfaces);

  vector<NodeType> connect(3);
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

           connect[0] = getNode(n0);
           connect[1] = getNode(n1);
           connect[2] = getNode(n2);
	   newface = new Face;
	   newface->setConnection(connect);
	   faces[i] = newface;
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
           connect[0] = getNode(n0);
           connect[1] = getNode(n1);
           connect[2] = getNode(n2);
           connect[3] = getNode(n3);
	   newface = new Face;
	   newface->setConnection(connect);
	   faces[i] = newface;
       }
  }
  
}

///////////////////////////////////////////////////////////////////////////////
void Mesh::readData( const string &fname)
{
   readNodes( fname );
   readEdges( fname );
   readFaces( fname );
   enumerate(2);
}

///////////////////////////////////////////////////////////////////////////////
int Jaal::readMeshData( iMesh_Instance &imesh, const string &fname)
{
    Mesh *m = new Mesh;
    m->readData(fname);
//  m = readOffData(fname);

    if( !m->isConsistentlyOriented() ) {
       cout << "Warning:Trying to make Triangle Mesh consistently oriented " << endl;
       m->makeConsistentlyOriented();
       if( m->isConsistentlyOriented() )
            cout << "Info: Triangle Mesh is now consistently oriented: Very good " << endl;
       else
            cout << "Alas ! Triangle Mesh is still inconsistently oriented: Check manually " << endl;
    }

    iBase_EntitySetHandle rootSet = 0;
    map<Vertex*, iBase_EntityHandle>  mapNodes;

    m->toMOAB(imesh, rootSet);

    delete m;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

