#include <meshkit/Mesh.hpp>

using namespace Jaal;

void
Mesh::save_off_format(const string &s)
{
    if (!isPruned())
    {
        prune();
        enumerate(0);
        enumerate(2);
    }
    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);

    if (!is_consistently_oriented())
    {
        cout << "Warning: Mesh is not conistently oriented " << endl;
    }

    size_t numnodes = nodes.size();
    size_t numfaces = faces.size();

    size_t nn = numnodes;
    ofile << "OFF" << endl;

    ofile << nn << " " << numfaces << " 0  " << endl;

    for (size_t i = 0; i < numnodes; i++)
    {
        Point3D p3d = nodes[i]->getXYZCoords();
        ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

    NodeSequence oldConnect, newConnect;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = faces[i];
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
}

/////////////////////////////////////////////////////////////////////////////////

void
Mesh::save_simple_format(const string &s)
{
    if (!isPruned())
    {
        prune();
        enumerate(0);
        enumerate(2);
    }
    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);

    if (!is_consistently_oriented())
    {
        cout << "Warning: Mesh is not conistently oriented " << endl;
    }

    size_t numnodes = nodes.size();
    size_t numfaces = faces.size();

    size_t nn = numnodes;

    ofile << nn << " " << numfaces << endl;

    for (size_t i = 0; i < numnodes; i++)
    {
        Point3D p3d = nodes[i]->getXYZCoords();
        ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

    vector<PNode> oldConnect, newConnect;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = faces[i];
        face->start_from_concave_corner();
        int nnodes = face->getSize(0);
        ofile << nnodes << " ";
        for (int j = 0; j < nnodes; j++)
        {
            Vertex *vtx = face->getNodeAt(j);
            size_t vid = vtx->getID();
            if (vid >= numnodes)
            {
                assert(!face->isRemoved());
                assert(!vtx->isRemoved());
                cout << face->getID() << endl;
                cout << face->isRemoved() << endl;
                cout << "Vertex indexing out of range " << vid << endl;
                exit(0);
            }
            ofile << vid << " ";
        }
        ofile << endl;
    }

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = nodes[i];
        ofile << vertex->getTag() << " ";
    }
    ofile << endl;

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = faces[i];
        ofile << face->getTag() << " ";
    }
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::read_simple_format_data( const string &fname)
{
  ifstream infile( fname.c_str(), ios::in);
  if( infile.fail() )  {
      cout << "Warning: cann't open node file " << fname << endl;
      return 1;
  }

  cout << " Reading Simple file " << endl;

  size_t id, numnodes, numfaces;

  infile >> numnodes >> numfaces;

  nodes.resize( numnodes );

  double x, y, z = 0.0;
  Point3D xyz;
  for( size_t i = 0; i < numnodes; i++)  {
       infile >> x >> y  >> z;
       xyz[0] = x;
       xyz[1] = y;
       xyz[2] = z;
       nodes[i] = Vertex::newObject();
       nodes[i]->setID(i);
       nodes[i]->setXYZCoords(xyz);
  }

  faces.resize( numfaces);

  int n0, n1, n2, n3, numelemnodes;
  NodeSequence connect;
  Face *newface;
  
  string restline;
  for( size_t i = 0; i < numfaces; i++) {
       infile >> numelemnodes;

       if( numelemnodes == 3 ) {
           connect.resize(3);
           infile >> n0 >> n1 >> n2;
           connect[0] = getNodeAt(n0);
           connect[1] = getNodeAt(n1);
           connect[2] = getNodeAt(n2);
	   newface = Face::newObject();
	   newface->setNodes(connect);
	   faces[i] = newface;
           faces[i]->setID(i);
       }

       if( numelemnodes == 4 ) {
           connect.resize(4);
           infile >> n0 >> n1 >> n2 >> n3;
           connect[0] = getNodeAt(n0);
           connect[1] = getNodeAt(n1);
           connect[2] = getNodeAt(n2);
           connect[3] = getNodeAt(n3);
	   newface = Face::newObject();
	   newface->setNodes(connect);
	   faces[i] = newface;
           faces[i]->setID(i);
       }
  }

  int itag;
  for( size_t i = 0; i < numnodes; i++)  {
      infile >> itag;
      nodes[i]->setTag(itag);
  }

  for( size_t i = 0; i < numfaces; i++)  {
      infile >> itag;
      faces[i]->setTag(itag);
  }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::read_vtk_format_data( const string &fname)
{
  ifstream infile( fname.c_str(), ios::in);
  if( infile.fail() )  {
      cout << "Warning: cann't open node file " << fname << endl;
      return 1;
  }

  size_t id, numnodes, numfaces;
  double x, y, z = 0.0;
  Point3D xyz;

  int n0, n1, n2, n3, numelemnodes;
  NodeSequence connect;
  Face *newface;

  string str;
  while( !infile.eof() ) {
        infile >> str;

	if( str == "POINTS") {
	    infile >> numnodes >> str;
            nodes.resize( numnodes );
           for( size_t i = 0; i < numnodes; i++)  {
                infile >> x >> y  >> z;
                xyz[0] = x;
                xyz[1] = y;
                xyz[2] = z;
                nodes[i] = Vertex::newObject();
                nodes[i]->setID(i);
                nodes[i]->setXYZCoords(xyz);
            }
        }

	if( str == "CELLS") {
	    infile >> numfaces >> str;
            faces.resize( numfaces);

            for( size_t i = 0; i < numfaces; i++) {
            infile >> numelemnodes;
            if( numelemnodes == 3 ) {
                connect.resize(3);
                infile >> n0 >> n1 >> n2;
                connect[0] = getNodeAt(n0);
                connect[1] = getNodeAt(n1);
                connect[2] = getNodeAt(n2);
	        newface = Face::newObject();
	        newface->setNodes(connect);
	        faces[i] = newface;
             }

            if( numelemnodes == 4 ) {
                connect.resize(4);
                infile >> n0 >> n1 >> n2 >> n3;
                connect[0] = getNodeAt(n0);
                connect[1] = getNodeAt(n1);
                connect[2] = getNodeAt(n2);
                connect[3] = getNodeAt(n3);
	        newface = Face::newObject();
	        newface->setNodes(connect);
	        faces[i] = newface;
            }
         }
       }
    }

  return 0;
}


//////////////////////////////////////////////////////////////////////////////////


void Mesh::readTriNodes( const string &fname)
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
}
///////////////////////////////////////////////////////////////////////////////
void Mesh::readTriEdges( const string &fname)
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

void Mesh::readTriFaces( const string &fname)
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

  NodeSequence connect(3);
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

           connect[0] = getNodeAt(n0);
           connect[1] = getNodeAt(n1);
           connect[2] = getNodeAt(n2);
	   newface = Face::newObject();
	   newface->setNodes(connect);
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
           connect[0] = getNodeAt(n0);
           connect[1] = getNodeAt(n1);
           connect[2] = getNodeAt(n2);
           connect[3] = getNodeAt(n3);
	   newface = Face::newObject();
	   newface->setNodes(connect);
	   faces[i] = newface;
       }
  }
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::readFromFile( const string &fname)
{
   cout << " Read File " << endl;

   int err = 1;
   if( fname.rfind(".vtk") != string::npos) {
       err = read_vtk_format_data( fname );
   }

   if( fname.rfind(".off") != string::npos)  {
       err = read_off_format_data( fname );
   }

   if( fname.rfind(".dat") != string::npos)  {
       err = read_simple_format_data( fname );
   }

   if( err )  {
       err = read_triangle_format_data(fname);
   }

   return err;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::read_triangle_format_data( const string &fname)
{
   readTriNodes( fname );
   readTriEdges( fname );
   readTriFaces( fname );
   enumerate(2);
   return 0;
}

///////////////////////////////////////////////////////////////////////////////

/*
#ifdef USE_MOAB
int Jaal::readMeshData( iMesh_Instance &imesh, const string &fname)
{
    Mesh *m = new Mesh; assert( m != NULL );
    m->readFromFile(fname);
// m = readOffData(fname);

    if( !m->isConsistentlyOriented() ) {
       cout << "Warning:Trying to make Triangle Mesh consistently oriented " << endl;
       m->makeConsistentlyOriented();
       if( m->isConsistentlyOriented() )
            cout << "Info: Triangle Mesh is now consistently oriented: Very good " << endl;
       else
            cout << "Alas ! Triangle Mesh is still inconsistently oriented: Check manually " << endl;
    }

    iBase_EntitySetHandle rootSet = 0;
//  map<Vertex*, iBase_EntityHandle>  mapNodes;

    m->toMOAB(imesh, rootSet);

    delete m;
    return 0;
}
#endif
*/
///////////////////////////////////////////////////////////////////////////////

