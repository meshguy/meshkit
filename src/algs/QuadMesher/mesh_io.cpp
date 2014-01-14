#include "meshkit/Mesh.hpp"

using namespace Jaal;

int 
MeshExporter::simple_file( Mesh *mesh, const string &s)
{
    if (!mesh->isPruned())
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
    }

    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) { 
        return 1;
    }

    if (!mesh->is_consistently_oriented())
    {
        cout << "Warning: Mesh is not conistently oriented " << endl;
    }

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    size_t nn = numnodes;

    ofile << nn << " " << numfaces << endl;

    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        const Point3D &p3d = v->getXYZCoords();
        ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

    vector<PNode> oldConnect, newConnect;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
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

/*
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        ofile << vertex->getTag() << " ";
    }
    ofile << endl;

    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        ofile << face->getTag() << " ";
    }
*/
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int MeshImporter::simple_file(const string &fname)
{
  ifstream infile( fname.c_str(), ios::in);
  if( infile.fail() )  {
      cout << "Warning: cann't open node file " << fname << endl;
      return 1;
  }

  size_t numnodes, numfaces;

  infile >> numnodes >> numfaces;

  mesh->reserve(numnodes,0);

  double x, y, z = 0.0;
  Point3D xyz;
  for( size_t i = 0; i < numnodes; i++)  {
       infile >> x >> y  >> z;
       xyz[0] = x;
       xyz[1] = y;
       xyz[2] = z;
       Vertex *v = Vertex::newObject();
       v->setID(i);
       v->setXYZCoords(xyz);
       mesh->addNode(v);
  }

  mesh->reserve( numfaces, 2);

  int n0, n1, n2, n3, numelemnodes;
  NodeSequence connect;
  Face *newface;
  
  string restline;
  for( size_t i = 0; i < numfaces; i++) {
       infile >> numelemnodes;

       if( numelemnodes == 3 ) {
           connect.resize(3);
           infile >> n0 >> n1 >> n2;
           connect[0] = mesh->getNodeAt(n0);
           connect[1] = mesh->getNodeAt(n1);
           connect[2] = mesh->getNodeAt(n2);
	   newface = Face::newObject();
	   newface->setNodes(connect);
           newface->setID(i);
           mesh->addFace( newface );
       }

       if( numelemnodes == 4 ) {
           connect.resize(4);
           infile >> n0 >> n1 >> n2 >> n3;
           connect[0] = mesh->getNodeAt(n0);
           connect[1] = mesh->getNodeAt(n1);
           connect[2] = mesh->getNodeAt(n2);
           connect[3] = mesh->getNodeAt(n3);
	   newface = Face::newObject();
	   newface->setNodes(connect);
           newface->setID(i);
           mesh->addFace( newface );
       }
  }

/*
  int itag;
  for( size_t i = 0; i < numnodes; i++)  {
      infile >> itag;
      Vertex *v = mesh->getNodeAt(i);
      v->setTag(itag);
  }

  for( size_t i = 0; i < numfaces; i++)  {
      infile >> itag;
      Face *f = mesh->getFaceAt(i);
      f->setTag(itag);
  }
*/

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Mesh::readFromFile( const string &fname)
{
   MeshImporter mimporter;
   Mesh *m = mimporter.load( fname, this);

   if( m ) return 0;

   return 1;
  
}

///////////////////////////////////////////////////////////////////////////////
