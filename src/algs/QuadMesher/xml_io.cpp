#include "meshkit/Mesh.hpp"
#include <iomanip>

using namespace Jaal;

int MeshImporter ::xml_file(const string &fname)
{
/*
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

  mesh->reserve( numNodes, 0);
  mesh->reserve( numFaces, 2);

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

       int err = face->setNodes( connect );
       if( !err ) 
          mesh->addFace(face);
       else {
          cout << "Fatal error:  Bad element " << endl;
          for( size_t j = 0; j < numNodes; j++) 
              cout << facevtx[j] << " ";
          exit(0);
          delete face;
       }
   }  

   int gid;
   infile >> str;
   if( str == "#NODE_GROUP")  {
       numNodes = mesh->getSize(0);
       for( size_t i = 0; i < numNodes; i++) {
            infile >> gid;
            Vertex *v = mesh->getNodeAt(i);
           v->setGroupID( gid );
//          v->setTag( gid );
       }
       infile >> str;
   }

   if( str == "#FACE_GROUP")  {
       numFaces = mesh->getSize(2);
       for( size_t i = 0; i < numFaces; i++) {
            infile >> gid;
            Face *f = mesh->getFaceAt(i);
//            f->setTag( gid );
            f->setGroupID( gid );
       }
   }
*/
   
   return 0;
}

//##############################################################################

int 
MeshExporter ::xml_file(Mesh *mesh, const string &s)
{
    size_t ncount;
    mesh->prune();
    mesh->enumerate(0);

    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) 
        return 1;

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    ofile << "<Mesh> " << endl;
    ofile << "<NumNodes> " <<  numnodes << " </NumNodes>" << endl;
    ofile << "<NumFaces> " <<  numfaces << " </NumFaces>" << endl;
    ofile << "<NodeCoordinates>" << endl;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        assert( v->isActive() );
        const Point3D &p3d = v->getXYZCoords();
        ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }
    ofile << "</NodeCoordinates>" << endl;

    // Write all edges ( mesh edges and feature edges ) 
    ncount = 0;
    if(ncount > 1 ) {
    ofile << "<Edges> " << endl;
    ofile << "<Count> " << ncount << " </Count" << endl;
    ofile << "</Edges> " << endl;
    }

    ncount = 0;
    if(ncount > 0) {
    ofile << "<FeatureEdges> " << endl;
    ofile << "<Count> " << ncount << " </Count" << endl;
    ofile << "</FeatureEdges> " << endl;
    }

    ofile << "<Faces>" << endl;

    // Write down all triangle faces 
    ncount = 0;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 3 ) ncount++;
    }

    if( ncount > 0) {
    ofile << "<Triangles> " << endl;
    ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 3 ) {
            for( int j = 0; j < 3; j++) 
                ofile << face->getNodeAt(j)->getID() << " ";
           ofile << endl;
        }
    }
    ofile << "</Triangles> " << endl;
    }

    ncount = 0;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 4 ) ncount++;

    }

    // Write down all quad  faces 
    if( ncount > 0) {
    ofile << "<Quads> " << endl;
    ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 4 ) {
            for( int j = 0; j < 4; j++) 
                ofile << face->getNodeAt(j)->getID() << " ";
           ofile << endl;
        }
    }
    ofile << "</Quads> " << endl;
    }

    ncount = 0;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) > 4 ) ncount++;
    }

    // Write down all polygonal faces 
    if( ncount > 0) {
    ofile << "<Polygons> " << endl;
    ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) > 4 ) {
            int nn = face->getSize(0);    
            ofile << nn << " ";
            for( int j = 0; j < nn; j++) 
                ofile << face->getNodeAt(j)->getID() << " ";
           ofile << endl;
        }
    }
    ofile << "</Polygons> " << endl;
    }

    ofile << "</Faces>" << endl;

    ncount = 0;
    if( ncount > 0) {
    ofile << "<Partition>" << endl;
    ofile << "<Count> " <<  ncount << " </Count" <<  endl;

    if( ncount > 1) {
    ofile << "<InterfaceEdges>" << endl;
    ofile << "</InterfaceEdges>" << endl;

    ofile << "<NodesPartID>" << endl;
    ofile << "</NodesPartID>" << endl;

    ofile << "<FacesPartID>" << endl;
    ofile << "</FacesPartID>" << endl;
    } 

    ofile << "</Partition>" << endl;
    }

    ofile << "<Attributes>" << endl;

    ofile << "<NodeAttributes>" << endl;
    ofile << "</NodeAttributes>" << endl;

    ofile << "<EdgeAttributes>" << endl;
    ofile << "</EdgeAttributes>" << endl;

    ofile << "<FaceAttributes>" << endl;
    ofile << "</FaceAttributes>" << endl;
 
    ofile << "</Attributes>" << endl;

/*
    size_t nn = numnodes;
    ofile << "OFF" << endl;

    NodeSequence oldConnect, newConnect;
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive()) {
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
                assert(!newConnect[j]->isRemoved());
                cout << face->getStatus() << endl;
                cout << "Vertex indexing out of range " << vid << " Total : " << numnodes << endl;
                exit(0);
            }
            ofile << vid << " ";
        }
        ofile << endl;
        }
    }

   ofile << "#NODE_GROUP ";
   size_t numNodes = mesh->getSize(0);
   for( size_t i = 0; i < numNodes; i++) {
       Vertex *v = mesh->getNodeAt(i);
//     ofile << v->getTag() << " ";
       ofile << v->getGroupID() << " ";
   }
   ofile << endl;

   ofile << "#FACE_GROUP ";
   size_t numFaces = mesh->getSize(2);
   for( size_t i = 0; i < numFaces; i++) {
        Face *f = mesh->getFaceAt(i);
//      ofile << f->getTag() << " ";
        ofile << f->getGroupID() << " ";
   }
*/
   ofile << "</Mesh>" << endl;
   
   return 0;
}
