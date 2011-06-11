#include "Mesh.hpp"

#include <iomanip>

using namespace Jaal;

int MeshExporter :: vtk_file(Mesh *mesh, const string &fname)
{
    ofstream ofile(fname.c_str(), ios::out);
    if( ofile.fail() ) return 1;

    size_t numNodes = mesh->getSize(0);

    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << " Jaal Mesh Generator " << endl;

    ofile << "ASCII " << endl;

    ofile << "DATASET UNSTRUCTURED_GRID " << endl;
    ofile << "POINTS " << numNodes << " float " << endl;
    ofile << std::setiosflags(ios::fixed);

    for( size_t i = 0; i < numNodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        const Point3D &xyz = vertex->getXYZCoords();
        ofile << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
    }

    size_t numFaces = mesh->getSize(2);
    size_t numConn = 0;

    for( size_t i = 0; i < numFaces; i++) {
        Face *face = mesh->getFaceAt(i);
        numConn += face->getSize(0);
    }

    ofile << "CELLS " << numFaces << "  " << numFaces + numConn << endl;

    for( size_t i = 0; i < numFaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nsize = face->getSize(0);
        ofile << nsize << " ";
        for(int i = 0; i < nsize; i++) {
            Vertex *v = face->getNodeAt(i);
            ofile << v->getID() << " ";
        }
        ofile << endl;
    }

    ofile << "CELL_TYPES " << numFaces << endl;
    for( size_t i = 0; i < numFaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nsize = face->getSize(0);
        int type  = 7;
        if( nsize == 3 ) type = 5;
        if( nsize == 4 ) type = 9;
        ofile << type << " ";
    }
    return 0;
}

int MeshImporter ::vtk_file( const string &fname)
{
    cout << " VTK File " << endl;
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << fname << endl;
        return 1;
    }

    size_t numnodes, numfaces;
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
            mesh->reserve( numnodes, 0);
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
        }

        if( str == "CELLS") {
            cout << " READ CELLS " << endl;
            infile >> numfaces >> str;
            mesh->reserve( numfaces, 2);

            for( size_t i = 0; i < numfaces; i++) {
                infile >> numelemnodes;
                cout <<  numelemnodes << endl;

                switch( numelemnodes ) {
                case 3:
                    connect.resize(3);
                    infile >> n0 >> n1 >> n2;
                    connect[0] = mesh->getNodeAt(n0);
                    connect[1] = mesh->getNodeAt(n1);
                    connect[2] = mesh->getNodeAt(n2);
                    newface = Face::newObject();
                    newface->setNodes(connect);
                    mesh->addFace( newface );
                    break;
                case 4:
                    connect.resize(4);
                    infile >> n0 >> n1 >> n2 >> n3;
                    connect[0] = mesh->getNodeAt(n0);
                    connect[1] = mesh->getNodeAt(n1);
                    connect[2] = mesh->getNodeAt(n2);
                    connect[3] = mesh->getNodeAt(n3);
                    newface = Face::newObject();
                    newface->setNodes(connect);
                    mesh->addFace( newface );
                    break;
                default:
                    for( int i = 0; i < numelemnodes; i++)
                        infile >> n0;
                }
            }
        }
    }

    return 0;
}



