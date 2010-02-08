#include <assert.h>
#include <fstream>

#include <Mesquite_all_headers.hpp>

#include <vector>
#include <map>
#include <set>

using namespace std;

vector<double>  vCoords;
vector<unsigned long >  vConnect;

map<int, vector<int> > relations02;  

vector<int>  vfixed;

std::map<int,int> global2local;

using namespace Mesquite;

///////////////////////////////////////////////////////////////////////////////

void readNodes( const string &fname)
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

  vCoords.resize( 3*numnodes );

  double x, y, z = 0.0;
  for( int i = 0; i < numnodes; i++)  {
       infile >> id >> x >> y ;
       if( boundflag ) infile >> bmark;
       global2local[id] = i;
       vCoords[3*i+0]  = x;
       vCoords[3*i+1]  = y;
       vCoords[3*i+2]  = z;
  }
}

///////////////////////////////////////////////////////////////////////////////

void saveNodes( const string &fname)
{
  string filename = fname + ".node";

  ofstream infile( filename.c_str(), ios::out);
  if( infile.fail() )  {
      cout << "Warning: cann't open node file " << filename << endl;
      return;
  }

  int id, numnodes, ndim, numattrib, boundflag, bmark;
  numnodes  = vCoords.size()/3;

  infile <<  numnodes << " 2 0 0 " << endl;

  for( int i = 0; i < numnodes; i++) 
       infile << i  << " " << vCoords[3*i] << " " << vCoords[3*i+1] << endl;
}

///////////////////////////////////////////////////////////////////////////////

void readFaces( const string &fname)
{
  string filename = fname + ".ele";
  ifstream infile( filename.c_str(), ios::in);
  if( infile.fail() ) {
      cout << "Warning: cann't open node file " << filename << endl;
      return;
  }

  cout << "Reading face file " << filename << endl;

  int id, numfaces, numelemnodes, boundflag, bmark;
  int n0, n1, n2, n3;
  
  infile >> numfaces >> numelemnodes >> boundflag;
  assert( numelemnodes == 4 );

  vConnect.resize( 4*numfaces);
  for( int i = 0; i < numfaces; i++)  {
       infile >> id >> n0 >> n1 >> n2 >> n3;
       if( boundflag )  infile >> bmark;
       n0 = global2local[n0];
       n1 = global2local[n1];
       n2 = global2local[n2];
       n3 = global2local[n3];
       vConnect[4*i+0] = n0;
       vConnect[4*i+1] = n1;
       vConnect[4*i+2] = n2;
       vConnect[4*i+3] = n3;
  }
}

////////////////////////////////////////////////////////////////////////////////

void saveFaces( const string &fname)
{
  string filename = fname + ".ele";
  ofstream infile( filename.c_str(), ios::out);
  if( infile.fail() ) {
      cout << "Warning: cann't open node file " << filename << endl;
      return;
  }

  int id, numfaces, numelemnodes, boundflag, bmark;
  int n0, n1, n2, n3;

  numfaces = vConnect.size()/4;
  
  infile << numfaces << " 4 0 " << endl;

  for( int i = 0; i < numfaces; i++)  {
       infile << i << " " 
              << vConnect[4*i+0] << "  "
              << vConnect[4*i+1] << "  "
              << vConnect[4*i+2] << "  "
              << vConnect[4*i+3] << endl;
  }
}

///////////////////////////////////////////////////////////////////////////////

void readData( const string &fname)
{
   readNodes( fname );
   readFaces( fname );
}

///////////////////////////////////////////////////////////////////////////////

void saveData()
{
   string fname = "smooth";
   saveNodes( fname );
   saveFaces( fname );
}

///////////////////////////////////////////////////////////////////////////////
bool isBoundary(int n0, int n1 )
{
   vector<int>  eCommon;

   set_intersection( relations02[n0].begin(), relations02[n0].end(),
                     relations02[n1].begin(), relations02[n1].end(),
                     inserter( eCommon, eCommon.begin() ) );
   if( eCommon.size() == 1) return 1;
   return 0;
}

///////////////////////////////////////////////////////////////////////////////

void filter_boundary_nodes()
{
   relations02.clear();

   size_t numQuads = vConnect.size()/4;

   for( int iface = 0; iface < numQuads; iface++) 
   {
        int n0 = vConnect[4*iface+0];
        int n1 = vConnect[4*iface+1];
        int n2 = vConnect[4*iface+2];
        int n3 = vConnect[4*iface+3];
        relations02[n0].push_back( iface );
        relations02[n1].push_back( iface );
        relations02[n2].push_back( iface );
        relations02[n3].push_back( iface );
   }

   set<int> boundNodes;
   for( int iface = 0; iface < numQuads; iface++) 
   {
        int n0 = vConnect[4*iface+0];
        int n1 = vConnect[4*iface+1];
        int n2 = vConnect[4*iface+2];
        int n3 = vConnect[4*iface+3];
        if( isBoundary( n0, n1 ) )
        {
            boundNodes.insert( n0 );
            boundNodes.insert( n1 );
        }

        if( isBoundary( n1, n2 ) )
        {
            boundNodes.insert( n1 );
            boundNodes.insert( n2 );
        }

        if( isBoundary( n2, n3 ) ) 
        {
            boundNodes.insert( n2 );
            boundNodes.insert( n3 );
        }

        if( isBoundary(n3, n0 )  )
        {
            boundNodes.insert( n3 );
            boundNodes.insert( n0 );
        }
    }
    cout << "Number of boundary nodes " << boundNodes.size() << endl;
    int numNodes = vCoords.size()/3;

    vfixed.resize( numNodes );
    for( int i = 0; i < numNodes; i++)
         vfixed[i] = 0;

    set<int>::const_iterator it;
    for( it = boundNodes.begin(); it != boundNodes.end(); ++it)
         vfixed[*it] = 1;
}

///////////////////////////////////////////////////////////////////////////////

void get_valency_info()
{
   map<int,set<int> > relations00;

   size_t numQuads = vConnect.size()/4;

   for( int i = 0; i < numQuads; i++) 
   {
        int n0 = vConnect[4*i+0];
        int n1 = vConnect[4*i+1];
        int n2 = vConnect[4*i+2];
        int n3 = vConnect[4*i+3];

        relations00[n0].insert(n1);
        relations00[n1].insert(n0);

        relations00[n1].insert(n2);
        relations00[n2].insert(n1);

        relations00[n2].insert(n3);
        relations00[n3].insert(n2);

        relations00[n3].insert(n0);
        relations00[n0].insert(n3);
   }
}

///////////////////////////////////////////////////////////////////////////////
void MeshQuiteOpt( )
{
    filter_boundary_nodes();

    unsigned long numNodes = vCoords.size()/3;
    unsigned long numQuads = vConnect.size()/4;

    // Create a mesh for the MeshQuite
    ArrayMesh mesh(3, numNodes, &vCoords[0], &vfixed[0], numQuads,
                   QUADRILATERAL, &vConnect[0]);

    PlanarDomain domain( PlanarDomain::XY );

    MsqError  err;
 
    LaplacianIQ laplacian_smoother;
    laplacian_smoother.run_instructions(&mesh, &domain, err);

    if (err) return;

/*
    ShapeImprovementWrapper  shape_wrapper(err);
    if( err ) {
        cout << "Shape wrapper error " << err << endl;
        exit(2);
    }
    shape_wrapper.run_instructions( &mesh, &domain, err);

    if( err ) 
    {
       cout << "Error smoothing mesh " << err << endl;
       return;
    }
*/

    saveData();
}
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char **argv )
{
    if( argc != 2 ) {
        cout << " Usage : executable <mesh file> " << endl;
        return 1;
    }

    readData( argv[1] );

    MeshQuiteOpt( );
 
    return 0;
}









