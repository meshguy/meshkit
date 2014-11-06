/** \file test_solidmesher.cpp \test
 *
 * Test the SolidMesher for a basic example.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/SolidSurfaceMesher.hpp"
#include "meshkit/SolidCurveMesher.hpp"
#include "meshkit/ModelEnt.hpp"


using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

#ifdef HAVE_ACIS
std::string extension = ".sat";
#elif HAVE_OCC
std::string extension = ".stp";
#endif

// Basic Tests
void read_cube_tris_test();
void read_cube_curves_test();
void read_cube_surfs_test();
void read_cube_vols_test();
void read_cube_vertex_pos_test();

//Connectivity Tests
void cube_verts_connectivity_test();
void cube_tris_connectivity_test();
void cube_tri_curve_coincidence_test();
void cube_edge_adjacencies_test(); 

//Other functions
int count_topo_for_dim( int dim, iMesh::EntityTopology topo);
void match_tri_edges_w_curves( std::vector<iMesh::EntityHandle> edges, MEntVector curves);

int main(int argc, char **argv) 
{
  
  // start up MK and load the geometry
  mk = new MKCore();

  std::string filename = "cube" ;

  filename = TestDir + "/" + filename + extension;

  mk->load_geometry(&filename[0]);

  MEntVector surfs;
  mk->get_entities_by_dimension(2,surfs);
  SolidSurfaceMesher *ssm;

  ssm = (SolidSurfaceMesher*) mk->construct_meshop("SolidSurfaceMesher", surfs);

  double facet_tol = 1e-04, geom_resabs = 1e-06;
  ssm->set_mesh_params(facet_tol, geom_resabs);

  mk->setup();
  mk->execute();


  //RUN TESTS
  int num_fail = 0;
  num_fail += RUN_TEST(read_cube_tris_test);
  num_fail += RUN_TEST(read_cube_curves_test);
  num_fail += RUN_TEST(read_cube_surfs_test);
  num_fail += RUN_TEST(read_cube_vols_test);
  num_fail += RUN_TEST(read_cube_vertex_pos_test);
  num_fail += RUN_TEST(cube_verts_connectivity_test);
  num_fail += RUN_TEST(cube_tris_connectivity_test);

#if HAVE_OCC
  return 0;
#else
  return num_fail;
#endif
}


//Tests
// NOTE: all tests should be performed using the iMesh interface as that is where our faceted
//       information lives. It can, however, be compared to the iGeom information if desired.


void read_cube_tris_test()
{
  
  // get the number of tris in the surface model ents
  int num_of_tris = count_topo_for_dim( 2, iMesh_TRIANGLE);

  //For a cube, there should be exactly 2 triangles per face, totaling 12 for the cube.
  CHECK_EQUAL(12, num_of_tris);

}

void read_cube_curves_test()
{
  
  MEntVector ents;
  
  // there should be 12 curve (dim==1) ModelEnts for a cube
  mk->get_entities_by_dimension(1,ents);

  int num_of_curves = ents.size();

  CHECK_EQUAL(12, num_of_curves);
}

void read_cube_surfs_test()
{
  
  MEntVector ents;
  
  // there should be 6 surf (dim==1) ModelEnts for a cube
  mk->get_entities_by_dimension(2,ents);

  int num_of_surfs = ents.size();

  CHECK_EQUAL(6, num_of_surfs);
}


void read_cube_vols_test()
{
  
  MEntVector ents;
  
  // there should be 1 vol (dim==1) ModelEnts for a cube
  mk->get_entities_by_dimension(3,ents);

  int num_of_surfs = ents.size();

  CHECK_EQUAL(1, num_of_surfs);
}

void read_cube_vertex_pos_test()
{

  MEntVector ents;

  mk->get_entities_by_dimension(0, ents);

  int num_of_verts = ents.size();

  // should be 8 vertex model ents for the cube
  CHECK_EQUAL(8, num_of_verts);

  // get the vertex coordinates
  double x[8];
  double y[8];
  double z[8];

  for(unsigned int i = 0; i < ents.size(); i++)
    {
      //get the iMesh EntityHandle for the vert
      iMesh::EntitySetHandle sh = IBSH(ents[i]->mesh_handle());

      std::vector<iMesh::EntityHandle> vert;
      mk->imesh_instance()->getEntities(sh, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES, vert);

      int vert_size = vert.size();
      // there should only be one entity in this set. If not, there is a problem.
      assert(1 == vert_size);
  
      //now get the vertex coords      
      mk->imesh_instance()->getVtxCoord(vert[0], x[i], y[i], z[i]);
    }

  //Check against known locations of the vertices

  std::vector<double> x_ref;
  std::vector<double> y_ref;
  std::vector<double> z_ref;

  // Vertex 1
  x_ref.push_back( 5 );
  y_ref.push_back( -5 );
  z_ref.push_back( 5 );

  // Vertex 2
  x_ref.push_back( 5 );
  y_ref.push_back( 5 );
  z_ref.push_back( 5 );

  // Vertex 3
  x_ref.push_back( -5 );
  y_ref.push_back( 5 );
  z_ref.push_back( 5 );

  // Vertex 4
  x_ref.push_back( -5 );
  y_ref.push_back( -5 );
  z_ref.push_back( 5 );

  // Vertex 5
  x_ref.push_back( 5 );
  y_ref.push_back( 5 );
  z_ref.push_back( -5 );

  // Vertex 6
  x_ref.push_back( 5 );
  y_ref.push_back( -5 );
  z_ref.push_back( -5 );

  // Vertex 7
  x_ref.push_back( -5 );
  y_ref.push_back( -5 );
  z_ref.push_back( -5 );
 
  // Vertex 8
  x_ref.push_back( -5 );
  y_ref.push_back( 5 );
  z_ref.push_back( -5 );

  
  for(unsigned int i=0; i<ents.size(); i++)
    {
      for(unsigned int j=0; j<x_ref.size(); j++)
	{
	  if( x[i]==x_ref[j] && y[i]==y_ref[j] && z[i]==z_ref[j] )
            {
              x_ref.erase( x_ref.begin()+j );
              y_ref.erase( y_ref.begin()+j );
              z_ref.erase( z_ref.begin()+j );
            }
	}
    }
  
  //After looping through each vertex loaded from the mesh
  //there should be no entities left in the reference vector
  int leftovers = x_ref.size();
  CHECK_EQUAL( 0, leftovers );

}

void cube_verts_connectivity_test()
{

  MEntVector verts;
  mk->get_entities_by_dimension(0, verts);

  MEntVector::iterator i;
  for( i = verts.begin(); i != verts.end(); i++)
    {
      std::vector<iMesh::EntityHandle> adj_tris;
      std::vector<int> dum;
      iMesh::EntitySetHandle sh = IBSH((*i)->mesh_handle());
      mk->imesh_instance()->getAdjEntities(sh, iBase_ALL_TYPES, iMesh_POINT, iBase_FACE, adj_tris, dum);

      int num_adj_tris = adj_tris.size();
      CHECK( num_adj_tris >=4 && num_adj_tris <=6);

    }

}


void cube_tris_connectivity_test()
{
  
  //get all triangles from the model ents
  MEntVector surfs; 
  mk->get_entities_by_dimension(2, surfs);
  
  int expected_num_of_adj_tris = 3;
  
  MEntVector::iterator i; 
  
  for( i = surfs.begin(); i != surfs.end(); i++) 
    {
      //get the triangles for each surface
      std::vector<int> dum;
      iMesh::EntitySetHandle sh = IBSH((*i)->mesh_handle()); 
      std::vector<iMesh::EntityHandle> surf_tris;
      mk->imesh_instance()->getEntities( sh, iBase_FACE, iMesh_TRIANGLE, surf_tris); 
      
      //each triangle should be adjacent to exactly 3 other triangles 
      for( std::vector<iMesh::EntityHandle>::iterator j = surf_tris.begin();
	   j != surf_tris.end(); j++)
	{
	  std::vector<iMesh::EntityHandle> adjacent_tris; 
	  mk->imesh_instance()->getEnt2ndAdj( *j, iBase_EDGE, iBase_FACE, adjacent_tris); 
	  

	  for(unsigned int k = 0; k < adjacent_tris.size(); k++) 
	    {

	      iMesh::EntityTopology topo; 
	      mk->imesh_instance()->getEntTopo(adjacent_tris[k], topo); 
	      std::cout << topo <<std::endl; 

	    }

	  CHECK(expected_num_of_adj_tris == int(adjacent_tris.size())); 
	}
    }

}

void cube_tri_curve_coincidence_test()
{
  //get all curves from the mesh 
  MEntVector curves; 
  mk->get_entities_by_dimension(1, curves); 
  
  //get all triangles from the mesh
  std::vector<iMesh::EntityHandle> tris; 
  mk->imesh_instance()->getEntities(0, iBase_FACE, iMesh_TRIANGLE, tris); 

  std::vector<iMesh::EntityHandle>::iterator i; 
  for( i = tris.begin(); i != tris.end(); i++)
    {

      //get the edges for this triangle 
      std::vector<iMesh::EntityHandle> tri_edges;
      mk->imesh_instance()->getEntAdj( *i, iBase_EDGE, tri_edges); 
      //make sure we've retrieved two edges for this triangle 
      CHECK( 2 == int(tri_edges.size()) ); 
      match_tri_edges_w_curves( tri_edges, curves ); 
      
    }
  


}

void match_tri_edges_w_curves( std::vector<iMesh::EntityHandle> edges, MEntVector curves)
{

  int match_counter = 0; 
  int num_of_tri_edges = edges.size(); 
  CHECK(num_of_tri_edges);

  for(std::vector<iMesh::EntityHandle>::iterator i = edges.begin(); 
      i != edges.end(); i++)
    {
      for(MEntVector::iterator j = curves.begin(); j != curves.end(); 
	  j++)
	{
	  //get the curve edges (there should be only one per curve for a cube)
	  iMesh::EntitySetHandle sh = IBSH((*j)->mesh_handle());
	  std::vector<iMesh::EntityHandle> cedges;
	  mk->imesh_instance()->getEntities( sh, iBase_EDGE, iMesh_LINE_SEGMENT, cedges); 

	  CHECK( 1 == int(cedges.size()) );

	  iMesh::EntityHandle curve_edge_handle = cedges[0];
	  iMesh::EntityHandle tri_edge_handle = *i;
	  if( tri_edge_handle == curve_edge_handle ) match_counter++;
	}

      
    }

  //make sure we found a match for each triangle edge to a curve edge
  CHECK( num_of_tri_edges == match_counter ); 
}

void cube_edge_adjacencies_test()
{

  //get all the curves of the cube
  MEntVector curves;
  mk->get_entities_by_dimension( 1, curves); 
  
  for(MEntVector::iterator i = curves.begin(); i != curves.end(); i++)
    {

      //get the iMesh set handle for the curve
      iMesh::EntitySetHandle sh = IBSH( (*i)->mesh_handle() ); 

      //get the curve edges
      std::vector<iMesh::EntityHandle> curve_edges; 
      mk->imesh_instance()->getEntities( sh, iBase_EDGE, iMesh_LINE_SEGMENT, curve_edges); 
      
      //for a cube there should only be one edge per curve
      CHECK( 1 == int(curve_edges.size()) );

      //check that each edge is adjacent to no more than 2 triangles 
      for(unsigned int i = 0; i < curve_edges.size(); i++) 
	{
	  
	  //get the adjacent triangles to the curve edge 
	  std::vector<iMesh::EntityHandle> adj_tris; 
	  mk->imesh_instance()->getEntAdj( curve_edges[i], iBase_FACE, adj_tris ); 
	  
	  //check that the entities returned are triangles 
	  for( unsigned int j = 0; j < adj_tris.size(); j++)
	    {
	      iMesh::EntityTopology topo; 
	      mk->imesh_instance()->getEntTopo( adj_tris[j], topo); 
	      CHECK( iMesh_TRIANGLE == topo ); 
	    }

	  //check that there are no more than two triangles adjacent to this edge
	  CHECK( 2 <= int(adj_tris.size()) ); 

	} //end curve edges loop
    } // end curves loop

}


int count_topo_for_dim( int dim, iMesh::EntityTopology topo)
{

  MEntVector ents;
  mk->get_entities_by_dimension(dim, ents);

  int counter = 0;

  for(unsigned int i = 0; i < ents.size(); i++)
    {
      int temp = 0;
      iMesh::EntitySetHandle set_handle = IBSH(ents[i]->mesh_handle());
      mk->imesh_instance()->getNumOfTopo(set_handle, topo, temp);

      counter += temp;
      temp = 0;
    }

  return counter;
}

