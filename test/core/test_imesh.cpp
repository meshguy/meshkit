//** \file test_imesh.cpp

#include "meshkit/MKCore.hpp"
#include <stdlib.h>

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk; 

void test_Adj(); 
void test_setAdjTable();


int main(int argc, char **argv) 
{


  mk = new MKCore(); 
  int num_fail = 0; 
  #ifdef HAVE_MOAB 
  num_fail+= RUN_TEST(test_Adj); 
  num_fail+= RUN_TEST(test_setAdjTable);
  #endif
  return num_fail; 
}

void test_Adj()
{

  //create a new tet in the mesh
  iMesh::EntityHandle v1, v2, v3, v4; 

  iMesh::AdjTableType adjTable =  mk->imesh_instance()->getAdjTable();

  //setup the adjacency table values such that intermediate dimension entities are created
  int table_vals[16] = { 1, 0, 0, 0, 
                         0, 1, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 0, 0};
  int table_size = 16;     
  mk->imesh_instance()->setAdjTable( table_vals, table_size ); 


  //create the vertices
  mk->imesh_instance()->createVtx( 0, 0, 0, v1); 
  mk->imesh_instance()->createVtx( 1, 0, 0, v2); 
  mk->imesh_instance()->createVtx( 0, 1, 0, v3);
  mk->imesh_instance()->createVtx( 0, 0, 1, v4); 

  std::vector<iMesh::EntityHandle> tet_tris(4);
  //generate the triangles using the proper vertices 
  iMesh::EntityHandle verts[] = {v1, v2, v3, 
                                 v2, v3, v4,
                                 v3, v4, v1,
                                 v1, v3, v4};
  mk->imesh_instance()->createEntArr( iMesh_TRIANGLE, verts, 12, &tet_tris[0]); 
  
 
  //get all of the triangles in the model 
  std::vector<iMesh::EntityHandle> tris;
  mk->imesh_instance()->getEntities( 0, iBase_FACE, iMesh_TRIANGLE, tris); 

  //make sure we get the correct number of triangles
  int expected_num_of_tris = 4; 
  CHECK( expected_num_of_tris == int(tris.size())); 

  //now start checking the adjacencies
  int expected_num_of_adj = 3; 

  //check first adj
  std::vector<iMesh::EntityHandle>::iterator i; 
  for( i = tris.begin(); i != tris.end(); i++)
    {
 
      //each triangle should be adjacent to three vertices and three edges
      std::vector<iMesh::EntityHandle> adj; 
      mk->imesh_instance()->getEntAdj( *i, iBase_VERTEX, adj); 
      CHECK( expected_num_of_adj == int(adj.size()) ); 
      
      adj.clear(); 
      mk->imesh_instance()->getEntAdj( *i, iBase_EDGE, adj); 
      CHECK( expected_num_of_adj == int(adj.size()) ); 
      
    }

  //check 2nd adj
  for( i = tris.begin(); i != tris.end(); i++)
    {
      std::vector<iMesh::EntityHandle> adj; 
      //each triangle should be adjacent to three other triangles via the edges 
      adj.clear(); 
      mk->imesh_instance()->getEnt2ndAdj( *i, iBase_EDGE, iBase_FACE, adj);

      //make sure these are all triangles
      std::vector<iMesh::EntityHandle>::iterator j; 
      for( j = adj.begin(); j != adj.end(); j++) 
	{
	  iMesh::EntityTopology topo; 
	  mk->imesh_instance()->getEntTopo( *j, topo); 
	  CHECK( iMesh_TRIANGLE == topo ); 
	}

      CHECK( expected_num_of_adj == int(adj.size()) );

    }

}

void test_setAdjTable()
{


  //setup the test adjacency table
  int table_vals[16] = { 1, 1, 1, 1, 
                         1, 1, 1, 1,
                         1, 1, 1, 1,
                         1, 1, 1, 1};
  int table_size = 16;     

  //set the adjacency table for the iMesh instance
  mk->imesh_instance()->setAdjTable( table_vals, table_size ); 
  
  //get the adjacancy table back and make sure there aren't any zero-values
  iMesh::AdjTableType adjTable = mk->imesh_instance()->getAdjTable();

  for( unsigned int i = 0; i < 4; i++)
    for(unsigned int j = 0; j < 4; j++)
      {
	{
	  if( i == j ) CHECK( iBase_AVAILABLE == adjTable[i][j] );
	  else CHECK( iBase_ALL_ORDER_1 == adjTable[i][j] );
	}
      }


}
  
  
