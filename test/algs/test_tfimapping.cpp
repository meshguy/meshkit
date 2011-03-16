/** \file test_edgemesh.cpp
 *
 * Test the EdgeMesher for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/EdgeMesher.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

void test_TFImapping();
void test_TFImappingcubit();

int main(int argc, char **argv) 
{
  
    // start up MK and load the geometry
  mk = new MKCore();

  int num_fail = 0;
  
  num_fail += RUN_TEST(test_TFImapping);
  num_fail += RUN_TEST(test_TFImappingcubit);
  
}

void test_TFImappingcubit()
{
	std::string file_name = TestDir + "/square.cub";
  	mk->load_geometry(file_name.c_str());

	mk->load_mesh(file_name.c_str());

	// populate mesh to relate geometry entities and mesh sets
  	mk->populate_mesh();

	MEntVector surfs, curves, loops;
	mk->get_entities_by_dimension(2, surfs);
	CHECK_EQUAL(1, (int)surfs.size());

	//now, do the TFIMapping
	TFIMapping *tm = (TFIMapping*)mk->construct_meshop("TFIMapping", surfs);
	mk->setup_and_execute();


	mk->clear_graph();

	

	
}

void test_TFImapping() 
{
  std::string file_name = TestDir + "/squaresurf.sat";
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  CHECK_EQUAL(1, (int)surfs.size());
  
    // make an edge mesher
  mk->get_entities_by_dimension(1, curves);
  //test there are 4 edges bounding the surface
  CHECK_EQUAL(4, (int)curves.size());
  EdgeMesher *em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves);

    // make a sizing function and set it on the surface
  SizingFunction esize(mk, 10, -1);
  surfs[0]->sizing_function_index(esize.core_index());
  
  // mesh the edges, by calling execute
  mk->setup_and_execute();

    // make sure we got the right number of edges
  moab::Range edges;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 1, edges);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(40, (int)edges.size());

  std::cout << "we are done with edge mesher" << std::endl;
  //ok, we are done with edge mesher
  //now, do the TFIMapping
  TFIMapping *tm = (TFIMapping*)mk->construct_meshop("TFIMapping", surfs);
  mk->setup_and_execute();
  
  mk->populate_mesh();

  //check whether we got the right number of quads after TFIMapping
  moab::Range faces;
  rval = mk->moab_instance()->get_entities_by_dimension(0, 2, faces);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(40, (int)faces.size());

  //output the mesh to vtk file
  mk->save_mesh("TFIMapping.vtk");

    // clean up
  delete em;
  delete tm;
  //delete mk->vertex_mesher();
  mk->clear_graph();
}
