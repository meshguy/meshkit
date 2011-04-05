#include "meshkit/MKCore.hpp"
#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/VertexMesher.hpp"


using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

void test_brick();


int main(int argc, char **argv) 
{
  
  // start up MK and load the geometry
  mk = new MKCore();

  int num_fail = 0;
  
  num_fail += RUN_TEST(test_brick);
  
}

void test_brick()
{
	std::string file_name = TestDir + "/BrickWithSrcMeshed1.cub";
	
	//load the geometry
	mk->load_geometry(file_name.c_str(), 0, 0, false);
	//load the mesh	
	mk->load_mesh(file_name.c_str());
	//populate the mesh and create the relation between the geometry and mesh
	mk->populate_mesh();

	// get the volumes
	MEntVector vols, surfs, curves, vertices;
	mk->get_entities_by_dimension(3, vols);
	
	std::cout << "Volume size = " << vols.size() << std::endl;
	
	ModelEnt *this_vol = (*vols.rbegin());

	// test getting surfaces
	this_vol->get_adjacencies(2, surfs);
	CHECK_EQUAL(6, (int)surfs.size());

	// test getting edges
	this_vol->get_adjacencies(1, curves);
	CHECK_EQUAL(12, (int)curves.size());

	// test getting vertices
	this_vol->get_adjacencies(0, vertices);
	CHECK_EQUAL(8, (int)vertices.size());

	//make a one-to-one sweeping
	OneToOneSwept *sw = (OneToOneSwept*) mk->construct_meshop("OneToOneSwept", vols);

	sw->SetSourceSurface(1);
	sw->SetTargetSurface(0);

	//set up the size
	SizingFunction swSize(mk, 10, -1);
	this_vol->sizing_function_index(swSize.core_index());

	//set up for the sweeping and sweep
	mk->setup_and_execute();

	//check the number of cells after OneToOneSwept
	moab::Range hex;
  	moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 3, hex);
  	CHECK_EQUAL(moab::MB_SUCCESS, rval);
  	CHECK_EQUAL(100, (int)hex.size());

	mk->save_mesh("OneToOneSwept.vtk");

	delete sw;
	delete mk->vertex_mesher();
	mk->clear_graph();
		
}

