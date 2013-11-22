/** \file test_SubMapping.cpp
 *
 * Test the SubMapping for a few examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/SubMapping.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/SubMapping.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

void test_SubMappingcubit();

int main(int argc, char **argv) 
{
  
    // start up MK and load the geometry
  mk = new MKCore();

  int num_fail = 0;
  
  //num_fail += RUN_TEST(test_TFImapping);
  
  num_fail += RUN_TEST(test_SubMappingcubit);


 #if HAVE_OCC
  return 0;
#else
  return num_fail;
#endif 
}

void test_SubMappingcubit()
{
    std::string file_name = TestDir + "/submapping_simple_example.cub";
    mk->load_geometry_mesh(file_name.c_str(), file_name.c_str());

	//check the number of geometrical edges
	MEntVector surfs, curves, loops;
  	mk->get_entities_by_dimension(2, surfs);
	ModelEnt *this_surf = (*surfs.rbegin());

	this_surf->get_adjacencies(1, curves);

	//CHECK_EQUAL(4, (int)curves.size());
	
	//check the number of mesh line segments
	moab::Range edges;
	moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 1, edges);
	CHECK_EQUAL(moab::MB_SUCCESS, rval);
	//CHECK_EQUAL(40, (int)edges.size());

	SizingFunction swSize(mk, -1, 0.5);
	this_surf->sizing_function_index(swSize.core_index());
	//now, do the SubMapping
	SubMapping *tm = (SubMapping*)mk->construct_meshop("SubMapping", surfs);
	//tm->SetupMeshSize(5.0);

	mk->setup_and_execute();

	//check the number of quads
	moab::Range faces;
	rval = mk->moab_instance()->get_entities_by_dimension(0, 2, faces);
	CHECK_EQUAL(moab::MB_SUCCESS, rval);
	//CHECK_EQUAL(100, (int)faces.size());

	mk->save_mesh("submapping.vtk");

	delete tm;
	mk->clear_graph();
	
}

  
