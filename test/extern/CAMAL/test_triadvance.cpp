/** \file test_triadvance.cpp
 *
 * Test the CAMALTriAdvance for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/CAMALTriAdvance.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;
bool save_mesh = false;

void holysurf_test();
void singleholesurf_test();
void singleholesurfimprinted_test();
void mesh_test(std::string filebase);

#ifdef HAVE_ACIS
std::string extension = ".sat";
#elif HAVE_OCC
std::string extension = ".stp";
#endif

int main(int argc, char **argv) 
{
  
    // start up MK and load the geometry
  mk = new MKCore;
  int num_fail = 0;
  
  if (argc == 2) save_mesh = true;
  
  num_fail += RUN_TEST(holysurf_test);
  num_fail += RUN_TEST(singleholesurf_test);
  num_fail += RUN_TEST(singleholesurfimprinted_test);
  return num_fail;
}

void holysurf_test() 
{
  mesh_test("holysurf");
}

void singleholesurf_test() 
{
  mesh_test("singleholesurf");
}

void singleholesurfimprinted_test()
{
  mesh_test("singleholesurfimprinted");
}

void mesh_test(std::string filebase)
{
  std::string file_name = TestDir + "/" + filebase + extension;
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector dum, surfs;
  mk->get_entities_by_dimension(2, dum);
  surfs.push_back(*dum.rbegin());
  CAMALTriAdvance *tm = (CAMALTriAdvance*) mk->construct_meshop("CAMALTriAdvance", surfs);

    // make a sizing function and set it on the surface
  SizingFunction esize(mk, -1, 0.25);
  surfs[0]->sizing_function_index(esize.core_index());
  
    // mesh the surface, by calling execute
  mk->setup_and_execute();

    // report the number of tris
  moab::Range tris;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() << " tris generated." << std::endl;

  if (save_mesh) {
        // output mesh
    std::string outfile = filebase + std::string(".vtk");
    moab::EntityHandle out_set = surfs[0]->mesh_handle();
    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL, &out_set, 1);
    MBERRCHK(rval, mk->moab_instance());
  }

  mk->clear_graph();
}

