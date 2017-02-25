/** \file test_ngtetmesher.cpp
 *
 * Test the NGTetMesher for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/NGTetMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;
bool save_mesh = false;

void threeholecube_test();
void simpletet_test();
void mesh_test(std::string filebase, double size);

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
  
  // Test 1                                                                                                                                                                            
  num_fail += RUN_TEST(simpletet_test);
  // clear mk, geometry and mesh                                                                                                                                                       
  mk->delete_all();
  // Test 2                                                                                                                                                                            
  num_fail += RUN_TEST(threeholecube_test);
  
  return num_fail;
}

void simpletet_test() 
{
  mesh_test("simpletet", 0.25);
}

void threeholecube_test() 
{
  mesh_test("threeholecube", 0.06);
}

void mesh_test(std::string filebase, double size)
{
  std::string file_name = TestDir + "/" + filebase + extension;
  mk->load_geometry(file_name.c_str());

    // get the volume
  MEntVector dum, vols;
  mk->get_entities_by_dimension(3, dum);
  vols.push_back(*dum.rbegin());
  /*NGTetMesher *tm = (NGTetMesher*)*/ mk->construct_meshop("NGTetMesher", vols);

    // make a sizing function and set it on the surface
  SizingFunction esize(mk, -1, size);
  vols[0]->sizing_function_index(esize.core_index());
  
    // mesh the surface, by calling execute
  mk->setup_and_execute();

    // report the number of tets
  moab::Range tets;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 3, tets);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tets.size() << " tets generated." << std::endl;

  if (save_mesh) {
        // output mesh
    std::string outfile = filebase + std::string(".vtk");
    moab::EntityHandle out_set = vols[0]->mesh_handle();
    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL, &out_set, 1);
    MBERRCHK(rval, mk->moab_instance());
  }

  mk->clear_graph();
}

