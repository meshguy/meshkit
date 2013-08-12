/*!
\example example_camaltriadvance.cpp

\section camaltriadvance_cpp_title <pretty-name-of-this-file>

\subsection camaltriadvance_cpp_in Input
\image html camaltriadvance.in.jpg
There is no input.

\subsection camaltriadvance_cpp_out Output
\image html camaltriadvance.out.jpg

\subsection camaltriadvance_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection camaltriadvance_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/CAMALTriAdvance.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;



MKCore *mk = NULL;
bool save_mesh = false;

void holysurf_test();
void singleholesurf_test();
void singleholesurfimprinted_test();
void square_test();
void mesh_test(std::string filebase);
void var_size_test();

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
  

  holysurf_test();
  singleholesurf_test();
  singleholesurfimprinted_test();
  square_test();
  var_size_test();
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

void square_test()
{
  mesh_test("squaresurf");
}
void mesh_test(std::string filebase)
{
  std::string file_name = TestDir + "/" + filebase + extension;
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector dum, surfs;
  mk->get_entities_by_dimension(2, dum);
  surfs.push_back(*dum.rbegin());
  mk->construct_meshop("CAMALTriAdvance", surfs);

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

void var_size_test()
{
  std::string file_name = TestDir + "/" + std::string("squaresurf") + extension;
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector dum, surfs;
  mk->get_entities_by_dimension(2, dum);
  surfs.push_back(*dum.rbegin());
  mk->construct_meshop("CAMALTriAdvance", surfs);

    // make a sizing function and set it on the surface
  // make a sizing function and set it on the surface
  SizingFunctionVar * svar = new SizingFunctionVar(mk, -1, 0.1);

  // these could be read from a file, or something
  double point0[3] = {0, 0, 0};
  double coeffs[4] = {0.05, 0.05, 0.05, 0.1};
  svar->set_linear_coeff(point0, coeffs);

  ModelEnt *this_surf = (*surfs.rbegin());
  this_surf->sizing_function_index(svar->core_index());
  // how to know if we use a var size during setup of the Edge Mesher?

    // mesh the surface, by calling execute
  mk->setup_and_execute();

    // report the number of tris
  moab::Range tris;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 2, tris);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  std::cout << tris.size() << " tris generated." << std::endl;

  if (save_mesh) {
        // output mesh
    std::string outfile = std::string("squaresurf") + std::string(".vtk");
    moab::EntityHandle out_set = surfs[0]->mesh_handle();
    rval = mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL, &out_set, 1);
    MBERRCHK(rval, mk->moab_instance());
  }

  mk->clear_graph();
}
