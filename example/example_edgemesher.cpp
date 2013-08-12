/*!
\example example_edgemesher.cpp

\section edgemesher_cpp_title <pretty-name-of-this-file>

\subsection edgemesher_cpp_in Input
\image html edgemesher.in.jpg
There is no input.

\subsection edgemesher_cpp_out Output
\image html edgemesher.out.jpg

\subsection edgemesher_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection edgemesher_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/edgemesherer.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/VertexMesher.hpp"

using namespace MeshKit;



MKCore *mk = NULL;

void edgemesher_hole();
void edgemesher_square();
void edgemesher_brick();
void edgemesher_var();
void edgemesher_dual();

#ifdef HAVE_ACIS
std::string extension = ".sat";
#elif HAVE_OCC
std::string extension = ".stp";
#endif

int main(int argc, char **argv) 
{
  
    // start up MK and load the geometry
  mk = new MKCore();

  int num_fail = 0;
  edgemesher_hole();
  edgemesher_square();
  edgemesher_brick();
  edgemesher_var();
  edgemesher_dual();
#if HAVE_OCC
  return 0;
#else
  return num_fail;
#endif
}

void edgemesher_hole() 
{
  std::string file_name = TestDir + "/holysurf" + extension;
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  CHECK_EQUAL(1, (int)surfs.size());

    // test getting loops
  std::vector<int> senses, loop_sizes;
  surfs[0]->boundary(0, loops, &senses, &loop_sizes);
  CHECK_EQUAL(4, (int)loop_sizes.size());
  
    // add up the loop sizes, should add to 8
  unsigned int l = 0;
  for (unsigned int i = 0; i < loop_sizes.size(); i++) l += loop_sizes[i];
  CHECK_EQUAL(8, (int)l);
  
    // make an edge mesher
  mk->get_entities_by_dimension(1, curves);
  edgemesherer *em = (edgemesherer*) mk->construct_meshop("edgemesherer", curves);

    // make a sizing function and set it on the surface
  SizingFunction esize(mk, -1, 0.25);
  surfs[0]->sizing_function_index(esize.core_index());
  
    // mesh the edges, by calling execute
  mk->setup_and_execute();

    // make sure we got the right number of edges
  moab::Range edges;
  moab::ErrorCode rval = mk->moab_instance()->get_entities_by_dimension(0, 1, edges);
  CHECK_EQUAL(moab::MB_SUCCESS, rval);
  CHECK_EQUAL(79, (int)edges.size());

    // clean up
  delete em;
  delete mk->vertex_mesher();
}

void edgemesher_square() 
{
  std::string file_name = TestDir + "/squaresurf" + extension;
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  ModelEnt *this_surf = (*surfs.rbegin());

    // test getting loops
  std::vector<int> senses, loop_sizes;
  this_surf->boundary(0, loops, &senses, &loop_sizes);
  CHECK_EQUAL(1, (int)loop_sizes.size());
  
    // make an edge mesher
  this_surf->get_adjacencies(1, curves);
  mk->construct_meshop("edgemesherer", curves);

    // make a sizing function and set it on the surface
  SizingFunction esize(mk, 1, 1.0);
  this_surf->sizing_function_index(esize.core_index());
  
    // mesh the edges, by calling execute
  mk->setup_and_execute();

  senses.clear();
  loop_sizes.clear();
  std::vector<moab::EntityHandle> mloops;
  this_surf->boundary(0, mloops, &senses, &loop_sizes);
  CHECK_EQUAL(4, (int)mloops.size());

    // print the vertex positions
  std::vector<double> coords(12);
  moab::ErrorCode rval = mk->moab_instance()->get_coords(&mloops[0], mloops.size(), &coords[0]);
  MBERRCHK(rval, mk->moab_instance());
//  std::cout << "Vertex positions:" << std::endl;
//  for (unsigned int i = 0; i < 4; i++)
//    std::cout << coords[3*i] << ", " << coords[3*i+1] << ", " << coords[3*i+2] << std::endl;

    // compute the cross product vector and print that
  double tmp[] = {0.0, 0.0, 1.0};
  Vector<3> p0(&coords[0]), p1(&coords[3]), p2(&coords[6]), unitz(tmp);
  p0 -= p1;
  p2 -= p1;
  p1 = vector_product(p2, p0);
  double dot = inner_product(unitz, p1);
  CHECK_REAL_EQUAL(1.0, dot, 1.0e-6);
  
  mk->clear_graph();
}

void edgemesher_brick() 
{
  std::string file_name = TestDir + "/brick" + extension;
  mk->load_geometry(file_name.c_str());

    // get the vol, surfs, curves
  MEntVector vols, surfs, curves, loops;
  mk->get_entities_by_dimension(3, vols);
  ModelEnt *this_vol = (*vols.rbegin());
  this_vol->get_adjacencies(2, surfs);
  this_vol->get_adjacencies(1, curves);

    // mesh all the edges
  mk->construct_meshop("edgemesherer", curves);
  SizingFunction esize(mk, 1, 1.0);
  this_vol->sizing_function_index(esize.core_index());
  mk->setup_and_execute();

    // test loop directionality for each surface
  for (MEntVector::iterator vit = surfs.begin(); vit != surfs.end(); vit++) {
    
    std::vector<int> senses, loop_sizes;
    std::vector<moab::EntityHandle> mloops;
    (*vit)->boundary(0, mloops, &senses, &loop_sizes);
    CHECK_EQUAL(4, (int)mloops.size());

      // get the vertex positions
    std::vector<double> coords(12);
    moab::ErrorCode rval = mk->moab_instance()->get_coords(&mloops[0], mloops.size(), &coords[0]);
    MBERRCHK(rval, mk->moab_instance());

      // compute the cross product vector and compare it to the surface normal
    Vector<3> p0(&coords[0]), p1(&coords[3]), p2(&coords[6]), p3;
    p0 -= p1;
    p2 -= p1;
    p1 = vector_product(p2, p0);
    (*vit)->evaluate(coords[0], coords[1], coords[2], NULL, p3.data());
    double dot = inner_product(p3, p1);
    CHECK_REAL_EQUAL(1.0, dot, 1.0e-6);
  }
  
  mk->clear_graph();
}
void edgemesher_var()
{
  // delete existing mesh
  mk->moab_instance()->delete_mesh();

  std::string file_name = TestDir + "/squaresurf" + extension;
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  ModelEnt *this_surf = (*surfs.rbegin());

    // test getting loops
  std::vector<int> senses, loop_sizes;
  this_surf->boundary(0, loops, &senses, &loop_sizes);
  CHECK_EQUAL(1, (int)loop_sizes.size());

    // make an edge mesher
  this_surf->get_adjacencies(1, curves);
  edgemesherer * emesher = (edgemesherer *)mk->construct_meshop("edgemesherer", curves);

  emesher->set_edge_scheme(edgemesherer::VARIABLE);

    // make a sizing function and set it on the surface
  SizingFunctionVar * svar = new SizingFunctionVar(mk, -1, 0.1);

  // these could be read from a file, or something
  double point0[3] = {0, 0, 0};
  double coeffs[4] = {0.05, 0.05, 0.05, 0.1};
  svar->set_linear_coeff(point0, coeffs);

  this_surf->sizing_function_index(svar->core_index());

  mk->setup_and_execute();

  //mk->moab_instance()->write_file("var1.vtk");
  mk->clear_graph();

}
void edgemesher_dual()
{
  // delete existing mesh
  mk->moab_instance()->delete_mesh();

  std::string file_name = TestDir + "/squaresurf" + extension;
  mk->load_geometry(file_name.c_str());

    // get the surface
  MEntVector surfs, curves, loops;
  mk->get_entities_by_dimension(2, surfs);
  ModelEnt *this_surf = (*surfs.rbegin());

    // test getting loops
  std::vector<int> senses, loop_sizes;
  this_surf->boundary(0, loops, &senses, &loop_sizes);
  CHECK_EQUAL(1, (int)loop_sizes.size());

    // make an edge mesher
  this_surf->get_adjacencies(1, curves);
  edgemesherer * emesher = (edgemesherer *)mk->construct_meshop("edgemesherer", curves);

  emesher->set_edge_scheme(edgemesherer::DUAL);

    // make a sizing function and set it on the surface
  SizingFunction * svar = new SizingFunction(mk, 10, -1);

  emesher->set_ratio(1.3);

  this_surf->sizing_function_index(svar->core_index());

  mk->setup_and_execute();

  //mk->moab_instance()->write_file("dual1.1.vtk");
  mk->clear_graph();

}
