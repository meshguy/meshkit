/** \file test_edgemesh.cpp
 *
 * Test the EdgeMesher for a few challenging examples.
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/VertexMesher.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk = NULL;

void edgemesh_hole();
void edgemesh_square();
void edgemesh_brick();

int main(int argc, char **argv) 
{
  
    // start up MK and load the geometry
  mk = new MKCore();

  int num_fail = 0;
  
  num_fail += RUN_TEST(edgemesh_hole);
  num_fail += RUN_TEST(edgemesh_square);
  num_fail += RUN_TEST(edgemesh_brick);
  
}

void edgemesh_hole() 
{
  std::string file_name = TestDir + "/holysurf.sat";
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
  EdgeMesher *em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves);

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

void edgemesh_square() 
{
  std::string file_name = TestDir + "/squaresurf.sat";
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
  EdgeMesher *em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves);

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

void edgemesh_brick() 
{
  std::string file_name = TestDir + "/brick.sat";
  mk->load_geometry(file_name.c_str());

    // get the vol, surfs, curves
  MEntVector vols, surfs, curves, loops;
  mk->get_entities_by_dimension(3, vols);
  ModelEnt *this_vol = (*vols.rbegin());
  this_vol->get_adjacencies(2, surfs);
  this_vol->get_adjacencies(1, curves);

    // mesh all the edges
  EdgeMesher *em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves);
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

