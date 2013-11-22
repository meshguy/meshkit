/*!
\example example_intassign.cpp

\section example_intassign_cpp_title Interval Assignment Example

\subsection example_intassign_cpp_in Input
Loads a simple quadface geometry from data directory.

\subsection example_intassign_cpp_out Output
Mesh file from TFI Mapping and Interval Assigment.

\warning Needs more work.

\subsection example_intassign_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"

#include "meshkit/IAInterface.hpp"
#include "meshkit/IAVariable.hpp"
#include "meshkit/TFIMapping.hpp"

#include <stdio.h>
#include <iostream>
#ifdef HAVE_ACIS
#define INTASSIGN_TEST_FILE "quadface.sat"
#elif defined(HAVE_OCC)
#define INTASSIGN_TEST_FILE "quadface.stp"
#endif

MeshKit::MKCore *mk;

MeshKit::IAInterface *new_ia_interface()
{
  return
      (MeshKit::IAInterface*) mk->construct_meshop("IntervalAssignment");
}

void mapping_test() 
{

}

int main(int argv, char* argc[])
{
  // currently unable to create more than one mk called IntervalAssignment
  mk = new MeshKit::MKCore();
  MeshKit::IAInterface *ia_interface = new_ia_interface();
  ia_interface->destroy_data();

  std::string file_name = (std::string) MESH_DIR +  "/" + INTASSIGN_TEST_FILE;
  printf("opening %s\n", file_name.c_str());
  mk->load_geometry_mesh(file_name.c_str(), NULL);

  //check the number of geometrical edges
  MeshKit::MEntVector surfs, loops;
  mk->get_entities_by_dimension(2, surfs);
  MeshKit::ModelEnt *this_surf = (*surfs.rbegin());

  MeshKit::MEntVector curves;
  std::vector<int> senses, loop_sizes;
  this_surf->boundary(1, curves, &senses, &loop_sizes);

  MeshKit::SizingFunction esize(mk, 3, -1);
  surfs[0]->sizing_function_index(esize.core_index());

  MeshKit::MEntVector side1, side2;
  side1.push_back(curves[0]); side2.push_back(curves[2]);
  ia_interface->constrain_sum_equal(ia_interface->make_constraint_group(side1),
                                    ia_interface->make_constraint_group(side2));
  side1.clear(); side2.clear();
  side1.push_back(curves[1]); side2.push_back(curves[3]);
  ia_interface->constrain_sum_equal(ia_interface->make_constraint_group(side1),
                                    ia_interface->make_constraint_group(side2));

  // if there are loops, and the loops have strictly less than 4 curves, then
  // ia_interface->constrain_sum_even( ia_interface->make_constraint_group(curves in loop) );

  //now, do the TFIMapping
  (MeshKit::TFIMapping*) mk->construct_meshop("TFIMapping", surfs);
  mk->setup_and_execute();
  mk->save_mesh("intassign.exo");
  return 0;
}
