/*! 
\example example_basic.cpp

\section example_basic_cpp_title basic Example

The goal of this is example is to start from no input and create a meshed geometry with proper Neumann and material sets.

\subsection example_basic_cpp_in Input
(none)

\subsection example_basic_cpp_out Output
\image html basic.out.jpg "(description of image)"

\subsection example_basic_cpp_inf Misc. Information
\author Brett Rhodes
\date 7-16-2013
\bug Currently requires multiple setup/execute cycles.
\warning This example is not currently complete.

\subsection example_basic_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpTemplate.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/OneToOneSwept.hpp"

using namespace MeshKit;

MKCore *mk;
bool save_mesh = true;

void test_mesh_op_template();

int main(int argc, char **argv)
{
  mk = new MKCore();
  int num_fail = 0;

  test_mesh_op_template();

  delete mk;
  return num_fail;
}

void test_mesh_op_template()
{
  MEntVector all_curves, all_surfs, the_vol, tfi_surface;
  int source_surface = 0; // These are
  int target_surface = 1; //   chosen by magic!
                                                                               // | Geometry Op,
  // Make the brick!                                                           // v Don't need this
  MeshOpTemplate *mot = (MeshOpTemplate*) mk->construct_meshop("MeshOpTemplate", MEntVector());
  printf("Got here fine\n");
  mot->set_name("MeshOpTemplate");
  //Vector<3> a;
  mk->setup_and_execute();
  mk->populate_model_ents();
  mk->clear_graph();

  // Get entity data
  mk->get_entities_by_dimension(1, all_curves);
  mk->get_entities_by_dimension(2, all_surfs);
  mk->get_entities_by_dimension(3, the_vol);

  // Size some things!
  ModelEnt *some_vol = (*the_vol.rbegin());
  SizingFunction wsize(mk, 4, -1);
  some_vol->sizing_function_index(wsize.core_index());

  // Mesh Surface
  (tfi_surface = all_surfs).resize(1);
  TFIMapping *tfi = (TFIMapping*) mk->construct_meshop("TFIMapping", tfi_surface);
  tfi->set_name("TFI");

  // Sweep Surface to Surface
  OneToOneSwept *oto = (OneToOneSwept*) mk->construct_meshop("OneToOneSwept", the_vol);
  oto->set_name("OTO");
  oto->SetSourceSurface(source_surface);
  oto->SetTargetSurface(target_surface);

  // Do work!
  printf("Pre-Setup\n");
  mk->setup();
  printf("Pre-Execute\n");
  mk->execute();
  printf("Pre-Clear\n");
  mk->clear_graph();

  // Some declarations for tags
  moab::Interface *mb = mk->moab_instance();
  moab::Tag matTag, neuTag;
  moab::Range matEnts, neuEnts;
  moab::EntityHandle matEntSet, neuEntSet;
  int matVal = 2, neuVal = 4; // these values are arbitrary

// Tag the material set
  // Collect the entites into a set
  mb->create_meshset(moab::MESHSET_SET, matEntSet);
  mb->get_entities_by_dimension(mb->get_root_set(), 3, matEnts, true);
  mb->add_entities(matEntSet, matEnts);
  
  // Grab a tag and attach it to the set
  mb->tag_get_handle("MATERIAL_SET", matTag);
  mb->tag_set_data(matTag, &matEntSet, 1, (void*) &matVal);

// The same thing for the Neumann set!
  mb->create_meshset(moab::MESHSET_SET, neuEntSet);
  mb->get_entities_by_dimension(mb->get_root_set(), 2, neuEnts, true);
  mb->add_entities(neuEntSet, neuEnts);

  mb->tag_get_handle("NEUMANN_SET", neuTag);
  mb->tag_set_data(neuTag, &neuEntSet, 1, (void*) &neuVal);

  // Possibly output (as files) our hard work!
  if(save_mesh) {
    #ifdef HAVE_ACIS
      mk->save_geometry("un_meshed_brick.sat");
      mk->save_mesh("meshed_brick.exo");
    #elif defined(HAVE_OCC)
      mk->save_geometry("un_meshed_brick.stp");
      mk->save_mesh("meshed_brick.exo");
    #endif
  }

  return;
}
