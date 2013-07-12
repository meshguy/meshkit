/*!
\example example_basic.cpp

\section title Basic Example

\subsection inf Misc. Information
\author Brett Rhodes
\date 7-12-2013
\bug Currently requires multiple setup/execute cycles.
\warning This example is not currently complete.

\subsection goal Goal
The goal of this is example is to start from no input and create a meshed geometry with proper Neumann and material sets.

\subsection cw Code Walkthrough
Code is subject to change, this will be written at a later point.

\subsection in Input
There is no input.

\subsection out Output
\image html example_basic.out.jpg

\subsection src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/MeshOpTemplate.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/EdgeMesher.hpp"
#include <algorithm>

//#include "meshkit/EBMesher.hpp"
//#include "meshkit/ModelEnt.hpp"
//#include "meshkit/VertexMesher.hpp"
//#include "meshkit/EdgeMesher.hpp"
//#include "meshkit/SizingFunctionVar.hpp"
//#include "meshkit/CAMALPaver.hpp"

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
  mot->set_name("MeshOpTemplate");
  Vector<3> a;
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

//  mk->imesh_instance()->createTag("TAG_SET_A", 1, iBase_INTEGER, mesh_tag);
/*
  // Tag things!
  int err;
  moab::Interface * mb = mk->moab_instance();
//  iBase_TagHandle mesh_tag;
  moab::Tag neu, mat;
  moab::Range tagged_ents;
  std::vector<moab::EntityHandle> ent_list;
  mb->tag_get_handle("MATERIAL_SET", mat);
  mb->get_entities_by_dimension(mb->get_root_set(), 3, ent_list, true);
  for (int i = 0; (unsigned int)i < ent_list.size(); i++) {
    iMesh_setEntSetIntData(mk->imesh_instance()->instance(), (iBase_EntitySetHandle)ent_list[i],
                            (iBase_TagHandle)mat, 1, &err);
    if (!err)
      printf("i = %d\n", i);
  }

  ent_list.resize(0);
  mb->tag_get_handle("NEUMANN_SET", neu);
  mb->get_entities_by_dimension(mb->get_root_set(), 2, ent_list, true);
  mb->tag_set_data(neu, ent_list, ent_list.size(), );
  for (int i = 0; (unsigned int)i < ent_list.size(); i++) {
    iMesh_setEntSetIntData(mk->imesh_instance()->instance(), iBase_EntitySetHandle(ent_list[i]),
                            (iBase_TagHandle)neu, 1, &err);
    if (!err)
      printf("i = %d\n", i);
  }

  mk->moab_instance()->tag_get_handle("NEUMANN_SET", &neu);
  mk->moab_instance()->get_entities_by_dimension(2, ent_list);
  mk->moab_instance()->iMesh_setEntSetIntData(mk->imesh_instance(), ent_list, neu, 1, new int);
  */

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
