/** \file test_meshoptemplate.cpp
 *
 * Test MeshOpTemplate
 *
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

