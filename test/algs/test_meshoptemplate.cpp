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

//#include "meshkit/EBMesher.hpp"
//#include "meshkit/ModelEnt.hpp"
//#include "meshkit/VertexMesher.hpp"
//#include "meshkit/EdgeMesher.hpp"
//#include "meshkit/SizingFunctionVar.hpp"
//#include "meshkit/CAMALPaver.hpp"

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk;
bool save_mesh = true;
void test_mesh_op_template();

int main(int argc, char **argv)
{
  mk = new MKCore();
  int num_fail = 0;

  num_fail += RUN_TEST(test_mesh_op_template);

  delete mk;
  return num_fail;
}

void test_mesh_op_template()
{
  MEntVector vols, surfs;

  // Make the brick!
  MeshOpTemplate *mot = (MeshOpTemplate*) mk->construct_meshop("MeshOpTemplate", surfs);
  mot->set_name("MeshOpTemplate");
  Vector<3> a;
  mk->setup_and_execute();
  mk->populate_model_ents();
  mk->clear_graph();

  // Mesh 1 Surface
  mk->get_entities_by_dimension(2, surfs);
  mk->get_entities_by_dimension(3, vols);
  surfs.resize(1);
  TFIMapping *tfi = (TFIMapping*) mk->construct_meshop("TFIMapping", surfs);
  tfi->set_name("TFIMapping");

  // Sweep
  OneToOneSwept *oto = (OneToOneSwept*) mk->construct_meshop("OneToOneSwept", vols);
  oto->set_name("OneToOneSwept");
  oto->SetSourceSurface(0);
  oto->SetTargetSurface(1);
  ModelEnt *some_vol = (*vols.rbegin());
  SizingFunction wsize(mk, 3, -1);
  some_vol->sizing_function_index(wsize.core_index());

  // Do work!
  mk->setup_and_execute();
  mk->clear_graph();

  if(save_mesh) {
    #ifdef HAVE_ACIS
    mk->save_geometry("un_meshed_brick.sat");
    mk->save_mesh("meshed_brick.exo");
    #elif defined(HAVE_OCC)
    mk->save_geometry("meshoptemplate.stp");
    mk->save_mesh("meshed_brick.exo");
//    mk->save_mesh("idunno.whatevs");
    #endif
  }

  return;
}

