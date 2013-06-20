/** \file test_meshoptemplate.cpp
 *
 * Test MeshOpTemplate
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/EBMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOpTemplate.hpp"
#include "meshkit/VertexMesher.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/CAMALPaver.hpp"
#include "meshkit/OneToOneSwept.hpp"
#include "meshkit/SizingFunctionVar.hpp"
#include <iostream>
using std::cout;
using std::endl;

using namespace MeshKit;

#include "TestUtil.hpp"

MKCore *mk;

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
  // Add mot to the graph
  MeshOpTemplate *mot = (MeshOpTemplate*) mk->construct_meshop("MeshOpTemplate", surfs);
  // Initialize if any?
  mot->set_name("MeshOpTemplate");
  // Run it
  mk->setup_and_execute();
  // Save Pre-liminary Geometry
  mk->save_geometry("mot1.sat");
  // Populate Entities
  mk->populate_model_ents();
  // Remove mot from the graph (delete or clear both work)
  delete mot;
  
//  mk->clear_graph();

// Mesh the brick!
  // GET THE SURFACE!
  mk->get_entities_by_dimension(2, surfs);
  /*for (int i = 0; i < 5; i++) {
    printf("Doing GEBD(%d)\n", i);
    mk->get_entities_by_dimension(i, surfs);
    if (surfs.size())
      printf("We have entities!\n");
  }*/

  // Surface Mesh 1 face

  printf("construct%d\n", (int)surfs.size());
  surfs.resize(1);
  TFIMapping *cp = (TFIMapping*) mk->construct_meshop("TFIMapping", surfs);
  SizingFunction esize(mk, 3, .5);
  surfs[0]->sizing_function_index(esize.core_index());
  printf("setup\n");
  mk->setup();
  printf("execute\n");
  mk->execute();
  printf("after\n");
  mk->clear_graph();
  mk->populate_model_ents();



  printf("construct\n");
  mk->get_entities_by_dimension(3, vols);
  ModelEnt *some_vol = (*vols.rbegin());
  OneToOneSwept *oto = (OneToOneSwept*) mk->construct_meshop("OneToOneSwept", vols);
  oto->SetSourceSurface(0);
  oto->SetTargetSurface(1);
  SizingFunction wsize(mk, 3, -1);
  some_vol->sizing_function_index(wsize.core_index());
  printf("setup\n");
  mk->setup();
  printf("execute\n");
  mk->execute();
  printf("after\n");

  


  if(1) {
    #ifdef HAVE_ACIS
    mk->save_geometry("mot2.sat");
    mk->save_mesh("motmesh.exo");
    #elif defined(HAVE_OCC)
    mk->save_geometry("meshoptemplate.stp");
    mk->save_mesh("idunno.whatevs");
    #endif
  }

  // delete the mot instance
  mk->clear_graph();
}

