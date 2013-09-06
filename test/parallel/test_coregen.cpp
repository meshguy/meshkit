/** \file test_coregen.cpp
 *
 * Test CoreGen
 *
 */

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/CoreGen.hpp"
#include "TestUtil.hpp"
#include "mpi.h"
using namespace MeshKit;


void test_coregen_default(int argc, char **argv);

int main(int argc, char *argv[])
{
  test_coregen_default(argc, argv);
  return 0;
}

void test_coregen_default(int argc, char **argv)
{
  // create a model entity vector for construting assygen meshop, note that NO model entities are required for assygen meshop.
  MEntVector volso;
  int nrank = 0, nprocs = 1;
  //Initialize MPI
#ifdef HAVE_PARALLEL_MOAB
  MPI::Init(argc, argv);
  nprocs = MPI::COMM_WORLD.Get_size();
  nrank = MPI::COMM_WORLD.Get_rank();
#endif
  MKCore *mk;
  mk = new MKCore();
  // construct the meshop and set name
  CoreGen *cg = (CoreGen*) mk->construct_meshop("CoreGen", volso);
  cg->set_name("coregen");

  // setup input/output files for creating the 'Reactor Core' model

  cg->prepareIO(argc, argv, nrank, nprocs, TestDir);
  mk->setup_and_execute();
  mk->save_mesh("coregen_t.exo");

#ifdef HAVE_PARALLEL_MOAB
  MPI::COMM_WORLD.Barrier();
  MPI::Finalize();
#endif
  delete cg;
  delete mk;
}



