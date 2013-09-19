/*!
 * \file test_coregen.cpp
 * \test This file contains 1 test. Algorithms tested include: CoreGen (Internally CopyMesh/CopyGeom).
 *
 * It's more a mesh application that uses other mesh algorithms.
 * This program passes command line arguments: options and CoreGen input filename
 * CoreGen MeshOp is initialized and setup and execute is called to create the core model
     * Steps in CoreGen MeshOp are:
     * Read input file (on all processors)
     * Read mesh or geometry files (load parallel if requested)
     * If parallel divide copy/move work
     * Loop through number of occurences and copy/move
     * Merge
     * Save
     * Write Makefile
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
    // serial process
    int nrank = 0, nprocs =1;

    //Initialize MPI
#ifdef HAVE_PARALLEL_MOAB
    MPI::Init(argc, argv);
    nprocs = MPI::COMM_WORLD.Get_size();
    nrank = MPI::COMM_WORLD.Get_rank();
#endif

    MKCore *mk;
    mk = new MKCore();
    // create a model entity vector for construting coregen meshop, note that NO model entities are required.
    MEntVector nullMEntVec;

    // construct the meshop and set name
    CoreGen *cg = (CoreGen*) mk->construct_meshop("CoreGen", nullMEntVec);
    cg->set_name("coregen");

    // setup input/output files for creating the 'Reactor Core' model
    cg->prepareIO(argc, argv, nrank, nprocs, TestDir);
    mk->setup_and_execute();
    mk->save_mesh("cgd.h5m");
#ifdef HAVE_PARALLEL_MOAB
    MPI::COMM_WORLD.Barrier();
    MPI::Finalize();
#endif

    delete mk;
}



