/*!
\example example_coregen.cpp

\section example_CoreGen_cpp_title Assembly Gen

Some text explaining objectives

\subsection example_CoreGen_cpp_in Input
\image html CoreGen.in.jpg "(description of image)"

\subsection example_CoreGen_cpp_out Output
\image html CoreGen.out.jpg "(description of image)"

\subsection example_CoreGen_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning Currently requires multiple setup/execute cycles.

\subsection example_CoreGen_cpp_src Source Code
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



