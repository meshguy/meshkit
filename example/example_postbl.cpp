/*!
\example example_postbl.cpp

\section example_postbl_cpp_title Post-Mesh Boundary Layer Generation Tool

\subsection example_postbl_cpp_in Input
Reads in a mesh file and boundary layer parameters.

\subsection example_postbl_cpp_out Output
Outputs a mesh file.

\subsection example_postbl_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/PostBL.hpp"

#include "meshkit/SizingFunction.hpp"
#include "meshkit/CopyGeom.hpp"
using namespace MeshKit;

MKCore *mk;

void test_PostBL_default(int argc, char **argv);

int main(int argc, char *argv[])
{
    mk = new MKCore();
    //! run the default test case
    test_PostBL_default(argc, argv);
    delete mk;
    return 0;
}

void test_PostBL_default(int argc, char **argv)
{
    //! create a model entity vector for construting PostBL meshop, note that model entities(mesh) input for PostBL meshop is read from a file.
    MEntVector volso;

    //! construct the meshop and set name
    PostBL *pbl = (PostBL*) mk->construct_meshop("PostBL", volso);
    pbl->set_name("PostBL");

    //!setup and execute PostBL graph node, point the executable to PostBL input file,
    pbl->PrepareIO(argc, argv, std::string(MESH_DIR));
    pbl->setup_this();
    pbl->execute_this();
    mk->save_mesh("out_postbl.exo");
    delete pbl;
}



