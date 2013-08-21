/*!
\example example_postbl.cpp

\section example_postbl_cpp_title <pretty-name-of-this-file>

\subsection example_postbl_cpp_in Input
\image html example_postbl.in.jpg "(description of image)"
There is no input.

\subsection example_postbl_cpp_out Output
\image html example_postbl.out.jpg "(description of image)"

\subsection example_postbl_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection example_postbl_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/PostBL.hpp"

#include "meshkit/CAMALTetMesher.hpp"
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
  PostBL *ag = (PostBL*) mk->construct_meshop("PostBL", volso);
   ag->set_name("PostBL");

  //!setup and execute PostBL graph node, point the executable to PostBL input file, 
  //ag->PrepareIO(argc, argv, TestDir);
  ag->setup_this();
  ag->execute_this();

  delete ag;
}



