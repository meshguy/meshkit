/*!
\example example_assygen.cpp

\section example_assygen_cpp_title Assembly Geometry and Mesh Script Generator
\subsection example_assygen_cpp_in Input
Keyword-based text input file describing the assembly geometry.\n
File can be found in data folder: assygen_default.inp
\subsection example_assygen_cpp_out Output
\image html assygen_out.jpg "Output toy assembly geometry"
\subsection example_assygen_cpp_inf Misc. Information
\author Rajeev Jain
\date 09-30-2013
\warning The old design is more stable. Can be found: MeshKit/rgg/assygen

\subsection example_assygen_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/AssyGen.hpp"
#include "meshkit/CopyGeom.hpp"

using namespace MeshKit;

MKCore *mk;
bool save_mesh = true;

int main(int argc, char *argv[])
{
  mk = new MKCore();
  // create a model entity vector for construting AssyGen meshop, note that NO model entities are required for AssyGen meshop.
  MEntVector volso;

  // construct the meshop and set name
  AssyGen *ag = (AssyGen*) mk->construct_meshop("AssyGen", volso);
  ag->set_name("AssyGen");

  // setup input/output AssyGen files for creating the 'Reactor Assembly' geometry
  ag->PrepareIO(argc, argv, MESH_DIR);
  mk->setup_and_execute();

  // AssyGen MeshOp internally saves the resulting geometry with the name specified in input file

  delete mk;
  return 0;
}

