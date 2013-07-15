/*!
\example example_meshoptemplate.cpp

\section example_meshoptemplate_cpp_title <pretty-name-of-this-file>

\subsection example_meshoptemplate_cpp_inf Misc. Information
\author <your-name-here>
\date 7-15-2013
\bug <placeholder>
\warning <placeholder>

\subsection example_meshoptemplate_cpp_goal Goal

\subsection example_meshoptemplate_cpp_cw Code Walkthrough

\subsection example_meshoptemplate_cpp_in Input
\image html example_meshoptemplate.in.jpg
There is no input.

\subsection example_meshoptemplate_cpp_out Output
\image html example_meshoptemplate.out.jpg

\subsection example_meshoptemplate_cpp_src Source Code
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

bool save_mesh = false;
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
  // Make the brick!                                                           // v Don't need this
  MeshOpTemplate *mot = (MeshOpTemplate*) mk->construct_meshop("MeshOpTemplate", MEntVector());
  mot->set_name("MeshOpTemplate");

  // Do work!
  mk->setup();
  mk->execute();
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

