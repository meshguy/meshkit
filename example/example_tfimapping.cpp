/*!
\example example_tfimapping.cpp

\section tfimapping_cpp_title TFI Mapping

\subsection tfimapping_cpp_in Input
\image html tfimapping_in.jpg

\subsection tfimapping_cpp_out Output
\image html tfimapping_out.jpg

\subsection tfimapping_cpp_inf Misc. Information
\author Brett Rhodes
\date 9-12-2013
\warning Possible changes when Interval assignment is integrated

\subsection tfimapping_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/TFIMapping.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "example_utils.hpp"

using namespace MeshKit;

const int NUM_INTERVALS = 3; // we want 3 intervals on each side (for no reason)
const int INTERVAL_SIZE = -1; // in sizing functions, -1 means not specified
const bool save_mesh = true;

#ifdef HAVE_ACIS
string extension = ".sat";
#elif HAVE_OCC
string extension = ".stp";
#endif

int main(int argc, char **argv)
{
  MKCore * mk;         // handle for the instance of MeshKit
  MEntVector curves;   // handle for the curve we need to retrieve, is a vector
  MEntVector surfaces; // handle for the surface we need to retrieve, is a vector
  EdgeMesher * em;     // handle for a MeshOp that helps fulfill a pre-req for TFIMapping
  TFIMapping * tfi;    // handle for our TFIMapping MeshOp

// Prepare MK
  mk = new MKCore(); // Start up MK
  mk->load_geometry( (example_dir + string("rectangle") + extension).c_str() ); // Load the geometry

// Prepare EdgeMesher (TFIMapping requires that 1 edge be meshed)
  mk->get_entities_by_dimension(1, curves); // get all 1D entities and store into "curves"
  curves.resize(1); // We need to mesh one curve for TFI to suceed
  em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves); // create an EdgeMesher to mesh 1 side

// Prepare TFIMapping
  mk->get_entities_by_dimension(2, surfaces); // get all 2D entities and store into "surfaces" (we only have 1)
  tfi = (TFIMapping*) mk->construct_meshop("TFIMapping", surfaces); // create the TFIMapping MeshOp instance, will operate on the entities stored in surfaces
  mk->get_graph().addArc(em->get_node(), tfi->get_node()); // TFIMapping depends on EdgeMesher (tfi needs a meshed edge)

// Specify Sizes
  SizingFunction sf(mk, NUM_INTERVALS, INTERVAL_SIZE); // create a sizing function
  surfaces[0]->sizing_function_index(sf.core_index()); // apply the same sizing function to them

// Execute
  mk->setup(); // calls setup_this() on all nodes in the graph
  mk->execute(); // calls execute_this() on all nodes in the graph

// Save
  if (save_mesh)
    mk->save_mesh("tfimapping_out.exo");

  return 0;
}
