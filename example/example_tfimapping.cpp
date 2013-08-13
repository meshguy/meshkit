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

const int NUM_INTERVALS = -1; // on our curve, we want 10 intervals
const int INTERVAL_SIZE = 0.5; // in sizing functions, -1 means not specified
const bool save_mesh = false;

#ifdef HAVE_ACIS
string extension = ".sat";
#elif HAVE_OCC
string extension = ".stp";
#endif

int main(int argc, char **argv)
{
  MKCore * mk;
  MEntVector curves;
  MEntVector surfaces;
  EdgeMesher * em;
  TFIMapping * tfi;

// Prepare MK
  mk = new MKCore(); // Start up MK
  mk->load_geometry( (example_dir + string("rectangle") + extension).c_str() ); // Load the geometry and mesh

// Prepare EdgeMesher (TFIMapping requires that 1 edge be meshed)
  mk->get_entities_by_dimension(1, curves); // get all 1D entites and store into "curves" (we need to mesh side of the rectangle in order for TFI to suceed)
  curves.resize(1); // We only need to mesh one curve
  em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves); // create an EdgeMesher to mesh 1 side

// Prepare TFIMapping
  mk->get_entities_by_dimension(2, surfaces); // get all 2D entities and store into "surfaces" (we only have 1)
  tfi = (TFIMapping*) mk->construct_meshop("TFIMapping", surfaces); // create the TFIMapping MeshOp instance, will operate on the entities stored in surfaces
  mk->get_graph().addArc(em->get_node(), tfi->get_node()); // TFIMapping depends on EdgeMesher (tfi needs a meshed edge)

// Specify Sizes
  SizingFunction sf(mk, NUM_INTERVALS, INTERVAL_SIZE); // create a sizing function
  surfaces[0]->sizing_function_index(sf.core_index()); // and apply it to our curve (we know curves[0] is our curve, because we only have 1)

// Execute
  mk->setup(); // calls setup_this() on all nodes in the graph
  mk->execute(); // calls execute_this() on all nodes in the graph

// Save
  if (save_mesh)
    mk->save_mesh("tfimapping_out.exo");

  return 0;
}
