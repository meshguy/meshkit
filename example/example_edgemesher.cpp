/*!
\example example_edgemesher.cpp

\section edgemesher_cpp_title Edge Mesher

\subsection edgemesher_cpp_in Input
\image html edgemesher_in.jpg "a 1D line"

\subsection edgemesher_cpp_out Output
\image html edgemesher_out.jpg "a meshed 1D line"

\subsection edgemesher_cpp_inf Misc. Information
\author Brett Rhodes
\date 9-12-2013

\subsection edgemesher_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/EdgeMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/Matrix.hpp"

using namespace MeshKit;

const int NUM_INTERVALS = 10; // on our curve, we want 10 intervals
const int INTERVAL_SIZE = -1; // in sizing functions, -1 means not specified
const bool save_mesh = true;

#ifdef HAVE_ACIS
std::string extension = ".sat";
#elif HAVE_OCC
std::string extension = ".stp";
#endif

int main(int argc, char **argv) 
{
  MKCore * mk;       // handle for the instance of MeshKit
  MEntVector curves; // handle for the curve we need to retrieve, is a vector
  EdgeMesher * em;   // handle for our MeshOp that we will create

// Prepare MK
  mk = new MKCore();  // Start up MK
  mk->load_geometry( (string(MESH_DIR) + "/" + string("spline") + extension).c_str() ); // Load the geometry

// Prepare EdgeMesher
  mk->get_entities_by_dimension(1, curves); // get all 1D entites and store into "curves" (we only have 1)
  em = (EdgeMesher*) mk->construct_meshop("EdgeMesher", curves); // create the EdgeMesher MeshOp instance, will operate on the entities stored in curves
  SizingFunction sf(mk, NUM_INTERVALS, INTERVAL_SIZE); // create a sizing function
  curves[0]->sizing_function_index(sf.core_index()); // and apply it to our curve (we know curves[0] is our curve, because we only have 1)

// Execute
  mk->setup(); // calls setup_this() on all nodes in the graph
  mk->execute(); // calls execute_this() on all nodes in the graph

// Save stuff
  if (save_mesh)
    mk->save_mesh("edgemesher_out.vtk"); // save the meshed file for the curve we just meshed

  return 0;
}
