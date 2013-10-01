/*!
\example example_ngtetmesher.cpp

\section example_ngtetmesher_cpp_title NetGen Mesh Operation

\subsection example_ngtetmesher_cpp_in Input
Loads a geometry file threeholecube
\subsection example_ngtetmesher_cpp_out Output
saves a vtk file the same name.
\subsection example_ngtetmesher_cpp_inf Misc. Information
It internally needs a tri mesher.

\subsection example_ngtetmesher_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/NGTetMesher.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/ModelEnt.hpp"

using namespace MeshKit;

MKCore *mk = NULL;
bool save_mesh = true;

#ifdef HAVE_ACIS
std::string extension = ".sat";
#elif HAVE_OCC
std::string extension = ".stp";
#endif

int main(int argc, char **argv) 
{

    // start up MK and load the geometry
    mk = new MKCore;
    std::string filebase = "threeholecube";
    std::string file_name = (std::string) MESH_DIR + "/" + filebase + extension;
    mk->load_geometry(file_name.c_str());

    // get the volume
    MEntVector dum, vols;
    mk->get_entities_by_dimension(3, dum);
    vols.push_back(*dum.rbegin());
    (NGTetMesher*) mk->construct_meshop("NGTetMesher", vols);

    // make a sizing function and set it on the surface
    SizingFunction esize(mk, -1, 0.25);
    vols[0]->sizing_function_index(esize.core_index());

    // mesh the surface, by calling execute
    mk->setup_and_execute();

    // report the number of tets
    moab::Range tets;
    mk->moab_instance()->get_entities_by_dimension(0, 3, tets);
    std::cout << tets.size() << " tets generated." << std::endl;

    if (save_mesh) {
        // output mesh
        std::string outfile = filebase + std::string(".vtk");
        moab::EntityHandle out_set = vols[0]->mesh_handle();
        mk->moab_instance()->write_file(outfile.c_str(), NULL, NULL, &out_set, 1);
    }

    return 0;
}

