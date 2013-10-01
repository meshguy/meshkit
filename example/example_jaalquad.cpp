/*!
\example example_jaalquad.cpp

\section example_jaalquad_cpp_title QuadMesher Algorithm in MeshKit

\subsection example_jaalquad_cpp_in Input
Reads in a quadface geometry from data folder.

\subsection example_jaalquad_cpp_out Output
QuadMesh is output. Algorithms needs a tri mesher (tri meshing node is created during setup of QuadMesher).
Currently works with only CAMALPaver.
Note sometimes linking CUBIT and CAMAL together might cause problems.
\warning Some tri's are left after quad mesher is finished.
\subsection example_jaalquad_cpp_src Source Code
*/
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/Matrix.hpp"
#include "meshkit/SizingFunction.hpp"
#include "meshkit/QuadMesh.hpp"

#ifdef HAVE_ACIS
#define TEST_QUADFACE "quadface.sat"
#elif defined(HAVE_OCC)
#define TEST_QUADFACE "quadface.stp"
#endif

using namespace MeshKit;

MKCore* core = 0;
int main()
{
    core = new MKCore(); // Start up MK
    moab::Interface* moab = core->moab_instance();
    std::string filename = (std::string) MESH_DIR + "/" + TEST_QUADFACE;
    core->load_geometry(filename.c_str());
    core->populate_model_ents(0, -1, -1);

    // get the tris
    MEntVector tris, dum;
    core->get_entities_by_dimension(2, dum);
    tris.push_back(*dum.rbegin());

    // run tri to quad mesher
    (MeshKit::QuadMesher*) core->construct_meshop( "QuadMesher", tris );
    double size = 1.1;
    SizingFunction esize(core, -1, size);
    tris[0]->sizing_function_index(esize.core_index());

    core->setup_and_execute();

    // removing tri's
    moab::Range tri;
    moab->get_entities_by_type( 0, moab::MBTRI, tri );
    if(tri.size() != 0)
        moab->delete_entities(tri);

    core->save_mesh("jaalquad.vtk");
    return 0;
}
