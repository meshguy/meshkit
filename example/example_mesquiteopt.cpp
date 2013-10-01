/*!
\example example_mesquiteopt.cpp

\section example_mesquiteopt_cpp_title MesquiteOpt: Mesh Quality Optimization

\subsection example_mesquiteopt_cpp_in Input
cropSE.h5m from data directory

\subsection example_mesquiteopt_cpp_out Output
VTK file o/p after using "smooth with fixed boundary option"

\subsection example_mesquiteopt_cpp_src Source Code
*/

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/MesquiteOpt.hpp"
#include "meshkit/iGeom.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/Matrix.hpp"


using namespace MeshKit;
MKCore* core;
int main()
{

    MKCore my_core;
    core = &my_core;

    MEntVector v;
    std::string file_name = (std::string) MESH_DIR + "/" + "cropSE.h5m";
    core->load_mesh(file_name.c_str());
    core->get_entities_by_dimension(3,v);
    MesquiteOpt op(core,v);
    //op.smooth_with_free_boundary();
    op.smooth_with_fixed_boundary();

    op.setup_this();
    op.execute_this();

    core->save_mesh("load_and_smooth.vtk");
    return 0;
}
