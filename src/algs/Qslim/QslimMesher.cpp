#include "meshkit/QslimMesher.hpp"
#include "meshkit/QslimOptions.hpp"
#include "QslimDecimation.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "meshkit/RegisterMeshOp.hpp"
#include "moab/ReadUtilIface.hpp"
#include <vector>
#include <math.h>

namespace MeshKit
{

//---------------------------------------------------------------------------//
//Entity Type initialization for edge meshing
moab::EntityType QslimMesher_tps[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBMAXTYPE};
const moab::EntityType* QslimMesher::output_types()
  { return QslimMesher_tps; }

//---------------------------------------------------------------------------//
// Construction Function for Qslim Mesher
QslimMesher::QslimMesher(MKCore *mk_core, const MEntVector &me_vec)
        : MeshScheme(mk_core, me_vec)
{

}
//---------------------------------------------------------------------------//
// setup function: set up the number of intervals for edge meshing through the 
// sizing function
void QslimMesher::setup_this()
{

  // make decimate mesher
  // will take only one input, a set from ModelEnts ?
  ModelEnt * me = mentSelection.begin()->first;
  moab::EntityHandle set=me->mesh_handle(); // this is the initial set
  _worker = new MeshKit::QslimDecimation(mk_core()->moab_instance(), set);
}

//---------------------------------------------------------------------------//
// execute function: Generate the mesh for decimation
void QslimMesher::execute_this()
{
  // put the result in the range, eventually
	_worker->decimate(_opts, mentSelection.begin()->second);
	// this will modify the initialSet !!
}

//---------------------------------------------------------------------------//
// Deconstruction function
QslimMesher::~QslimMesher()
{
  delete _worker;
}


} // MeshKit namespace

