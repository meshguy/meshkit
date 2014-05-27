

#include "meshkit/MKCore.hpp"
#include "meshkit/SolidSurfaceMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "moab/GeomTopoTool.hpp"


namespace MeshKit
{

// Output mesh types for this class
moab::EntityType SolidSurfaceMesher_types[] = {moab::MBVERTEX, moab::MBEDGE, moab::MBTRI, moab::MBMAXTYPE};
const moab::EntityType* SolidSurfaceMesher::output_types()
{ return SolidSurfaceMesher_types; }

SolidSurfaceMesher::SolidSurfaceMesher(MKCore *mk_core, const MEntVector &ments)
  : MeshScheme(mk_core, ments)
{
}

SolidSurfaceMesher::~SolidSurfaceMesher()
{
}

void SolidSurfaceMesher::setup_this()
{
}

void SolidSurfaceMesher::execute_this()
{
}

}
