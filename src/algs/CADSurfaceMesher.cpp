

#include "meshkit/MKCore.hpp"
#include "meshkit/CADSurfaceMesher.hpp"
#include "meshkit/ModelEnt.hpp"
#include <iostream>
#include <iGeom.h>

#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "moab/GeomTopoTool.hpp"


namespace MeshKit
{

  CADSurfaceMesher::CADSurfaceMesher(MKCore *mk_core, const MEntVector &ments)
    : MeshScheme(mk_core, ments)
  {
  }

}
