

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

  SolidSurfaceMesher::SolidSurfaceMesher(MKCore *mk_core, const MEntVector &ments)
    : MeshScheme(mk_core, ments)
  {
  }

}
