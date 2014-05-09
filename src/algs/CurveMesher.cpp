
#include "meshkit/MKCore.hpp"
#include "meshkit/CurveMesher.hpp"
#include <iostream>
#include <iGeom.h>

namespace MeshKit
{
// Construction Function for CurveMesher



CurveMesher::CurveMesher(MKCore *mk_core, const MEntVector &ments, double faceting_tolerance, double geom_resabs)
{
  facet_tol = faceting_tolerance;
}

// Destructor Function for Curvemesher
CurveMesher::~CurveMesher()
{
}



}
