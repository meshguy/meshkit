

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include <iGeom.h>
#include <iMesh.h>
#include <set>
#include <iRel.h>
#include <vector>
#include "meshkit/MKCore.hpp"

namespace MeshKit
{
 
class CurveMesher
{
public: 
CurveMesher(MKCore *mk, const MEntVector &ments, double faceting_tolerance, double geom_resabs);

  ~CurveMesher();
private:
  double facet_tol;


};
  
}
