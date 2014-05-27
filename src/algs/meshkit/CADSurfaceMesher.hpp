

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
#include "meshkit/MeshScheme.hpp"


namespace MeshKit
{

  using namespace std; 

  class CADSurfaceMesher : public MeshScheme
  {
   public:
     CADSurfaceMesher(MKCore *mk, const MEntVector &ments);
   
  };

}
