#ifndef __MESHIMPROVE_HPP
#define __MESHIMPROVE_HPP

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits.h>

#include "MKVersion.h"
#include "meshkit/MKCore.hpp"
#include "meshkit/iMesh.hpp"
#include "meshkit/MeshOp.hpp"
#include "meshkit/ModelEnt.hpp"
#include <set>
#include <vector>
#include <map>

//#include <MsqError.hpp>
//#include <ShapeImprovementWrapper.hpp>
//#include <MsqIMesh.hpp>
//#include <MsqIGeom.hpp>

using namespace std;

namespace MeshKit {

class MeshImprove
{	
public:
	//public function
	MeshImprove(MKCore* core, bool isLaplacian = false, bool isUntangle = true, bool isShapeImprove = true, bool isSizeAdapt = true, iGeom * geom_inst = NULL);
	~MeshImprove();
	void SurfMeshImprove(iBase_EntityHandle surface, iBase_EntitySetHandle surfMesh, iBase_EntityType entity_type);
	void VolumeMeshImprove(iBase_EntitySetHandle volMesh, iBase_EntityType entity_type);
private:
	//member variables
	MKCore* mk_core;
//	double max;
//	iBase_EntitySetHandle mesh_root_set;
//	iBase_TagHandle mesh_id_tag;
	bool IsLaplacian, IsUntangle, IsShapeImprove, IsSizeAdapt;
	iGeom * igeom_inst;
};

}
#endif

