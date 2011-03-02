/*
 * primitives.h
 *
 *  Created on: Mar 22, 2010
 *      Author: iulian
 */

#ifndef PRIMITIVES_H_
#define PRIMITIVES_H_

#include "std.h"
#include "3D.h"
#include "Mat4.h"
#include "moab/Interface.hpp"
// this is for moab::ErrorCode ?

#include "meshkit/QslimOptions.hpp"
#include <vector>

extern int validFaceCount; // defined in QslimDecimation
extern moab::Tag validTag;   // defined in QslimDecimation
// the plane data will be stored for each triangle, and recomputed for each triangle
// it will be needed only for the -m option active (so do not recompute if you do not need it)
extern moab::Tag planeDataTag;  // defined in QslimDecimation
extern int ehIsValid(moab::EntityHandle v);  // defined in QslimDecimation; maybe we should pass Interface too
extern QslimOptions opts;

extern void filterValid(moab::Interface * mb, std::vector<moab::EntityHandle> & io);
extern int classifyVertex(moab::Interface * mb, moab::EntityHandle v);

extern moab::ErrorCode contractionRegion(moab::Interface * mb, moab::EntityHandle v1, moab::EntityHandle v2,
		std::vector<moab::EntityHandle> & changed);

extern Vec3 getVec3FromMBVertex(moab::Interface * mb, moab::EntityHandle v);

// this will just retrieve it from the tag data
extern Plane trianglePlane(moab::Interface * mb, moab::EntityHandle tri);

extern void computeTrianglePlane (moab::Interface * mb, moab::EntityHandle tri);

extern moab::ErrorCode contract (moab::Interface * mb, moab::EntityHandle v0, moab::EntityHandle v1, Vec3 & vnew, std::vector<moab::EntityHandle> & changed );

#endif /* PRIMITIVES_H_ */
